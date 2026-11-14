(ns fpnum.optim-test
  "Pruebas de fpnum.optim: gradiente por AD frente al analítico, descenso de
  gradiente con paso fijo y con backtracking, y mínimos cuadrados frente a las
  ecuaciones normales."
  (:require [clojure.test :refer [deftest is testing]]
            [clojure.test.check.clojure-test :refer [defspec]]
            [clojure.test.check.generators :as gen]
            [clojure.test.check.properties :as prop]
            [fpnum.ad.forward :as ad]
            [fpnum.linear :as la]
            [fpnum.optim :as opt]
            [fpnum.svd.core :as svd]))

(defn- approx-num= [a b tol]
  (<= (Math/abs (- (double a) (double b)))
      (* tol (+ 1.0 (Math/abs (double a)) (Math/abs (double b))))))

(defn- approx-vec= [x y tol]
  (and (= (count x) (count y))
       (every? true? (map #(approx-num= %1 %2 tol) x y))))

(def ^:private gen-num
  (gen/double* {:min -5.0 :max 5.0 :NaN? false :infinite? false}))

(defn- gen-problema
  "Problema de mínimos cuadrados aleatorio: [X y w] con X de m×n (m ≥ 3)."
  []
  (gen/bind (gen/tuple (gen/choose 3 8) (gen/choose 1 4))
            (fn [[m n]]
              (gen/tuple (gen/vector (gen/vector gen-num n) m)
                         (gen/vector gen-num m)
                         (gen/vector gen-num n)))))

;; ---------------------------------------------------------------------------
;; Funciones objetivo y gradientes
;; ---------------------------------------------------------------------------

(deftest rosenbrock-en-el-minimo
  (is (= 0.0 (:a (ad/lift (opt/rosenbrock (mapv ad/lift [1.0 1.0]))))))
  (is (= [0.0 0.0] (ad/grad opt/rosenbrock [1.0 1.0])))
  (is (= [0.0 0.0] (opt/rosenbrock-gradient [1.0 1.0]))))

(defspec rosenbrock-ad-coincide-con-gradiente-analitico 200
  (prop/for-all [x (gen/double* {:min -2.0 :max 2.0 :NaN? false :infinite? false})
                 y (gen/double* {:min -2.0 :max 2.0 :NaN? false :infinite? false})]
    (approx-vec= (ad/grad opt/rosenbrock [x y])
                 (opt/rosenbrock-gradient [x y])
                 1e-9)))

(defspec objetivo-minimos-cuadrados-es-norma-al-cuadrado 200
  (prop/for-all [[x y w] (gen-problema)]
    (let [f (opt/least-squares-objective x y)
          residuo (la/v- (la/mat-vec x w) y)]
      (approx-num= (:a (ad/lift (f (mapv ad/lift w))))
                   (la/dot residuo residuo)
                   1e-9))))

(defspec gradiente-minimos-cuadrados-ad-coincide-con-analitico 200
  (prop/for-all [[x y w] (gen-problema)]
    (approx-vec= (ad/grad (opt/least-squares-objective x y) w)
                 (opt/least-squares-gradient x y w)
                 1e-9)))

;; ---------------------------------------------------------------------------
;; Datos ruidosos y ecuaciones normales
;; ---------------------------------------------------------------------------

(deftest datos-reproducibles
  (let [params {:n 20 :w-true [1.0 -2.0] :sigma 0.5 :seed 7}
        a (opt/noisy-regression-data params)
        b (opt/noisy-regression-data params)
        c (opt/noisy-regression-data (assoc params :seed 8))]
    (testing "misma semilla, mismos datos"
      (is (= a b)))
    (testing "otra semilla, otros datos"
      (is (not= (:x a) (:x c))))
    (testing "dimensiones"
      (is (= [20 2] (la/dims (:x a))))
      (is (= 20 (count (:y a))))
      (is (= [1.0 -2.0] (:w-true a))))))

(deftest sin-ruido-las-ecuaciones-normales-recuperan-w
  (let [{:keys [x y w-true]} (opt/noisy-regression-data
                              {:n 30 :w-true [2.0 -1.0 0.5] :sigma 0.0 :seed 1})]
    (is (approx-vec= w-true (opt/normal-equations x y) 1e-8))))

(deftest ecuaciones-normales-con-columnas-dependientes
  (testing "columnas duplicadas → XᵀX singular → error claro"
    (let [x [[1.0 1.0] [2.0 2.0] [3.0 3.0]]
          y [1.0 2.0 3.0]]
      (is (thrown-with-msg? clojure.lang.ExceptionInfo #"singular"
                            (opt/normal-equations x y))))))

(deftest paso-seguro
  (let [{:keys [x]} (opt/noisy-regression-data
                     {:n 50 :w-true [1.0 2.0 3.0] :sigma 0.1 :seed 3})
        nf (la/norm-fro x)
        alpha (opt/safe-step x)]
    (is (pos? alpha))
    (is (< alpha (/ 1.0 (* nf nf))) "α < 1/‖X‖²_F")
    (testing "es menor que 2/L, con L = 2·λmax(XᵀX) ≤ 2‖X‖²_F"
      (is (< alpha (/ 2.0 (* 2.0 nf nf)))))))

;; ---------------------------------------------------------------------------
;; Descenso con paso fijo
;; ---------------------------------------------------------------------------

(deftest paso-fijo-en-cuadratica-simple
  (let [f (fn [[x]] (ad/d* (ad/d- x 3) (ad/d- x 3)))   ; (x − 3)², mínimo en 3
        r (opt/gradient-descent f [0.0] {:alpha 0.1 :max-iter 1000 :tol 1e-10})]
    (is (:converged? r))
    (is (= :converged (:reason r)))
    (is (approx-num= 3.0 (first (:x r)) 1e-8))))

(deftest razones-de-parada
  (let [f (fn [[x]] (ad/d* x x))]
    (testing "max-iter"
      (let [r (opt/gradient-descent f [10.0] {:alpha 0.001 :max-iter 5})]
        (is (= :max-iter (:reason r)))
        (is (not (:converged? r)))
        (is (= 5 (:iterations r)))))
    (testing "paso demasiado grande: diverge (α = 1.5 duplica |x| en cada paso)"
      (let [r (opt/gradient-descent f [1.0] {:alpha 1.5 :max-iter 100000})]
        (is (= :diverged (:reason r)))
        (is (not (:converged? r)))))
    (testing "punto inicial ya óptimo: converge en 0 iteraciones"
      (let [r (opt/gradient-descent f [0.0] {:alpha 0.1})]
        (is (= :converged (:reason r)))
        (is (= 0 (:iterations r)))))))

(deftest trayectoria-bien-formada
  (let [f (fn [[x]] (ad/d* x x))
        r (opt/gradient-descent f [4.0] {:alpha 0.25 :max-iter 20 :tol 1e-12})
        tray (:trajectory r)]
    (testing "una fila por punto visitado"
      (is (= (inc (:iterations r)) (count tray)))
      (is (= (range (count tray)) (map :iter tray))))
    (testing "la primera fila es el punto inicial"
      (is (= [4.0] (:x (first tray)))))
    (testing "α registrado en cada fila, y nil en la última"
      (is (every? #(= 0.25 (:alpha %)) (butlast tray)))
      (is (nil? (:alpha (last tray)))))
    (testing "la última fila coincide con el resultado"
      (is (= (:x r) (:x (last tray))))
      (is (= (:f r) (:f (last tray)))))))

;; ---------------------------------------------------------------------------
;; Mínimos cuadrados con paso seguro
;; ---------------------------------------------------------------------------

(deftest minimos-cuadrados-por-descenso-coincide-con-ecuaciones-normales
  (let [{:keys [x y]} (opt/noisy-regression-data
                       {:n 100 :w-true [2.0 -1.0 0.5] :sigma 0.1 :seed 42})
        alpha (opt/safe-step x)
        r (opt/gradient-descent (opt/least-squares-objective x y) (la/zeros 3)
                                {:alpha alpha :max-iter 5000 :tol 1e-8})
        w-ne (opt/normal-equations x y)
        tray (:trajectory r)]
    (testing "converge"
      (is (:converged? r)))
    (testing "coincide con las ecuaciones normales"
      (is (approx-vec= w-ne (:x r) 1e-6)))
    (testing "con el paso seguro f decrece en cada iteración (salvo redondeo)"
      ;; Cerca del mínimo f ≈ 0.99 y los cambios son de ~1e-16, del orden del
      ;; error de redondeo, así que f puede subir unas pocas unidades de eps.
      (is (every? true? (map (fn [a b] (<= (:f b) (+ (:f a) 1e-12))) tray (rest tray))))
      (testing "y mientras el gradiente no es diminuto, baja de verdad"
        (is (->> (map vector tray (rest tray))
                 (filter (fn [[a _]] (> (:grad-norm a) 1e-4)))
                 (every? (fn [[a b]] (< (:f b) (:f a))))))))
    (testing "el α registrado es el paso seguro"
      (is (every? #(= alpha (:alpha %)) (butlast tray))))
    (testing "el gradiente en la solución es casi cero"
      (is (< (la/norm2 (opt/least-squares-gradient x y (:x r))) 1e-6)))))

;; ---------------------------------------------------------------------------
;; Backtracking
;; ---------------------------------------------------------------------------

(def ^:private rosenbrock-result
  (delay (opt/gradient-descent opt/rosenbrock [-1.2 1.0]
                               {:backtracking? true :max-iter 100000 :tol 1e-6})))

(deftest backtracking-converge-cerca-de-1-1
  (let [r @rosenbrock-result]
    (is (:converged? r))
    (is (= :converged (:reason r)))
    (is (approx-vec= [1.0 1.0] (:x r) 1e-4))
    (is (< (:f r) 1e-8))))

(deftest backtracking-f-baja-en-cada-paso-aceptado
  (let [tray (:trajectory @rosenbrock-result)
        c 1e-4]
    (testing "f no crece en ningún paso"
      (is (every? true? (map (fn [a b] (<= (:f b) (:f a))) tray (rest tray)))))
    (testing "cada paso cumple la condición de Armijo con el α registrado"
      (is (every? true?
                  (map (fn [a b]
                         (<= (:f b)
                             (+ (- (:f a) (* c (:alpha a) (let [g (:grad-norm a)] (* g g))))
                                1e-12)))
                       tray (rest tray)))))
    (testing "el α registrado es positivo y a lo más α0 = 1"
      (is (every? #(< 0.0 (:alpha %) 1.000001) (butlast tray))))
    (testing "la última fila no tiene α"
      (is (nil? (:alpha (last tray)))))))

(deftest backtracking-reduce-alpha-si-hace-falta
  (let [f (fn [[x]] (ad/d* x x))
        r (opt/gradient-descent f [1.0] {:backtracking? true :alpha0 100.0
                                         :max-iter 50 :tol 1e-10})]
    (testing "α0 = 100 no sirve en x²; se reduce hasta uno aceptable"
      (is (< (:alpha (first (:trajectory r))) 100.0)))
    (is (:converged? r))
    (is (approx-num= 0.0 (first (:x r)) 1e-9))))

(deftest backtracking-sin-alpha-aceptable
  (let [f (fn [[x]] (ad/d* x x))
        r (opt/gradient-descent f [1.0] {:backtracking? true :alpha0 1.0
                                         :min-alpha 2.0})]
    (is (= :line-search-failed (:reason r)))
    (is (not (:converged? r)))))

;; ---------------------------------------------------------------------------
;; Mínimos cuadrados por pseudoinversa (SVD)
;; ---------------------------------------------------------------------------

(deftest least-squares-svd-coincide-con-ecuaciones-normales
  (let [{:keys [x y]} (opt/noisy-regression-data
                       {:n 100 :w-true [2.0 -1.0 0.5] :sigma 0.1 :seed 42})]
    (is (approx-vec= (opt/normal-equations x y) (opt/least-squares-svd x y) 1e-8))))

(deftest least-squares-svd-recupera-la-solucion-exacta
  (let [{:keys [x y w-true]} (opt/noisy-regression-data
                              {:n 30 :w-true [2.0 -1.0 0.5] :sigma 0.0 :seed 1})]
    (is (approx-vec= w-true (opt/least-squares-svd x y) 1e-9))))

(deftest least-squares-svd-es-pinv-por-y
  (let [{:keys [x y]} (opt/noisy-regression-data
                       {:n 25 :w-true [1.0 2.0 3.0 4.0] :sigma 0.3 :seed 9})]
    (is (approx-vec= (la/mat-vec (svd/pinv x) y) (opt/least-squares-svd x y) 1e-9))))

(deftest least-squares-svd-gradiente-nulo-en-la-solucion
  (let [{:keys [x y]} (opt/noisy-regression-data
                       {:n 40 :w-true [1.0 -2.0 0.5] :sigma 0.2 :seed 5})
        w (opt/least-squares-svd x y)]
    (is (< (la/norm2 (opt/least-squares-gradient x y w)) 1e-9))))

(deftest least-squares-svd-rango-deficiente-da-la-solucion-de-minima-norma
  (let [g (vec (for [i (range 12)] [(Math/sin (* 1.3 i)) (Math/cos (* 0.7 i)) (/ (inc i) 12.0)]))
        x (mapv #(conj % (first %)) g)          ; 4ª columna = 1ª → rango 3
        y (mapv #(+ (* 2.0 (first %)) (second %)) g)
        w (opt/least-squares-svd x y)
        nulo [1.0 0.0 0.0 -1.0]]
    (testing "las ecuaciones normales fallan"
      (is (thrown-with-msg? clojure.lang.ExceptionInfo #"singular" (opt/normal-equations x y))))
    (testing "w es solución de mínimos cuadrados (gradiente nulo)"
      (is (< (la/norm2 (opt/least-squares-gradient x y w)) 1e-9)))
    (testing "w es ortogonal al espacio nulo → mínima norma"
      (is (< (Math/abs (la/dot w nulo)) 1e-10)))
    (testing "cualquier otra solución w + t·nulo tiene mayor norma"
      (doseq [t [-1.0 0.5 2.0]]
        (let [w2 (la/v+ w (la/scale t nulo))]
          (is (< (la/norm2 (opt/least-squares-gradient x y w2)) 1e-9) "sigue siendo solución")
          (is (> (la/norm2 w2) (la/norm2 w))))))))

(deftest least-squares-svd-dimensiones-incompatibles
  (is (thrown? clojure.lang.ExceptionInfo
               (opt/least-squares-svd [[1.0 2.0] [3.0 4.0] [5.0 6.0]] [1.0 2.0]))))

;; ---------------------------------------------------------------------------
;; :stop-when
;; ---------------------------------------------------------------------------

(deftest stop-when-detiene-con-su-propio-criterio
  (let [f (fn [[x]] (ad/d* x x))
        r (opt/gradient-descent f [10.0] {:alpha 0.1 :tol 0.0 :max-iter 1000
                                          :stop-when (fn [{:keys [x]}] (< (Math/abs (first x)) 1.0))})]
    (is (= :stopped (:reason r)))
    (is (not (:converged? r)) "no es convergencia por gradiente")
    (is (< (Math/abs (first (:x r))) 1.0))
    (testing "se detiene en la primera iteración que cumple el criterio"
      (is (every? #(>= (Math/abs (first (:x %))) 1.0) (butlast (:trajectory r)))))
    (is (= (inc (:iterations r)) (count (:trajectory r))))))

(deftest stop-when-recibe-la-fila-completa
  (let [f (fn [[x]] (ad/d* x x))
        vistas (atom [])
        _ (opt/gradient-descent f [3.0] {:alpha 0.1 :max-iter 3 :tol 0.0
                                         :stop-when (fn [fila] (swap! vistas conj fila) false)})]
    (is (= [0 1 2 3] (map :iter @vistas)))
    (is (every? #(= #{:iter :x :f :grad-norm :alpha} (set (keys %))) @vistas))
    (is (= [3.0] (:x (first @vistas))))))

(deftest stop-when-que-nunca-se-cumple-no-interfiere
  (let [f (fn [[x]] (ad/d* x x))
        r (opt/gradient-descent f [10.0] {:alpha 0.1 :max-iter 5 :tol 0.0 :stop-when (constantly false)})]
    (is (= :max-iter (:reason r)))
    (is (= 5 (:iterations r)))))

(deftest convergencia-por-gradiente-tiene-prioridad-sobre-stop-when
  (let [f (fn [[x]] (ad/d* x x))
        r (opt/gradient-descent f [0.0] {:alpha 0.1 :stop-when (constantly true)})]
    (is (= :converged (:reason r)))))
