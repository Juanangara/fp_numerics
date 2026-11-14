(ns fpnum.ad.reverse-test
  "Pruebas de fpnum.ad.reverse (reverse mode) y de las operaciones genéricas.

  La referencia principal es el modo forward (ya validado contra derivadas
  conocidas y diferencias finitas): una misma función, escrita con las
  operaciones genéricas, debe dar el mismo gradiente en los dos modos. Además
  se compara con los gradientes analíticos de Rosenbrock y de mínimos cuadrados
  (fpnum.optim)."
  (:require [clojure.test :refer [deftest is testing]]
            [clojure.test.check.clojure-test :refer [defspec]]
            [clojure.test.check.generators :as gen]
            [clojure.test.check.properties :as prop]
            [fpnum.ad.forward :as fwd]
            [fpnum.ad.ops :as ad]
            [fpnum.ad.reverse :as rev]
            [fpnum.linear :as la]
            [fpnum.optim :as opt]))

;; ---------------------------------------------------------------------------
;; Utilidades
;; ---------------------------------------------------------------------------

(defn- nan? [x] (Double/isNaN (double x)))

(defn- igual?
  "Igualdad numérica con tolerancia relativa; dos NaN son iguales y dos
  infinitos del mismo signo también (los dos modos deben coincidir en los
  casos límite)."
  [a b tol]
  (let [a (double a) b (double b)]
    (cond (and (nan? a) (nan? b)) true
          (or (nan? a) (nan? b)) false
          (or (Double/isInfinite a) (Double/isInfinite b)) (= a b)
          :else (<= (Math/abs (- a b)) (* tol (+ 1.0 (Math/abs a) (Math/abs b)))))))

(defn- vec-igual? [x y tol]
  (and (= (count x) (count y)) (every? true? (map #(igual? %1 %2 tol) x y))))

(def ^:private gen-x
  (gen/double* {:min 0.5 :max 3.0 :NaN? false :infinite? false}))

;; ---------------------------------------------------------------------------
;; Operaciones genéricas (fachada ops)
;; ---------------------------------------------------------------------------

(deftest la-fachada-reexporta-las-operaciones-de-forward
  (is (identical? ad/d+ fwd/d+))
  (is (identical? ad/d* fwd/d*))
  (is (identical? ad/dsin fwd/dsin))
  (is (identical? ad/dpow fwd/dpow)))

(deftest con-numeros-y-duales-todo-sigue-igual
  (is (= (fwd/dual 3 0) (ad/d+ 1 2)))
  (is (= (fwd/dual 6 5) (ad/d* (fwd/dual 3 1) (fwd/dual 2 1))))
  (is (= (fwd/dual 0 0) (ad/d+)))
  (is (= (fwd/dual 1 0) (ad/d*))))

(deftest nodos-y-numeros-en-cualquier-orden
  (testing "número primero, nodo después"
    (is (= {:value -3.0 :grad [-1.0]} (rev/value-and-grad (fn [[x]] (ad/d- 2.0 x)) [5.0])))
    (is (= {:value 15.0 :grad [3.0]} (rev/value-and-grad (fn [[x]] (ad/d* 3 x)) [5.0])))
    (is (= {:value 0.2 :grad [-0.04]} (rev/value-and-grad (fn [[x]] (ad/d-div 1 x)) [5.0]))))
  (testing "nodo primero, número después"
    (is (= {:value 3.0 :grad [1.0]} (rev/value-and-grad (fn [[x]] (ad/d- x 2.0)) [5.0])))
    (is (= {:value 25.0 :grad [10.0]} (rev/value-and-grad (fn [[x]] (ad/dpow x 2)) [5.0]))))
  (testing "con un solo argumento"
    (is (= {:value -5.0 :grad [-1.0]} (rev/value-and-grad (fn [[x]] (ad/d- x)) [5.0])))
    (is (= {:value 5.0 :grad [1.0]} (rev/value-and-grad (fn [[x]] (ad/d+ x)) [5.0])))
    (is (= {:value 5.0 :grad [1.0]} (rev/value-and-grad (fn [[x]] (ad/d* x)) [5.0])))))

(deftest no-se-pueden-mezclar-duales-con-derivada-y-nodos
  (let [capturar (fn [f] (try (f) :sin-error (catch clojure.lang.ExceptionInfo e (ex-message e))))
        nodo (volatile! nil)
        _ (rev/grad (fn [[x]] (vreset! nodo x) x) [1.0])
        d (fwd/dual 1 1)]
    (is (re-find #"mezclar" (capturar #(ad/d+ d @nodo))))
    (is (re-find #"mezclar" (capturar #(ad/d* @nodo d))))
    (is (re-find #"mezclar" (capturar #(ad/dpow d @nodo))))))

(deftest los-duales-constantes-se-tratan-como-constantes-en-reverse
  ;; (d* 2.0 3.0) o (d+) devuelven un Dual constante (derivada 0): combinado con
  ;; un nodo es la constante que es. Esto cubre el plegado de constantes dentro
  ;; de una función que se diferencia con reverse.
  (testing "aritmética con constantes dentro de la función"
    (is (= {:value 30.0 :grad [6.0]} (rev/value-and-grad (fn [[x]] (ad/d* x (ad/d* 2.0 3.0))) [5.0])))
    (is (= {:value 30.0 :grad [6.0]} (rev/value-and-grad (fn [[x]] (ad/d* (ad/d* 2.0 3.0) x)) [5.0])))
    (is (= {:value 3.0 :grad [-1.0]} (rev/value-and-grad (fn [[x]] (ad/d- (ad/d+ 1 7) x)) [5.0])))
    (is (= {:value 0.5 :grad [-0.1]} (rev/value-and-grad (fn [[x]] (ad/d-div (ad/d* 2.5 1) x)) [5.0])) "2.5/x")
    (is (igual? 25.0 (:value (rev/value-and-grad (fn [[x]] (ad/dpow x (ad/d+ 1 1))) [5.0])) 1e-12)
        "un exponente que es un Dual constante"))
  (testing "sumas y productos vacíos"
    (is (= {:value 5.0 :grad [1.0]} (rev/value-and-grad (fn [[x]] (ad/d+ x (ad/d+))) [5.0])))
    (is (= {:value 5.0 :grad [1.0]} (rev/value-and-grad (fn [[x]] (ad/d* x (ad/d*))) [5.0])))
    (is (= {:value 5.0 :grad [1.0]} (rev/value-and-grad (fn [[x]] (apply ad/d+ x [])) [5.0]))))
  (testing "función cuyo resultado es solo un Dual constante"
    (is (= [0.0] (rev/grad (fn [_] (ad/d* 2.0 3.0)) [5.0])))))

;; ---------------------------------------------------------------------------
;; Límite de la regla heurística «Dual con derivada 0 = constante»
;; ---------------------------------------------------------------------------

(deftest la-regla-del-dual-con-derivada-cero-es-heuristica
  ;; Ver docs/decisiones.md (D2). La regla comprueba la derivada y no el origen
  ;; del Dual, así que ESTA PRUEBA FIJA SU LÍMITE conocido, no un comportamiento
  ;; deseable.
  (let [capturar (fn [f] (try (f) (catch clojure.lang.ExceptionInfo e (ex-message e))))]
    (testing "un Dual con derivada exactamente 0 se acepta como constante"
      (is (= {:value 8.0 :grad [1.0]}
             (rev/value-and-grad (fn [[x]] (ad/d+ x (fwd/dual 3 0))) [5.0]))))
    (testing "un Dual con derivada distinta de 0 es una mezcla de modos y falla"
      (is (re-find #"mezclar" (capturar #(rev/grad (fn [[x]] (ad/d+ x (fwd/dual 3 1))) [5.0])))))
    (testing "LÍMITE CONOCIDO: un Dual de forward cuya derivada casualmente se anula (x² en 0)
              también se acepta en silencio como constante, sin aviso"
      (let [dual-de-forward (fwd/dual 0.0 0.0)]    ; p. ej. el valor de x² en x = 0 al derivar en forward
        (is (= {:value 5.0 :grad [1.0]}
               (rev/value-and-grad (fn [[x]] (ad/d+ x dual-de-forward)) [5.0])))))))

;; ---------------------------------------------------------------------------
;; Derivadas conocidas
;; ---------------------------------------------------------------------------

(deftest derivadas-conocidas
  (testing "polinomios"
    (is (igual? 12.0 (rev/diff (fn [x] (ad/dpow x 3)) 2.0) 1e-12))
    (is (igual? 6.0 (rev/diff (fn [x] (ad/d+ (ad/d* x x x) (ad/d* 3 x) 2)) 1.0) 1e-12)))
  (testing "funciones elementales"
    (is (igual? 1.0 (rev/diff ad/dsin 0.0) 1e-12))
    (is (igual? (- (Math/sin 1.0)) (rev/diff ad/dcos 1.0) 1e-12))
    (is (igual? Math/E (rev/diff ad/dexp 1.0) 1e-12))
    (is (igual? 0.5 (rev/diff ad/dlog 2.0) 1e-12)))
  (testing "producto, cociente y cadena"
    (let [x 1.3]
      (is (igual? (+ (Math/sin x) (* x (Math/cos x)))
                  (rev/diff (fn [t] (ad/d* t (ad/dsin t))) x) 1e-12)))
    (let [x 2.0]
      (is (igual? (/ (- 1.0 (* x x)) (Math/pow (+ 1.0 (* x x)) 2))
                  (rev/diff (fn [t] (ad/d-div t (ad/d+ 1 (ad/d* t t)))) x) 1e-12)))
    (let [x 0.9
          s (Math/sin (* x x))]
      (is (igual? (* (Math/exp s) (Math/cos (* x x)) 2.0 x)
                  (rev/diff (fn [t] (ad/dexp (ad/dsin (ad/d* t t)))) x) 1e-12)))))

(deftest gradiente-de-dos-variables
  (let [g (fn [[x y]] (ad/d+ (ad/d* x x y) (ad/dsin y)))
        [gx gy] (rev/grad g [1.0 2.0])]
    (is (igual? 4.0 gx 1e-12))
    (is (igual? (+ 1.0 (Math/cos 2.0)) gy 1e-12))))

(deftest value-and-grad-da-el-valor
  (let [r (rev/value-and-grad (fn [[x y]] (ad/d+ (ad/d* x x y) (ad/dsin y))) [1.0 2.0])]
    (is (igual? (+ 2.0 (Math/sin 2.0)) (:value r) 1e-12))
    (is (= 2 (count (:grad r))))))

;; ---------------------------------------------------------------------------
;; Casos límite: mismo comportamiento que forward
;; ---------------------------------------------------------------------------

(deftest casos-limite-iguales-a-forward
  (doseq [[nombre f x] [["x² en 0" (fn [t] (ad/dpow t 2)) 0.0]
                        ["x³ en -1" (fn [t] (ad/dpow t 3)) -1.0]
                        ["x⁰ en 5" (fn [t] (ad/dpow t 0)) 5.0]
                        ["x⁰ en 0" (fn [t] (ad/dpow t 0)) 0.0]
                        ["x² en -2" (fn [t] (ad/dpow t 2)) -2.0]
                        ["2ˣ en 3" (fn [t] (ad/dpow 2 t)) 3.0]
                        ["xˣ en 1" (fn [t] (ad/dpow t t)) 1.0]
                        ["xˣ en 2" (fn [t] (ad/dpow t t)) 2.0]
                        ["log en 0" ad/dlog 0.0]
                        ["log en -1" ad/dlog -1.0]
                        ["1/x en 0" (fn [t] (ad/d-div 1 t)) 0.0]
                        ["1/x² en 0" (fn [t] (ad/d-div 1 (ad/d* t t))) 0.0]]]
    (testing nombre
      (is (igual? (fwd/diff f x) (rev/diff f x) 1e-12)
          (str "forward " (fwd/diff f x) " frente a reverse " (rev/diff f x))))))

(deftest funcion-constante-da-gradiente-cero
  (is (= [0.0 0.0] (rev/grad (fn [_] 7.0) [1.0 2.0])))
  (is (= [0.0 0.0] (rev/grad (fn [_] (ad/d+ 1 2)) [1.0 2.0])) "un dual constante también es constante")
  (is (= 0.0 (rev/diff (fn [_] 5.0) 3.0)))
  (testing "una variable que no se usa tiene derivada cero"
    (is (= [2.0 0.0] (rev/grad (fn [[x _]] (ad/d* x 2)) [3.0 4.0])))))

(deftest f-debe-devolver-un-nodo-o-un-numero
  (is (thrown? clojure.lang.ExceptionInfo (rev/grad (fn [_] "no es un número") [1.0])))
  (is (thrown? clojure.lang.ExceptionInfo (rev/grad (fn [_] nil) [1.0]))))

;; ---------------------------------------------------------------------------
;; Nodos compartidos
;; ---------------------------------------------------------------------------

(deftest un-nodo-usado-en-varios-sitios-suma-sus-adjuntos
  (testing "f(x) = x + x + ... (100 veces) → 100"
    (is (= [100.0] (rev/grad (fn [[x]] (apply ad/d+ (repeat 100 x))) [3.0]))))
  (testing "f(x) = Σ x² (50 veces) → 100x"
    (is (igual? 300.0 (first (rev/grad (fn [[x]] (apply ad/d+ (repeat 50 (ad/d* x x)))) [3.0])) 1e-12)))
  (testing "diamante: y = (x+1)(x+2) → 2x + 3"
    (is (igual? 13.0 (rev/diff (fn [x] (ad/d* (ad/d+ x 1) (ad/d+ x 2))) 5.0) 1e-12)))
  (testing "subexpresión compartida: t = sin x; f = t·t + t → 2t·cos x + cos x"
    (let [x 0.7 t (Math/sin x)]
      (is (igual? (+ (* 2.0 t (Math/cos x)) (Math/cos x))
                  (rev/diff (fn [v] (let [s (ad/dsin v)] (ad/d+ (ad/d* s s) s))) x) 1e-12)))))

;; ---------------------------------------------------------------------------
;; Mismo resultado que forward
;; ---------------------------------------------------------------------------

(def ^:private funciones
  [(fn [x] (ad/d* x x (ad/dsin x)))
   (fn [x] (ad/d-div (ad/dexp x) (ad/d+ 1 (ad/d* x x))))
   (fn [x] (ad/dlog (ad/d+ 1 (ad/d* x x))))
   (fn [x] (ad/dsin (ad/dcos x)))
   (fn [x] (ad/dpow x 1.5))
   (fn [x] (ad/dpow x x))])

(defspec reverse-coincide-con-forward-en-una-variable 300
  (prop/for-all [i (gen/choose 0 (dec (count funciones)))
                 x gen-x]
    (let [f (nth funciones i)]
      (igual? (fwd/diff f x) (rev/diff f x) 1e-12))))

;; Árboles de expresión aleatorios. Los denominadores y los argumentos del
;; logaritmo se protegen con 1 + e² (≥ 1) para evitar dominios inválidos.

(def ^:private gen-arbol
  (gen/recursive-gen
   (fn [interno]
     (gen/one-of
      [(gen/fmap (fn [[op a b]] [op a b]) (gen/tuple (gen/elements [:+ :- :*]) interno interno))
       (gen/fmap (fn [[a b]] [:div a b]) (gen/tuple interno interno))
       (gen/fmap (fn [[op a]] [op a]) (gen/tuple (gen/elements [:sin :cos :exp :log]) interno))
       (gen/fmap (fn [[n a]] [:pow n a]) (gen/tuple (gen/elements [0 1 2 3]) interno))]))
   (gen/one-of [(gen/elements [[:var 0] [:var 1] [:var 2]])
                (gen/fmap (fn [c] [:const c]) (gen/double* {:min -2.0 :max 2.0 :NaN? false :infinite? false}))])))

(defn- evaluar-arbol [arbol xs]
  (let [ev #(evaluar-arbol % xs)]
    (case (first arbol)
      :var (xs (second arbol))
      :const (second arbol)
      :+ (ad/d+ (ev (nth arbol 1)) (ev (nth arbol 2)))
      :- (ad/d- (ev (nth arbol 1)) (ev (nth arbol 2)))
      :* (ad/d* (ev (nth arbol 1)) (ev (nth arbol 2)))
      :div (let [b (ev (nth arbol 2))] (ad/d-div (ev (nth arbol 1)) (ad/d+ 1.0 (ad/d* b b))))
      :sin (ad/dsin (ev (nth arbol 1)))
      :cos (ad/dcos (ev (nth arbol 1)))
      :exp (ad/dexp (ad/dsin (ev (nth arbol 1))))
      :log (let [a (ev (nth arbol 1))] (ad/dlog (ad/d+ 1.0 (ad/d* a a))))
      :pow (ad/dpow (ev (nth arbol 2)) (second arbol)))))

(defspec reverse-coincide-con-forward-en-expresiones-aleatorias 500
  (prop/for-all [arbol gen-arbol
                 p (gen/vector gen-x 3)]
    (let [f (fn [xs] (evaluar-arbol arbol xs))]
      (vec-igual? (fwd/grad f p) (rev/grad f p) 1e-9))))

(defspec el-valor-de-reverse-coincide-con-el-de-forward 300
  (prop/for-all [arbol gen-arbol
                 p (gen/vector gen-x 3)]
    (let [f (fn [xs] (evaluar-arbol arbol xs))
          en-forward (let [s (f (mapv #(fwd/->Dual % 0.0) p))] (if (number? s) s (:a s)))]
      (igual? en-forward (:value (rev/value-and-grad f p)) 1e-9))))

;; ---------------------------------------------------------------------------
;; Frente a los gradientes analíticos de Rosenbrock y de mínimos cuadrados
;; ---------------------------------------------------------------------------

(defspec rosenbrock-reverse-coincide-con-el-analitico 200
  (prop/for-all [x (gen/double* {:min -2.0 :max 2.0 :NaN? false :infinite? false})
                 y (gen/double* {:min -2.0 :max 2.0 :NaN? false :infinite? false})]
    (vec-igual? (opt/rosenbrock-gradient [x y]) (rev/grad opt/rosenbrock [x y]) 1e-9)))

(def ^:private gen-problema
  (let [num (gen/double* {:min -5.0 :max 5.0 :NaN? false :infinite? false})]
    (gen/bind (gen/tuple (gen/choose 3 8) (gen/choose 1 4))
              (fn [[m n]] (gen/tuple (gen/vector (gen/vector num n) m)
                                     (gen/vector num m)
                                     (gen/vector num n))))))

(defspec minimos-cuadrados-reverse-coincide-con-el-analitico-y-con-forward 200
  (prop/for-all [[x y w] gen-problema]
    (let [f (opt/least-squares-objective x y)
          analitico (opt/least-squares-gradient x y w)]
      (and (vec-igual? analitico (rev/grad f w) 1e-9)
           (vec-igual? (fwd/grad f w) (rev/grad f w) 1e-9)))))

;; ---------------------------------------------------------------------------
;; Inmutabilidad, repetibilidad y presentación
;; ---------------------------------------------------------------------------

(deftest no-hay-estado-compartido
  (let [f (fn [[x y]] (ad/d+ (ad/d* x x y) (ad/dsin y)))
        p [1.0 2.0]
        copia (vec p)]
    (is (= (rev/grad f p) (rev/grad f p)) "dos llamadas dan exactamente lo mismo")
    (is (= copia p) "la entrada no se modifica")
    (is (= (rev/grad f [1.0 2.0]) (rev/grad f [1.0 2.0])))
    (testing "llamadas anidadas en un mismo hilo con funciones distintas no se interfieren"
      (let [a (rev/grad (fn [[x]] (ad/d* x x)) [3.0])
            b (rev/grad (fn [[x]] (ad/d* x x x)) [3.0])
            a2 (rev/grad (fn [[x]] (ad/d* x x)) [3.0])]
        (is (= [6.0] a a2))
        (is (= [27.0] b))))))

(deftest los-nodos-se-imprimen-sin-recorrer-el-grafo
  ;; 1000 términos: con identidad estructural de ids (enfoque 1a) las cadenas de
  ;; ~5000 o más desbordan la pila en la igualdad de ids anidados (ver
  ;; fpnum.ad.reverse-pending y scripts/measure_reverse_sharing.clj).
  (let [capturado (volatile! nil)
        _ (rev/value-and-grad
           (fn [[x]]
             (let [cadena (reduce (fn [acc _] (ad/d+ (ad/d* acc 0.999) x)) x (range 1000))]
               (vreset! capturado cadena)
               cadena))
           [1.0])]
    (is (rev/node? @capturado))
    (is (< (count (pr-str @capturado)) 60) "se muestra solo la operación y el valor")
    (is (re-find #"^#Node\[add " (pr-str @capturado)))))
