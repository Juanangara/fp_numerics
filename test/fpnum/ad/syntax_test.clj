(ns fpnum.ad.syntax-test
  "Pruebas de la macro ad-fn (incremento 1: reescritura básica con macroexpand,
  sin manejo de sombreado). Incluye pruebas que FIJAN LOS LÍMITES documentados:
  describen lo que la macro NO hace hoy, no un comportamiento deseable."
  (:require [clojure.test :refer [deftest is testing]]
            [fpnum.ad.forward :as fwd]
            [fpnum.ad.ops :as ad]
            [fpnum.ad.reverse :as rev]
            [fpnum.ad.syntax :as syn :refer [ad-fn]]
            [fpnum.optim :as opt]))

(defn- igual? [a b tol]
  (<= (Math/abs (- (double a) (double b))) (* tol (+ 1.0 (Math/abs (double a)) (Math/abs (double b))))))

(defn- vec-igual? [x y tol]
  (and (= (count x) (count y)) (every? true? (map #(igual? %1 %2 tol) x y))))

;; ---------------------------------------------------------------------------
;; La reescritura (sin evaluar)
;; ---------------------------------------------------------------------------

(deftest reescritura-de-operadores
  (is (= `(ad/d+ (ad/d* ~'x ~'x) 1) (syn/ad-rewrite '(+ (* x x) 1))))
  (is (= `(ad/d- ~'a ~'b) (syn/ad-rewrite '(- a b))))
  (is (= `(ad/d- ~'a) (syn/ad-rewrite '(- a))) "resta unaria")
  (is (= `(ad/d-div ~'a ~'b) (syn/ad-rewrite '(/ a b))))
  (is (= `(ad/d-div ~'a) (syn/ad-rewrite '(/ a))) "inverso")
  (testing "las formas calificadas de clojure.core"
    (is (= `(ad/d+ ~'x 1) (syn/ad-rewrite '(clojure.core/+ x 1))))
    (is (= `(ad/d* ~'x ~'x) (syn/ad-rewrite '(clojure.core/* x x))))
    (is (= `(ad/d- ~'x 1) (syn/ad-rewrite '(clojure.core/- x 1))))
    (is (= `(ad/d-div ~'x 2) (syn/ad-rewrite '(clojure.core// x 2))))))

(deftest reescritura-de-funciones-de-math
  (is (= `(ad/dsin ~'x) (syn/ad-rewrite '(Math/sin x))))
  (is (= `(ad/dcos ~'x) (syn/ad-rewrite '(Math/cos x))))
  (is (= `(ad/dexp ~'x) (syn/ad-rewrite '(Math/exp x))))
  (is (= `(ad/dlog ~'x) (syn/ad-rewrite '(Math/log x))))
  (is (= `(ad/dpow ~'x 3) (syn/ad-rewrite '(Math/pow x 3))))
  (is (= `(syn/ad-sqrt ~'x) (syn/ad-rewrite '(Math/sqrt x))))
  (testing "con el nombre completo de la clase"
    (is (= `(ad/dsin ~'x) (syn/ad-rewrite '(java.lang.Math/sin x))))
    (is (= `(ad/dpow ~'x 2) (syn/ad-rewrite '(java.lang.Math/pow x 2)))))
  (testing "la forma de interoperabilidad explícita"
    (is (= `(ad/dsin ~'x) (syn/ad-rewrite '(. Math sin x))))
    (is (= `(ad/dpow ~'x 2) (syn/ad-rewrite '(. java.lang.Math pow x 2)))))
  (testing "otras funciones de Math no se tocan"
    (is (= '(Math/abs x) (syn/ad-rewrite '(Math/abs x))))
    (is (= 'Math/PI (syn/ad-rewrite 'Math/PI)))))

(deftest los-operadores-como-valor-tambien-se-reescriben
  ;; (list 'reduce ...) y no la comilla sintáctica, que calificaría reduce como
  ;; clojure.core/reduce y la forma reescrita conserva el símbolo tal cual.
  (is (= (list 'reduce `ad/d+ 0 'xs) (syn/ad-rewrite '(reduce + 0 xs))))
  (is (= (list 'apply `ad/d* 'xs) (syn/ad-rewrite '(apply * xs))))
  (is (= (list 'map `ad/d* 'a 'b) (syn/ad-rewrite '(map * a b)))))

(deftest las-macros-se-expanden-antes-de-reescribir
  (testing "-> y ->>"
    (is (= `(ad/d+ (ad/d* ~'x 2) 1) (syn/ad-rewrite '(-> x (* 2) (+ 1)))))
    (is (= `(ad/d- 10 (ad/d* 2 ~'x)) (syn/ad-rewrite '(->> x (* 2) (- 10))))))
  (testing "let se convierte en let* y se reescribe por dentro"
    (is (= `(let* [~'t (ad/d* ~'x ~'x)] (ad/d+ ~'t 1))
           (syn/ad-rewrite '(let [t (* x x)] (+ t 1))))))
  (testing "las funciones anónimas de dentro del cuerpo también se reescriben"
    (let [r (syn/ad-rewrite '(apply + (map #(* % %) w)))]
      (is (= `ad/d+ (second r)))
      (is (= `ad/d* (first (nth (nth (nth r 2) 1) 2))) "el * de #(* % %)"))))

(deftest la-macro-expande-a-fn
  (is (= `(fn [~'x] (ad/d+ ~'x 1))
         (macroexpand-1 '(fpnum.ad.syntax/ad-fn [x] (+ x 1)))))
  (testing "admite varias formas en el cuerpo y destructuración"
    (is (= `(fn [[~'x ~'y]] (ad/d+ ~'x ~'y) (ad/d* ~'x ~'y))
           (macroexpand-1 '(fpnum.ad.syntax/ad-fn [[x y]] (+ x y) (* x y)))))))

;; ---------------------------------------------------------------------------
;; Evaluación: la misma función natural en los dos modos
;; ---------------------------------------------------------------------------

(deftest una-variable-en-los-dos-modos
  (let [f (ad-fn [x] (+ (* x x) (Math/sin x)))
        x 1.3
        exacta (+ (* 2.0 x) (Math/cos x))]
    (is (igual? exacta (fwd/diff f x) 1e-12))
    (is (igual? exacta (rev/diff f x) 1e-12))
    (testing "con números normales devuelve un Dual constante, no un número (límite 7)"
      (is (fwd/dual? (f x)))
      (is (igual? (+ (* x x) (Math/sin x)) (:a (f x)) 1e-12))
      (is (= 0.0 (:b (f x)))))))

(deftest varias-variables-con-destructuracion
  (let [g (ad-fn [[x y]] (+ (* x x y) (Math/sin y)))
        exacto [4.0 (+ 1.0 (Math/cos 2.0))]]
    (is (vec-igual? exacto (fwd/grad g [1.0 2.0]) 1e-12))
    (is (vec-igual? exacto (rev/grad g [1.0 2.0]) 1e-12))))

(deftest todas-las-funciones-de-la-tabla
  (doseq [[nombre f derivada x] [["exp" (ad-fn [t] (Math/exp t)) #(Math/exp %) 0.7]
                                 ["log" (ad-fn [t] (Math/log t)) #(/ 1.0 %) 2.0]
                                 ["cos" (ad-fn [t] (Math/cos t)) #(- (Math/sin %)) 0.9]
                                 ["pow" (ad-fn [t] (Math/pow t 3)) #(* 3.0 % %) 1.5]
                                 ["sqrt" (ad-fn [t] (Math/sqrt t)) #(/ 0.5 (Math/sqrt %)) 4.0]
                                 ["inverso" (ad-fn [t] (/ t)) #(- (/ 1.0 (* % %))) 2.0]
                                 ["opuesto" (ad-fn [t] (- t)) (constantly -1.0) 3.0]
                                 ["cociente de 3" (ad-fn [t] (/ t 2 3)) (constantly (/ 1.0 6.0)) 5.0]]]
    (testing nombre
      (is (igual? (derivada x) (fwd/diff f x) 1e-12))
      (is (igual? (derivada x) (rev/diff f x) 1e-12)))))

(deftest bucles-sobre-datos-dentro-del-cuerpo
  (let [w [1.0 -2.0 3.0]
        esperado (mapv #(* 2.0 %) w)]
    (testing "apply + con una función anónima"
      (let [f (ad-fn [v] (apply + (map #(* % %) v)))]
        (is (vec-igual? esperado (fwd/grad f w) 1e-12))
        (is (vec-igual? esperado (rev/grad f w) 1e-12))))
    (testing "reduce + con for"
      (let [f (ad-fn [v] (reduce + (for [a v] (* a a))))]
        (is (vec-igual? esperado (fwd/grad f w) 1e-12))
        (is (vec-igual? esperado (rev/grad f w) 1e-12))))
    (testing "producto escalar con map * y apply +"
      (let [c [2.0 3.0 4.0]
            f (ad-fn [v] (apply + (map * c v)))]
        (is (vec-igual? c (fwd/grad f w) 1e-12))
        (is (vec-igual? c (rev/grad f w) 1e-12))))))

(deftest macros-de-hilado-y-let-dentro-del-cuerpo
  (let [f (ad-fn [x] (let [t (* x x)] (-> t (+ 1) (* 2))))   ; 2(x² + 1)
        x 3.0]
    (is (igual? (* 4.0 x) (fwd/diff f x) 1e-12))
    (is (igual? (* 4.0 x) (rev/diff f x) 1e-12))
    (is (igual? 20.0 (:a (f x)) 1e-12))))

(deftest funciones-constantes
  (is (= [0.0] (rev/grad (ad-fn [x] 5.0) [1.0])))
  (is (= [0.0] (rev/grad (ad-fn [x] (* 2 3)) [1.0])))
  (is (= 0.0 (fwd/diff (ad-fn [x] (* 2 3)) 1.0))))

(deftest rosenbrock-con-sintaxis-natural
  (let [f (ad-fn [[x y]] (+ (* (- 1 x) (- 1 x)) (* 100 (- y (* x x)) (- y (* x x)))))]
    (doseq [p [[-1.2 1.0] [0.5 0.5] [1.0 1.0] [2.0 -1.0]]]
      (is (vec-igual? (opt/rosenbrock-gradient p) (fwd/grad f p) 1e-9) (str "forward " p))
      (is (vec-igual? (opt/rosenbrock-gradient p) (rev/grad f p) 1e-9) (str "reverse " p)))))

(deftest sirve-con-el-descenso-de-gradiente-en-los-dos-modos
  (let [f (ad-fn [[x y]] (+ (* (- x 3) (- x 3)) (* (+ y 1) (+ y 1))))   ; mínimo en (3, -1)
        opts {:alpha 0.1 :max-iter 500 :tol 1e-9}
        con-forward (opt/gradient-descent f [0.0 0.0] opts)
        con-reverse (opt/gradient-descent f [0.0 0.0] (assoc opts :grad-fn rev/grad))]
    (is (:converged? con-forward))
    (is (vec-igual? [3.0 -1.0] (:x con-forward) 1e-6))
    (is (vec-igual? [3.0 -1.0] (:x con-reverse) 1e-6))
    (is (= (:iterations con-forward) (:iterations con-reverse)))))

;; ---------------------------------------------------------------------------
;; LÍMITES documentados (describen lo que la macro NO hace en este incremento)
;; ---------------------------------------------------------------------------

(defn- cuadrado-de-clojure-core
  "Función definida FUERA de ad-fn: usa el * de clojure.core."
  [v]
  (* v v))

(deftest limite-1-no-reescribe-las-funciones-definidas-fuera
  (let [f (ad-fn [x] (+ 1 (cuadrado-de-clojure-core x)))]
    (is (igual? 5.0 (:a (f 2.0)) 1e-12) "con números normales funciona (el resultado es un Dual constante)")
    (is (thrown? ClassCastException (fwd/diff f 2.0)) "con un Dual falla")
    (is (thrown? ClassCastException (rev/diff f 2.0)) "con un nodo falla")))

(deftest limite-2-reescribe-tambien-la-aritmetica-de-indices
  (testing "(+ 0 1) pasa a (d+ 0 1), que devuelve un Dual y no sirve de índice, ni siquiera con números normales"
    (is (thrown? ClassCastException ((ad-fn [x] (nth [1.0 2.0 3.0] (+ 0 1))) 5.0))))
  (testing "con - * / pasa lo mismo"
    (is (thrown? ClassCastException ((ad-fn [x] (nth [1.0 2.0 3.0] (* 1 1))) 5.0))))
  (testing "inc y dec NO se reescriben: sirven de contador"
    (is (= 3.0 ((ad-fn [x] (nth [1.0 2.0 3.0] (inc (inc 0)))) 5.0)))
    (is (= '(inc n) (syn/ad-rewrite '(inc n))))
    (is (= '(dec n) (syn/ad-rewrite '(dec n))))))

(deftest limite-5-las-comparaciones-no-se-reescriben
  (let [f (ad-fn [x] (if (< x 3) (* x x) 0))]
    (is (= 4.0 (:a (f 2.0))) "con números normales funciona (el resultado es un Dual constante)")
    (is (thrown? Exception (fwd/diff f 2.0)) "con un Dual, `<` falla")
    (is (thrown? Exception (rev/diff f 2.0)) "con un nodo, `<` falla")))

(deftest limite-3-sin-manejo-de-sombreado
  (testing "un local llamado + se sustituye igualmente y el código no compila"
    (is (thrown? Exception
                 (eval '(fpnum.ad.syntax/ad-fn [x] (let [+ (fn [a b] a)] (+ x 1)))))))
  (testing "la reescritura no mira si el + es el de clojure.core"
    (is (= `(ad/d+ ~'a ~'b) (syn/ad-rewrite '(+ a b))))))

(deftest limite-4-las-formas-citadas-tambien-se-reescriben
  (is (= `(quote (ad/d+ 1 2)) (syn/ad-rewrite '(quote (+ 1 2))))
      "dentro de un quote la reescritura no debería entrar, pero lo hace en este incremento"))
