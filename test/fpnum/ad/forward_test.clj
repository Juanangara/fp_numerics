(ns fpnum.ad.forward-test
  "Pruebas de fpnum.ad.forward: derivadas conocidas, regla de la cadena y
  propiedades contra diferencias finitas centrales."
  (:require [clojure.test :refer [deftest is testing]]
            [clojure.test.check.clojure-test :refer [defspec]]
            [clojure.test.check.generators :as gen]
            [clojure.test.check.properties :as prop]
            [fpnum.ad.forward :as ad]))

(defn- approx= [a b tol]
  (or (= a b)
      (<= (Math/abs (- (double a) (double b)))
          (* tol (+ 1.0 (Math/abs (double a)) (Math/abs (double b)))))))

;; ---------------------------------------------------------------------------
;; Tipo Dual
;; ---------------------------------------------------------------------------

(deftest tipo-dual
  (is (= (ad/->Dual 1.0 2.0) (ad/dual 1 2)))
  (is (ad/dual? (ad/dual 1 2)))
  (is (not (ad/dual? 3.0)))
  (is (= (ad/->Dual 3.0 0.0) (ad/lift 3)))
  (let [d (ad/dual 1 2)] (is (identical? d (ad/lift d)))))

;; ---------------------------------------------------------------------------
;; Aritmética básica
;; ---------------------------------------------------------------------------

(deftest aritmetica
  (let [x (ad/dual 3 1)
        y (ad/dual 2 1)]
    (is (= (ad/dual 5 2) (ad/d+ x y)))
    (is (= (ad/dual 1 0) (ad/d- x y)))
    (is (= (ad/dual 6 5) (ad/d* x y)))            ; (3+ε)(2+ε) = 6 + 5ε
    (is (= (ad/dual 1.5 (/ (- 2.0 3.0) 4.0)) (ad/d-div x y)))
    (is (= (ad/dual -3 -1) (ad/d- x)))
    (is (= (ad/dual 0.5 -0.25) (ad/d-div (ad/dual 2 1))))))

(deftest aridades-variadicas
  (is (= (ad/dual 0 0) (ad/d+)))
  (is (= (ad/dual 1 0) (ad/d*)))
  (is (= (ad/dual 6 0) (ad/d+ 1 2 3)))
  (is (= (ad/dual 24 0) (ad/d* 2 3 4)))
  (is (= (ad/dual 4 0) (ad/d- 10 3 2 1))))

(deftest numeros-se-tratan-como-constantes
  (is (= (ad/dual 6 2) (ad/d* 2.0 (ad/dual 3 1))))
  (is (= (ad/dual 4 1) (ad/d+ (ad/dual 3 1) 1))))

;; ---------------------------------------------------------------------------
;; Derivadas conocidas
;; ---------------------------------------------------------------------------

(deftest derivadas-conocidas
  (testing "polinomios"
    (is (approx= 12.0 (ad/diff (fn [x] (ad/dpow x 3)) 2.0) 1e-12))
    (is (approx= 6.0 (ad/diff (fn [x] (ad/d+ (ad/d* x x x) (ad/d* 3 x) 2)) 1.0) 1e-12)))
  (testing "funciones elementales"
    (is (approx= 1.0 (ad/diff ad/dsin 0.0) 1e-12))
    (is (approx= (- (Math/sin 1.0)) (ad/diff ad/dcos 1.0) 1e-12))
    (is (approx= Math/E (ad/diff ad/dexp 1.0) 1e-12))
    (is (approx= 0.5 (ad/diff ad/dlog 2.0) 1e-12)))
  (testing "función constante: derivada 0"
    (is (= 0.0 (ad/diff (fn [_] 5.0) 3.0)))))

(deftest regla-del-producto-cociente-y-cadena
  (testing "producto: (x·sin x)' = sin x + x cos x"
    (let [x 1.3]
      (is (approx= (+ (Math/sin x) (* x (Math/cos x)))
                   (ad/diff (fn [t] (ad/d* t (ad/dsin t))) x) 1e-12))))
  (testing "cociente: (x/(1+x²))' = (1 − x²)/(1+x²)²"
    (let [x 2.0]
      (is (approx= (/ (- 1.0 (* x x)) (Math/pow (+ 1.0 (* x x)) 2))
                   (ad/diff (fn [t] (ad/d-div t (ad/d+ 1 (ad/d* t t)))) x) 1e-12))))
  (testing "cadena: (sin x²)' = 2x cos x²"
    (let [x 0.7]
      (is (approx= (* 2.0 x (Math/cos (* x x)))
                   (ad/diff (fn [t] (ad/dsin (ad/d* t t))) x) 1e-12))))
  (testing "cadena triple: (exp(sin(x²)))' = exp(sin x²)·cos x²·2x"
    (let [x 0.9
          s (Math/sin (* x x))]
      (is (approx= (* (Math/exp s) (Math/cos (* x x)) 2.0 x)
                   (ad/diff (fn [t] (ad/dexp (ad/dsin (ad/d* t t)))) x) 1e-12)))))

;; ---------------------------------------------------------------------------
;; Potencias
;; ---------------------------------------------------------------------------

(deftest potencias
  (is (= 0.0 (ad/diff (fn [x] (ad/dpow x 2)) 0.0)) "x² en 0 → 0")
  (is (approx= 3.0 (ad/diff (fn [x] (ad/dpow x 3)) -1.0) 1e-12))
  (is (approx= -4.0 (ad/diff (fn [x] (ad/dpow x 2)) -2.0) 1e-12))
  (is (= 0.0 (ad/diff (fn [x] (ad/dpow x 0)) 5.0)) "x⁰ = 1 es constante")
  (is (= 0.0 (ad/diff (fn [x] (ad/dpow x 0)) 0.0)) "sin NaN en x = 0")
  (testing "base constante, exponente variable: (2ˣ)' = ln 2 · 2ˣ"
    (is (approx= (* (Math/log 2.0) 8.0)
                 (ad/diff (fn [x] (ad/dpow 2 x)) 3.0) 1e-12)))
  (testing "(xˣ)' = xˣ (ln x + 1); en x = 1 vale 1"
    (is (approx= 1.0 (ad/diff (fn [x] (ad/dpow x x)) 1.0) 1e-12))
    (is (approx= (* (Math/pow 2.0 2.0) (+ (Math/log 2.0) 1.0))
                 (ad/diff (fn [x] (ad/dpow x x)) 2.0) 1e-12))))

;; ---------------------------------------------------------------------------
;; Dominios
;; ---------------------------------------------------------------------------

(deftest dominios
  (is (Double/isNaN (:a (ad/dlog (ad/dual -1 1)))))
  (is (Double/isInfinite (:a (ad/dlog (ad/dual 0 1)))))
  (is (Double/isInfinite (:a (ad/d-div 1 (ad/dual 0 1))))))

;; ---------------------------------------------------------------------------
;; Gradiente
;; ---------------------------------------------------------------------------

(deftest gradiente
  (testing "g(x,y) = x²y + sin y → ∇g = [2xy, x² + cos y]"
    (let [g (fn [[x y]] (ad/d+ (ad/d* x x y) (ad/dsin y)))
          [gx gy] (ad/grad g [1.0 2.0])]
      (is (approx= 4.0 gx 1e-12))
      (is (approx= (+ 1.0 (Math/cos 2.0)) gy 1e-12))))
  (testing "norma al cuadrado: ∇(Σ xᵢ²) = 2x"
    (let [f (fn [xs] (apply ad/d+ (map #(ad/d* % %) xs)))]
      (is (= [2.0 -4.0 6.0] (ad/grad f [1.0 -2.0 3.0])))))
  (testing "función constante: gradiente cero"
    (is (= [0.0 0.0] (ad/grad (fn [_] 7.0) [1.0 2.0])))))

;; ---------------------------------------------------------------------------
;; Propiedades contra diferencias finitas centrales
;; ---------------------------------------------------------------------------

(defn- fd-central [f x]
  (let [h 1e-5]
    (/ (- (f (+ x h)) (f (- x h))) (* 2.0 h))))

;; Pares [función con duales, misma función con doubles]
(def ^:private funciones
  [[(fn [x] (ad/d* x x (ad/dsin x)))
    (fn [x] (* x x (Math/sin x)))]
   [(fn [x] (ad/d-div (ad/dexp x) (ad/d+ 1 (ad/d* x x))))
    (fn [x] (/ (Math/exp x) (+ 1.0 (* x x))))]
   [(fn [x] (ad/dlog (ad/d+ 1 (ad/d* x x))))
    (fn [x] (Math/log (+ 1.0 (* x x))))]
   [(fn [x] (ad/dsin (ad/dcos x)))
    (fn [x] (Math/sin (Math/cos x)))]
   [(fn [x] (ad/dpow x 1.5))
    (fn [x] (Math/pow x 1.5))]
   [(fn [x] (ad/dpow x x))
    (fn [x] (Math/pow x x))]])

(defspec ad-coincide-con-diferencias-finitas 300
  (prop/for-all [i (gen/choose 0 (dec (count funciones)))
                 x (gen/double* {:min 0.5 :max 3.0 :NaN? false :infinite? false})]
    (let [[f-dual f-num] (nth funciones i)]
      (approx= (ad/diff f-dual x) (fd-central f-num x) 1e-6))))

(defspec polinomio-cuadratico-es-exacto 200
  (prop/for-all [a (gen/double* {:min -10.0 :max 10.0 :NaN? false :infinite? false})
                 b (gen/double* {:min -10.0 :max 10.0 :NaN? false :infinite? false})
                 c (gen/double* {:min -10.0 :max 10.0 :NaN? false :infinite? false})
                 x (gen/double* {:min -10.0 :max 10.0 :NaN? false :infinite? false})]
    (approx= (+ (* 2.0 a x) b)
             (ad/diff (fn [t] (ad/d+ (ad/d* a t t) (ad/d* b t) c)) x)
             1e-12)))

(defspec regla-de-la-cadena-sin-exp 200
  (prop/for-all [x (gen/double* {:min -2.0 :max 2.0 :NaN? false :infinite? false})]
    (approx= (* (Math/cos (Math/exp x)) (Math/exp x))
             (ad/diff (fn [t] (ad/dsin (ad/dexp t))) x)
             1e-12)))

(defspec linealidad-de-la-derivada 200
  (prop/for-all [k (gen/double* {:min -5.0 :max 5.0 :NaN? false :infinite? false})
                 x (gen/double* {:min -3.0 :max 3.0 :NaN? false :infinite? false})]
    (approx= (+ (* k (ad/diff ad/dsin x)) (ad/diff ad/dcos x))
             (ad/diff (fn [t] (ad/d+ (ad/d* k (ad/dsin t)) (ad/dcos t))) x)
             1e-12)))

(defspec gradiente-coincide-con-diferencias-finitas 200
  (prop/for-all [x (gen/double* {:min 0.5 :max 2.0 :NaN? false :infinite? false})
                 y (gen/double* {:min 0.5 :max 2.0 :NaN? false :infinite? false})]
    (let [g-dual (fn [[a b]] (ad/d+ (ad/d* a b) (ad/dsin b) (ad/dlog a)))
          g-num (fn [a b] (+ (* a b) (Math/sin b) (Math/log a)))
          [gx gy] (ad/grad g-dual [x y])]
      (and (approx= gx (fd-central #(g-num % y) x) 1e-6)
           (approx= gy (fd-central #(g-num x %) y) 1e-6)))))
