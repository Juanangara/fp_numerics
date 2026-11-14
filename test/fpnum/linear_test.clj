(ns fpnum.linear-test
  "Pruebas de fpnum.linear: casos concretos hechos a mano y propiedades
  algebraicas con test.check.

  Con doubles las igualdades algebraicas (como la asociatividad) no son
  exactas por el redondeo, así que las propiedades se comparan con una
  tolerancia relativa (`approx=`) y no con `=`."
  (:require [clojure.test :refer [deftest is testing]]
            [clojure.test.check.clojure-test :refer [defspec]]
            [clojure.test.check.generators :as gen]
            [clojure.test.check.properties :as prop]
            [fpnum.linear :as la]))

;; ---------------------------------------------------------------------------
;; Comparación con tolerancia
;; ---------------------------------------------------------------------------

(def ^:private eps 1e-9)

(defn- approx-num= [a b]
  (<= (Math/abs (- (double a) (double b)))
      (* eps (+ 1.0 (Math/abs (double a)) (Math/abs (double b))))))

(defn- approx-vec= [x y]
  (and (= (count x) (count y))
       (every? true? (map approx-num= x y))))

(defn- approx-mat= [a b]
  (and (= (count a) (count b))
       (every? true? (map approx-vec= a b))))

;; ---------------------------------------------------------------------------
;; Generadores (valores acotados para que el redondeo no domine)
;; ---------------------------------------------------------------------------

(def ^:private gen-num
  (gen/double* {:min -100.0 :max 100.0 :NaN? false :infinite? false}))

(def ^:private gen-dim (gen/choose 1 5))

(defn- gen-vec [n] (gen/vector gen-num n))

(defn- gen-mat [m n] (gen/vector (gen-vec n) m))

;; ---------------------------------------------------------------------------
;; Casos concretos
;; ---------------------------------------------------------------------------

(deftest constructores
  (is (= [0.0 0.0 0.0] (la/zeros 3)))
  (is (= [[0.0 0.0 0.0] [0.0 0.0 0.0]] (la/zero-matrix 2 3)))
  (is (= [[1.0 0.0] [0.0 1.0]] (la/identity-matrix 2)))
  (is (= [2 3] (la/dims [[1.0 2.0 3.0] [4.0 5.0 6.0]]))))

(deftest operaciones-de-vectores
  (is (= [5.0 7.0 9.0] (la/v+ [1.0 2.0 3.0] [4.0 5.0 6.0])))
  (is (= [-3.0 -3.0 -3.0] (la/v- [1.0 2.0 3.0] [4.0 5.0 6.0])))
  (is (= [2.0 4.0 6.0] (la/scale 2.0 [1.0 2.0 3.0])))
  (is (= 32.0 (la/dot [1.0 2.0 3.0] [4.0 5.0 6.0])))
  (is (approx-num= 5.0 (la/norm2 [3.0 4.0]))))

(deftest operaciones-de-matrices
  (let [a [[1.0 2.0] [3.0 4.0]]
        b [[5.0 6.0] [7.0 8.0]]]
    (is (= [[1.0 3.0] [2.0 4.0]] (la/transpose a)))
    (is (= [[1.0 4.0] [2.0 5.0] [3.0 6.0]]
           (la/transpose [[1.0 2.0 3.0] [4.0 5.0 6.0]])))
    (is (= [17.0 39.0] (la/mat-vec a [5.0 6.0])))
    (is (= [[19.0 22.0] [43.0 50.0]] (la/mat-mat a b)))
    (is (approx-num= (Math/sqrt 30.0) (la/norm-fro a)))))

(deftest matriz-rectangular
  (testing "(2x3)·(3x2) da una matriz de 2x2"
    (is (= [[58.0 64.0] [139.0 154.0]]
           (la/mat-mat [[1.0 2.0 3.0] [4.0 5.0 6.0]]
                       [[7.0 8.0] [9.0 10.0] [11.0 12.0]])))))

(deftest dimensiones-incompatibles
  (is (thrown? clojure.lang.ExceptionInfo (la/v+ [1.0 2.0] [1.0 2.0 3.0])))
  (is (thrown? clojure.lang.ExceptionInfo (la/v- [1.0 2.0] [1.0])))
  (is (thrown? clojure.lang.ExceptionInfo (la/dot [1.0] [1.0 2.0])))
  (is (thrown? clojure.lang.ExceptionInfo (la/mat-vec [[1.0 2.0]] [1.0 2.0 3.0])))
  (is (thrown? clojure.lang.ExceptionInfo (la/mat-mat [[1.0 2.0]] [[1.0 2.0]])))
  (testing "el mapa de datos informa las dimensiones"
    (is (= {:op 'mat-mat :a [1 2] :b [1 2]}
           (try (la/mat-mat [[1.0 2.0]] [[1.0 2.0]])
                (catch clojure.lang.ExceptionInfo e (ex-data e)))))))

(deftest inmutabilidad
  (let [x [1.0 2.0 3.0]
        a [[1.0 2.0] [3.0 4.0]]]
    (la/v+ x x)
    (la/scale 2.0 x)
    (la/transpose a)
    (la/mat-mat a a)
    (is (= [1.0 2.0 3.0] x))
    (is (= [[1.0 2.0] [3.0 4.0]] a))))

;; ---------------------------------------------------------------------------
;; Propiedades: vectores
;; ---------------------------------------------------------------------------

(defspec suma-conmutativa 200
  (prop/for-all [[x y] (gen/bind gen-dim (fn [n] (gen/tuple (gen-vec n) (gen-vec n))))]
    (approx-vec= (la/v+ x y) (la/v+ y x))))

(defspec suma-asociativa 200
  (prop/for-all [[x y z] (gen/bind gen-dim
                                   (fn [n] (gen/tuple (gen-vec n) (gen-vec n) (gen-vec n))))]
    (approx-vec= (la/v+ (la/v+ x y) z) (la/v+ x (la/v+ y z)))))

(defspec suma-neutro-cero 200
  (prop/for-all [x (gen/bind gen-dim gen-vec)]
    (approx-vec= x (la/v+ x (la/zeros (count x))))))

(defspec resta-es-suma-del-opuesto 200
  (prop/for-all [[x y] (gen/bind gen-dim (fn [n] (gen/tuple (gen-vec n) (gen-vec n))))]
    (approx-vec= (la/v- x y) (la/v+ x (la/scale -1.0 y)))))

(defspec dot-simetrico 200
  (prop/for-all [[x y] (gen/bind gen-dim (fn [n] (gen/tuple (gen-vec n) (gen-vec n))))]
    (approx-num= (la/dot x y) (la/dot y x))))

(defspec dot-lineal 200
  (prop/for-all [[x y z] (gen/bind gen-dim
                                   (fn [n] (gen/tuple (gen-vec n) (gen-vec n) (gen-vec n))))
                 k gen-num]
    (approx-num= (la/dot (la/v+ (la/scale k x) y) z)
                 (+ (* k (la/dot x z)) (la/dot y z)))))

(defspec norma2-no-negativa-y-homogenea 200
  (prop/for-all [x (gen/bind gen-dim gen-vec)
                 k gen-num]
    (and (>= (la/norm2 x) 0.0)
         (approx-num= (la/norm2 (la/scale k x))
                      (* (Math/abs k) (la/norm2 x))))))

(defspec cauchy-schwarz 200
  (prop/for-all [[x y] (gen/bind gen-dim (fn [n] (gen/tuple (gen-vec n) (gen-vec n))))]
    (<= (Math/abs (la/dot x y))
        (+ (* (la/norm2 x) (la/norm2 y)) (* eps (+ 1.0 (* (la/norm2 x) (la/norm2 y))))))))

(defspec desigualdad-triangular 200
  (prop/for-all [[x y] (gen/bind gen-dim (fn [n] (gen/tuple (gen-vec n) (gen-vec n))))]
    (let [lado-izq (la/norm2 (la/v+ x y))
          lado-der (+ (la/norm2 x) (la/norm2 y))]
      (<= lado-izq (+ lado-der (* eps (+ 1.0 lado-der)))))))

;; ---------------------------------------------------------------------------
;; Propiedades: matrices
;; ---------------------------------------------------------------------------

(defspec transpuesta-involutiva 200
  (prop/for-all [a (gen/bind (gen/tuple gen-dim gen-dim) (fn [[m n]] (gen-mat m n)))]
    (= a (la/transpose (la/transpose a)))))

(defspec transpuesta-cambia-dimensiones 200
  (prop/for-all [[m n] (gen/tuple gen-dim gen-dim)]
    (let [mat (la/zero-matrix m n)]
      (= [n m] (la/dims (la/transpose mat))))))

(defspec transpuesta-del-producto 200
  (prop/for-all [[a b] (gen/bind (gen/tuple gen-dim gen-dim gen-dim)
                                 (fn [[m n p]] (gen/tuple (gen-mat m n) (gen-mat n p))))]
    (approx-mat= (la/transpose (la/mat-mat a b))
                 (la/mat-mat (la/transpose b) (la/transpose a)))))

(defspec producto-asociativo 200
  (prop/for-all [[a b c] (gen/bind (gen/tuple gen-dim gen-dim gen-dim gen-dim)
                                   (fn [[m n p q]]
                                     (gen/tuple (gen-mat m n) (gen-mat n p) (gen-mat p q))))]
    (approx-mat= (la/mat-mat (la/mat-mat a b) c)
                 (la/mat-mat a (la/mat-mat b c)))))

(defspec identidad-es-neutro 200
  (prop/for-all [a (gen/bind (gen/tuple gen-dim gen-dim) (fn [[m n]] (gen-mat m n)))]
    (let [[m n] (la/dims a)]
      (and (approx-mat= a (la/mat-mat (la/identity-matrix m) a))
           (approx-mat= a (la/mat-mat a (la/identity-matrix n)))))))

(defspec mat-vec-distribuye-sobre-suma 200
  (prop/for-all [[a x y] (gen/bind (gen/tuple gen-dim gen-dim)
                                   (fn [[m n]] (gen/tuple (gen-mat m n) (gen-vec n) (gen-vec n))))]
    (approx-vec= (la/mat-vec a (la/v+ x y))
                 (la/v+ (la/mat-vec a x) (la/mat-vec a y)))))

(defspec mat-vec-es-mat-mat-con-columna 200
  (prop/for-all [[a x] (gen/bind (gen/tuple gen-dim gen-dim)
                                 (fn [[m n]] (gen/tuple (gen-mat m n) (gen-vec n))))]
    (let [columna (mapv vector x)]
      (approx-vec= (la/mat-vec a x)
                   (mapv first (la/mat-mat a columna))))))

(defspec frobenius-es-norma2-de-las-filas 200
  (prop/for-all [a (gen/bind (gen/tuple gen-dim gen-dim) (fn [[m n]] (gen-mat m n)))]
    (let [por-filas (Math/sqrt (reduce + 0.0 (map #(let [r (la/norm2 %)] (* r r)) a)))]
      (approx-num= (la/norm-fro a) por-filas))))

(defspec frobenius-invariante-por-transpuesta 200
  (prop/for-all [a (gen/bind (gen/tuple gen-dim gen-dim) (fn [[m n]] (gen-mat m n)))]
    (approx-num= (la/norm-fro a) (la/norm-fro (la/transpose a)))))

;; ---------------------------------------------------------------------------
;; solve (eliminación gaussiana con pivoteo parcial)
;; ---------------------------------------------------------------------------

(deftest solve-casos-concretos
  (testing "sistema 2x2: 2x + y = 5, x + 3y = 10 → x = 1, y = 3"
    (is (approx-vec= [1.0 3.0] (la/solve [[2.0 1.0] [1.0 3.0]] [5.0 10.0]))))
  (testing "sistema 3x3"
    (is (approx-vec= [2.0 3.0 -1.0]
                     (la/solve [[2.0 1.0 -1.0] [-3.0 -1.0 2.0] [-2.0 1.0 2.0]]
                               [8.0 -11.0 -3.0]))))
  (testing "necesita intercambio de filas (pivote cero en la diagonal)"
    (is (approx-vec= [3.0 2.0] (la/solve [[0.0 1.0] [1.0 0.0]] [2.0 3.0]))))
  (testing "identidad"
    (is (approx-vec= [4.0 5.0 6.0]
                     (la/solve (la/identity-matrix 3) [4.0 5.0 6.0])))))

(deftest solve-matriz-singular
  (testing "filas proporcionales"
    (is (thrown-with-msg? clojure.lang.ExceptionInfo #"singular"
                          (la/solve [[1.0 2.0] [2.0 4.0]] [1.0 2.0]))))
  (testing "matriz de ceros"
    (is (thrown-with-msg? clojure.lang.ExceptionInfo #"singular"
                          (la/solve (la/zero-matrix 3 3) [1.0 2.0 3.0]))))
  (testing "3x3 de rango 2 (la tercera fila es combinación de las otras)"
    (is (thrown-with-msg? clojure.lang.ExceptionInfo #"singular"
                          (la/solve [[1.0 2.0 3.0] [4.0 5.0 6.0] [7.0 8.0 9.0]]
                                    [1.0 2.0 3.0]))))
  (testing "el mapa de datos indica la operación"
    (is (= 'solve (:op (try (la/solve [[1.0 2.0] [2.0 4.0]] [1.0 2.0])
                            (catch clojure.lang.ExceptionInfo e (ex-data e))))))))

(deftest solve-dimensiones-invalidas
  (is (thrown? clojure.lang.ExceptionInfo (la/solve [[1.0 2.0 3.0] [4.0 5.0 6.0]] [1.0 2.0])))
  (is (thrown? clojure.lang.ExceptionInfo (la/solve [[1.0 0.0] [0.0 1.0]] [1.0 2.0 3.0])))
  (is (thrown? clojure.lang.ExceptionInfo (la/solve [] []))))

(defn- gen-sistema-bien-condicionado
  "Matriz n×n con diagonal dominante (siempre invertible) y vector b."
  []
  (gen/bind gen-dim
            (fn [n]
              (gen/fmap (fn [[a b]]
                          [(vec (map-indexed
                                 (fn [i fila] (update fila i + (* 1000.0 (inc n))))
                                 a))
                           b])
                        (gen/tuple (gen-mat n n) (gen-vec n))))))

(defspec solve-residuo-pequeno 200
  (prop/for-all [[a b] (gen-sistema-bien-condicionado)]
    (approx-vec= b (la/mat-vec a (la/solve a b)))))

(defspec solve-recupera-x-conocido 200
  (prop/for-all [[a x] (gen/bind gen-dim
                                 (fn [n]
                                   (gen/tuple (gen/fmap (fn [m] (vec (map-indexed
                                                                      (fn [i f] (update f i + (* 1000.0 (inc n))))
                                                                      m)))
                                                        (gen-mat n n))
                                              (gen-vec n))))]
    (approx-vec= x (la/solve a (la/mat-vec a x)))))
