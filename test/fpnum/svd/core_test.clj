(ns fpnum.svd.core-test
  "Pruebas de fpnum.svd.core.

  * Valores singulares contra NumPy (fixture generado con
    scripts/generate_svd_fixture.py; src/ no depende de NumPy).
  * Reconstrucción U·Σ·Vᵀ = A y ortonormalidad de U y V.
  * Matriz de Hilbert, mal condicionada, con una referencia independiente:
    el determinante exacto (con racionales de Clojure).
  * Por qué se evita AᵀA: pasar por AᵀA destruye los valores singulares
    pequeños.
  * Propiedades con test.check."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :refer [deftest is testing]]
            [clojure.test.check.clojure-test :refer [defspec]]
            [clojure.test.check.generators :as gen]
            [clojure.test.check.properties :as prop]
            [fpnum.linear :as la]
            [fpnum.svd.core :as svd]))

;; ---------------------------------------------------------------------------
;; Utilidades
;; ---------------------------------------------------------------------------

(defn- error-reconstruccion
  "‖A − U Σ Vᵀ‖_F / ‖A‖_F (o el valor absoluto si A es cero)."
  [a r]
  (let [num (la/norm-fro (mapv la/v- a (svd/reconstruct r)))
        den (la/norm-fro a)]
    (if (zero? den) num (/ num den))))

(defn- error-ortonormal
  "Mayor |(QᵀQ − I)ᵢⱼ|: mide cuánto se aleja Q de tener columnas ortonormales."
  [q]
  (let [g (la/mat-mat (la/transpose q) q)
        n (count g)]
    (apply max (for [i (range n) j (range n)]
                 (Math/abs (- (get-in g [i j]) (if (= i j) 1.0 0.0)))))))

(defn- max-dif-abs [x y]
  (apply max 0.0 (map #(Math/abs (- (double %1) (double %2))) x y)))

(defn- error-relativo [x ref]
  (/ (Math/abs (- x ref)) (Math/abs ref)))

(defn- bien-formada?
  "Comprueba dimensiones, orden de s y ortonormalidad de U y V."
  [a r]
  (let [[m n] (la/dims a)
        k (min m n)
        s (:s r)]
    (and (= [m k] (la/dims (:u r)))
         (= [n k] (la/dims (:v r)))
         (= k (count s))
         (every? #(>= % 0.0) s)
         (apply >= s)
         (< (error-ortonormal (:u r)) 1e-12)
         (< (error-ortonormal (:v r)) 1e-12))))

;; ---------------------------------------------------------------------------
;; Fixture de NumPy
;; ---------------------------------------------------------------------------

(def ^:private fixture
  (delay (edn/read-string (slurp (io/resource "fpnum/svd/svd_reference.edn")))))

(defn- caso [nombre]
  (first (filter #(= nombre (:name %)) @fixture)))

(deftest fixture-cargado
  (is (>= (count @fixture) 12))
  (is (every? #(and (:name %) (:matrix %) (:s %)) @fixture)))

(deftest valores-singulares-coinciden-con-numpy
  (doseq [{:keys [name matrix s]} @fixture]
    (testing name
      (let [r (svd/svd matrix)
            s1 (first s)]
        (is (= (count s) (count (:s r))))
        ;; LAPACK solo garantiza error absoluto ~ eps·σ₁; se compara con esa
        ;; escala (1e-13·σ₁, unas 450 veces eps).
        (is (< (max-dif-abs (:s r) s) (* 1e-13 s1))
            (str "máxima diferencia con NumPy: " (max-dif-abs (:s r) s)))))))

(deftest reconstruccion-y-ortonormalidad-en-todo-el-fixture
  (doseq [{:keys [name matrix]} @fixture]
    (testing name
      (let [r (svd/svd matrix)]
        (is (:converged? r))
        (is (< (error-reconstruccion matrix r) 1e-13))
        (is (bien-formada? matrix r))))))

;; ---------------------------------------------------------------------------
;; Matriz de Hilbert (mal condicionada)
;; ---------------------------------------------------------------------------

(defn- det-hilbert
  "Determinante exacto de la matriz de Hilbert n×n:
  det Hₙ = cₙ⁴ / c₂ₙ con cₙ = Π_{i=1}^{n−1} i!. Se calcula con racionales."
  [n]
  (let [fact (fn [k] (reduce * 1N (range 1 (inc k))))
        c (fn [k] (reduce * 1N (map fact (range 1 k))))]
    (double (/ (* (c n) (c n) (c n) (c n)) (c (* 2 n))))))

(deftest hilbert-reconstruccion-y-ortonormalidad
  (doseq [n [5 8 10]]
    (let [{:keys [matrix]} (caso (str "hilbert-" n))
          r (svd/svd matrix)]
      (testing (str "Hilbert " n)
        (is (:converged? r))
        (is (< (error-reconstruccion matrix r) 1e-13))
        (is (< (error-ortonormal (:u r)) 1e-12))
        (is (< (error-ortonormal (:v r)) 1e-12))))))

(deftest hilbert-numero-de-condicion
  (testing "κ(H₈) ≈ 1.5e10 coincide con la referencia de NumPy"
    (let [{:keys [matrix cond]} (caso "hilbert-8")
          s (:s (svd/svd matrix))]
      (is (< (error-relativo (/ (first s) (peek s)) cond) 1e-3)))))

(deftest hilbert-producto-de-valores-singulares-es-el-determinante
  ;; |det A| = Π σᵢ. Referencia independiente de NumPy: el determinante exacto.
  ;; La matriz de Hilbert guardada en doubles difiere de la exacta en ~1e-16 por
  ;; entrada, lo que limita la precisión alcanzable en los σ pequeños: por eso
  ;; la tolerancia crece con κ.
  (doseq [[n tol] [[5 1e-9] [8 1e-5] [10 1e-3]]]
    (testing (str "Hilbert " n)
      (let [{:keys [matrix]} (caso (str "hilbert-" n))
            producto (reduce * (:s (svd/svd matrix)))]
        (is (< (error-relativo producto (det-hilbert n)) tol))))))

;; ---------------------------------------------------------------------------
;; Por qué se evita formar AᵀA
;; ---------------------------------------------------------------------------

(defn- sigmas-via-ata
  "Valores singulares calculados como √(autovalores de AᵀA), con la propia
  SVD de AᵀA. Es el camino que se evita: κ(AᵀA) = κ(A)²."
  [a]
  (mapv #(Math/sqrt (max 0.0 %))
        (:s (svd/svd (la/mat-mat (la/transpose a) a)))))

(deftest evitar-ata-conserva-los-valores-singulares-pequenos
  (doseq [nombre ["decaimiento-geometrico-12x12" "hilbert-8" "hilbert-10"]]
    (testing nombre
      (let [{:keys [matrix s]} (caso nombre)
            s-min (peek s)
            directo (peek (:s (svd/svd matrix)))
            via-ata (peek (sigmas-via-ata matrix))]
        (testing "Jacobi directo: σ_min con error relativo pequeño"
          (is (< (error-relativo directo s-min) 1e-3)))
        (testing "vía AᵀA: κ² supera la precisión del double y σ_min es basura"
          (is (> (error-relativo via-ata s-min) 10.0)))))))

;; ---------------------------------------------------------------------------
;; Casos especiales
;; ---------------------------------------------------------------------------

(deftest matriz-identidad
  (let [r (svd/svd (la/identity-matrix 4))]
    (is (every? #(< (Math/abs (- % 1.0)) 1e-14) (:s r)))
    (is (< (error-reconstruccion (la/identity-matrix 4) r) 1e-14))))

(deftest diagonal-con-signos-y-desordenada
  (let [a [[-1.0 0.0 0.0] [0.0 3.0 0.0] [0.0 0.0 -2.0]]
        r (svd/svd a)]
    (is (= [3.0 2.0 1.0] (:s r)) "no negativos y en orden decreciente")
    (is (< (error-reconstruccion a r) 1e-14))))

(deftest matriz-de-ceros
  (let [a (la/zero-matrix 4 3)
        r (svd/svd a)]
    (is (= [0.0 0.0 0.0] (:s r)))
    (testing "U sigue siendo ortonormal (columnas completadas)"
      (is (< (error-ortonormal (:u r)) 1e-14)))
    (is (< (error-ortonormal (:v r)) 1e-14))
    (is (= (la/zero-matrix 4 3) (svd/reconstruct r)))))

(deftest columna-de-ceros
  (let [a [[1.0 0.0 2.0] [2.0 0.0 1.0] [3.0 0.0 0.0] [0.0 0.0 1.0]]
        r (svd/svd a)]
    (is (< (peek (:s r)) 1e-14))
    (is (bien-formada? a r))
    (is (< (error-reconstruccion a r) 1e-14))))

(deftest rango-deficiente
  (testing "rango 1: un solo valor singular distinto de cero"
    (let [a (mapv (fn [x] (mapv #(* x %) [1.0 -1.0 2.0 0.5])) [1.0 2.0 3.0 4.0])
          r (svd/svd a)]
      (is (:converged? r))
      (is (< (/ (second (:s r)) (first (:s r))) 1e-14))
      (is (bien-formada? a r))
      (is (< (error-reconstruccion a r) 1e-14)))))

(deftest vectores-columna-y-fila
  (testing "columna de 5×1: σ = ‖x‖"
    (let [x [3.0 4.0 0.0 0.0 12.0]
          a (mapv vector x)
          r (svd/svd a)]
      (is (< (Math/abs (- 13.0 (first (:s r)))) 1e-14))
      (is (bien-formada? a r))
      (is (< (error-reconstruccion a r) 1e-14))))
  (testing "fila de 1×5"
    (let [a [[3.0 4.0 0.0 0.0 12.0]]
          r (svd/svd a)]
      (is (< (Math/abs (- 13.0 (first (:s r)))) 1e-14))
      (is (bien-formada? a r))
      (is (< (error-reconstruccion a r) 1e-14)))))

(deftest matriz-1x1
  (let [r (svd/svd [[-3.0]])]
    (is (= [3.0] (:s r)))
    (is (< (error-reconstruccion [[-3.0]] r) 1e-15))))

(deftest escalas-extremas
  ;; Sin el escalado por potencias de 2, α·β se desborda por abajo con entradas
  ;; ≈1e-76 y α se desborda con entradas ≈1e200. Aquí se compara siempre con la
  ;; SVD de la matriz base, dividiendo por k (no se usa norm-fro, que
  ;; desbordaría).
  (let [base [[1.0 2.0] [3.0 4.0] [5.0 6.0]]
        s-base (:s (svd/svd base))]
    (doseq [k [1e-300 1e-200 1e-76 1e-20 1e20 1e150 1e300]]
      (testing (str "k = " k)
        (let [a (mapv #(la/scale k %) base)
              r (svd/svd a)
              recon (mapv #(la/scale (/ 1.0 k) %) (svd/reconstruct r))]
          (is (:converged? r))
          (is (< (max-dif-abs (map #(/ % k) (:s r)) s-base) 1e-13))
          (is (< (error-ortonormal (:u r)) 1e-13))
          (is (< (error-ortonormal (:v r)) 1e-13))
          (is (< (max-dif-abs (flatten recon) (flatten base)) 1e-12)))))))

(deftest entrada-invalida
  (is (thrown? clojure.lang.ExceptionInfo (svd/svd [])))
  (is (thrown? clojure.lang.ExceptionInfo (svd/svd [[]])))
  (is (thrown? clojure.lang.ExceptionInfo (svd/svd [[1.0 2.0] [3.0]]))))

(deftest max-sweeps-limita-las-iteraciones
  (let [{:keys [matrix]} (caso "aleatoria-10x10")
        r (svd/svd matrix {:max-sweeps 1})]
    (is (= 1 (:sweeps r)))
    (is (not (:converged? r)))
    (testing "aun sin converger, A = U Σ Vᵀ se cumple (las rotaciones son exactas)"
      (is (< (error-reconstruccion matrix r) 1e-13)))))

(deftest inmutabilidad
  (let [a [[1.0 2.0] [3.0 4.0] [5.0 6.0]]
        copia (mapv vec a)]
    (svd/svd a)
    (is (= copia a))))

;; ---------------------------------------------------------------------------
;; Propiedades
;; ---------------------------------------------------------------------------

(def ^:private gen-num
  (gen/double* {:min -10.0 :max 10.0 :NaN? false :infinite? false}))

(def ^:private gen-matriz
  (gen/bind (gen/tuple (gen/choose 1 6) (gen/choose 1 6))
            (fn [[m n]] (gen/vector (gen/vector gen-num n) m))))

(defspec reconstruye-y-es-ortonormal 200
  (prop/for-all [a gen-matriz]
    (let [r (svd/svd a)]
      (and (:converged? r)
           (bien-formada? a r)
           (< (error-reconstruccion a r) 1e-12)))))

(defspec identidad-de-frobenius 200
  ;; ‖A‖_F² = Σ σᵢ²
  (prop/for-all [a gen-matriz]
    (let [f (la/norm-fro a)
          suma (reduce + (map #(* % %) (:s (svd/svd a))))]
      (<= (Math/abs (- (* f f) suma)) (* 1e-12 (+ 1.0 (* f f)))))))

(defspec transpuesta-tiene-los-mismos-valores-singulares 200
  (prop/for-all [a gen-matriz]
    (let [s (:s (svd/svd a))
          st (:s (svd/svd (la/transpose a)))]
      (and (= (count s) (count st))
           (<= (max-dif-abs s st) (* 1e-12 (+ 1.0 (first s))))))))

(defspec escalar-multiplica-los-valores-singulares 200
  (prop/for-all [a gen-matriz
                 k gen-num]
    (let [s (:s (svd/svd a))
          sk (:s (svd/svd (mapv #(la/scale k %) a)))]
      (<= (max-dif-abs sk (map #(* (Math/abs k) %) s))
          (* 1e-12 (+ 1.0 (first s) (Math/abs k)))))))

(defspec rango-deficiente-aleatorio 200
  ;; A = B·C con B m×r y C r×n tiene rango ≤ r: σ_{r+1..} ≈ 0 y U sigue siendo
  ;; ortonormal gracias a la completación de columnas nulas.
  (prop/for-all [[b c r] (gen/bind (gen/tuple (gen/choose 3 6) (gen/choose 3 6)
                                              (gen/choose 1 2))
                                   (fn [[m n r]]
                                     (gen/tuple (gen/vector (gen/vector gen-num r) m)
                                                (gen/vector (gen/vector gen-num n) r)
                                                (gen/return r))))]
    (let [a (la/mat-mat b c)
          res (svd/svd a)
          s (:s res)]
      (and (bien-formada? a res)
           (< (error-reconstruccion a res) 1e-12)
           (every? #(<= % (* 1e-12 (+ 1.0 (first s)))) (drop r s))))))

;; ---------------------------------------------------------------------------
;; Pruebas diferenciales: svd (rápida, arrays) frente a svd-reference (inmutable)
;; ---------------------------------------------------------------------------
;; Ambas implementaciones hacen las mismas operaciones en el mismo orden, así
;; que deben coincidir en σ, reconstrucción y ortonormalidad (tolerancia), y de
;; hecho coinciden exactamente.

(defn- hilbert [n]
  (vec (for [i (range 1 (inc n))]
         (vec (for [j (range 1 (inc n))] (/ 1.0 (+ i j -1)))))))

(defn- coinciden?
  "Verdadero si `svd` y `svd-reference` coinciden sobre la matriz `a` en σ,
  reconstrucción, ortonormalidad, barridos y convergencia."
  [a]
  (let [rapida (svd/svd a)
        ref (svd/svd-reference a)
        s1 (max (first (:s ref)) 1e-300)]
    (and (= (count (:s rapida)) (count (:s ref)))
         (< (max-dif-abs (:s rapida) (:s ref)) (* 1e-13 s1))
         (< (error-reconstruccion a rapida) 1e-13)
         (< (error-reconstruccion a ref) 1e-13)
         (< (error-ortonormal (:u rapida)) 1e-12)
         (< (error-ortonormal (:v rapida)) 1e-12)
         (< (error-ortonormal (:u ref)) 1e-12)
         (< (error-ortonormal (:v ref)) 1e-12)
         (= (:sweeps rapida) (:sweeps ref))
         (= (:converged? rapida) (:converged? ref)))))

(deftest diferencial-sobre-el-fixture
  (doseq [{:keys [name matrix]} @fixture]
    (testing name
      (is (coinciden? matrix)))))

(deftest diferencial-resultados-identicos-bit-a-bit
  (doseq [{:keys [name matrix]} @fixture]
    (testing name
      (let [rapida (svd/svd matrix)
            ref (svd/svd-reference matrix)]
        (is (= (:s rapida) (:s ref)))
        (is (= (:u rapida) (:u ref)))
        (is (= (:v rapida) (:v ref)))))))

(deftest diferencial-hilbert-numericamente-singular
  ;; Desde n = 13 el σ menor de Hilbert (≈1e-16) cae bajo el umbral de ruido
  ;; de redondeo: ejercita la ruta de columnas nulas y completación de U.
  (doseq [n [13 15]]
    (testing (str "Hilbert " n)
      (let [a (hilbert n)]
        (is (coinciden? a))
        (is (= (:s (svd/svd a)) (:s (svd/svd-reference a))))))))

(deftest diferencial-rango-deficiente-y-columnas-nulas
  (doseq [[nombre a] {"ceros 4x3" (la/zero-matrix 4 3)
                      "columna de ceros" [[1.0 0.0 2.0] [2.0 0.0 1.0] [3.0 0.0 0.0] [0.0 0.0 1.0]]
                      "rango 1" (mapv (fn [x] (mapv #(* x %) [1.0 -1.0 2.0 0.5])) [1.0 2.0 3.0 4.0])
                      "columnas repetidas" [[1.0 1.0 1.0] [2.0 2.0 2.0] [3.0 3.0 3.0]]
                      "ancha de rango 1" [[1.0 2.0 3.0 4.0] [2.0 4.0 6.0 8.0]]
                      "1x1" [[-3.0]]}]
    (testing nombre
      (is (coinciden? a)))))

(deftest diferencial-escalas-extremas
  (let [base [[1.0 2.0] [3.0 4.0] [5.0 6.0]]]
    (doseq [k [1e-300 1e-76 1e20 1e300]]
      (testing (str "k = " k)
        (let [a (mapv #(la/scale k %) base)
              rapida (svd/svd a)
              ref (svd/svd-reference a)]
          (is (= (:s rapida) (:s ref)))
          (is (= (:sweeps rapida) (:sweeps ref))))))))

(deftest diferencial-opciones-y-errores
  (let [{:keys [matrix]} (caso "aleatoria-10x10")]
    (testing "max-sweeps: ambas se detienen igual, sin converger"
      (let [rapida (svd/svd matrix {:max-sweeps 1})
            ref (svd/svd-reference matrix {:max-sweeps 1})]
        (is (= 1 (:sweeps rapida) (:sweeps ref)))
        (is (not (:converged? rapida)))
        (is (= (:s rapida) (:s ref)))))
    (testing "tol más laxa: ambas hacen los mismos barridos"
      (is (= (:sweeps (svd/svd matrix {:tol 1e-6}))
             (:sweeps (svd/svd-reference matrix {:tol 1e-6}))))))
  (testing "entrada inválida: ambas lanzan ex-info"
    (doseq [f [svd/svd svd/svd-reference]]
      (is (thrown? clojure.lang.ExceptionInfo (f [])))
      (is (thrown? clojure.lang.ExceptionInfo (f [[1.0 2.0] [3.0]]))))))

(deftest diferencial-no-modifica-la-entrada
  (let [a [[1.0 2.0] [3.0 4.0] [5.0 6.0]]
        copia (mapv vec a)]
    (svd/svd a)
    (svd/svd-reference a)
    (is (= copia a))
    (testing "llamadas repetidas dan el mismo resultado (sin estado compartido)"
      (is (= (svd/svd a) (svd/svd a))))))

(defspec diferencial-matrices-aleatorias 300
  (prop/for-all [a gen-matriz]
    (coinciden? a)))

(defspec diferencial-rango-deficiente-aleatorio 200
  (prop/for-all [a (gen/bind (gen/tuple (gen/choose 3 6) (gen/choose 3 6) (gen/choose 1 2))
                             (fn [[m n r]]
                               (gen/fmap (fn [[b c]] (la/mat-mat b c))
                                         (gen/tuple (gen/vector (gen/vector gen-num r) m)
                                                    (gen/vector (gen/vector gen-num n) r)))))]
    (coinciden? a)))

;; ---------------------------------------------------------------------------
;; truncate (aproximación de rango k)
;; ---------------------------------------------------------------------------

(deftest truncate-dimensiones-y-contenido
  (let [{:keys [matrix]} (caso "aleatoria-20x8")
        r (svd/svd matrix)
        t (svd/truncate r 3)]
    (is (= [20 3] (la/dims (:u t))))
    (is (= [8 3] (la/dims (:v t))))
    (is (= (subvec (:s r) 0 3) (:s t)))
    (testing "conserva las primeras columnas, sin alterarlas"
      (is (= (mapv #(subvec % 0 3) (:u r)) (:u t)))
      (is (= (mapv #(subvec % 0 3) (:v r)) (:v t))))
    (testing "las columnas truncadas siguen siendo ortonormales"
      (is (< (error-ortonormal (:u t)) 1e-12))
      (is (< (error-ortonormal (:v t)) 1e-12)))))

(deftest truncate-con-k-completo-reconstruye-todo
  (let [{:keys [matrix]} (caso "aleatoria-10x10")
        r (svd/svd matrix)]
    (is (= (svd/reconstruct r) (svd/reconstruct (svd/truncate r 10))))
    (is (< (error-reconstruccion matrix (svd/truncate r 10)) 1e-13))))

(deftest truncate-k-fuera-de-rango
  (let [r (svd/svd [[1.0 2.0] [3.0 4.0] [5.0 6.0]])]
    (is (thrown? clojure.lang.ExceptionInfo (svd/truncate r 0)))
    (is (thrown? clojure.lang.ExceptionInfo (svd/truncate r 3)) "solo hay 2 valores singulares")
    (is (thrown? clojure.lang.ExceptionInfo (svd/truncate r -1)))
    (is (thrown? clojure.lang.ExceptionInfo (svd/truncate r 1.5)))
    (is (= 1 (count (:s (svd/truncate r 1)))))
    (is (= 2 (count (:s (svd/truncate r 2)))))))

(deftest truncate-rango-1-de-una-matriz-de-rango-1-es-exacta
  (let [a (mapv (fn [x] (mapv #(* x %) [1.0 -1.0 2.0 0.5])) [1.0 2.0 3.0 4.0])
        r (svd/svd a)]
    (is (< (error-reconstruccion a (svd/truncate r 1)) 1e-14))))

;; ---------------------------------------------------------------------------
;; Número de condición, rango numérico y pseudoinversa
;; ---------------------------------------------------------------------------

(defn- frob-dif [a b] (la/norm-fro (mapv la/v- a b)))

(defn- cumple-moore-penrose?
  "Las cuatro condiciones de Moore–Penrose entre A y su pseudoinversa A⁺, con una
  tolerancia que crece con κ_eff = ‖A‖_F‖A⁺‖_F (el error de redondeo esperado
  es del orden de ε·κ_eff)."
  [a ap]
  (let [na (la/norm-fro a)
        nap (la/norm-fro ap)
        k-eff (* na nap)
        tol (* 1e-12 (+ 1.0 k-eff) (+ 1.0 na nap))
        aap (la/mat-mat a ap)
        apa (la/mat-mat ap a)]
    (and (< (frob-dif (la/mat-mat aap a) a) tol)
         (< (frob-dif (la/mat-mat apa ap) ap) tol)
         (< (frob-dif aap (la/transpose aap)) tol)
         (< (frob-dif apa (la/transpose apa)) tol))))

(deftest pinv-moore-penrose-en-todo-el-fixture
  (doseq [{:keys [name matrix]} @fixture]
    (testing name
      (let [ap (svd/pinv matrix)]
        (is (= [(second (la/dims matrix)) (first (la/dims matrix))] (la/dims ap)))
        (is (cumple-moore-penrose? matrix ap))))))

;; Las propiedades de pinv usan entradas redondeadas a milésimas: con entradas
;; arbitrariamente pequeñas (p. ej. 1e-162) pinv devuelve entradas ≈1e161,
;; correctas, pero `la/norm-fro` desborda (eleva al cuadrado) al medirlas.
(def ^:private gen-num-milesimas
  (gen/fmap #(/ (Math/rint (* % 1000.0)) 1000.0) gen-num))

(def ^:private gen-matriz-milesimas
  (gen/bind (gen/tuple (gen/choose 1 6) (gen/choose 1 6))
            (fn [[m n]] (gen/vector (gen/vector gen-num-milesimas n) m))))

(defspec pinv-moore-penrose-aleatorias 200
  (prop/for-all [a gen-matriz-milesimas]
    (cumple-moore-penrose? a (svd/pinv a))))

(defspec pinv-moore-penrose-rango-deficiente 200
  (prop/for-all [a (gen/bind (gen/tuple (gen/choose 3 6) (gen/choose 3 6) (gen/choose 1 2))
                             (fn [[m n r]]
                               (gen/fmap (fn [[b c]] (la/mat-mat b c))
                                         (gen/tuple (gen/vector (gen/vector gen-num-milesimas r) m)
                                                    (gen/vector (gen/vector gen-num-milesimas n) r)))))]
    (cumple-moore-penrose? a (svd/pinv a))))

(deftest pinv-casos-concretos
  (testing "diagonal: invierte las entradas distintas de cero"
    (is (every? true? (map #(< (max-dif-abs %1 %2) 1e-14)
                           (svd/pinv [[2.0 0.0] [0.0 4.0]])
                           [[0.5 0.0] [0.0 0.25]])))
    (is (every? true? (map #(< (max-dif-abs %1 %2) 1e-14)
                           (svd/pinv [[2.0 0.0] [0.0 0.0]])
                           [[0.5 0.0] [0.0 0.0]]))))
  (testing "matriz invertible: A⁺ = A⁻¹"
    (let [a [[2.0 1.0] [1.0 3.0]]
          prod (la/mat-mat a (svd/pinv a))]
      (is (< (frob-dif prod (la/identity-matrix 2)) 1e-14))))
  (testing "alta de columnas independientes: A⁺ es inversa por la izquierda"
    (let [a [[1.0 2.0] [3.0 4.0] [5.0 6.0]]]
      (is (< (frob-dif (la/mat-mat (svd/pinv a) a) (la/identity-matrix 2)) 1e-13))))
  (testing "ancha de filas independientes: A⁺ es inversa por la derecha"
    (let [a [[1.0 2.0 3.0] [4.0 5.0 6.0]]]
      (is (< (frob-dif (la/mat-mat a (svd/pinv a)) (la/identity-matrix 2)) 1e-13))))
  (testing "matriz de ceros: A⁺ = 0 con las dimensiones traspuestas"
    (is (= (la/zero-matrix 3 4) (svd/pinv (la/zero-matrix 4 3)))))
  (testing "acepta una matriz o el resultado de svd, con el mismo resultado"
    (let [a [[1.0 2.0] [3.0 4.0] [5.0 6.0]]]
      (is (= (svd/pinv a) (svd/pinv (svd/svd a)))))))

(deftest pinv-rcond-regulariza
  (let [a [[1.0 0.0] [0.0 1e-8]]]
    (testing "por defecto se invierte 1e-8 (está por encima de max(m,n)·ε)"
      (is (< (max-dif-abs (flatten (svd/pinv a)) (flatten [[1.0 0.0] [0.0 1e8]])) 1.0)))
    (testing "con rcond = 1e-6 se descarta el valor singular pequeño"
      (is (= [[1.0 0.0] [0.0 0.0]] (mapv #(mapv (fn [x] (+ 0.0 x)) %) (svd/pinv a 1e-6)))))))

(deftest inverse-singular-values-pone-ceros-bajo-la-tolerancia
  ;; σ = [4 2 0] (orden decreciente) → inversos [1/4 1/2 0]
  (is (= [0.25 0.5 0.0]
         (let [r (svd/svd [[4.0 0.0 0.0] [0.0 2.0 0.0] [0.0 0.0 0.0]])]
           (svd/inverse-singular-values r))))
  (is (= [1.0 0.0] (mapv #(+ 0.0 %) (svd/inverse-singular-values [[1.0 0.0] [0.0 1e-8]] 1e-6)))))

(deftest numero-de-condicion
  (is (< (Math/abs (- 4.0 (svd/condition-number [[4.0 0.0 0.0] [0.0 -2.0 0.0] [0.0 0.0 1.0]]))) 1e-13))
  (is (< (Math/abs (- 1.0 (svd/condition-number (la/identity-matrix 5)))) 1e-14))
  (testing "matriz singular: κ enorme o infinito"
    (is (> (svd/condition-number [[1.0 2.0] [2.0 4.0]]) 1e15))
    (is (= Double/POSITIVE_INFINITY (svd/condition-number (la/zero-matrix 3 3)))))
  (testing "Hilbert 8 coincide con la referencia de NumPy"
    (let [{:keys [matrix cond]} (caso "hilbert-8")]
      (is (< (error-relativo (svd/condition-number matrix) cond) 1e-3))))
  (testing "acepta un resultado de svd"
    (let [a [[3.0 0.0] [0.0 1.0]]]
      (is (= (svd/condition-number a) (svd/condition-number (svd/svd a)))))))

(deftest rango-numerico
  (is (= 5 (svd/numerical-rank (la/identity-matrix 5))))
  (is (= 0 (svd/numerical-rank (la/zero-matrix 3 4))))
  (is (= 1 (svd/numerical-rank (:matrix (caso "rango-1-4x4")))))
  (is (= 2 (svd/numerical-rank (:matrix (caso "rango-2-6x5")))))
  (is (= 3 (svd/numerical-rank [[1.0 2.0 3.0 4.0] [2.0 4.0 6.0 8.0] [0.0 1.0 0.0 1.0] [1.0 0.0 0.0 0.0]]))
      "una fila es múltiplo de otra → rango 3")
  (is (= 8 (svd/numerical-rank (:matrix (caso "hilbert-8")))))
  (testing "Hilbert: rango completo hasta n = 12; desde n = 13 el σ menor cae bajo el ruido"
    (is (= 10 (svd/numerical-rank (hilbert 10))))
    (is (= 12 (svd/numerical-rank (hilbert 12))))
    (doseq [n [13 14 15]]
      (is (< (svd/numerical-rank (hilbert n)) n) (str "Hilbert " n))))
  (testing "rcond: con una tolerancia mayor se cuentan menos direcciones"
    (let [a [[1.0 0.0] [0.0 1e-8]]]
      (is (= 2 (svd/numerical-rank a)))
      (is (= 1 (svd/numerical-rank a 1e-6)))))
  (testing "matrices anchas"
    (is (= 2 (svd/numerical-rank [[1.0 2.0 3.0] [4.0 5.0 6.0]])))))
