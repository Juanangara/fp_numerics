;; Genera los CSV de precisión de la SVD por Jacobi de un solo lado.
;;
;; Respalda en archivos las cifras que el README cita sobre la exactitud de
;; `fpnum.svd.core/svd`:
;;   docs/data/svd-accuracy.csv             los 12 casos del fixture de NumPy:
;;                                          σ frente a NumPy, error de
;;                                          reconstrucción y ortonormalidad
;;   docs/data/svd-ata-vs-directo.csv       σ mínimo por Jacobi directo frente a
;;                                          pasar por AᵀA (por qué no se forma AᵀA)
;;   docs/data/svd-hilbert-determinante.csv el producto de los σ de Hilbert frente
;;                                          a su determinante exacto
;;
;; Necesita el fixture test/fpnum/svd/svd_reference.edn (valores singulares de
;; NumPy; se regenera con scripts/generate_svd_fixture.py). El código de src/
;; no usa NumPy.
;;
;; Uso, desde la raíz del proyecto:
;;   clojure -M -i scripts/generate_accuracy_csvs.clj
;; Es determinista: dos ejecuciones dan archivos idénticos.

(require '[clojure.edn :as edn]
         '[clojure.java.io :as io]
         '[clojure.string :as str]
         '[fpnum.linear :as la]
         '[fpnum.svd.core :as svd])

(def fixture-path "test/fpnum/svd/svd_reference.edn")
(def data-dir "docs/data")

(def casos (edn/read-string (slurp fixture-path)))
(def por-nombre (into {} (map (juxt :name identity)) casos))

;; ---------------------------------------------------------------------------
;; Medidas
;; ---------------------------------------------------------------------------

(defn error-ortonormal
  "Mayor |(QᵀQ − I)ᵢⱼ|: cuánto se aleja Q de tener columnas ortonormales."
  [q]
  (let [g (la/mat-mat (la/transpose q) q)
        n (count g)]
    (apply max (for [i (range n) j (range n)]
                 (Math/abs (- (get-in g [i j]) (if (= i j) 1.0 0.0)))))))

(defn error-reconstruccion
  "‖A − U Σ Vᵀ‖_F / ‖A‖_F."
  [a r]
  (let [num (la/norm-fro (mapv la/v- a (svd/reconstruct r)))
        den (la/norm-fro a)]
    (if (zero? den) num (/ num den))))

(defn max-dif-abs [x y]
  (apply max 0.0 (map #(Math/abs (- (double %1) (double %2))) x y)))

(defn error-relativo [x ref]
  (/ (Math/abs (- x ref)) (Math/abs ref)))

(defn sigmas-via-ata
  "Valores singulares como √(valores singulares de AᵀA): el camino que se evita,
  porque κ(AᵀA) = κ(A)²."
  [a]
  (mapv #(Math/sqrt (max 0.0 %))
        (:s (svd/svd (la/mat-mat (la/transpose a) a)))))

(defn det-hilbert
  "Determinante exacto de la matriz de Hilbert n×n: det Hₙ = cₙ⁴ / c₂ₙ con
  cₙ = Π_{i=1}^{n−1} i!, calculado con racionales."
  [n]
  (let [fact (fn [k] (reduce * 1N (range 1 (inc k))))
        c (fn [k] (reduce * 1N (map fact (range 1 k))))]
    (double (/ (* (c n) (c n) (c n) (c n)) (c (* 2 n))))))

;; ---------------------------------------------------------------------------
;; Tablas
;; ---------------------------------------------------------------------------

(defn fila-accuracy [{:keys [name matrix s cond]}]
  (let [[m n] (la/dims matrix)
        r (svd/svd matrix)]
    [name m n (or cond "")
     (/ (max-dif-abs (:s r) s) (first s))
     (error-reconstruccion matrix r)
     (error-ortonormal (:u r))
     (error-ortonormal (:v r))
     (:sweeps r)
     (:converged? r)]))

(def casos-ata ["cuadrada-3x3" "aleatoria-10x10" "decaimiento-geometrico-12x12"
                "hilbert-5" "hilbert-8" "hilbert-10"])

(defn fila-ata [nombre]
  (let [{:keys [matrix s cond]} (por-nombre nombre)
        s-min (peek s)
        directo (peek (:s (svd/svd matrix)))
        via-ata (peek (sigmas-via-ata matrix))]
    [nombre cond s-min directo (error-relativo directo s-min) via-ata (error-relativo via-ata s-min)]))

(defn fila-determinante [n]
  (let [{:keys [matrix cond]} (por-nombre (str "hilbert-" n))
        producto (reduce * (:s (svd/svd matrix)))
        exacto (det-hilbert n)]
    [n cond exacto producto (error-relativo producto exacto)]))

(defn escribir-csv! [nombre cabecera filas]
  (let [path (str data-dir "/" nombre)]
    (io/make-parents path)
    (spit path (str/join "\n" (cons cabecera (map #(str/join "," %) filas))))
    (println "escrito" path (str "(" (count filas) " filas)"))))

;; ---------------------------------------------------------------------------

(println "== Precisión de la SVD por Jacobi: CSV de respaldo ==")
(let [accuracy (mapv fila-accuracy casos)
      ata (mapv fila-ata casos-ata)
      det (mapv fila-determinante [5 8 10])]
  (escribir-csv! "svd-accuracy.csv"
                 "case,m,n,cond_numpy,max_abs_diff_over_sigma1,rel_recon_error,orth_u,orth_v,sweeps,converged"
                 accuracy)
  (escribir-csv! "svd-ata-vs-directo.csv"
                 "case,cond_numpy,sigma_min_numpy,sigma_min_direct,rel_err_direct,sigma_min_via_ata,rel_err_via_ata"
                 ata)
  (escribir-csv! "svd-hilbert-determinante.csv"
                 "n,cond_numpy,det_exact,product_of_sigmas,rel_error"
                 det)
  (println)
  (println (format "Frente a NumPy: máxima diferencia de σ / σ₁ en los %d casos = %.1e"
                   (count accuracy) (apply max (map #(nth % 4) accuracy))))
  (println (format "Reconstrucción: error relativo máximo = %.1e | ortonormalidad máxima = %.1e"
                   (apply max (map #(nth % 5) accuracy)) (apply max (mapcat #(subvec % 6 8) accuracy))))
  (doseq [[nombre _ _ _ ed _ ea] ata]
    (println (format "  %-30s σ mínimo, error relativo: directo %.1e | vía AᵀA %.1e" nombre ed ea)))
  (doseq [[n _ _ _ e] det]
    (println (format "  Hilbert %2d: producto de σ frente al determinante exacto, error relativo %.1e" n e))))
(shutdown-agents)
