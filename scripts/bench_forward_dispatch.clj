;; Medición A/B: ¿se volvió más lento el modo forward al añadir la delegación al
;; protocolo ADValue (refactor a operaciones genéricas)?
;;
;; Compara, en la MISMA JVM y con las rondas intercaladas (orden ABBA para
;; cancelar derivas de la CPU), dos versiones del namespace forward:
;;   before: scripts/baseline/forward_baseline.clj  (copia congelada del forward original)
;;   after : src/fpnum/ad/forward.clj             (con la rama de delegación)
;;
;; Criterio (acordado): solo es regresión si los rangos [min, max] de las dos
;; versiones NO se solapan Y la mediana de `after` supera en más de un 10 % a la
;; de `before`.
;;
;; Uso, desde la raíz del proyecto:
;;   clojure -M -i scripts/bench_forward_dispatch.clj
;; Escribe docs/data/ad-forward-dispatch.csv

(load-file "scripts/baseline/forward_baseline.clj")

(require '[clojure.java.io :as io]
         '[clojure.string :as str]
         '[fpnum.ad.forward :as after]
         '[fpnum.ad.forward-baseline :as before]
         '[fpnum.linear :as la])

(def rondas 15)
(def csv-path "docs/data/ad-forward-dispatch.csv")

(defn ops-de [ns-sym]
  (into {} (for [s '[d+ d- d* d-div dsin dexp dpow grad diff]]
             [(keyword (name s)) @(ns-resolve ns-sym s)])))

(def variantes {:before (ops-de 'fpnum.ad.forward-baseline)
                :after  (ops-de 'fpnum.ad.forward)})

(defn gauss-matriz [m n seed]
  (let [r (java.util.Random. (long seed))]
    (vec (repeatedly m #(vec (repeatedly n (fn [] (.nextGaussian r))))))))

(def X (gauss-matriz 40 4 1))
(def y (mapv first (gauss-matriz 40 1 2)))
(def w0 [0.1 0.2 0.3 0.4])

(defn objetivo-ls [{:keys [d+ d- d*]}]
  (fn [w]
    (apply d+ (map (fn [fila yi]
                     (let [r (d- (apply d+ (map d* fila w)) yi)]
                       (d* r r)))
                   X y))))

(defn rosenbrock [{:keys [d+ d- d*]}]
  (fn [[x y]]
    (let [u (d- 1 x)
          v (d- y (d* x x))]
      (d+ (d* u u) (d* 100 v v)))))

;; Cada caso devuelve un thunk que hace una cantidad fija de trabajo y un
;; resultado (suma de control) que debe ser idéntico en las dos versiones.
(def casos
  [["grad-minimos-cuadrados-40x4"
    (fn [{:keys [grad] :as o}]
      (let [f (objetivo-ls o)]
        #(reduce (fn [acc i] (+ acc (first (grad f (assoc w0 0 (* 0.001 i))))))
                 0.0 (range 2000))))]
   ["diff-funcion-compuesta"
    (fn [{:keys [d+ d* d-div dsin dexp diff]}]
      (let [f (fn [x] (d-div (dexp (dsin x)) (d+ 1.0 (d* x x))))]
        #(reduce (fn [acc i] (+ acc (diff f (* 0.00001 i)))) 0.0 (range 300000))))]
   ["descenso-3000-iteraciones-ls"
    ;; paso 0.004 < 1/||X||_F^2 ~ 0.0056: converge (con 0.05 diverge a NaN y mide aritmetica no representativa)
    (fn [{:keys [grad] :as o}]
      (let [f (objetivo-ls o)]
        #(first (reduce (fn [w _] (mapv (fn [wi gi] (- wi (* 0.004 gi))) w (grad f w)))
                        w0 (range 3000)))))]
   ["grad-rosenbrock-100000"
    (fn [{:keys [grad] :as o}]
      (let [f (rosenbrock o)]
        #(reduce (fn [acc i] (+ acc (first (grad f [(* 0.00001 i) 1.0])))) 0.0 (range 100000))))]])

(defn medir [thunk]
  (let [t0 (System/nanoTime)
        r (thunk)]
    [(/ (- (System/nanoTime) t0) 1e6) r]))

(defn mediana [xs]
  (let [s (vec (sort xs)) n (count s)]
    (if (odd? n) (s (quot n 2)) (/ (+ (s (dec (quot n 2))) (s (quot n 2))) 2.0))))

(println "== A/B del despacho en forward ==")
(println "Calentando la JVM con las dos versiones...")
(doseq [[_ mk] casos, v [:before :after]]
  (let [t (mk (variantes v))] (dotimes [_ 3] (t))))

(def resultados
  (into {}
        (for [[nombre mk] casos]
          (let [thunks {:before (mk (variantes :before)) :after (mk (variantes :after))}
                muestras (atom {:before [] :after []})
                checks (atom {})]
            (doseq [ronda (range rondas)]
              (doseq [v (if (even? ronda) [:before :after] [:after :before])]
                (let [[ms r] (medir (thunks v))]
                  (swap! muestras update v conj ms)
                  (swap! checks assoc v r))))
            [nombre {:muestras @muestras :checks @checks}]))))

(defn estadisticas [xs] {:n (count xs) :mediana (mediana xs) :min (apply min xs) :max (apply max xs)})

(doseq [[nombre {:keys [checks]}] resultados, [v r] checks]
  (when-not (and (number? r) (Double/isFinite (double r)))
    (throw (ex-info "resultado no finito: el caso no mide trabajo representativo" {:caso nombre :variante v :resultado r}))))

(def filas
  (vec (for [[nombre {:keys [muestras checks]}] resultados
             :let [b (estadisticas (:before muestras))
                   a (estadisticas (:after muestras))
                   razon (/ (:mediana a) (:mediana b))
                   solapan? (and (<= (:min a) (:max b)) (<= (:min b) (:max a)))
                   veredicto (cond (and (not solapan?) (> razon 1.10)) "REGRESION"
                                   (and (not solapan?) (< razon 0.90)) "mas rapido"
                                   :else "sin diferencia (rangos solapados o <10%)")]]
         {:caso nombre :before b :after a :razon razon :solapan? solapan?
          :veredicto veredicto :identicos? (= (:before checks) (:after checks))})))

(doseq [{:keys [caso before after razon solapan? veredicto identicos?]} filas]
  (println (format "\n%s" caso))
  (println (format "  before: mediana %8.1f ms  rango [%.1f, %.1f]" (:mediana before) (:min before) (:max before)))
  (println (format "  after : mediana %8.1f ms  rango [%.1f, %.1f]" (:mediana after) (:min after) (:max after)))
  (println (format "  after/before = %.3f | rangos solapados: %s | resultados identicos: %s | => %s"
                   razon solapan? identicos? veredicto)))

(io/make-parents csv-path)
(spit csv-path
      (str/join "\n"
                (cons "case,variant,runs,median_ms,min_ms,max_ms,median_ratio_after_over_before,ranges_overlap,verdict"
                      (mapcat (fn [{:keys [caso before after razon solapan? veredicto]}]
                                [(str/join "," [caso "before" (:n before) (:mediana before) (:min before) (:max before) "" "" ""])
                                 (str/join "," [caso "after" (:n after) (:mediana after) (:min after) (:max after) razon solapan? veredicto])])
                              filas))))
(println "\nCSV escrito en" csv-path)
(println (if (some #(= "REGRESION" (:veredicto %)) filas)
           "\nRESULTADO: HAY REGRESION en al menos un caso"
           "\nRESULTADO: sin regresion (ningun caso con rangos disjuntos y >10% mas lento)"))
(shutdown-agents)
