(ns fpnum.demo.reverse
  "Demo ejecutable: diferenciación automática forward frente a reverse.

  Teoría
  ------
  Forward obtiene una derivada parcial por pasada, así que un gradiente de n
  variables cuesta n evaluaciones de f: O(n · costo(f)). Reverse obtiene todo el
  gradiente con una pasada hacia adelante y una hacia atrás: O(costo(f)),
  sin depender de n. La razón forward/reverse debe crecer con n.

  Experimento
  -----------
  Tiempo de un gradiente con `fpnum.ad.forward/grad` y `fpnum.ad.reverse/grad`,
  para n = 10, 50, 100 y 200 variables, con dos funciones escritas UNA sola
  vez con las operaciones genéricas (sirven en los dos modos):

  * mínimos cuadrados ‖X·w − y‖² con m = 100 filas: costo O(m·n), así que
    forward cuesta O(m·n²) y reverse O(m·n).
  * suma de cuadrados Σ wᵢ²: costo O(n), así que forward cuesta O(n²) y
    reverse O(n). Es la forma más limpia de ver la ley de escala.

  Cada muestra repite el gradiente las veces necesarias para durar al menos
  ~20 ms (los gradientes pequeños duran microsegundos y el reloj no los
  resolvería), y se reporta el tiempo POR gradiente: mediana y rango de varias
  muestras, tras calentar la JVM. Las muestras de forward y de reverse se
  INTERCALAN, y la razón se calcula por rondas pareadas (ver `time-pair`):
  el rendimiento de la máquina cambia entre corridas (hasta 2.5 veces) y
  dentro de una corrida, y solo las razones pareadas son comparables. También
  se comprueba que los dos modos dan el mismo gradiente.

  Salida: docs/data/ad-forward-vs-reverse.csv

  Uso: clojure -M -m fpnum.demo.reverse"
  (:require [clojure.java.io :as io]
            [clojure.string :as str]
            [fpnum.ad.forward :as fwd]
            [fpnum.ad.ops :as ad]
            [fpnum.ad.reverse :as rev]
            [fpnum.optim :as opt]))

(def default-ns
  "Números de variables n del experimento."
  [10 50 100 200])

(def default-m
  "Filas m del problema de mínimos cuadrados (fijas mientras varía n)."
  100)

(def csv-path
  "CSV con el tiempo de un gradiente por forward y por reverse, para cada n."
  "docs/data/ad-forward-vs-reverse.csv")

;; ---------------------------------------------------------------------------
;; Problemas (escritos una vez con las operaciones genéricas)
;; ---------------------------------------------------------------------------

(defn gaussian-matrix
  "Matriz m×n con entradas gaussianas estándar, reproducible con `seed`."
  [m n seed]
  (let [rnd (java.util.Random. (long seed))]
    (vec (repeatedly m #(vec (repeatedly n (fn [] (.nextGaussian rnd))))))))

(defn least-squares-case
  "Objetivo de mínimos cuadrados ‖X·w − y‖² con X de m×n. Devuelve
  {:name :f :point :n :m}; el punto es w = (0.1, ..., 0.1)."
  [m n]
  (let [x (gaussian-matrix m n 1)
        y (mapv first (gaussian-matrix m 1 2))]
    {:name "minimos-cuadrados" :n n :m m
     :f (opt/least-squares-objective x y)
     :point (vec (repeat n 0.1))}))

(defn sum-of-squares-case
  "Suma de cuadrados Σ wᵢ² (gradiente 2w), en el punto wᵢ = 0.01·i."
  [n]
  {:name "suma-de-cuadrados" :n n :m nil
   :f (fn [w] (apply ad/d+ (map (fn [wi] (ad/d* wi wi)) w)))
   :point (mapv #(* 0.01 (inc %)) (range n))})

;; ---------------------------------------------------------------------------
;; Medición
;; ---------------------------------------------------------------------------

(defn median
  "Mediana de una colección de números (promedio de los dos centrales si hay
  una cantidad par)."
  [xs]
  (let [s (vec (sort xs))
        n (count s)]
    (if (odd? n)
      (s (quot n 2))
      (/ (+ (s (dec (quot n 2))) (s (quot n 2))) 2.0))))

(defn- ms-since [t0] (/ (- (System/nanoTime) t0) 1e6))

(defn time-per-call
  "Milisegundos POR llamada de `(thunk)`: calibra cuántas repeticiones hacen
  falta para que cada muestra dure al menos `min-sample-ms`, calienta con
  `warmup` muestras y devuelve {:median :min :max :reps} de `runs` muestras."
  [thunk {:keys [runs warmup min-sample-ms] :or {runs 5 warmup 2 min-sample-ms 20.0}}]
  (let [t0 (System/nanoTime)
        _ (thunk)
        una (max 1e-3 (ms-since t0))
        reps (max 1 (long (Math/ceil (/ min-sample-ms una))))
        muestra (fn [] (let [t (System/nanoTime)]
                         (dotimes [_ reps] (thunk))
                         (/ (ms-since t) reps)))]
    (dotimes [_ warmup] (muestra))
    (let [ms (vec (repeatedly runs muestra))]
      {:median (median ms) :min (apply min ms) :max (apply max ms) :reps reps})))

(defn- calibrar
  "Repeticiones necesarias para que una muestra de `thunk` dure al menos
  `min-sample-ms`."
  [thunk min-sample-ms]
  (let [t0 (System/nanoTime)
        _ (thunk)
        una (max 1e-3 (ms-since t0))]
    (max 1 (long (Math/ceil (/ min-sample-ms una))))))

(defn- una-muestra
  "Milisegundos por llamada de una muestra de `reps` repeticiones."
  [thunk reps]
  (let [t (System/nanoTime)]
    (dotimes [_ reps] (thunk))
    (/ (ms-since t) reps)))

(defn- estadisticas [ms reps]
  {:median (median ms) :min (apply min ms) :max (apply max ms) :reps reps})

(defn time-pair
  "Mide dos funciones sin argumentos, `thunk-a` y `thunk-b`, con las muestras
  INTERCALADAS (ronda par: A y luego B; ronda impar: B y luego A).

  Es necesario porque en este tipo de máquina (un portátil) el rendimiento
  cambia con la carga y la temperatura, hasta 2.5 veces entre una corrida y
  otra y también dentro de una corrida: si se midiera todo A y luego todo B,
  cada uno vería un estado distinto. Con las muestras intercaladas, las dos
  de una misma ronda comparten el mismo estado, y la razón se calcula POR
  RONDA.

  Devuelve {:a {:median :min :max :reps} :b {...} :ratios [a/b por ronda]}, con
  los tiempos en ms por llamada."
  [thunk-a thunk-b {:keys [runs warmup min-sample-ms] :or {runs 7 warmup 2 min-sample-ms 20.0}}]
  (let [ra (calibrar thunk-a min-sample-ms)
        rb (calibrar thunk-b min-sample-ms)]
    (dotimes [_ warmup] (una-muestra thunk-a ra) (una-muestra thunk-b rb))
    (let [rondas (vec (for [i (range runs)]
                        (if (even? i)
                          (let [a (una-muestra thunk-a ra) b (una-muestra thunk-b rb)] [a b])
                          (let [b (una-muestra thunk-b rb) a (una-muestra thunk-a ra)] [a b]))))
          as (mapv first rondas)
          bs (mapv second rondas)]
      {:a (estadisticas as ra)
       :b (estadisticas bs rb)
       :ratios (mapv (fn [[a b]] (/ a b)) rondas)})))

(defn warm-up!
  "Calienta la JVM: ejecuta durante `ms` milisegundos gradientes de tamaño medio
  en los dos modos, para que HotSpot haya compilado ambos caminos antes de medir.
  Sin esto, las primeras filas de la tabla miden código aún sin optimizar."
  [ms]
  (doseq [{:keys [f point]} [(least-squares-case 100 50) (sum-of-squares-case 100)]
          grad [fwd/grad rev/grad]]
    (let [fin (+ (System/nanoTime) (* ms 1000000))]
      (while (< (System/nanoTime) fin)
        (grad f point)))))

(defn max-rel-diff
  "Mayor diferencia relativa |a − b|/(1 + |a|) entre dos gradientes."
  [a b]
  (apply max 0.0 (map (fn [u v] (/ (Math/abs (- (double u) (double v))) (+ 1.0 (Math/abs (double u))))) a b)))

(defn scaling-row
  "Mide forward y reverse sobre un caso ({:name :f :point :n :m}) con muestras
  intercaladas (ver `time-pair`) y devuelve la fila
  {:function :n :m :forward :reverse :ratio :ratio-range :max-rel-diff}.

  :forward y :reverse son {:median :min :max :reps} en ms por gradiente.
  :ratio es la MEDIANA de las razones forward/reverse de cada ronda (no la razón
  de las medianas), y :ratio-range es [mínima máxima] de esas razones."
  [{:keys [name f point n m]} opts]
  (let [{:keys [a b ratios]} (time-pair #(fwd/grad f point) #(rev/grad f point) opts)]
    {:function name :n n :m m
     :forward a :reverse b
     :ratio (median ratios)
     :ratio-range [(apply min ratios) (apply max ratios)]
     :max-rel-diff (max-rel-diff (fwd/grad f point) (rev/grad f point))}))
;; ---------------------------------------------------------------------------
;; Salida
;; ---------------------------------------------------------------------------

(defn scaling-csv-lines
  "Líneas CSV de las filas de `scaling-row`."
  [rows]
  (cons "function,n,m,forward_ms,forward_min_ms,forward_max_ms,reverse_ms,reverse_min_ms,reverse_max_ms,ratio_forward_over_reverse,ratio_min,ratio_max,max_rel_diff"
        (for [{:keys [function n m forward reverse ratio ratio-range max-rel-diff]} rows]
          (str/join "," [function n (or m "")
                         (:median forward) (:min forward) (:max forward)
                         (:median reverse) (:min reverse) (:max reverse)
                         ratio (first ratio-range) (second ratio-range) max-rel-diff]))))

(defn- demo-misma-funcion []
  (println "-- La misma función en los dos modos --")
  (let [g (fn [[x y]] (ad/d+ (ad/d* x x y) (ad/dsin y)))]
    (println "g(x, y) = x²y + sin y en (1, 2)")
    (println "  forward:" (fwd/grad g [1.0 2.0]))
    (println "  reverse:" (rev/grad g [1.0 2.0]))
    (println "  exacto: " [(* 2.0 1.0 2.0) (+ 1.0 (Math/cos 2.0))])))

(defn run-experiment!
  "Ejecuta el experimento contra n, imprime la tabla y escribe el CSV.
  Opciones: :ns (tamaños), :m (filas de mínimos cuadrados), :runs, :warmup,
  :min-sample-ms, :warm-up-ms (calentamiento global previo, 1500 por defecto) y
  :csv-path. Devuelve las filas."
  ([] (run-experiment! {}))
  ([{:keys [ns m csv-path] :or {ns default-ns m default-m} :as opts}]
   (let [csv-path (or csv-path fpnum.demo.reverse/csv-path)
         casos (concat (for [n ns] (least-squares-case m n))
                       (for [n ns] (sum-of-squares-case n)))]
     (println "Calentando la JVM...")
     (warm-up! (get opts :warm-up-ms 1500))
     (println "-- Tiempo de un gradiente: forward frente a reverse --")
     (println "(ms por gradiente, mediana [mín–máx]; razón = mediana de forward/reverse por ronda [mín–máx])")
     (println "función            |   n |    m | forward                  | reverse                  | razón            | dif. máx.")
     (let [rows (vec (for [caso casos]
                       (let [r (scaling-row caso opts)
                             f (:forward r)
                             v (:reverse r)]
                         (println (format "%-18s | %3d | %4s | %8.3f [%.3f–%.3f] | %8.3f [%.3f–%.3f] | %5.2f [%.2f–%.2f] | %.1e"
                                          (:function r) (:n r) (or (:m r) "-")
                                          (:median f) (:min f) (:max f)
                                          (:median v) (:min v) (:max v)
                                          (:ratio r) (first (:ratio-range r)) (second (:ratio-range r))
                                          (:max-rel-diff r)))
                         r)))]
       (io/make-parents csv-path)
       (spit csv-path (str/join "\n" (scaling-csv-lines rows)))
       (println "\nCSV escrito en" csv-path)
       rows))))

(defn -main
  "Compara forward y reverse: la misma función en los dos modos y el tiempo de
  un gradiente para n = 10, 50, 100 y 200."
  [& _args]
  (println "== AD forward frente a reverse ==\n")
  (demo-misma-funcion)
  (println)
  (run-experiment!))
