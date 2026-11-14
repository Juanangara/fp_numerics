;; Medición de reverse mode con subexpresiones muy compartidas y cadenas largas.
;;
;; Se escribió para evaluar la identidad estructural de nodos (enfoque 1a) y se
;; conserva para la identidad por contador (1b, la vigente). Resultados:
;;   docs/data/ad-reverse-sharing-1a.csv  (1a: StackOverflow desde ~5000 y S2 exponencial)
;;   docs/data/ad-reverse-sharing-1b.csv  (1b: lineal en los tres casos)
;; Ver docs/decisiones.md.
;;
;; Casos:
;;   L1  cadena lineal larga:  acc <- acc*0.999 + x,  n términos, construida UNA vez.
;;       Detecta StackOverflow en hash o igualdad de ids anidados.
;;   S1  cadena de cuadrados ((x²)²)²... de profundidad d, construida UNA vez.
;;       Como árbol tendría 2^d caminos; con nodos compartidos debe ser lineal.
;;   S2  la misma cadena de cuadrados construida DOS veces y sumada. Los dos
;;       nodos finales son iguales pero distintos como objetos: obligan a
;;       comparar ids anidados iguales.
;;
;; Uso, desde la raíz del proyecto:
;;   clojure -M -i scripts/measure_reverse_sharing.clj
;; Escribe docs/data/ad-reverse-sharing-1b.csv (cambiar `csv-path` para otra estrategia)

(require '[clojure.java.io :as io]
         '[clojure.string :as str]
         '[fpnum.ad.ops :as ad]
         '[fpnum.ad.reverse :as rev])

(def csv-path "docs/data/ad-reverse-sharing-1b.csv")

(defn medir
  "Ejecuta (thunk) y devuelve {:resultado :ok|:stack-overflow|:error :segundos :valor}."
  [thunk]
  (let [t0 (System/nanoTime)]
    (try
      (let [v (thunk)]
        {:resultado :ok :segundos (/ (- (System/nanoTime) t0) 1e9) :valor v})
      (catch StackOverflowError _
        {:resultado :stack-overflow :segundos (/ (- (System/nanoTime) t0) 1e9)})
      (catch Throwable e
        {:resultado :error :segundos (/ (- (System/nanoTime) t0) 1e9) :valor (.getName (class e))}))))

;; --- L1: cadena lineal -------------------------------------------------------
(defn gradiente-cadena-lineal [n]
  (first (rev/grad (fn [[x]] (reduce (fn [acc _] (ad/d+ (ad/d* acc 0.999) x)) x (range n))) [1.0])))

(defn esperado-cadena-lineal [n]
  ;; d/dx de f_n con f_0 = x, f_{k+1} = 0.999 f_k + x:  Σ_{k=0}^{n} 0.999^k
  (/ (- 1.0 (Math/pow 0.999 (inc n))) (- 1.0 0.999)))

;; --- S1 / S2: cadenas de cuadrados ------------------------------------------
(defn cadena-de-cuadrados [x d]
  (reduce (fn [a _] (ad/d* a a)) x (range d)))

(defn gradiente-s1 [d]
  (first (rev/grad (fn [[x]] (cadena-de-cuadrados x d)) [1.0])))

(defn gradiente-s2 [d]
  (first (rev/grad (fn [[x]] (ad/d+ (cadena-de-cuadrados x d) (cadena-de-cuadrados x d))) [1.0])))

(def filas (atom []))

(defn registrar! [caso n r esperado-ok?]
  (swap! filas conj [caso n (name (:resultado r)) (format "%.4f" (:segundos r))
                     (if (= :ok (:resultado r)) (str esperado-ok?) "")])
  (println (format "%-3s n/d = %6d | %-14s | %9.4f s | gradiente correcto: %s"
                   caso n (name (:resultado r)) (:segundos r)
                   (if (= :ok (:resultado r)) esperado-ok? "-"))))

(println "== Identidad estructural de ids: cadenas largas y subexpresiones compartidas ==")

(println "\n-- L1: cadena lineal, construida una vez --")
(doseq [n [500 1000 2000 3000 5000 10000 20000 50000 100000]]
  (let [r (medir #(gradiente-cadena-lineal n))]
    (registrar! "L1" n r (and (= :ok (:resultado r))
                              (< (Math/abs (- (:valor r) (esperado-cadena-lineal n)))
                                 (* 1e-9 (esperado-cadena-lineal n)))))))

(println "\n-- S1: cadena de cuadrados, construida una vez (x = 1: el gradiente es 2^d) --")
(doseq [d [10 50 100 200 300 500 1000]]
  (let [r (medir #(gradiente-s1 d))]
    (registrar! "S1" d r (and (= :ok (:resultado r))
                              (< (Math/abs (- (:valor r) (Math/pow 2.0 d))) (* 1e-9 (Math/pow 2.0 d)))))))

(println "\n-- S2: la misma cadena construida dos veces y sumada (gradiente 2·2^d) --")
(loop [[d & resto] [8 10 12 14 16 18 20 22 24 26 28 30 100 300 1000]
       ultimo 0.0]
  (when d
    (if (> (* 2 ultimo) 20.0)
      (println (format "(se omiten d >= %d: el paso anterior ya tardó %.1f s y cada +2 de profundidad lo cuadruplica)" d ultimo))
      (let [r (medir #(gradiente-s2 d))]
        (registrar! "S2" d r (and (= :ok (:resultado r))
                                  (< (Math/abs (- (:valor r) (* 2.0 (Math/pow 2.0 d)))) (* 1e-9 (* 2.0 (Math/pow 2.0 d))))))
        (recur resto (:segundos r))))))

(io/make-parents csv-path)
(spit csv-path (str/join "\n" (cons "case,n_or_depth,result,seconds,gradient_correct"
                                    (map #(str/join "," %) @filas))))
(println "\nCSV escrito en" csv-path)
(shutdown-agents)
