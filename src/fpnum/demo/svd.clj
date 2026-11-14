(ns fpnum.demo.svd
  "Demo ejecutable de la SVD. Tiene dos modos.

  1. Compresión de imágenes (por defecto). Carga una imagen, la reduce a unos
     192 px de lado mayor, calcula la SVD de cada canal R, G y B y reconstruye
     la imagen con distintos k (valores singulares conservados). Escribe:
       docs/demo-svd/index.html     página con un deslizador sobre k
       docs/demo-svd/img/           PNG original y uno por cada k
       docs/data/svd-image-*.csv    valores singulares y métricas
     Si no encuentra la imagen, usa una sintética y lo avisa.

       clojure -M -m fpnum.demo.svd                      ; usa docs/img/original.jpg
       clojure -M -m fpnum.demo.svd ruta/a/imagen.png

  2. Benchmark (`benchmark`). Mide el costo de la SVD por Jacobi sobre
     matrices aleatorias cuadradas (entradas gaussianas, semilla fija),
     comparando las dos implementaciones de `fpnum.svd.core`:
       * `svd-reference`: inmutable, legible, lenta.
       * `svd`: rápida, con arrays double[] y normas de columna en caché.
     Para cada tamaño informa tiempo, barridos, error de reconstrucción y
     ortonormalidad de U y V, y guarda un CSV por implementación en docs/data/.
     El benchmark vive aquí y no en las pruebas porque medir tiempo es E/S del
     entorno, no una propiedad del algoritmo.

       clojure -M -m fpnum.demo.svd benchmark           ; tamaños 64, 128 y 256
       clojure -M -m fpnum.demo.svd benchmark 32 64     ; tamaños a elección"
  (:require [clojure.java.io :as io]
            [clojure.string :as str]
            [fpnum.linear :as la]
            [fpnum.svd.core :as svd]
            [fpnum.svd.image :as image]))

(def reference-csv-path
  "CSV del benchmark de `svd-reference` (la versión inmutable)."
  "docs/data/svd-benchmark-reference.csv")

(def fast-csv-path
  "CSV del benchmark de `svd` (la versión rápida, con arrays)."
  "docs/data/svd-benchmark-fast.csv")

(def default-sizes
  "Tamaños n de las matrices del benchmark."
  [64 128 256])

(defn random-matrix
  "Matriz n×n con entradas gaussianas estándar, reproducible con `seed`."
  [n seed]
  (let [rnd (java.util.Random. (long seed))]
    (vec (repeatedly n #(vec (repeatedly n (fn [] (.nextGaussian rnd))))))))

(defn- error-ortonormal
  "Mayor |(QᵀQ − I)ᵢⱼ|."
  [q]
  (let [g (la/mat-mat (la/transpose q) q)
        n (count g)]
    (apply max (for [i (range n) j (range n)]
                 (Math/abs (- (get-in g [i j]) (if (= i j) 1.0 0.0)))))))

(defn benchmark-row
  "Calcula la SVD con `svd-fn` (por ejemplo `svd/svd` o `svd/svd-reference`)
  de una matriz aleatoria n×n y devuelve un mapa con el tiempo en segundos,
  los barridos y los errores."
  [svd-fn n]
  (let [a (random-matrix n n)
        t0 (System/nanoTime)
        r (svd-fn a)
        segundos (/ (- (System/nanoTime) t0) 1e9)]
    {:n n
     :seconds segundos
     :sweeps (:sweeps r)
     :converged? (:converged? r)
     :recon-error (/ (la/norm-fro (mapv la/v- a (svd/reconstruct r)))
                     (la/norm-fro a))
     :orth-u (error-ortonormal (:u r))
     :orth-v (error-ortonormal (:v r))}))

(defn write-csv!
  "Escribe las filas del benchmark en `path` (crea las carpetas necesarias)."
  [rows path]
  (io/make-parents path)
  (spit path
        (str/join "\n"
                  (cons "n,seconds,sweeps,rel_recon_error,orth_u,orth_v"
                        (map (fn [{:keys [n seconds sweeps recon-error orth-u orth-v]}]
                               (str/join "," [n seconds sweeps recon-error orth-u orth-v]))
                             rows)))))

(defn- imprimir-fila [r]
  (println (format "n = %3d | %8.2f s | %2d barridos | error recon. %.1e | ortog. U %.1e V %.1e"
                   (:n r) (:seconds r) (:sweeps r)
                   (:recon-error r) (:orth-u r) (:orth-v r))))

(defn- medir
  "Mide `svd-fn` para cada tamaño (tras calentar la JVM), imprime cada fila y
  guarda el CSV en `path`. Devuelve las filas."
  [titulo svd-fn sizes path]
  (println (str "\n-- " titulo " --"))
  (benchmark-row svd-fn 32)             ; calentamiento
  (let [rows (vec (for [n sizes]
                    (doto (benchmark-row svd-fn n) imprimir-fila)))]
    (write-csv! rows path)
    (println "CSV escrito en" path)
    rows))

(defn run-benchmark!
  "Mide las dos implementaciones para los tamaños dados (o 64, 128 y 256),
  imprime una comparación y guarda un CSV por implementación."
  [sizes]
  (let [sizes (if (seq sizes) (mapv #(Long/parseLong %) sizes) default-sizes)]
    (println "== SVD por Jacobi de un solo lado: costo ==")
    (let [ref-rows (medir "svd-reference (inmutable)" svd/svd-reference sizes reference-csv-path)
          fast-rows (medir "svd (rápida, arrays double[])" svd/svd sizes fast-csv-path)]
      (println "\n-- Aceleración --")
      (doseq [[r f] (map vector ref-rows fast-rows)]
        (println (format "n = %3d | %8.2f s -> %6.2f s | x%.1f"
                         (:n r) (:seconds r) (:seconds f)
                         (/ (:seconds r) (:seconds f))))))))

;; ---------------------------------------------------------------------------
;; Compresión de imágenes
;; ---------------------------------------------------------------------------

(def default-image
  "Imagen que comprime la demo. Si el archivo no existe se usa una imagen
  sintética reproducible."
  "docs/img/original.jpg")

(def default-html-dir
  "Carpeta de la página interactiva y de los PNG de cada k."
  "docs/demo-svd")

(def default-data-dir
  "Carpeta de los CSV de valores singulares y métricas."
  "docs/data")

(def default-max-side
  "Lado mayor, en píxeles, al que se reduce la imagen antes de la SVD."
  192)

(defn- cargar-canales
  "Canales de la imagen en `image-path`, redimensionada a `max-side`. Si el
  archivo no existe usa una imagen sintética y lo avisa."
  [image-path max-side]
  (if (.isFile (io/file image-path))
    (-> (image/load-image image-path)
        (image/resize max-side)
        image/image->channels)
    (do (println (str "AVISO: no existe " image-path "; se usa una imagen sintética."))
        (image/synthetic-channels max-side (quot (* max-side 2) 3) 42))))

(defn run-compression!
  "Comprime la imagen y escribe la página, los PNG y los CSV (ver el docstring
  del namespace). Opciones: :max-side, :html-dir y :data-dir. Devuelve el
  resultado de `fpnum.svd.image/compress` junto con las rutas escritas."
  ([] (run-compression! default-image {}))
  ([image-path] (run-compression! image-path {}))
  ([image-path {:keys [max-side html-dir data-dir]
                :or {max-side default-max-side
                     html-dir default-html-dir
                     data-dir default-data-dir}}]
   (println "== Compresión de imágenes con SVD ==")
   (let [canales (cargar-canales image-path max-side)
         h (count (:r canales))
         w (count (first (:r canales)))
         t0 (System/nanoTime)
         result (image/compress canales (image/default-ks (min h w)))
         segundos (/ (- (System/nanoTime) t0) 1e9)
         salidas (image/write-outputs! result {:html-dir html-dir :data-dir data-dir})]
     (println (format "Imagen de %d×%d; 3 SVD y %d reconstrucciones en %.1f s"
                      w h (count (:steps result)) segundos))
     (println "\n   k | energía  | error rel. | PSNR (dB) | compresión")
     (doseq [{:keys [k metrics]} (:steps result)]
       (println (format "%4d | %6.2f %% | %10.4f | %9.2f | %9.2fx"
                        k (* 100 (:energy metrics)) (:rel-error metrics)
                        (:psnr metrics) (:ratio metrics))))
     (println "\nPágina:" (:html salidas))
     (println "Imágenes:" (:images salidas))
     (println "CSV:" (:singular-values-csv salidas) "y" (:metrics-csv salidas))
     (assoc result :outputs salidas))))

(defn -main
  "Sin argumentos o con la ruta de una imagen: compresión de imágenes. Con
  `benchmark [tamaños...]`: mide el costo de la SVD."
  [& args]
  (cond
    (= "benchmark" (first args)) (run-benchmark! (rest args))
    (seq args) (run-compression! (first args))
    :else (run-compression!)))
