(ns fpnum.svd.image-test
  "Pruebas de fpnum.svd.image.

  * Teoría: Eckart–Young–Mirsky (‖A − A_k‖₂ = σ_{k+1} y ‖A − A_k‖_F² = Σ_{i>k} σᵢ²).
  * La reconstrucción incremental coincide con la directa.
  * Métricas con valores conocidos (PSNR, energía, razón de compresión).
  * Entrada y salida con imágenes sintéticas en archivos temporales (no hay
    imágenes binarias en test/)."
  (:require [clojure.java.io :as io]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]
            [fpnum.linear :as la]
            [fpnum.svd.core :as svd]
            [fpnum.svd.image :as image])
  (:import (java.awt.image BufferedImage)
           (java.io File)))

;; ---------------------------------------------------------------------------
;; Utilidades
;; ---------------------------------------------------------------------------

(defn- approx= [a b tol]
  (<= (Math/abs (- (double a) (double b)))
      (* tol (+ 1.0 (Math/abs (double a)) (Math/abs (double b))))))

(defn- max-dif [a b]
  (apply max 0.0 (map #(Math/abs (- (double %1) (double %2))) (flatten a) (flatten b))))

(defn- temp-dir ^File []
  (let [f (File/createTempFile "svd-img-test" "")]
    (.delete f)
    (.mkdirs f)
    f))

(defn- borrar-recursivo [^File d]
  (doseq [^File f (reverse (file-seq d))] (.delete f)))

(defmacro ^:private con-temp-dir [[sym] & cuerpo]
  `(let [~sym (temp-dir)]
     (try ~@cuerpo
          (finally (borrar-recursivo ~sym)))))

(def ^:private canales-24x16 (delay (image/synthetic-channels 24 16 7)))

;; ---------------------------------------------------------------------------
;; Imagen sintética
;; ---------------------------------------------------------------------------

(deftest imagen-sintetica
  (let [c (image/synthetic-channels 30 20 1)]
    (testing "dimensiones: h filas por w columnas, en los tres canales"
      (doseq [k [:r :g :b]]
        (is (= [20 30] (la/dims (c k))))))
    (testing "valores enteros en [0, 255]"
      (is (every? #(and (<= 0.0 % 255.0) (== % (Math/rint %)))
                  (flatten (vals c)))))
    (testing "reproducible con la misma semilla, distinta con otra"
      (is (= c (image/synthetic-channels 30 20 1)))
      (is (not= c (image/synthetic-channels 30 20 2))))))

;; ---------------------------------------------------------------------------
;; Teoría: Eckart–Young–Mirsky
;; ---------------------------------------------------------------------------

(deftest eckart-young-mirsky
  (doseq [canal [:r :g :b]]
    (let [a (@canales-24x16 canal)
          r (svd/svd a)
          s (:s r)]
      (doseq [k [1 3 5 10]]
        (testing (str "canal " (name canal) ", k = " k)
          (let [ak (svd/reconstruct (svd/truncate r k))
                dif (mapv la/v- a ak)]
            (testing "norma de Frobenius: ‖A − A_k‖_F² = Σ_{i>k} σᵢ²"
              (let [f (la/norm-fro dif)]
                (is (approx= (* f f) (reduce + (map #(* % %) (drop k s))) 1e-9))))
            (testing "norma 2: ‖A − A_k‖₂ = σ_{k+1}"
              (is (< (Math/abs (- (first (:s (svd/svd dif))) (nth s k)))
                     (* 1e-9 (first s)))))))))))

;; ---------------------------------------------------------------------------
;; compress
;; ---------------------------------------------------------------------------

(def ^:private ks [1 2 4 8 16])

(def ^:private resultado (delay (image/compress @canales-24x16 ks)))

(deftest estructura-de-compress
  (let [r @resultado]
    (is (= [16 24] (:dims r)))
    (is (= 16 (:rank r)))
    (is (= ks (mapv :k (:steps r))))
    (is (= [:r :g :b] (vec (sort-by {:r 0 :g 1 :b 2} (keys (:sigmas r))))))
    (is (every? #(= 16 (count %)) (vals (:sigmas r))))
    (testing "valores singulares no negativos y decrecientes"
      (is (every? (fn [s] (and (every? #(>= % 0.0) s) (apply >= s))) (vals (:sigmas r)))))
    (testing "cada reconstrucción tiene el tamaño de la imagen"
      (doseq [{:keys [channels]} (:steps r) c [:r :g :b]]
        (is (= [16 24] (la/dims (channels c))))))))

(deftest reconstruccion-incremental-coincide-con-la-directa
  (let [r @resultado]
    (doseq [{:keys [k channels]} (:steps r) c [:r :g :b]]
      (testing (str "k = " k ", canal " (name c))
        (let [directa (svd/reconstruct (svd/truncate (svd/svd (@canales-24x16 c)) k))]
          (is (< (max-dif (channels c) directa) 1e-9)))))))

(deftest metricas-monotonas-y-coherentes-con-la-teoria
  (let [{:keys [steps sigmas]} @resultado
        total (reduce + (for [c [:r :g :b] x (sigmas c)] (* x x)))]
    (testing "la energía sube y el error baja al crecer k"
      (is (apply < (map #(get-in % [:metrics :energy]) steps)))
      (is (apply > (map #(get-in % [:metrics :rel-error]) steps)))
      (is (apply <= (map #(get-in % [:metrics :psnr]) steps))))
    (testing "error relativo = √(Σ_{i>k} σᵢ² / Σ σᵢ²), sumando los tres canales"
      (doseq [{:keys [k metrics]} steps]
        (let [cola (reduce + (for [c [:r :g :b] x (drop k (sigmas c))] (* x x)))]
          (is (approx= (:rel-error metrics) (Math/sqrt (/ cola total)) 1e-8)))))
    (testing "energía + (energía no retenida) = 1"
      (doseq [{:keys [k metrics]} steps]
        (let [cola (reduce + (for [c [:r :g :b] x (drop k (sigmas c))] (* x x)))]
          (is (approx= (:energy metrics) (- 1.0 (/ cola total)) 1e-9)))))
    (testing "razón de compresión m·n / (k(m+n+1))"
      (doseq [{:keys [k metrics]} steps]
        (is (approx= (:ratio metrics) (/ (* 16.0 24.0) (* k (+ 16 24 1))) 1e-12))))
    (testing "con k = rango la reconstrucción es exacta"
      (let [{:keys [metrics]} (last steps)]
        (is (approx= 1.0 (:energy metrics) 1e-12))
        (is (< (:rel-error metrics) 1e-12))
        (is (> (:psnr metrics) 100.0) "tras redondear a 8 bits la imagen es idéntica (PSNR infinito)")))))

(deftest ks-invalidos
  (let [c @canales-24x16]
    (is (thrown? clojure.lang.ExceptionInfo (image/compress c [])))
    (is (thrown? clojure.lang.ExceptionInfo (image/compress c [3 2])))
    (is (thrown? clojure.lang.ExceptionInfo (image/compress c [2 2])))
    (is (thrown? clojure.lang.ExceptionInfo (image/compress c [0 1])))
    (is (thrown? clojure.lang.ExceptionInfo (image/compress c [1 17])) "el rango máximo es min(16, 24) = 16")))

(deftest default-ks-incluye-el-rango
  (is (= [1 2 3 5 8 10 15 20 30 40 50 75 100 128] (image/default-ks 128)))
  (is (= [1 2 3 5 7] (image/default-ks 7)))
  (is (= [1] (image/default-ks 1))))

;; ---------------------------------------------------------------------------
;; Métricas con valores conocidos
;; ---------------------------------------------------------------------------

(defn- gris [v] {:r [[v v v] [v v v]] :g [[v v v] [v v v]] :b [[v v v] [v v v]]})

(deftest psnr-valores-conocidos
  (testing "diferencia constante de 1: MSE = 1, PSNR = 10·log₁₀(255²) ≈ 48.1308 dB"
    (is (approx= 48.13080361 (image/psnr (gris 100.0) (gris 101.0)) 1e-9)))
  (testing "diferencia constante de 10: MSE = 100"
    (is (approx= (* 10.0 (Math/log10 (/ 65025.0 100.0))) (image/psnr (gris 100.0) (gris 110.0)) 1e-12)))
  (testing "imágenes idénticas: MSE = 0 y PSNR = +Infinity"
    (is (= Double/POSITIVE_INFINITY (image/psnr (gris 100.0) (gris 100.0))))))

(deftest mse-promedia-sobre-todos-los-canales
  (let [a (gris 100.0)
        b (assoc (gris 100.0) :r [[102.0 102.0 102.0] [102.0 102.0 102.0]])]
    (is (approx= (/ 4.0 3.0) (image/mse a b) 1e-12) "solo un canal de tres difiere en 2 → (4+0+0)/3")))

(deftest error-relativo
  (is (approx= 1.0 (image/relative-error [[3.0 4.0]] [[0.0 0.0]]) 1e-12))
  (is (approx= 0.0 (image/relative-error [[3.0 4.0]] [[3.0 4.0]]) 1e-12))
  (is (approx= 5.0 (image/relative-error [[0.0 0.0]] [[3.0 4.0]]) 1e-12) "A = 0: norma de la diferencia"))

(deftest energia-retenida
  (is (approx= 0.36 (image/energy-retained [3.0 4.0] 1) 1e-12))
  (is (approx= 1.0 (image/energy-retained [3.0 4.0] 2) 1e-12))
  (is (approx= 1.0 (image/energy-retained [0.0 0.0] 1) 1e-12) "energía total cero"))

(deftest razon-de-compresion
  (is (approx= (/ 10000.0 2010.0) (image/compression-ratio 100 100 10) 1e-12))
  (testing "con k grande el almacenamiento supera al original (razón < 1)"
    (is (< (image/compression-ratio 100 100 60) 1.0)))
  (testing "k(m+n+1) = m·n da razón 1"
    (is (approx= 1.0 (image/compression-ratio 10 10 (/ 100.0 21.0)) 1e-12))))

(deftest cuantizar
  (is (= [[0.0 0.0 1.0 255.0 255.0]]
         (image/quantize [[-3.2 0.4 0.5 254.6 300.0]]))))

;; ---------------------------------------------------------------------------
;; Imágenes: canales, PNG, redimensionado
;; ---------------------------------------------------------------------------

(deftest canales-ida-y-vuelta
  (let [c @canales-24x16
        img (image/channels->image c)]
    (is (instance? BufferedImage img))
    (is (= [24 16] [(.getWidth img) (.getHeight img)]))
    (is (= c (image/image->channels img)))))

(deftest channels->image-recorta-y-redondea
  (let [c {:r [[300.0 -5.0 127.5]] :g [[0.4 0.6 255.0]] :b [[10.0 20.0 30.0]]}
        out (image/image->channels (image/channels->image c))]
    (is (= [[255.0 0.0 128.0]] (:r out)))
    (is (= [[0.0 1.0 255.0]] (:g out)))
    (is (= [[10.0 20.0 30.0]] (:b out)))))

(deftest png-sin-perdida
  (con-temp-dir [d]
    (let [c @canales-24x16
          f (io/file d "sub" "x.png")]
      (image/save-png! (image/channels->image c) f)
      (is (.isFile f) "crea las carpetas intermedias")
      (is (= c (image/image->channels (image/load-image f)))))))

(deftest carga-con-errores
  (is (thrown-with-msg? clojure.lang.ExceptionInfo #"No existe"
                        (image/load-image "no/existe/imagen.png")))
  (con-temp-dir [d]
    (let [f (io/file d "texto.png")]
      (spit f "esto no es una imagen")
      (is (thrown-with-msg? clojure.lang.ExceptionInfo #"Formato"
                            (image/load-image f))))))

(deftest redimensionado
  (let [grande (image/channels->image (image/synthetic-channels 405 270 3))]
    (testing "conserva la proporción: 405×270 con lado máximo 192 → 192×128"
      (let [r (image/resize grande 192)]
        (is (= [192 128] [(.getWidth r) (.getHeight r)]))
        (is (= BufferedImage/TYPE_INT_RGB (.getType r)))))
    (testing "una imagen alta se limita por el alto"
      (let [alta (image/channels->image (image/synthetic-channels 60 120 3))
            r (image/resize alta 40)]
        (is (= [20 40] [(.getWidth r) (.getHeight r)]))))
    (testing "nunca agranda"
      (let [chica (image/channels->image @canales-24x16)
            r (image/resize chica 192)]
        (is (= [24 16] [(.getWidth r) (.getHeight r)]))
        (is (= @canales-24x16 (image/image->channels r)) "sin cambios de píxeles")))
    (testing "la reducción conserva el brillo medio (sin aliasing grosero)"
      (let [media (fn [c] (/ (reduce + (flatten (vals c))) (* 3.0 (count (flatten (:r c))))))
            antes (media (image/image->channels grande))
            despues (media (image/image->channels (image/resize grande 48)))]
        (is (< (Math/abs (- antes despues)) 3.0))))))

;; ---------------------------------------------------------------------------
;; Salidas: CSV, JSON, HTML y archivos
;; ---------------------------------------------------------------------------

(deftest csv-de-valores-singulares
  (let [lineas (image/singular-values-csv-lines @resultado)]
    (is (= "i,r,g,b" (first lineas)))
    (is (= 17 (count lineas)) "cabecera + 16 valores singulares")
    (is (every? #(= 4 (count (str/split % #","))) lineas))
    (is (= "1" (first (str/split (second lineas) #","))))))

(deftest csv-de-metricas
  (let [lineas (image/metrics-csv-lines @resultado)]
    (is (= "k,energy,rel_error,psnr_db,compression_ratio" (first lineas)))
    (is (= (inc (count ks)) (count lineas)))
    (is (every? #(= 5 (count (str/split % #","))) lineas))
    (is (= (map str ks) (map #(first (str/split % #",")) (rest lineas))))
    (is (str/includes? (last lineas) "Infinity") "k = rango: PSNR infinito")))

(deftest datos-y-html-del-visor
  (let [datos (image/viewer-data @resultado)
        html (image/render-viewer-html @resultado)]
    (is (= [16 24 16] [(:m datos) (:n datos) (:rank datos)]))
    (is (= (mapv #(str "img/" (image/image-file-name %)) ks) (:images datos)))
    (testing "el HTML incrusta los datos y no deja el marcador"
      (is (not (str/includes? html "__DATA__")))
      (is (str/includes? html "\"relError\":"))
      (is (str/includes? html "\"images\":[\"img/k-001.png\""))
      (is (str/includes? html "\"psnr\":null") "el PSNR infinito se serializa como null"))
    (testing "es un documento HTML completo y sin recursos externos"
      (is (str/starts-with? html "<!doctype html>"))
      (is (not (re-find #"https?://" (str/replace html "http://www.w3.org/2000/svg" "")))))))

(deftest nombre-de-archivo-por-k
  (is (= "k-001.png" (image/image-file-name 1)))
  (is (= "k-128.png" (image/image-file-name 128))))

(deftest escritura-de-resultados
  (con-temp-dir [d]
    (let [html-dir (io/file d "demo")
          data-dir (io/file d "data")
          salidas (image/write-outputs! @resultado {:html-dir html-dir :data-dir data-dir})]
      (is (.isFile (io/file (:html salidas))))
      (is (.isFile (io/file html-dir "img" "original.png")))
      (doseq [k ks]
        (is (.isFile (io/file html-dir "img" (image/image-file-name k))) (str "k = " k)))
      (is (.isFile (io/file data-dir "svd-image-singular-values.csv")))
      (is (.isFile (io/file data-dir "svd-image-metrics.csv")))
      (testing "con k = rango el PNG es idéntico al original (sin pérdida)"
        (is (= (image/image->channels (image/load-image (io/file html-dir "img" "original.png")))
               (image/image->channels (image/load-image (io/file html-dir "img" "k-016.png"))))))
      (testing "el CSV escrito coincide con las líneas generadas"
        (is (= (str/join "\n" (image/metrics-csv-lines @resultado))
               (slurp (io/file data-dir "svd-image-metrics.csv"))))))))
