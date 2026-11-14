(ns fpnum.demo.svd-test
  "Pruebas de fpnum.demo.svd: el benchmark devuelve filas bien formadas para
  las dos implementaciones y el CSV se escribe correctamente (con tamaños
  pequeños, para que sea rápido)."
  (:require [clojure.java.io :as io]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]
            [fpnum.demo.svd :as demo]
            [fpnum.svd.core :as svd]
            [fpnum.svd.image :as image]))

(deftest matriz-aleatoria-reproducible
  (is (= (demo/random-matrix 5 7) (demo/random-matrix 5 7)))
  (is (not= (demo/random-matrix 5 7) (demo/random-matrix 5 8)))
  (is (= 5 (count (demo/random-matrix 5 7))))
  (is (every? #(= 5 (count %)) (demo/random-matrix 5 7))))

(deftest fila-de-benchmark
  (doseq [[nombre f] [["svd" svd/svd] ["svd-reference" svd/svd-reference]]]
    (testing nombre
      (let [r (demo/benchmark-row f 12)]
        (is (= 12 (:n r)))
        (is (pos? (:seconds r)))
        (is (:converged? r))
        (is (< (:recon-error r) 1e-13))
        (is (< (:orth-u r) 1e-12))
        (is (< (:orth-v r) 1e-12))))))

(defmacro ^:private con-temp-dir [[sym] & cuerpo]
  `(let [~sym (doto (java.io.File/createTempFile "demo-svd-test" "") (.delete) (.mkdirs))]
     (try ~@cuerpo
          (finally (doseq [^java.io.File f# (reverse (file-seq ~sym))] (.delete f#))))))

(deftest compresion-con-una-imagen-real
  (con-temp-dir [d]
    (let [png (io/file d "entrada.png")]
      (image/save-png! (image/channels->image (image/synthetic-channels 80 60 5)) png)
      (let [res (volatile! nil)
            salida (with-out-str
                     (vreset! res (demo/run-compression! (str png)
                                                         {:max-side 32
                                                          :html-dir (str (io/file d "demo"))
                                                          :data-dir (str (io/file d "data"))})))
            r @res]
        (testing "la imagen se reduce a 32 px de lado mayor y se comprimen sus 3 canales"
          (is (= [24 32] (:dims r)))
          (is (= 24 (:rank r)))
          (is (= 24 (:k (last (:steps r))))))
        (testing "imprime la tabla y las rutas"
          (is (str/includes? salida "PSNR"))
          (is (str/includes? salida "index.html"))
          (is (not (str/includes? salida "AVISO"))))
        (testing "escribe la página, el original, un PNG por k y los CSV"
          (is (.isFile (io/file d "demo" "index.html")))
          (is (.isFile (io/file d "demo" "img" "original.png")))
          (is (= (count (:steps r))
                 (count (filter #(re-matches #"k-\d+\.png" (.getName ^java.io.File %))
                                (.listFiles (io/file d "demo" "img"))))))
          (is (.isFile (io/file d "data" "svd-image-metrics.csv")))
          (is (.isFile (io/file d "data" "svd-image-singular-values.csv"))))))))

(deftest compresion-sin-imagen-usa-la-sintetica
  (con-temp-dir [d]
    (let [res (volatile! nil)
          salida (with-out-str
                   (vreset! res (demo/run-compression! (str (io/file d "no-existe.jpg"))
                                                       {:max-side 24
                                                        :html-dir (str (io/file d "demo"))
                                                        :data-dir (str (io/file d "data"))})))]
      (is (str/includes? salida "AVISO"))
      (is (= [16 24] (:dims @res)))
      (is (.isFile (io/file d "demo" "index.html"))))))

(deftest escritura-csv
  (let [tmp (java.io.File/createTempFile "svd-bench" ".csv")]
    (try
      (demo/write-csv! [(demo/benchmark-row svd/svd 6) (demo/benchmark-row svd/svd 8)]
                       (.getPath tmp))
      (let [lineas (str/split-lines (slurp tmp))]
        (testing "cabecera y una fila por tamaño"
          (is (= "n,seconds,sweeps,rel_recon_error,orth_u,orth_v" (first lineas)))
          (is (= 3 (count lineas)))
          (is (= ["6" "8"] (map #(first (str/split % #",")) (rest lineas))))
          (is (every? #(= 6 (count (str/split % #","))) lineas))))
      (finally (io/delete-file tmp true)))))
