(ns fpnum.demo.ad-test
  "Pruebas de fpnum.demo.ad: la tabla de error frente a h tiene forma de U y
  el CSV se escribe correctamente."
  (:require [clojure.java.io :as io]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]
            [fpnum.demo.ad :as demo]))

(deftest tabla-de-error
  (let [rows (demo/error-table)
        mejor (demo/best-row rows)]
    (testing "tiene 31 filas con h decreciente"
      (is (= 31 (count rows)))
      (is (apply > (map :h rows))))
    (testing "curva en U: el mínimo está cerca de 1e-8"
      (is (< 1e-9 (:h mejor) 1e-7)))
    (testing "el error es mayor en ambos extremos que en el mínimo"
      (is (> (:fd-error (first rows)) (* 100 (:fd-error mejor))))
      (is (> (:fd-error (last rows)) (* 100 (:fd-error mejor)))))
    (testing "la AD es mucho más precisa que las diferencias finitas"
      (is (< (:ad-error (first rows)) 1e-15))
      (is (< (:ad-error (first rows)) (:fd-error mejor))))))

(deftest escritura-del-csv
  (let [tmp (java.io.File/createTempFile "error-vs-h" ".csv")
        path (.getPath tmp)]
    (try
      (demo/write-csv! (demo/error-table) path)
      (let [lineas (str/split-lines (slurp path))]
        (is (= "h,fd_error,ad_error" (first lineas)))
        (is (= 32 (count lineas)))
        (is (every? #(= 3 (count (str/split % #","))) lineas)))
      (finally (io/delete-file tmp true)))))
