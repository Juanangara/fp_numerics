(ns fpnum.demo.descent-test
  "Pruebas de fpnum.demo.descent: experimento de mínimos cuadrados y salida
  CSV de las trayectorias."
  (:require [clojure.java.io :as io]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]
            [fpnum.demo.descent :as demo]))

(deftest experimento-de-minimos-cuadrados
  (let [r (demo/least-squares-experiment)]
    (is (:converged? r))
    (is (< (:distance-to-normal-equations r) 1e-6))
    (is (pos? (:alpha r)))
    (testing "el resultado es reproducible (semilla fija)"
      (is (= (:x r) (:x (demo/least-squares-experiment)))))
    (testing "w se acerca a los coeficientes verdaderos (ruido σ = 0.1)"
      (is (every? #(< % 0.1)
                  (map (fn [a b] (Math/abs (- a b))) (:x r) (:w-true r)))))))

(def ^:private trayectoria-2d
  [{:iter 0 :x [1.0 2.0] :f 5.0 :grad-norm 4.0 :alpha 0.5}
   {:iter 1 :x [0.5 1.0] :f 1.25 :grad-norm 2.0 :alpha 0.5}
   {:iter 2 :x [0.25 0.5] :f 0.3 :grad-norm 1.0 :alpha 0.5}
   {:iter 3 :x [0.1 0.2] :f 0.05 :grad-norm 0.4 :alpha nil}])

(deftest resumen-de-los-experimentos
  (let [ls (demo/least-squares-experiment)
        rb {:iterations 13756 :reason :converged :f 6.1e-13 :x [0.9999992 0.9999984]}
        filas (demo/summary-rows ls rb)
        por-metrica (into {} (map (fn [[e m v]] [[e m] v]) filas))
        lineas (demo/summary-csv-lines filas)]
    (testing "mínimos cuadrados: el problema, el paso y la distancia a las ecuaciones normales"
      (is (= 100 (get por-metrica ["minimos-cuadrados" "observaciones"])))
      (is (= 42 (get por-metrica ["minimos-cuadrados" "semilla"])))
      (is (= (:alpha ls) (get por-metrica ["minimos-cuadrados" "alpha_seguro"])))
      (is (= (:iterations ls) (get por-metrica ["minimos-cuadrados" "iteraciones"])))
      (is (< (get por-metrica ["minimos-cuadrados" "distancia_a_ecuaciones_normales"]) 1e-6)))
    (testing "tres coeficientes de cada tipo"
      (doseq [pref ["w_descenso_" "w_ecuaciones_normales_" "w_verdadero_"] i [1 2 3]]
        (is (contains? por-metrica ["minimos-cuadrados" (str pref i)]) (str pref i))))
    (testing "Rosenbrock"
      (is (= 13756 (get por-metrica ["rosenbrock" "iteraciones"])))
      (is (= [-1.2 1.0] [(get por-metrica ["rosenbrock" "inicio_1"]) (get por-metrica ["rosenbrock" "inicio_2"])]))
      (is (= 0.9999992 (get por-metrica ["rosenbrock" "x_final_1"]))))
    (testing "formato CSV: cabecera y tres columnas"
      (is (= "experiment,metric,value" (first lineas)))
      (is (= (inc (count filas)) (count lineas)))
      (is (every? #(= 3 (count (str/split % #"," -1))) lineas)))))

(deftest formato-csv
  (let [lineas (demo/trajectory->csv-lines trayectoria-2d)]
    (is (= "iter,x1,x2,f,grad_norm,alpha" (first lineas)))
    (is (= 5 (count lineas)))
    (testing "el α de la última fila queda vacío"
      (is (str/ends-with? (last lineas) ",")))
    (testing "todas las filas tienen 6 columnas"
      (is (every? #(= 6 (count (str/split % #"," -1))) lineas)))))

(deftest csv-con-stride
  (let [lineas (demo/trajectory->csv-lines trayectoria-2d 2)]
    (testing "conserva iter 0 y 2 (múltiplos de 2) y siempre la última (3)"
      (is (= ["0" "2" "3"] (map #(first (str/split % #",")) (rest lineas)))))))

(deftest escritura-csv
  (let [tmp (java.io.File/createTempFile "descent" ".csv")]
    (try
      (demo/write-trajectory-csv! trayectoria-2d (.getPath tmp))
      (is (= 5 (count (str/split-lines (slurp tmp)))))
      (finally (io/delete-file tmp true)))))
