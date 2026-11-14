(ns fpnum.demo.reverse-test
  "Pruebas de fpnum.demo.reverse: los casos del experimento, la medición de
  tiempos y la salida CSV. Los experimentos usan tamaños pequeños para ser
  rápidos; no se afirma nada sobre qué modo es más rápido, porque eso depende
  del tamaño y de la máquina."
  (:require [clojure.java.io :as io]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]
            [fpnum.ad.forward :as fwd]
            [fpnum.ad.reverse :as rev]
            [fpnum.demo.reverse :as demo]
            [fpnum.optim :as opt]))

(def ^:private rapido {:runs 2 :warmup 1 :min-sample-ms 1.0})

(defmacro ^:private con-temp-dir [[sym] & cuerpo]
  `(let [~sym (doto (java.io.File/createTempFile "demo-reverse-test" "") (.delete) (.mkdirs))]
     (try ~@cuerpo
          (finally (doseq [^java.io.File f# (reverse (file-seq ~sym))] (.delete f#))))))

(deftest casos-del-experimento
  (testing "mínimos cuadrados: dimensiones y gradiente igual en los dos modos"
    (let [{:keys [n m point f] :as caso} (demo/least-squares-case 20 5)]
      (is (= [5 20] [n m]))
      (is (= 5 (count point)))
      (is (= "minimos-cuadrados" (:name caso)))
      (is (< (demo/max-rel-diff (fwd/grad f point) (rev/grad f point)) 1e-9))))
  (testing "suma de cuadrados: el gradiente es 2w en los dos modos"
    (let [{:keys [point f n]} (demo/sum-of-squares-case 6)
          esperado (mapv #(* 2.0 %) point)]
      (is (= 6 n))
      (is (< (demo/max-rel-diff esperado (fwd/grad f point)) 1e-12))
      (is (< (demo/max-rel-diff esperado (rev/grad f point)) 1e-12))))
  (testing "el problema es reproducible"
    (let [a (demo/least-squares-case 10 3)
          b (demo/least-squares-case 10 3)]
      (is (= (fwd/grad (:f a) (:point a)) (fwd/grad (:f b) (:point b)))))))

(deftest mediana
  (is (= 3 (demo/median [5 1 3])))
  (is (= 2.5 (demo/median [4 1 2 3])))
  (is (= 7 (demo/median [7]))))

(deftest diferencia-relativa-maxima
  (is (= 0.0 (demo/max-rel-diff [1.0 2.0] [1.0 2.0])))
  (is (< (Math/abs (- 0.5 (demo/max-rel-diff [1.0 5.0] [2.0 5.0]))) 1e-12)
      "|1 − 2|/(1 + 1) = 0.5")
  (is (= 0.0 (demo/max-rel-diff [] []))))

(deftest tiempo-por-llamada
  (let [r (demo/time-per-call #(reduce + (range 1000)) rapido)]
    (is (every? #(contains? r %) [:median :min :max :reps]))
    (is (pos? (:median r)))
    (is (<= (:min r) (:median r) (:max r)))
    (is (>= (:reps r) 1))))

(deftest una-fila-del-experimento
  (let [r (demo/scaling-row (demo/least-squares-case 20 5) rapido)]
    (is (= "minimos-cuadrados" (:function r)))
    (is (= [5 20] [(:n r) (:m r)]))
    (is (pos? (:ratio r)))
    (is (< (:max-rel-diff r) 1e-9) "forward y reverse dan el mismo gradiente")
    (testing "la razón es la mediana de las razones por ronda, dentro de su rango"
      (let [[rmin rmax] (:ratio-range r)]
        (is (<= rmin (:ratio r) rmax))
        (is (<= (/ (get-in r [:forward :min]) (get-in r [:reverse :max])) rmin) "cota inferior")
        (is (>= (/ (get-in r [:forward :max]) (get-in r [:reverse :min])) rmax) "cota superior")))))

(deftest tiempo-de-un-par-con-muestras-intercaladas
  (let [orden (atom [])
        a (fn [] (swap! orden conj :a) (reduce + (range 500)))
        b (fn [] (swap! orden conj :b) (reduce + (range 500)))
        r (demo/time-pair a b {:runs 4 :warmup 0 :min-sample-ms 0.0})]
    (is (= 4 (count (:ratios r))))
    (is (every? pos? (:ratios r)))
    (is (every? #(<= (:min %) (:median %) (:max %)) [(:a r) (:b r)]))
    (testing "tras la calibración (a, b), las rondas alternan el orden: A B, B A, A B, B A"
      ;; con min-sample-ms 0 cada muestra es una sola llamada
      (is (= [:a :b] (take 2 @orden)) "calibración")
      (is (= [:a :b :b :a :a :b :b :a] (drop 2 @orden))))))

(deftest lineas-csv
  (let [filas [{:function "minimos-cuadrados" :n 10 :m 100
                :forward {:median 2.0 :min 1.9 :max 2.1} :reverse {:median 4.0 :min 3.9 :max 4.2}
                :ratio 0.5 :ratio-range [0.4 0.6] :max-rel-diff 0.0}
               {:function "suma-de-cuadrados" :n 10 :m nil
                :forward {:median 0.2 :min 0.19 :max 0.21} :reverse {:median 0.1 :min 0.09 :max 0.12}
                :ratio 2.0 :ratio-range [1.8 2.2] :max-rel-diff 1e-16}]
        lineas (demo/scaling-csv-lines filas)]
    (is (= 3 (count lineas)))
    (is (= 13 (count (str/split (first lineas) #","))))
    (is (every? #(= 13 (count (str/split % #"," -1))) lineas))
    (testing "sin m (suma de cuadrados) la celda queda vacía"
      (is (str/starts-with? (nth lineas 2) "suma-de-cuadrados,10,,")))))

(deftest ejecucion-completa-con-tamanos-pequenos
  (con-temp-dir [d]
    (let [csv (str (io/file d "ad.csv"))
          res (volatile! nil)
          salida (with-out-str
                   (vreset! res (demo/run-experiment!
                                 (assoc rapido :ns [3 5] :m 10 :warm-up-ms 10 :csv-path csv))))
          filas @res]
      (testing "una fila por función y tamaño"
        (is (= 4 (count filas)))
        (is (= [["minimos-cuadrados" 3] ["minimos-cuadrados" 5]
                ["suma-de-cuadrados" 3] ["suma-de-cuadrados" 5]]
               (mapv (juxt :function :n) filas))))
      (testing "todos los gradientes coinciden"
        (is (every? #(< (:max-rel-diff %) 1e-9) filas)))
      (testing "imprime la tabla y escribe el CSV"
        (is (str/includes? salida "forward frente a reverse"))
        (is (str/includes? salida "CSV escrito"))
        (is (= 5 (count (str/split-lines (slurp csv)))))))))
