(ns fpnum.ad.reverse-sharing-test
  "Subexpresiones muy compartidas, cadenas largas y estado global (contador de
  ids) en reverse mode.

  Estos casos fallaban con la identidad estructural de ids (enfoque 1a):
  desbordaban la pila con sumas de ~4000 términos o más, y construir dos veces
  la misma expresión muy compartida era exponencial. Con ids de contador (1b)
  pasan todos. La historia completa está en docs/decisiones.md y la medición,
  con tiempos, en scripts/measure_reverse_sharing.clj
  (docs/data/ad-reverse-sharing-1a.csv y ad-reverse-sharing-1b.csv)."
  (:require [clojure.test :refer [deftest is testing]]
            [fpnum.ad.ops :as ad]
            [fpnum.ad.reverse :as rev])
  (:import (fpnum.ad.reverse Node)))

(defn con-limite
  "Ejecuta `thunk` en un hilo daemon y espera `ms` milisegundos. Devuelve
  {:ok valor}, {:stack-overflow true}, {:error excepción} o {:timeout true}.
  Un hilo daemon no impide que la JVM termine si el cálculo no acaba (por
  ejemplo, si fuera exponencial)."
  [ms thunk]
  (let [res (promise)]
    (doto (Thread. ^Runnable (fn []
                               (deliver res (try {:ok (thunk)}
                                                 (catch StackOverflowError _ {:stack-overflow true})
                                                 (catch Throwable e {:error e})))))
      (.setDaemon true)
      (.start))
    (deref res ms {:timeout true})))

(defn cadena-de-cuadrados [x d]
  (reduce (fn [a _] (ad/d* a a)) x (range d)))

(defn cadena-lineal [x n]
  (reduce (fn [acc _] (ad/d+ (ad/d* acc 0.999) x)) x (range n)))

(defn esperado-cadena-lineal
  "Σ_{k=0}^{n} 0.999^k: derivada de f_n con f_0 = x, f_{k+1} = 0.999·f_k + x."
  [n]
  (/ (- 1.0 (Math/pow 0.999 (inc n))) (- 1.0 0.999)))

(defn cerca? [a b tol] (<= (Math/abs (- (double a) (double b))) (* tol (Math/abs (double b)))))

;; ---------------------------------------------------------------------------
;; Compartición dentro de un mismo grafo
;; ---------------------------------------------------------------------------

(deftest cadena-de-cuadrados-compartida-es-lineal
  ;; ((x²)²)²… de profundidad d tiene 2^d caminos como árbol. Con x = 1 el
  ;; gradiente es 2^d, que sigue siendo representable hasta d = 1023.
  (doseq [d [10 50 100 300 1000]]
    (testing (str "profundidad " d)
      (let [r (con-limite 5000 #(first (rev/grad (fn [[x]] (cadena-de-cuadrados x d)) [1.0])))]
        (is (contains? r :ok) (str "resultado: " (keys r)))
        (is (cerca? (:ok r) (Math/pow 2.0 d) 1e-9))))))

(deftest la-misma-expresion-construida-dos-veces-es-lineal
  ;; Con 1a esto era EXPONENCIAL (d = 28 tardaba 20 s; d = 30 más de un minuto).
  ;; Con ids de contador los dos nodos finales son distintos y no hay
  ;; comparaciones profundas: milisegundos.
  (doseq [d [4 8 12 30 100 300 1000]]
    (testing (str "profundidad " d)
      (let [r (con-limite 5000
                          #(first (rev/grad (fn [[x]] (ad/d+ (cadena-de-cuadrados x d)
                                                              (cadena-de-cuadrados x d)))
                                            [1.0])))]
        (is (not (:timeout r)) "más de 5 s: la construcción duplicada no debe ser exponencial")
        (is (contains? r :ok) (str "resultado: " (keys r)))
        (is (cerca? (:ok r) (* 2.0 (Math/pow 2.0 d)) 1e-9))))))

(deftest tres-copias-de-una-expresion-compartida
  (let [r (con-limite 5000
                      #(first (rev/grad (fn [[x]] (ad/d+ (cadena-de-cuadrados x 200)
                                                          (cadena-de-cuadrados x 200)
                                                          (cadena-de-cuadrados x 200)))
                                        [1.0])))]
    (is (cerca? (:ok r) (* 3.0 (Math/pow 2.0 200)) 1e-9))))

;; ---------------------------------------------------------------------------
;; Cadenas y sumas largas (antes: StackOverflowError desde ~4000 niveles)
;; ---------------------------------------------------------------------------

(deftest cadena-lineal-larga
  (doseq [n [500 1000 5000 20000 50000]]
    (testing (str n " términos")
      (let [r (con-limite 30000 #(first (rev/grad (fn [[x]] (cadena-lineal x n)) [1.0])))]
        (is (not (:stack-overflow r)) "StackOverflowError")
        (is (contains? r :ok) (str "resultado: " (keys r)))
        (is (cerca? (:ok r) (esperado-cadena-lineal n) 1e-9))))))

(deftest una-variable-usada-muchisimas-veces
  ;; (apply d+ ...) pliega a la izquierda: n usos son n sumas anidadas.
  (doseq [n [1000 4000 5000 20000]]
    (testing (str n " usos")
      (is (= {:ok (double n)}
             (con-limite 30000 #(first (rev/grad (fn [[x]] (apply ad/d+ (repeat n x))) [3.0]))))))))

(deftest suma-de-veinte-mil-productos
  ;; Es lo que hace un objetivo de mínimos cuadrados con 20 000 observaciones.
  (let [r (con-limite 30000
                      #(first (rev/grad (fn [[x]] (apply ad/d+ (map (fn [i] (ad/d* x (double i))) (range 1 20001))))
                                        [1.0])))]
    (is (contains? r :ok) (str "resultado: " (keys r)))
    (is (cerca? (:ok r) (/ (* 20000.0 20001.0) 2.0) 1e-9))))

(deftest objetivo-de-minimos-cuadrados-con-muchas-observaciones
  ;; 5000 filas × 3 columnas: antes desbordaba la pila.
  (let [rnd (java.util.Random. 5)
        filas (vec (repeatedly 5000 #(vec (repeatedly 3 (fn [] (.nextGaussian rnd))))))
        y (vec (repeatedly 5000 #(.nextGaussian rnd)))
        w [0.3 -0.2 0.5]
        objetivo (fn [[a b c]]
                   (apply ad/d+ (map (fn [[x1 x2 x3] yi]
                                       (let [r (ad/d- (ad/d+ (ad/d* a x1) (ad/d* b x2) (ad/d* c x3)) yi)]
                                         (ad/d* r r)))
                                     filas y)))
        analitico (let [res (mapv (fn [[x1 x2 x3] yi] (- (+ (* 0.3 x1) (* -0.2 x2) (* 0.5 x3)) yi)) filas y)]
                    (mapv (fn [j] (* 2.0 (reduce + (map (fn [fila r] (* (fila j) r)) filas res)))) (range 3)))
        r (con-limite 30000 #(rev/grad objetivo w))]
    (is (contains? r :ok) (str "resultado: " (keys r)))
    (is (every? true? (map #(cerca? %1 %2 1e-9) (:ok r) analitico)))))

;; ---------------------------------------------------------------------------
;; Estado global (el contador de ids)
;; ---------------------------------------------------------------------------

(deftest llamadas-concurrentes-no-se-interfieren
  ;; El contador de ids es global. Cada nodo toma un id único aunque varios
  ;; hilos construyan grafos a la vez, y el gradiente no depende de los ids.
  (let [resultados (doall (pmap (fn [k]
                                  (first (rev/grad (fn [[x]] (apply ad/d+ (repeat k x))) [1.0])))
                                (range 1 301)))]
    (is (= (mapv double (range 1 301)) (vec resultados))))
  (let [f (fn [[x y]] (ad/d+ (ad/d* x x y) (ad/dsin y)))
        resultados (doall (pmap (fn [_] (rev/grad f [1.0 2.0])) (range 200)))]
    (is (apply = resultados) "el mismo gradiente en todos los hilos")))

(deftest el-resultado-no-depende-de-los-ids
  ;; Los ids solo distinguen nodos: llamar muchas veces (ids cada vez mayores)
  ;; da siempre exactamente el mismo gradiente.
  (let [f (fn [[x y]] (ad/d+ (ad/d* x x y) (ad/dsin y) (ad/dexp (ad/d* x y))))
        primero (rev/grad f [1.0 2.0])]
    (dotimes [_ 50]
      (is (= primero (rev/grad f [1.0 2.0]))))))

(defn- ids-crecientes?
  "Recorre el grafo desde `n` y comprueba que el id de cada nodo es mayor que el
  de todos sus argumentos (los ids crecen en el orden de creación y un nodo se
  crea después que sus argumentos)."
  [^Node n]
  (every? (fn [^Node a] (and (> (.-id n) (.-id a)) (ids-crecientes? a)))
          (.-args n)))

(deftest el-id-de-un-nodo-es-mayor-que-el-de-sus-argumentos
  (let [capturado (volatile! nil)]
    (rev/grad (fn [[x y]]
                (let [s (ad/d+ (ad/d* x y) (ad/dsin x) 2.0)
                      r (ad/d-div (ad/dexp s) (ad/d+ 1 (ad/d* s s)))]
                  (vreset! capturado r)
                  r))
              [1.0 2.0])
    (is (ids-crecientes? @capturado))
    (testing "un mismo id nunca designa a dos nodos distintos"
      ;; Un nodo compartido aparece varias veces al recorrer el grafo, pero
      ;; siempre es el mismo objeto; dos objetos distintos con el mismo id
      ;; serían un error de identidad.
      (let [vistos (volatile! {})
            conflictos (volatile! 0)
            recorrer (fn recorrer [^Node n]
                       (let [previo (get @vistos (.-id n))]
                         (cond (nil? previo) (do (vswap! vistos assoc (.-id n) n)
                                                 (doseq [a (.-args n)] (recorrer a)))
                               (not (identical? previo n)) (vswap! conflictos inc))))]
        (recorrer @capturado)
        (is (= 0 @conflictos))
        (is (> (count @vistos) 5) "se recorrió un grafo con varios nodos distintos")))))
