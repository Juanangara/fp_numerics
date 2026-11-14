(ns fpnum.demo.descent
  "Demo ejecutable del descenso de gradiente con diferenciación automática.

  Dos experimentos:
  1. Mínimos cuadrados sobre datos ruidosos reproducibles (semilla fija), con
     el paso seguro α < 1/‖X‖²_F, comparado con las ecuaciones normales.
  2. Función de Rosenbrock desde (−1.2, 1) con backtracking (Armijo).

  Las trayectorias (iteración, coordenadas, f, ‖∇f‖ y α usado) se guardan en
  CSV dentro de docs/data/ para graficarlas.

  Uso: clojure -M -m fpnum.demo.descent"
  (:require [clojure.java.io :as io]
            [clojure.string :as str]
            [fpnum.linear :as la]
            [fpnum.optim :as opt]))

(def least-squares-csv
  "CSV con la trayectoria del descenso de gradiente en mínimos cuadrados."
  "docs/data/descent-least-squares.csv")
(def rosenbrock-csv
  "CSV con la trayectoria (una de cada 100 iteraciones) del descenso sobre Rosenbrock."
  "docs/data/descent-rosenbrock.csv")
(def summary-csv
  "CSV con las cifras de resumen de los dos experimentos (formato largo:
  experimento, métrica, valor)."
  "docs/data/descent-resumen.csv")

;; ---------------------------------------------------------------------------
;; Experimentos (puros)
;; ---------------------------------------------------------------------------

(defn least-squares-experiment
  "Descenso de gradiente con paso seguro sobre datos ruidosos (semilla 42).
  Devuelve el resultado del descenso junto con la solución de referencia por
  ecuaciones normales y la solución verdadera."
  []
  (let [{:keys [x y w-true]} (opt/noisy-regression-data
                              {:n 100 :w-true [2.0 -1.0 0.5] :sigma 0.1 :seed 42})
        alpha (opt/safe-step x)
        res (opt/gradient-descent (opt/least-squares-objective x y)
                                  (la/zeros (count w-true))
                                  {:alpha alpha :max-iter 5000 :tol 1e-8})
        w-ne (opt/normal-equations x y)]
    (assoc res
           :problem {:observations 100 :sigma 0.1 :seed 42}
           :alpha alpha
           :w-normal-equations w-ne
           :w-true w-true
           :distance-to-normal-equations (la/norm2 (la/v- (:x res) w-ne)))))

(def ^:private rosenbrock-start [-1.2 1.0])
(def ^:private rosenbrock-tol 1e-6)

(defn rosenbrock-experiment
  "Descenso con backtracking sobre Rosenbrock desde (−1.2, 1)."
  []
  (opt/gradient-descent opt/rosenbrock rosenbrock-start
                        {:backtracking? true :max-iter 100000 :tol rosenbrock-tol}))

;; ---------------------------------------------------------------------------
;; Salida CSV
;; ---------------------------------------------------------------------------

(defn trajectory->csv-lines
  "Líneas CSV de una trayectoria: cabecera `iter,x1,...,xn,f,grad_norm,alpha`
  y una fila por punto. Con `stride` > 1 se conserva una de cada `stride`
  iteraciones (y siempre la última), para trayectorias largas."
  ([trajectory] (trajectory->csv-lines trajectory 1))
  ([trajectory stride]
   (let [n (count (:x (first trajectory)))
         ultima (:iter (peek trajectory))
         cabecera (str/join "," (concat ["iter"]
                                        (map #(str "x" (inc %)) (range n))
                                        ["f" "grad_norm" "alpha"]))
         fila->linea (fn [{:keys [iter x f grad-norm alpha]}]
                       (str/join "," (concat [iter] x [f grad-norm (or alpha "")])))]
     (cons cabecera
           (map fila->linea
                (filter #(or (zero? (mod (:iter %) stride)) (= (:iter %) ultima))
                        trajectory))))))

(defn write-trajectory-csv!
  "Escribe la trayectoria en `path` (crea las carpetas necesarias)."
  ([trajectory path] (write-trajectory-csv! trajectory path 1))
  ([trajectory path stride]
   (io/make-parents path)
   (spit path (str/join "\n" (trajectory->csv-lines trajectory stride)))))

(defn summary-rows
  "Filas [experimento métrica valor] con las cifras de los dos experimentos, a
  partir de los resultados de `least-squares-experiment` y de
  `rosenbrock-experiment`. Es lo que respalda en un archivo las cifras que
  también se imprimen por consola (iteraciones, α, distancia a las ecuaciones
  normales, punto final)."
  [ls rb]
  (let [coord (fn [prefijo experimento v]
                (map-indexed (fn [i x] [experimento (str prefijo (inc i)) x]) v))]
    (vec (concat
          [["minimos-cuadrados" "observaciones" (get-in ls [:problem :observations])]
           ["minimos-cuadrados" "sigma_ruido" (get-in ls [:problem :sigma])]
           ["minimos-cuadrados" "semilla" (get-in ls [:problem :seed])]
           ["minimos-cuadrados" "alpha_seguro" (:alpha ls)]
           ["minimos-cuadrados" "iteraciones" (:iterations ls)]
           ["minimos-cuadrados" "motivo" (name (:reason ls))]
           ["minimos-cuadrados" "distancia_a_ecuaciones_normales" (:distance-to-normal-equations ls)]
           ["minimos-cuadrados" "f_final" (:f ls)]]
          (coord "w_descenso_" "minimos-cuadrados" (:x ls))
          (coord "w_ecuaciones_normales_" "minimos-cuadrados" (:w-normal-equations ls))
          (coord "w_verdadero_" "minimos-cuadrados" (:w-true ls))
          [["rosenbrock" "tolerancia_gradiente" rosenbrock-tol]
           ["rosenbrock" "iteraciones" (:iterations rb)]
           ["rosenbrock" "motivo" (name (:reason rb))]
           ["rosenbrock" "f_final" (:f rb)]]
          (coord "inicio_" "rosenbrock" rosenbrock-start)
          (coord "x_final_" "rosenbrock" (:x rb))))))

(defn summary-csv-lines
  "Líneas CSV `experiment,metric,value` de las filas de `summary-rows`."
  [rows]
  (cons "experiment,metric,value"
        (map #(str/join "," %) rows)))

;; ---------------------------------------------------------------------------
;; Punto de entrada
;; ---------------------------------------------------------------------------

(defn -main
  "Ejecuta los dos experimentos, imprime un resumen y escribe los CSV
  (trayectorias y resumen)."
  [& _args]
  (println "== Descenso de gradiente con AD ==")
  (let [r (least-squares-experiment)
        rb (rosenbrock-experiment)]
    (println "\n-- Mínimos cuadrados (100 observaciones, ruido σ = 0.1, semilla 42) --")
    (println "α seguro (0.9/‖X‖²_F):" (:alpha r))
    (println "Iteraciones:" (:iterations r) "| motivo:" (:reason r))
    (println "w por descenso:       " (:x r))
    (println "w por ec. normales:   " (:w-normal-equations r))
    (println "w verdadero:          " (:w-true r))
    (println "‖w_descenso − w_normales‖ =" (:distance-to-normal-equations r))
    (write-trajectory-csv! (:trajectory r) least-squares-csv)
    (println "CSV escrito en" least-squares-csv)
    (println "\n-- Rosenbrock desde (-1.2, 1) con backtracking --")
    (println "Iteraciones:" (:iterations rb) "| motivo:" (:reason rb))
    (println "x final:" (:x rb) "| f final:" (:f rb))
    (write-trajectory-csv! (:trajectory rb) rosenbrock-csv 100)
    (println "CSV escrito en" rosenbrock-csv "(una de cada 100 iteraciones)")
    (io/make-parents summary-csv)
    (spit summary-csv (str/join "\n" (summary-csv-lines (summary-rows r rb))))
    (println "CSV de resumen escrito en" summary-csv)))
