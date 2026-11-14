(ns fpnum.demo.conditioning
  "Demo ejecutable: unión de los dos proyectos. Mínimos cuadrados min ‖Xw − y‖²
  resueltos de dos formas, y cómo se comportan según el número de condición
  κ(X) = σmax/σmin.

  Métodos
  -------
  * Pseudoinversa por SVD (directo): w = X⁺y. Su costo no depende de κ y su
    error es del orden de κ·ε (en las mediciones, como mucho ~6 veces κ·ε
    hasta κ = 1e10; con κ = 1e12 llega a ~24 veces).
  * Descenso de gradiente con el gradiente por AD (iterativo): en cada
    iteración el error en la dirección singular vᵢ se multiplica por
    (1 − 2α·σᵢ²). El modo más lento es el de σmin:

        ρ = 1 − 2α·σmin²

    Con el paso seguro α = 0.9/‖X‖²_F ≈ 0.9/σmax² resulta ρ ≈ 1 − 1.8/κ², así
    que las iteraciones necesarias crecen como κ²: si κ se multiplica por 10,
    las iteraciones se multiplican por 100.

  Experimentos
  ------------
  A. Descenso frente a SVD: iteraciones medidas frente a la predicción
     ln(tol·√n)/ln|ρ|, para varios κ. También muestra que parar cuando ‖∇f‖ es
     pequeño es engañoso: ∇f = 2XᵀX(w − w*), que en la dirección lenta vale
     2σmin²·error, así que con κ grande un gradiente diminuto convive con un
     error grande.
  B. Precisión frente a κ: SVD (error del orden de κ·ε) frente a las
     ecuaciones normales (del orden de κ²·ε, y `solve` las declara singulares
     entre κ = 1e6 y 1e7). Con κ ≈ 1e3 los dos métodos dan errores
     comparables; la diferencia clara empieza hacia κ ≈ 1e4. Es la lección de
     la SVD (por qué no formar AᵀA) aplicada a mínimos cuadrados.
  C. Rango deficiente (una columna duplicada): las ecuaciones normales fallan,
     la SVD da la solución de mínima norma, y el descenso desde w = 0 converge
     a esa misma solución, porque el gradiente siempre está en el espacio fila
     de X y el descenso nunca se sale de él.

  Datos reproducibles: X = U·diag(σ)·Vᵀ con U y V ortonormales obtenidos de la
  SVD de matrices gaussianas con semilla fija, σ geométricos de 1 a 1/κ, y
  w_verdadero = V·(1, ..., 1)/√n (excita todos los modos por igual) e
  y = X·w_verdadero.

  Salidas en docs/data/:
    conditioning-gd.csv            experimento A (una fila por κ)
    conditioning-convergence.csv   distancia a la solución frente a iteración
    conditioning-accuracy.csv      experimento B
    conditioning-resumen.csv       exponente medido (iteraciones ∝ κ^p) y experimento C

  Uso:
    clojure -M -m fpnum.demo.conditioning              ; κ = 1, 3, 10, 30, 100, 1000
    clojure -M -m fpnum.demo.conditioning 1 10 100     ; κ para el experimento A"
  (:require [clojure.java.io :as io]
            [clojure.string :as str]
            [fpnum.linear :as la]
            [fpnum.optim :as opt]
            [fpnum.svd.core :as svd]))

(def ^:private eps 2.220446049250313E-16)

(def default-gd-kappas
  "Números de condición del experimento A (descenso frente a SVD). κ = 1000 no
  converge en el tope de iteraciones: es el contraste."
  [1.0 3.0 10.0 30.0 100.0 1000.0])

(def default-accuracy-kappas
  "Números de condición del experimento B (precisión): de 1 a 1e12, de 10 en 10."
  (mapv #(Math/pow 10.0 %) (range 0 13)))

(def default-data-dir
  "Carpeta donde se escriben los CSV del experimento."
  "docs/data")

;; ---------------------------------------------------------------------------
;; Datos del problema (puros)
;; ---------------------------------------------------------------------------

(defn gaussian-matrix
  "Matriz m×n con entradas gaussianas estándar, reproducible con `seed`."
  [m n seed]
  (let [rnd (java.util.Random. (long seed))]
    (vec (repeatedly m #(vec (repeatedly n (fn [] (.nextGaussian rnd))))))))

(defn geometric-sigmas
  "n valores singulares que decrecen geométricamente de 1 a 1/kappa:
  σᵢ = κ^(−i/(n−1)). Con κ = 1 son todos 1."
  [n kappa]
  (if (= n 1)
    [1.0]
    (mapv #(Math/pow kappa (- (/ (double %) (dec n)))) (range n))))

(defn matrix-with-singular-values
  "Matriz X de m×n (n = cantidad de `sigmas`, m ≥ n) con exactamente esos
  valores singulares: X = U·diag(σ)·Vᵀ, con U (m×n) y V (n×n) de columnas
  ortonormales, obtenidas de la SVD de matrices gaussianas con semilla fija.
  Devuelve {:x X :u U :v V}."
  [m sigmas seed]
  (let [n (count sigmas)
        u (:u (svd/svd (gaussian-matrix m n seed)))
        v (:v (svd/svd (gaussian-matrix n n (inc seed))))
        x (la/mat-mat (mapv (fn [fila] (mapv * fila sigmas)) u)
                      (la/transpose v))]
    {:x x :u u :v v}))

(defn make-problem
  "Problema de mínimos cuadrados de m×n con κ(X) = kappa.
  Devuelve {:x X :y y :w-true w :sigmas σ}, con w = V·(1, ..., 1)/√n (norma 1)
  e y = X·w (sistema consistente: la solución exacta es w)."
  [m n kappa seed]
  (let [sigmas (geometric-sigmas n kappa)
        {:keys [x v]} (matrix-with-singular-values m sigmas seed)
        w-true (la/scale (/ 1.0 (Math/sqrt n)) (la/mat-vec v (vec (repeat n 1.0))))]
    {:x x :y (la/mat-vec x w-true) :w-true w-true :sigmas sigmas}))

;; ---------------------------------------------------------------------------
;; Teoría del descenso
;; ---------------------------------------------------------------------------

(defn contraction-factor
  "ρ = 1 − 2·α·σmin²: factor por iteración del error en la dirección singular
  más lenta, para el descenso con paso α sobre ‖Xw − y‖²."
  [alpha sigma-min]
  (- 1.0 (* 2.0 alpha sigma-min sigma-min)))

(defn predicted-iterations
  "Iteraciones predichas para que el error relativo baje de `tol`, partiendo de
  w = 0 con w_verdadero = V·(1, ..., 1)/√n: la componente inicial en la
  dirección lenta vale 1/√n del error total, y decae como |ρ|ᵏ, de modo que

      k = ln(tol·√n) / ln|ρ|

  Es Infinity si |ρ| ≥ 1 (no converge). Es una aproximación: supone que el
  modo lento domina el error final, lo cual es cierto para κ moderado o
  grande."
  [alpha sigma-min n tol]
  (let [rho (Math/abs (double (contraction-factor alpha sigma-min)))]
    (cond (>= rho 1.0) Double/POSITIVE_INFINITY
          (zero? rho) 1.0
          :else (/ (Math/log (* tol (Math/sqrt n))) (Math/log rho)))))

(defn loglog-slope
  "Pendiente de ln(iteraciones) frente a ln(κ), por mínimos cuadrados, usando
  las filas del experimento A con κ ≥ `kappa-min` (por defecto 10) que
  alcanzaron la tolerancia. Para el descenso con paso seguro es ≈ 2
  (iteraciones ∝ κ²). Devuelve nil si hay menos de 2 puntos."
  ([rows] (loglog-slope rows 10.0))
  ([rows kappa-min]
   (let [pts (for [{:keys [kappa iters-dist]} rows
                   :when (and iters-dist (>= kappa kappa-min) (pos? iters-dist))]
               [(Math/log kappa) (Math/log iters-dist)])]
     (when (>= (count pts) 2)
       (let [n (count pts)
             mx (/ (reduce + (map first pts)) n)
             my (/ (reduce + (map second pts)) n)
             sxx (reduce + (map #(Math/pow (- (first %) mx) 2) pts))
             sxy (reduce + (map #(* (- (first %) mx) (- (second %) my)) pts))]
         (/ sxy sxx))))))

;; ---------------------------------------------------------------------------
;; Experimento A: descenso de gradiente frente a SVD
;; ---------------------------------------------------------------------------

(defn- decimar
  "Una de cada `stride` filas, y siempre la última."
  [tray stride]
  (let [ultima (:iter (peek tray))]
    (filter #(or (zero? (mod (:iter %) stride)) (= (:iter %) ultima)) tray)))

(defn gd-experiment
  "Descenso de gradiente (AD, paso seguro, desde w = 0) frente a la solución por
  SVD, para un problema con κ(X) = kappa.

  Opciones: :m (40), :n (4), :seed (42), :tol (1e-6, distancia relativa a la
  solución de la SVD), :grad-tol (1e-8), :max-iter (300000) y :points (400,
  puntos de la curva de convergencia).

  Se detiene cuando la distancia relativa ≤ tol Y ‖∇f‖ ≤ grad-tol, y registra
  en qué iteración se alcanzó cada criterio por separado.

  Devuelve la fila {:kappa :cond-measured :sigma-min :alpha :rho :iters-pred
  :iters-dist :iters-grad :dist-at-grad :final-dist :iterations :reason}
  y, en :convergence, los puntos [{:iter :dist}] de la curva."
  [kappa {:keys [m n seed tol grad-tol max-iter points]
          :or {m 40 n 4 seed 42 tol 1e-6 grad-tol 1e-8 max-iter 300000 points 400}}]
  (let [{:keys [x y]} (make-problem m n kappa seed)
        r (svd/svd x)
        w-star (opt/least-squares-svd x y)
        w-norm (la/norm2 w-star)
        dist (fn [w] (/ (la/norm2 (la/v- w w-star)) w-norm))
        alpha (opt/safe-step x)
        sigma-min (peek (:s r))
        res (opt/gradient-descent (opt/least-squares-objective x y) (la/zeros n)
                                  {:alpha alpha :max-iter max-iter :tol 0.0
                                   :stop-when (fn [{:keys [x grad-norm]}]
                                                (and (<= (dist x) tol) (<= grad-norm grad-tol)))})
        tray (:trajectory res)
        primera (fn [pred] (some #(when (pred %) (:iter %)) tray))
        it-dist (primera #(<= (dist (:x %)) tol))
        it-grad (primera #(<= (:grad-norm %) grad-tol))
        stride (max 1 (quot (count tray) points))]
    {:kappa kappa
     :cond-measured (svd/condition-number r)
     :sigma-min sigma-min
     :alpha alpha
     :rho (contraction-factor alpha sigma-min)
     :iters-pred (predicted-iterations alpha sigma-min n tol)
     :iters-dist it-dist
     :iters-grad it-grad
     :dist-at-grad (when it-grad (dist (:x (nth tray it-grad))))
     :final-dist (dist (:x res))
     :iterations (:iterations res)
     :reason (:reason res)
     :convergence (mapv (fn [{:keys [iter x]}] {:iter iter :dist (dist x)})
                        (decimar tray stride))}))

;; ---------------------------------------------------------------------------
;; Experimento B: precisión frente a κ
;; ---------------------------------------------------------------------------

(defn accuracy-experiment
  "Error relativo ‖w − w_verdadero‖/‖w_verdadero‖ de la solución por SVD y por
  ecuaciones normales, para un sistema consistente con κ(X) = kappa.

  Opciones: :m (30), :n (6) y :seed (11). Si las ecuaciones normales fallan
  por singulares, :err-normal es nil y :normal-status es \"singular\".

  Referencias: κ·ε (lo que se espera de la SVD) y κ²·ε (de las ecuaciones
  normales)."
  [kappa {:keys [m n seed] :or {m 30 n 6 seed 11}}]
  (let [{:keys [x y w-true]} (make-problem m n kappa seed)
        error (fn [w] (/ (la/norm2 (la/v- w w-true)) (la/norm2 w-true)))
        w-ne (try (opt/normal-equations x y)
                  (catch clojure.lang.ExceptionInfo e
                    (if (re-find #"singular" (ex-message e)) nil (throw e))))]
    {:kappa kappa
     :err-svd (error (opt/least-squares-svd x y))
     :err-normal (when w-ne (error w-ne))
     :normal-status (if w-ne "ok" "singular")
     :kappa-eps (* kappa eps)
     :kappa2-eps (* kappa kappa eps)}))

;; ---------------------------------------------------------------------------
;; Experimento C: rango deficiente
;; ---------------------------------------------------------------------------

(defn rank-deficient-experiment
  "X de m×4 con la cuarta columna igual a la primera (rango 3). Compara:
  ecuaciones normales (fallan), SVD (solución de mínima norma) y descenso de
  gradiente desde w = 0.

  Devuelve {:rank :normal-status :w-svd :w-gd :gd-iterations :gd-reason
  :dist-gd-svd :null-dot :norm-svd}. :null-dot es w_svd·(1, 0, 0, −1)/√2:
  la solución de mínima norma es ortogonal al espacio nulo de X, así que debe
  ser ≈ 0."
  ([] (rank-deficient-experiment {}))
  ([{:keys [m seed] :or {m 20 seed 3}}]
   (let [g (gaussian-matrix m 3 seed)
         x (mapv #(conj % (first %)) g)
         y (mapv first (gaussian-matrix m 1 (inc seed)))
         w-ne (try (opt/normal-equations x y)
                   (catch clojure.lang.ExceptionInfo e
                     (if (re-find #"singular" (ex-message e)) nil (throw e))))
         w-svd (opt/least-squares-svd x y)
         res (opt/gradient-descent (opt/least-squares-objective x y) (la/zeros 4)
                                   {:alpha (opt/safe-step x) :max-iter 100000 :tol 1e-12})
         nulo (la/scale (/ 1.0 (Math/sqrt 2.0)) [1.0 0.0 0.0 -1.0])]
     {:rank (svd/numerical-rank x)
      :normal-status (if w-ne "ok" "singular")
      :w-svd w-svd
      :w-gd (:x res)
      :gd-iterations (:iterations res)
      :gd-reason (:reason res)
      :dist-gd-svd (la/norm2 (la/v- (:x res) w-svd))
      :null-dot (la/dot w-svd nulo)
      :norm-svd (la/norm2 w-svd)})))

;; ---------------------------------------------------------------------------
;; Salida CSV
;; ---------------------------------------------------------------------------

(defn- celda [x] (if (nil? x) "" (str x)))

(defn gd-csv-lines
  "Líneas CSV del experimento A (una fila por κ). Las celdas sin valor
  (tolerancia no alcanzada) quedan vacías. `iterations` es el número de
  iteraciones realizadas (el tope, si no se alcanzó la tolerancia)."
  [rows]
  (cons "kappa,cond_measured,sigma_min,alpha,rho,iters_pred,iters_dist,iters_grad,dist_at_grad,final_dist,iterations,reason"
        (for [{:keys [kappa cond-measured sigma-min alpha rho iters-pred iters-dist
                      iters-grad dist-at-grad final-dist iterations reason]} rows]
          (str/join "," (map celda [kappa cond-measured sigma-min alpha rho iters-pred
                                    iters-dist iters-grad dist-at-grad final-dist
                                    iterations (some-> reason name)])))))

(defn convergence-csv-lines
  "Líneas CSV kappa,iter,dist con las curvas de convergencia de todas las filas."
  [rows]
  (cons "kappa,iter,dist"
        (for [{:keys [kappa convergence]} rows
              {:keys [iter dist]} convergence]
          (str/join "," [kappa iter dist]))))

(defn accuracy-csv-lines
  "Líneas CSV del experimento B. err_normal queda vacío si fue singular."
  [rows]
  (cons "kappa,err_svd,err_normal,normal_status,kappa_eps,kappa2_eps"
        (for [{:keys [kappa err-svd err-normal normal-status kappa-eps kappa2-eps]} rows]
          (str/join "," (map celda [kappa err-svd err-normal normal-status kappa-eps kappa2-eps])))))

(defn summary-rows
  "Filas [experimento métrica valor] que respaldan en un archivo las cifras de
  los experimentos A y C que también se imprimen por consola: el exponente
  medido de iteraciones ∝ κ^p (ajuste log-log con κ ≥ 10, ver `loglog-slope`)
  y los resultados del caso de rango deficiente. `gd` son las filas de
  `gd-experiment` y `rd` el resultado de `rank-deficient-experiment`."
  [gd rd]
  (vec (concat
        [["descenso-frente-a-kappa" "exponente_loglog_iteraciones_vs_kappa" (loglog-slope gd)]
         ["descenso-frente-a-kappa" "kappa_minimo_del_ajuste" 10.0]
         ["descenso-frente-a-kappa" "kappas_sin_alcanzar_tolerancia"
          (str/join " " (map #(long (:kappa %)) (remove :iters-dist gd)))]
         ["rango-deficiente" "rango_numerico" (:rank rd)]
         ["rango-deficiente" "ecuaciones_normales" (:normal-status rd)]
         ["rango-deficiente" "iteraciones_descenso" (:gd-iterations rd)]
         ["rango-deficiente" "motivo_descenso" (name (:gd-reason rd))]
         ["rango-deficiente" "distancia_descenso_a_svd" (:dist-gd-svd rd)]
         ["rango-deficiente" "producto_con_espacio_nulo" (:null-dot rd)]
         ["rango-deficiente" "norma_solucion_svd" (:norm-svd rd)]]
        (map-indexed (fn [i x] ["rango-deficiente" (str "w_svd_" (inc i)) x]) (:w-svd rd))
        (map-indexed (fn [i x] ["rango-deficiente" (str "w_descenso_" (inc i)) x]) (:w-gd rd)))))

(defn summary-csv-lines
  "Líneas CSV `experiment,metric,value` de las filas de `summary-rows`. Una
  celda sin valor (por ejemplo, un ajuste imposible) queda vacía."
  [rows]
  (cons "experiment,metric,value"
        (map (fn [[e m v]] (str/join "," [e m (celda v)])) rows)))

(defn- escribir! [lineas path]
  (io/make-parents path)
  (spit path (str/join "\n" lineas)))

;; ---------------------------------------------------------------------------
;; Punto de entrada
;; ---------------------------------------------------------------------------

(defn- fmt-num [x] (if (nil? x) "—" (format "%.3g" (double x))))

(defn run-experiments!
  "Ejecuta los tres experimentos, imprime los resultados y escribe los CSV.
  Opciones: :gd-kappas, :accuracy-kappas, :data-dir y las opciones de
  `gd-experiment` (:m, :n, :max-iter...). Devuelve {:gd :accuracy :rank-deficient}."
  ([] (run-experiments! {}))
  ([{:keys [gd-kappas accuracy-kappas data-dir]
     :or {gd-kappas default-gd-kappas
          accuracy-kappas default-accuracy-kappas
          data-dir default-data-dir}
     :as opts}]
   (println "== Mínimos cuadrados: SVD frente a descenso de gradiente (AD) ==")
   (println "\n-- A. Descenso de gradiente: iteraciones frente a κ --")
   (println "(tolerancia: distancia relativa a la solución de la SVD ≤ 1e-6)")
   (println "     κ | κ medido |   α    |    ρ     | iter. predichas | iter. medidas | iter. con ‖∇f‖≤1e-8 | dist. en ese punto | motivo")
   (let [gd (vec (for [k gd-kappas]
                   (let [r (gd-experiment k opts)]
                     (println (format "%6.0f | %8.3g | %6.4f | %8.6f | %15s | %13s | %19s | %18s | %s"
                                      (:kappa r) (:cond-measured r) (:alpha r) (:rho r)
                                      (fmt-num (:iters-pred r))
                                      (if (:iters-dist r) (str (:iters-dist r)) (str ">" (:iterations r)))
                                      (if (:iters-grad r) (str (:iters-grad r)) "—")
                                      (fmt-num (:dist-at-grad r))
                                      (name (:reason r))))
                     (when-not (:iters-dist r)
                       (println (format "       (no alcanzó la tolerancia en %d iteraciones; distancia relativa final: %.3g)"
                                        (:iterations r) (:final-dist r))))
                     r)))
         slope (loglog-slope gd)]
     (when slope
       (println (format "\nExponente medido: iteraciones ∝ κ^%.2f (predicción: 2, por ρ ≈ 1 − 1.8/κ²)" slope)))
     (println "\n-- B. Precisión frente a κ (sistema consistente, solución conocida) --")
     (println "         κ |  error SVD | error ec. normales | κ·ε       | κ²·ε")
     (let [acc (vec (for [k accuracy-kappas]
                      (let [r (accuracy-experiment k opts)]
                        (println (format "%10.0e | %10.2e | %18s | %9.2e | %9.2e"
                                         (:kappa r) (:err-svd r)
                                         (if (:err-normal r) (format "%.2e" (:err-normal r)) "singular")
                                         (:kappa-eps r) (:kappa2-eps r)))
                        r)))
           rd (rank-deficient-experiment)]
       (println "\n-- C. Rango deficiente (cuarta columna = primera) --")
       (println "rango numérico:" (:rank rd) "| ecuaciones normales:" (:normal-status rd))
       (println "w (SVD, mínima norma):" (mapv #(format "%.4f" %) (:w-svd rd)))
       (println "w (descenso desde 0): " (mapv #(format "%.4f" %) (:w-gd rd))
                (str "(" (:gd-iterations rd) " iteraciones, " (name (:gd-reason rd)) ")"))
       (println (format "‖w_descenso − w_SVD‖ = %.2e | w_SVD·(espacio nulo) = %.2e" (:dist-gd-svd rd) (:null-dot rd)))
       (escribir! (gd-csv-lines gd) (str data-dir "/conditioning-gd.csv"))
       (escribir! (convergence-csv-lines gd) (str data-dir "/conditioning-convergence.csv"))
       (escribir! (accuracy-csv-lines acc) (str data-dir "/conditioning-accuracy.csv"))
       (escribir! (summary-csv-lines (summary-rows gd rd)) (str data-dir "/conditioning-resumen.csv"))
       (println "\nCSV escritos en" data-dir "(conditioning-gd, conditioning-convergence, conditioning-accuracy, conditioning-resumen)")
       {:gd gd :accuracy acc :rank-deficient rd}))))

(defn -main
  "Sin argumentos usa κ = 1, 3, 10, 30, 100 y 1000 en el experimento A. Con
  argumentos, los toma como los κ de ese experimento."
  [& args]
  (if (seq args)
    (run-experiments! {:gd-kappas (mapv #(Double/parseDouble %) args)})
    (run-experiments!)))
