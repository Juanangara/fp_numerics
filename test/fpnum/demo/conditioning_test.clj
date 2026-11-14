(ns fpnum.demo.conditioning-test
  "Pruebas de fpnum.demo.conditioning: construcción de problemas con κ
  conocido, predicción teórica de las iteraciones del descenso, experimentos A,
  B y C, y salida CSV. Los experimentos usan dimensiones pequeñas para ser
  rápidos."
  (:require [clojure.java.io :as io]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]
            [fpnum.demo.conditioning :as demo]
            [fpnum.linear :as la]
            [fpnum.svd.core :as svd]))

(defn- approx= [a b tol]
  (<= (Math/abs (- (double a) (double b)))
      (* tol (+ 1.0 (Math/abs (double a)) (Math/abs (double b))))))

(defn- error-ortonormal [q]
  (let [g (la/mat-mat (la/transpose q) q)
        n (count g)]
    (apply max (for [i (range n) j (range n)]
                 (Math/abs (- (get-in g [i j]) (if (= i j) 1.0 0.0)))))))

(defmacro ^:private con-temp-dir [[sym] & cuerpo]
  `(let [~sym (doto (java.io.File/createTempFile "conditioning-test" "") (.delete) (.mkdirs))]
     (try ~@cuerpo
          (finally (doseq [^java.io.File f# (reverse (file-seq ~sym))] (.delete f#))))))

;; ---------------------------------------------------------------------------
;; Construcción de problemas
;; ---------------------------------------------------------------------------

(deftest valores-singulares-geometricos
  (is (= [1.0] (demo/geometric-sigmas 1 100.0)))
  (is (= [1.0 1.0 1.0] (demo/geometric-sigmas 3 1.0)))
  (let [s (demo/geometric-sigmas 3 100.0)]
    (is (every? true? (map #(approx= %1 %2 1e-12) [1.0 0.1 0.01] s))))
  (testing "κ = σ₁/σₙ"
    (let [s (demo/geometric-sigmas 5 1e6)]
      (is (approx= 1e6 (/ (first s) (peek s)) 1e-9)))))

(deftest matriz-con-valores-singulares-dados
  (let [sigmas [1.0 0.3 0.05 0.001]
        {:keys [x u v]} (demo/matrix-with-singular-values 12 sigmas 7)]
    (is (= [12 4] (la/dims x)))
    (testing "U y V tienen columnas ortonormales"
      (is (< (error-ortonormal u) 1e-12))
      (is (< (error-ortonormal v) 1e-12)))
    (testing "los valores singulares de X son los pedidos"
      (is (every? true? (map #(approx= %1 %2 1e-9) sigmas (:s (svd/svd x))))))))

(deftest problema-reproducible-y-consistente
  (let [p (demo/make-problem 20 3 100.0 5)]
    (is (= p (demo/make-problem 20 3 100.0 5)))
    (is (not= (:x p) (:x (demo/make-problem 20 3 100.0 6))))
    (testing "w_verdadero tiene norma 1 e y = X·w"
      (is (approx= 1.0 (la/norm2 (:w-true p)) 1e-12))
      (is (= (:y p) (la/mat-vec (:x p) (:w-true p)))))
    (testing "κ(X) es el pedido"
      (is (approx= 100.0 (svd/condition-number (:x p)) 1e-8)))))

;; ---------------------------------------------------------------------------
;; Teoría
;; ---------------------------------------------------------------------------

(deftest factor-de-contraccion
  (is (approx= 0.982 (demo/contraction-factor 0.9 0.1) 1e-12))
  (is (approx= 0.55 (demo/contraction-factor 0.225 1.0) 1e-12)))

(deftest iteraciones-predichas
  (testing "α = 0.9, σmin = 0.1, n = 4: ρ = 0.982, k = ln(2e-6)/ln(0.982) ≈ 722.5"
    (is (approx= 722.5 (demo/predicted-iterations 0.9 0.1 4 1e-6) 1e-3)))
  (testing "si |ρ| ≥ 1 no converge"
    (is (= Double/POSITIVE_INFINITY (demo/predicted-iterations 3.0 1.0 4 1e-6))))
  (testing "crece como κ²: 10 veces más κ ⇒ ~100 veces más iteraciones"
    (let [alpha 0.9
          it (fn [kappa] (demo/predicted-iterations alpha (/ 1.0 kappa) 4 1e-6))]
      (is (< 90.0 (/ (it 1000.0) (it 100.0)) 110.0)))))

(deftest pendiente-log-log
  (testing "iteraciones = c·κ² da pendiente 2"
    (let [rows (for [k [10.0 30.0 100.0]] {:kappa k :iters-dist (int (* 7.7 k k))})]
      (is (approx= 2.0 (demo/loglog-slope rows) 1e-3))))
  (testing "ignora κ < 10 y las filas que no alcanzaron la tolerancia"
    (let [rows [{:kappa 1.0 :iters-dist 24} {:kappa 10.0 :iters-dist 770}
                {:kappa 100.0 :iters-dist 77000} {:kappa 1000.0 :iters-dist nil}]]
      (is (approx= 2.0 (demo/loglog-slope rows) 1e-3))))
  (testing "con menos de 2 puntos no hay pendiente"
    (is (nil? (demo/loglog-slope [{:kappa 10.0 :iters-dist 700}])))
    (is (nil? (demo/loglog-slope [])))))

;; ---------------------------------------------------------------------------
;; Experimento A: descenso frente a SVD
;; ---------------------------------------------------------------------------

(def ^:private chico {:m 20 :n 3})

(deftest descenso-bien-condicionado
  (let [r (demo/gd-experiment 1.0 chico)]
    (is (approx= 1.0 (:cond-measured r) 1e-9))
    (is (some? (:iters-dist r)))
    (is (< (:iters-dist r) 60))
    (is (<= (:final-dist r) (* 1.0001 1e-6)))
    (is (= :stopped (:reason r)))))

(deftest las-iteraciones-medidas-siguen-a-la-prediccion
  (let [r (demo/gd-experiment 10.0 chico)
        medidas (:iters-dist r)
        predichas (:iters-pred r)]
    (is (some? medidas))
    (is (< 0.8 (/ medidas predichas) 1.25)
        (str "medidas " medidas ", predichas " predichas))
    (testing "el factor ρ medido en la curva coincide con la teoría"
      (let [c (:convergence r)
            [a b] [(nth c (quot (count c) 2)) (peek c)]
            rho-medido (Math/pow (/ (:dist b) (:dist a)) (/ 1.0 (- (:iter b) (:iter a))))]
        (is (approx= (:rho r) rho-medido 1e-3))))))

(deftest las-iteraciones-crecen-como-kappa-al-cuadrado
  (let [rows (mapv #(demo/gd-experiment % chico) [10.0 20.0 40.0])
        its (mapv :iters-dist rows)]
    (is (every? some? its))
    (is (apply < its))
    (testing "duplicar κ multiplica las iteraciones por ≈ 4"
      (is (< 3.0 (/ (second its) (first its)) 5.0))
      (is (< 3.0 (/ (nth its 2) (second its)) 5.0)))
    (testing "el exponente medido es ≈ 2"
      (is (< 1.8 (demo/loglog-slope rows) 2.2)))))

(deftest un-gradiente-pequeno-no-implica-estar-cerca
  (let [r (demo/gd-experiment 30.0 chico)]
    (is (some? (:iters-grad r)))
    (is (> (:dist-at-grad r) 1e-6)
        "al llegar a ‖∇f‖ ≤ 1e-8 el error real todavía supera la tolerancia de 1e-6")
    (is (> (:dist-at-grad r) 2e-6))))

(deftest descenso-que-no-llega-a-la-tolerancia
  (let [r (demo/gd-experiment 100.0 (assoc chico :max-iter 200))]
    (is (nil? (:iters-dist r)))
    (is (= :max-iter (:reason r)))
    (is (= 200 (:iterations r)))
    (is (> (:final-dist r) 0.1) "con κ = 100 en 200 iteraciones casi no se acercó")))

(deftest curva-de-convergencia
  (let [r (demo/gd-experiment 10.0 (assoc chico :points 50))
        c (:convergence r)]
    (is (<= (count c) 100))
    (is (= 0 (:iter (first c))))
    (is (= (:iterations r) (:iter (peek c))) "incluye siempre el último punto")
    (is (approx= 1.0 (:dist (first c)) 1e-9) "parte de w = 0: distancia relativa 1")
    (is (apply >= (map :dist c)) "la distancia no sube con el paso seguro")))

;; ---------------------------------------------------------------------------
;; Experimento B: precisión frente a κ
;; ---------------------------------------------------------------------------

(deftest precision-bien-condicionado
  (let [r (demo/accuracy-experiment 1.0 {})]
    (is (= "ok" (:normal-status r)))
    (is (< (:err-svd r) 1e-13))
    (is (< (:err-normal r) 1e-13))))

(deftest las-ecuaciones-normales-pierden-el-doble-de-digitos
  ;; Con κ ≈ 1e3 ambos métodos dan errores comparables (~1e-12); la diferencia
  ;; clara empieza hacia κ ≈ 1e4 y a κ = 1e5 es de varios órdenes de magnitud.
  (let [r (demo/accuracy-experiment 1e5 {})]
    (is (= "ok" (:normal-status r)))
    (is (< (:err-svd r) (* 1000.0 (:kappa-eps r))) "SVD: error del orden de κ·ε")
    (is (> (:err-normal r) (* 100.0 (:err-svd r))) "ec. normales: error del orden de κ²·ε, mucho peor")
    (is (< (:err-normal r) (* 10.0 (:kappa2-eps r))) "y no peor que κ²·ε")))

(deftest con-kappa-enorme-las-ecuaciones-normales-fallan-y-la-svd-no
  (let [r (demo/accuracy-experiment 1e10 {})]
    (is (= "singular" (:normal-status r)))
    (is (nil? (:err-normal r)))
    (is (< (:err-svd r) 1e-3) "la SVD aún da la solución con ~6 cifras")))

(deftest error-de-la-svd-crece-como-kappa-por-epsilon
  (doseq [k [1.0 1e3 1e6 1e9]]
    (let [r (demo/accuracy-experiment k {})]
      (is (< (:err-svd r) (+ 1e-15 (* 1000.0 (:kappa-eps r)))) (str "κ = " k)))))

;; ---------------------------------------------------------------------------
;; Experimento C: rango deficiente
;; ---------------------------------------------------------------------------

(deftest rango-deficiente
  (let [r (demo/rank-deficient-experiment)]
    (is (= 3 (:rank r)))
    (is (= "singular" (:normal-status r)) "las ecuaciones normales fallan")
    (testing "el descenso desde 0 converge a la solución de mínima norma de la SVD"
      (is (= :converged (:gd-reason r)))
      (is (< (:dist-gd-svd r) 1e-8)))
    (testing "w_SVD es ortogonal al espacio nulo (1, 0, 0, −1)/√2"
      (is (< (Math/abs (:null-dot r)) 1e-10)))
    (testing "la columna duplicada recibe el mismo peso: w₁ = w₄"
      (is (approx= (first (:w-svd r)) (peek (:w-svd r)) 1e-9)))))

;; ---------------------------------------------------------------------------
;; Salida CSV y ejecución completa
;; ---------------------------------------------------------------------------

(deftest lineas-csv
  (let [fila {:kappa 10.0 :cond-measured 10.0 :sigma-min 0.1 :alpha 0.9 :rho 0.982
              :iters-pred 722.5 :iters-dist nil :iters-grad 700 :dist-at-grad nil
              :final-dist 0.5 :iterations 300000 :reason :max-iter
              :convergence [{:iter 0 :dist 1.0} {:iter 5 :dist 0.9}]}
        gd (demo/gd-csv-lines [fila])
        conv (demo/convergence-csv-lines [fila])
        acc (demo/accuracy-csv-lines [{:kappa 1e10 :err-svd 1e-6 :err-normal nil
                                       :normal-status "singular" :kappa-eps 1e-6 :kappa2-eps 1e4}])]
    (is (= "kappa,cond_measured,sigma_min,alpha,rho,iters_pred,iters_dist,iters_grad,dist_at_grad,final_dist,iterations,reason"
           (first gd)))
    (testing "las celdas sin valor quedan vacías, y los motivos se escriben sin dos puntos"
      (is (= 12 (count (str/split (second gd) #"," -1))))
      (is (str/includes? (second gd) ",0.5,300000,max-iter") "final_dist, iterations y reason al final")
      (is (str/includes? (second gd) ",,"))
      (is (str/ends-with? (second gd) ",max-iter")))
    (is (= ["kappa,iter,dist" "10.0,0,1.0" "10.0,5,0.9"] conv))
    (is (= 6 (count (str/split (second acc) #"," -1))))
    (is (str/includes? (second acc) ",,singular,"))))

(deftest resumen-del-experimento
  (let [gd [{:kappa 10.0 :iters-dist 921} {:kappa 20.0 :iters-dist 3700}
            {:kappa 40.0 :iters-dist 14800} {:kappa 1000.0 :iters-dist nil}]
        rd {:rank 3 :normal-status "singular" :gd-iterations 81 :gd-reason :converged
            :dist-gd-svd 2.8e-14 :null-dot 1.4e-17 :norm-svd 0.2
            :w-svd [0.1 0.2 0.3 0.1] :w-gd [0.1 0.2 0.3 0.1]}
        filas (demo/summary-rows gd rd)
        lineas (demo/summary-csv-lines filas)
        por-metrica (into {} (map (fn [[e m v]] [[e m] v]) filas))]
    (testing "el exponente es el ajuste log-log de las filas con κ ≥ 10 que alcanzaron la tolerancia"
      (is (= (demo/loglog-slope gd) (get por-metrica ["descenso-frente-a-kappa" "exponente_loglog_iteraciones_vs_kappa"])))
      (is (< 1.9 (get por-metrica ["descenso-frente-a-kappa" "exponente_loglog_iteraciones_vs_kappa"]) 2.1)))
    (testing "los κ que no alcanzaron la tolerancia se listan"
      (is (= "1000" (get por-metrica ["descenso-frente-a-kappa" "kappas_sin_alcanzar_tolerancia"]))))
    (testing "las cifras del experimento C"
      (is (= 3 (get por-metrica ["rango-deficiente" "rango_numerico"])))
      (is (= 81 (get por-metrica ["rango-deficiente" "iteraciones_descenso"])))
      (is (= "converged" (get por-metrica ["rango-deficiente" "motivo_descenso"]))))
    (testing "formato CSV"
      (is (= "experiment,metric,value" (first lineas)))
      (is (every? #(= 3 (count (str/split % #"," -1))) lineas))
      (is (= (inc (count filas)) (count lineas))))
    (testing "sin ajuste posible (menos de 2 puntos) la celda queda vacía"
      (let [l (demo/summary-csv-lines (demo/summary-rows [{:kappa 10.0 :iters-dist 900}] rd))]
        (is (some #(= "descenso-frente-a-kappa,exponente_loglog_iteraciones_vs_kappa," %) l))))))

(deftest ejecucion-completa-con-tamanos-pequenos
  (con-temp-dir [d]
    (let [res (volatile! nil)
          salida (with-out-str
                   (vreset! res (demo/run-experiments!
                                 {:gd-kappas [10.0 20.0 40.0]
                                  :accuracy-kappas [1.0 1e3 1e10]
                                  :data-dir (str d)
                                  :m 20 :n 3})))
          {:keys [gd accuracy rank-deficient]} @res]
      (is (= 3 (count gd)))
      (is (= 3 (count accuracy)))
      (is (= 3 (:rank rank-deficient)))
      (testing "imprime los tres experimentos y el exponente medido"
        (is (str/includes? salida "A. Descenso de gradiente"))
        (is (str/includes? salida "B. Precisión frente a κ"))
        (is (str/includes? salida "C. Rango deficiente"))
        (is (re-find #"iteraciones ∝ κ\^[12]\.\d\d" salida)))
      (testing "escribe los tres CSV"
        (is (= 4 (count (str/split-lines (slurp (io/file d "conditioning-gd.csv"))))))
        (is (= 4 (count (str/split-lines (slurp (io/file d "conditioning-accuracy.csv"))))))
        (is (str/starts-with? (slurp (io/file d "conditioning-convergence.csv")) "kappa,iter,dist"))
        (is (str/starts-with? (slurp (io/file d "conditioning-resumen.csv")) "experiment,metric,value"))))))
