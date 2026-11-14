(ns fpnum.demo.ad
  "Demo ejecutable de diferenciación automática.

  Compara, para f(x) = sin x en x = 1, el error de la derivada por diferencias
  finitas hacia adelante, (f(x+h) − f(x))/h, con el de la diferenciación
  automática, y guarda la tabla en un CSV para graficar la curva en U.

  La curva en U aparece por dos efectos opuestos:
  * h grande: error de truncamiento, proporcional a h.
  * h pequeño: error de cancelación, porque f(x+h) − f(x) resta números casi
    iguales y se pierden cifras; el error crece como ε_máquina / h.
  El mínimo queda cerca de h ≈ √ε_máquina ≈ 1e-8. La AD no depende de h.

  Uso: clojure -M -m fpnum.demo.ad"
  (:require [clojure.java.io :as io]
            [clojure.string :as str]
            [fpnum.ad.forward :as ad]))

(def csv-path
  "CSV con el error de la derivada por diferencias finitas y por AD frente a h."
  "docs/data/error-vs-h.csv")

(defn finite-difference
  "Derivada de `f` en `x` por diferencias finitas hacia adelante con paso `h`."
  [f x h]
  (/ (- (f (+ x h)) (f x)) h))

(defn error-table
  "Filas {:h :fd-error :ad-error} para f = sin en x = 1, con h = 10^e para
  cada exponente e entre -1 y -16 en pasos de 0.5. El valor exacto de
  referencia es cos(1)."
  []
  (let [x 1.0
        exacta (Math/cos x)
        ad-error (Math/abs (- (ad/diff ad/dsin x) exacta))]
    (vec (for [i (range 0 31)
               :let [h (Math/pow 10.0 (- -1.0 (* 0.5 i)))]]
           {:h h
            :fd-error (Math/abs (- (finite-difference #(Math/sin %) x h) exacta))
            :ad-error ad-error}))))

(defn best-row
  "Fila de la tabla con menor error por diferencias finitas."
  [rows]
  (apply min-key :fd-error rows))

(defn write-csv!
  "Escribe la tabla en `path` (crea las carpetas necesarias)."
  [rows path]
  (io/make-parents path)
  (spit path
        (str/join "\n"
                  (cons "h,fd_error,ad_error"
                        (map (fn [{:keys [h fd-error ad-error]}]
                               (str/join "," [h fd-error ad-error]))
                             rows)))))

(defn -main
  "Muestra derivadas de ejemplo y genera el CSV de error frente a h."
  [& _args]
  (println "== Diferenciación automática (modo forward) ==")
  (let [f (fn [x] (ad/d+ (ad/d* x x) (ad/dsin x)))]
    (println "f(x) = x² + sin x;  f'(1) por AD =" (ad/diff f 1.0)
             "; exacta =" (+ 2.0 (Math/cos 1.0))))
  (let [g (fn [[x y]] (ad/d+ (ad/d* x x y) (ad/dsin y)))]
    (println "g(x,y) = x²y + sin y;  ∇g(1,2) por AD =" (ad/grad g [1.0 2.0])
             "; exacto =" [(* 2.0 1.0 2.0) (+ 1.0 (Math/cos 2.0))]))
  (let [rows (error-table)
        mejor (best-row rows)]
    (write-csv! rows csv-path)
    (println "Mejor h por diferencias finitas:" (:h mejor)
             "con error" (:fd-error mejor))
    (println "Error de la AD:" (:ad-error (first rows)))
    (println "CSV escrito en" csv-path)))
