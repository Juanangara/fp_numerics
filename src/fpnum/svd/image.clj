(ns fpnum.svd.image
  "Compresión de imágenes con SVD de bajo rango.

  Idea
  ----
  Una imagen en color son tres matrices (canales R, G y B) con valores de 0 a
  255. Para cada canal A se calcula la SVD una sola vez, A = U Σ Vᵀ, y para
  cada k se conserva la aproximación de rango k (teorema de
  Eckart–Young–Mirsky: es la mejor posible):

      A_k = Σ_{i=1..k} σᵢ uᵢ vᵢᵀ,    ‖A − A_k‖₂ = σ_{k+1},    ‖A − A_k‖_F² = Σ_{i>k} σᵢ²

  Como cada término σᵢ uᵢ vᵢᵀ solo se suma una vez, pasar de un k al siguiente
  cuesta únicamente los términos nuevos (ver `compress`).

  El almacenamiento pasa de m·n a k(m+n+1) números por canal (k columnas de U,
  k de V y k valores singulares). Si k es grande, k(m+n+1) supera a m·n y
  \"comprimir\" ya no ahorra nada: la razón de compresión baja de 1.

  Métricas
  --------
  * Error relativo ‖A − A_k‖_F / ‖A‖_F (los tres canales juntos). Se calcula
    sobre la reconstrucción en punto flotante, de modo que coincide con la
    fórmula de σ.
  * Energía retenida Σ_{i≤k} σᵢ² / Σ σᵢ² (la matriz no se centra, así que k = 1
    ya retiene casi todo: captura el brillo medio).
  * PSNR = 10·log₁₀(255² / MSE), en dB. Se calcula sobre la imagen ya
    redondeada y recortada a [0, 255] (lo que realmente se vería y guardaría).
    Si la imagen es idéntica, MSE = 0 y el PSNR es +Infinity.
  * Razón de compresión m·n / (k(m+n+1)).

  Entrada y salida
  ----------------
  Este namespace, junto con fpnum.demo, es donde se permite la E/S: carga y
  guardado de imágenes (javax.imageio) y escritura de CSV y HTML. Los cálculos
  (`compress` y las métricas) son puros y trabajan sobre matrices de doubles
  como las de fpnum.linear.

  Los canales son matrices de h filas por w columnas (la fila es la coordenada
  y de la imagen)."
  (:require [clojure.java.io :as io]
            [clojure.string :as str]
            [fpnum.linear :as la]
            [fpnum.svd.core :as svd])
  (:import (java.awt Graphics2D RenderingHints)
           (java.awt.image BufferedImage)
           (java.io File)
           (javax.imageio ImageIO)))

(set! *warn-on-reflection* true)

;; ---------------------------------------------------------------------------
;; Imágenes: carga, redimensionado, canales y guardado (E/S)
;; ---------------------------------------------------------------------------

(defn load-image
  "Carga una imagen (PNG, JPEG, GIF o BMP) con ImageIO. Lanza `ex-info` si el
  archivo no existe o el formato no se reconoce."
  ^BufferedImage [path]
  (let [f (io/file path)]
    (when-not (.isFile f)
      (throw (ex-info "No existe la imagen" {:path (str path)})))
    (or (ImageIO/read f)
        (throw (ex-info "Formato de imagen no reconocido" {:path (str path)})))))

(defn- dibujar-escalada
  "Dibuja `src` escalada a w×h en una imagen RGB nueva, con la interpolación dada."
  ^BufferedImage [^BufferedImage src w h interpolacion]
  (let [out (BufferedImage. (int w) (int h) BufferedImage/TYPE_INT_RGB)
        ^Graphics2D g (.createGraphics out)]
    (.setRenderingHint g RenderingHints/KEY_INTERPOLATION interpolacion)
    (.drawImage g src 0 0 (int w) (int h) nil)
    (.dispose g)
    out))

(defn resize
  "Redimensiona `img` para que su lado mayor mida a lo sumo `max-side` píxeles,
  conservando la proporción. Devuelve una imagen RGB. Nunca agranda.

  Al reducir mucho (más de la mitad) se hace por pasos, mitad a mitad con
  interpolación bilineal, y un último ajuste bicúbico: reducir de golpe un
  factor grande salta píxeles y produce aliasing."
  ^BufferedImage [^BufferedImage img max-side]
  (let [w (.getWidth img)
        h (.getHeight img)
        escala (min 1.0 (/ (double max-side) (double (max w h))))
        nw (max 1 (int (Math/round (* w escala))))
        nh (max 1 (int (Math/round (* h escala))))]
    (loop [^BufferedImage cur img]
      (let [cw (.getWidth cur)
            ch (.getHeight cur)]
        (if (and (<= cw (* 2 nw)) (<= ch (* 2 nh)))
          (dibujar-escalada cur nw nh RenderingHints/VALUE_INTERPOLATION_BICUBIC)
          (recur (dibujar-escalada cur (max nw (quot cw 2)) (max nh (quot ch 2))
                                   RenderingHints/VALUE_INTERPOLATION_BILINEAR)))))))

(defn image->channels
  "Separa una imagen en sus canales: {:r :g :b}, cada uno una matriz de h filas
  por w columnas con valores 0.0–255.0."
  [^BufferedImage img]
  (let [w (.getWidth img)
        h (.getHeight img)
        ^ints px (.getRGB img 0 0 w h nil 0 w)
        canal (fn [desplazamiento]
                (vec (for [y (range h)]
                       (vec (for [x (range w)]
                              (double (bit-and (bit-shift-right (aget px (+ (* y w) x))
                                                                (long desplazamiento))
                                               255)))))))]
    {:r (canal 16) :g (canal 8) :b (canal 0)}))

(defn clamp-round
  "Redondea al entero más cercano y recorta a [0, 255]."
  [x]
  (let [v (Math/round (double x))]
    (cond (< v 0) 0 (> v 255) 255 :else v)))

(defn quantize
  "Matriz con cada valor redondeado y recortado a [0, 255] (lo que quedaría al
  guardar la imagen en 8 bits por canal)."
  [mat]
  (mapv (fn [fila] (mapv #(double (clamp-round %)) fila)) mat))

(defn channels->image
  "Arma una imagen RGB a partir de {:r :g :b}, redondeando y recortando cada
  valor a [0, 255]."
  ^BufferedImage [{:keys [r g b]}]
  (let [h (count r)
        w (count (first r))
        img (BufferedImage. (int w) (int h) BufferedImage/TYPE_INT_RGB)
        px (int-array (for [y (range h) x (range w)]
                        (unchecked-int
                         (bit-or (bit-shift-left (clamp-round (get-in r [y x])) 16)
                                 (bit-shift-left (clamp-round (get-in g [y x])) 8)
                                 (clamp-round (get-in b [y x]))))))]
    (.setRGB img 0 0 (int w) (int h) ^ints px 0 (int w))
    img))

(defn save-png!
  "Guarda `img` como PNG en `path` (crea las carpetas necesarias)."
  [^BufferedImage img path]
  (let [f (io/file path)]
    (io/make-parents f)
    (when-not (ImageIO/write img "png" f)
      (throw (ex-info "No se pudo escribir el PNG" {:path (str path)})))
    f))

(defn synthetic-channels
  "Imagen sintética reproducible de w×h (degradado de cielo, un disco, franjas
  diagonales y algo de ruido, con semilla fija). Sirve para pruebas y como
  respaldo cuando no hay una foto. Devuelve {:r :g :b}."
  [w h seed]
  (let [rnd (java.util.Random. (long seed))
        cx (* 0.65 w)
        cy (* 0.4 h)
        radio (* 0.22 (min w h))
        pixel (fn [x y]
                (let [t (/ (double y) h)
                      disco? (< (+ (Math/pow (- x cx) 2) (Math/pow (- y cy) 2)) (* radio radio))
                      franja? (and (> y (* 0.6 h)) (even? (quot (+ x y) 6)))
                      ruido (* 6.0 (.nextGaussian rnd))
                      [r g b] (cond disco? [235.0 200.0 60.0]
                                    franja? [40.0 120.0 60.0]
                                    :else [(+ 40.0 (* 80.0 t)) (+ 90.0 (* 90.0 t)) (- 230.0 (* 60.0 t))])]
                  [(+ r ruido) (+ g ruido) (+ b ruido)]))
        px (vec (for [y (range h)] (vec (for [x (range w)] (pixel x y)))))
        canal (fn [i] (mapv (fn [fila] (mapv #(double (clamp-round (nth % i))) fila)) px))]
    {:r (canal 0) :g (canal 1) :b (canal 2)}))

;; ---------------------------------------------------------------------------
;; Métricas (puras)
;; ---------------------------------------------------------------------------

(defn- sumar-cuadrados-matriz [a]
  (let [f (la/norm-fro a)] (* f f)))

(defn relative-error
  "‖A − B‖_F / ‖A‖_F, donde B aproxima a A. Si A es cero devuelve la norma de
  la diferencia."
  [a b]
  (let [num (la/norm-fro (mapv la/v- a b))
        den (la/norm-fro a)]
    (if (zero? den) num (/ num den))))

(defn mse
  "Error cuadrático medio entre dos conjuntos de canales {:r :g :b}, promediado
  sobre todos los píxeles y canales."
  [orig approx]
  (let [canales [:r :g :b]
        suma (reduce + (for [c canales]
                         (sumar-cuadrados-matriz (mapv la/v- (orig c) (approx c)))))
        n (reduce + (for [c canales] (* (count (orig c)) (count (first (orig c))))))]
    (/ suma n)))

(defn psnr
  "PSNR en dB entre dos conjuntos de canales {:r :g :b}: 10·log₁₀(255²/MSE).
  Si son idénticos (MSE = 0) devuelve +Infinity."
  [orig approx]
  (let [e (mse orig approx)]
    (if (zero? e)
      Double/POSITIVE_INFINITY
      (* 10.0 (Math/log10 (/ (* 255.0 255.0) e))))))

(defn energy-retained
  "Fracción de energía Σ_{i≤k} σᵢ² / Σ σᵢ² que retienen los k primeros valores
  singulares de `s`. Es 1.0 si la energía total es cero."
  [s k]
  (let [total (reduce + (map #(* % %) s))]
    (if (zero? total)
      1.0
      (/ (reduce + (map #(* % %) (take k s))) total))))

(defn compression-ratio
  "Razón de compresión m·n / (k(m+n+1)) de una matriz m×n con rango k."
  [m n k]
  (/ (double (* m n)) (double (* k (+ m n 1)))))

;; ---------------------------------------------------------------------------
;; Compresión (pura)
;; ---------------------------------------------------------------------------

(defn- sumar-terminos
  "A_k previa + los términos σᵢ uᵢ vᵢᵀ con i en [desde, hasta): actualiza una
  aproximación de rango `desde` a una de rango `hasta`."
  [acc {:keys [u s v]} desde hasta]
  (mapv (fn [fila-u fila-acc]
          (let [w (mapv #(* (s %) (fila-u %)) (range desde hasta))]
            (mapv (fn [fila-v a] (+ a (la/dot w (subvec fila-v desde hasta))))
                  v fila-acc)))
        u acc))

(defn default-ks
  "Valores de k razonables para mostrar, hasta `rank` (siempre incluye rank)."
  [rank]
  (vec (distinct (concat (filter #(<= % rank) [1 2 3 5 8 10 15 20 30 40 50 75 100 150 200 256])
                         [rank]))))

(defn compress
  "Calcula la SVD de cada canal de `channels` ({:r :g :b}) una sola vez y las
  aproximaciones de rango k para cada k de `ks`.

  `ks` debe ser estrictamente creciente y sus valores estar entre 1 y el rango
  máximo min(h, w); si no, lanza `ex-info`. Cada aproximación se obtiene
  sumando solo los términos nuevos a la anterior.

  Devuelve:
    :dims     [h w]
    :rank     min(h, w)
    :sigmas   {:r [...] :g [...] :b [...]} valores singulares de cada canal
    :original los canales de entrada
    :steps    un vector, uno por cada k, de mapas
                {:k k
                 :channels {:r A_k ...}      reconstrucción sin redondear
                 :metrics {:energy :rel-error :psnr :ratio}}"
  [channels ks]
  (let [{:keys [r]} channels
        h (count r)
        w (count (first r))
        rank (min h w)]
    (when (or (empty? ks)
              (not (apply < ks))
              (not (<= 1 (first ks)))
              (not (<= (last ks) rank)))
      (throw (ex-info "ks debe ser creciente y estar entre 1 y min(alto, ancho)"
                      {:op 'compress :ks ks :rank rank})))
    (let [svds (update-vals (select-keys channels [:r :g :b]) svd/svd)
          sigmas (update-vals svds :s)
          total (reduce + (for [c [:r :g :b] x (sigmas c)] (* x x)))
          ceros (la/zero-matrix h w)
          original-total (reduce + (for [c [:r :g :b]] (sumar-cuadrados-matriz (channels c))))
          inicial {:prev 0 :acc (zipmap [:r :g :b] (repeat ceros)) :steps []}
          {:keys [steps]}
          (reduce
           (fn [{:keys [prev acc steps]} k]
             (let [acc' (into {} (for [c [:r :g :b]]
                                   [c (sumar-terminos (acc c) (svds c) prev k)]))
                   dif (reduce + (for [c [:r :g :b]]
                                   (sumar-cuadrados-matriz (mapv la/v- (channels c) (acc' c)))))
                   energia (/ (reduce + (for [c [:r :g :b] x (take k (sigmas c))] (* x x)))
                              (if (zero? total) 1.0 total))
                   cuant (into {} (for [c [:r :g :b]] [c (quantize (acc' c))]))]
               {:prev k
                :acc acc'
                :steps (conj steps
                             {:k k
                              :channels acc'
                              :metrics {:energy energia
                                        :rel-error (if (zero? original-total)
                                                     (Math/sqrt dif)
                                                     (Math/sqrt (/ dif original-total)))
                                        :psnr (psnr channels cuant)
                                        :ratio (compression-ratio h w k)}})}))
           inicial
           ks)]
      {:dims [h w]
       :rank rank
       :sigmas sigmas
       :original channels
       :steps steps})))

;; ---------------------------------------------------------------------------
;; Salidas: CSV, JSON y HTML
;; ---------------------------------------------------------------------------

(defn singular-values-csv-lines
  "Líneas CSV `i,r,g,b` con los valores singulares de cada canal."
  [{:keys [sigmas rank]}]
  (cons "i,r,g,b"
        (for [i (range rank)]
          (str/join "," [(inc i) (get-in sigmas [:r i]) (get-in sigmas [:g i]) (get-in sigmas [:b i])]))))

(defn metrics-csv-lines
  "Líneas CSV `k,energy,rel_error,psnr_db,compression_ratio`, una por k."
  [{:keys [steps]}]
  (cons "k,energy,rel_error,psnr_db,compression_ratio"
        (for [{:keys [k metrics]} steps]
          (str/join "," [k (:energy metrics) (:rel-error metrics) (:psnr metrics) (:ratio metrics)]))))

(defn- json
  "Serializa a JSON datos simples: nil, booleanos, números (los no finitos
  salen como null), cadenas, keywords, mapas y secuencias."
  [x]
  (cond
    (nil? x) "null"
    (boolean? x) (str x)
    (number? x) (if (and (double? x) (or (Double/isNaN x) (Double/isInfinite x))) "null" (str x))
    (string? x) (str \" (str/escape x {\" "\\\"" \\ "\\\\" \newline "\\n"}) \")
    (keyword? x) (json (name x))
    (map? x) (str "{" (str/join "," (map (fn [[k v]] (str (json (name k)) ":" (json v))) x)) "}")
    (sequential? x) (str "[" (str/join "," (map json x)) "]")
    :else (throw (ex-info "Tipo no serializable a JSON" {:value x}))))

(defn image-file-name
  "Nombre del PNG de la reconstrucción con k valores singulares."
  [k]
  (format "k-%03d.png" k))

(defn viewer-data
  "Datos que consume la página HTML (mapa listo para serializar a JSON)."
  [{:keys [dims rank sigmas steps]}]
  (let [[h w] dims]
    {:m h :n w :rank rank
     :original "img/original.png"
     :sigmas sigmas
     :images (mapv #(str "img/" (image-file-name (:k %))) steps)
     :metrics (mapv (fn [{:keys [k metrics]}]
                      {:k k
                       :energy (:energy metrics)
                       :relError (:rel-error metrics)
                       :psnr (:psnr metrics)
                       :ratio (:ratio metrics)})
                    steps)}))

(defn render-viewer-html
  "Página HTML autocontenida (sin librerías externas) con un deslizador de k.
  Los datos van incrustados en JSON; las imágenes se referencian con rutas
  relativas (img/...)."
  [result]
  (let [plantilla (slurp (io/resource "fpnum/svd/viewer.html"))
        marca "__DATA__"
        datos (json (viewer-data result))]
    ;; str/replace con cadena literal (no regex), para que ningún carácter de los
    ;; datos se interprete como patrón.
    (str/replace plantilla marca (java.util.regex.Matcher/quoteReplacement datos))))

(defn write-outputs!
  "Escribe todos los resultados de `compress`:
    <html-dir>/index.html          página con el deslizador
    <html-dir>/img/original.png    imagen original (ya redimensionada)
    <html-dir>/img/k-NNN.png       una reconstrucción por cada k
    <data-dir>/svd-image-singular-values.csv
    <data-dir>/svd-image-metrics.csv
  Devuelve un mapa con las rutas escritas."
  [{:keys [original steps] :as result} {:keys [html-dir data-dir]}]
  (let [img-dir (io/file html-dir "img")
        sv-csv (io/file data-dir "svd-image-singular-values.csv")
        m-csv (io/file data-dir "svd-image-metrics.csv")
        html (io/file html-dir "index.html")]
    (save-png! (channels->image original) (io/file img-dir "original.png"))
    (doseq [{:keys [k channels]} steps]
      (save-png! (channels->image channels) (io/file img-dir (image-file-name k))))
    (doseq [[f lineas] [[sv-csv (singular-values-csv-lines result)]
                        [m-csv (metrics-csv-lines result)]]]
      (io/make-parents f)
      (spit f (str/join "\n" lineas)))
    (io/make-parents html)
    (spit html (render-viewer-html result))
    {:html (str html) :images (str img-dir) :singular-values-csv (str sv-csv) :metrics-csv (str m-csv)}))
