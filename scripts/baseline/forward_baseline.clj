;; COPIA CONGELADA de src/fpnum/ad/forward.clj tal como era antes del refactor a
;; operaciones genéricas por protocolo (forward original, solo con duales y números).
;; Solo sirve como línea base para la medición A/B de scripts/bench_forward_dispatch.clj.
;; NO está en el classpath del proyecto. No editar.
(ns fpnum.ad.forward-baseline
  "Diferenciación automática en modo forward con números duales.

  Teoría
  ------
  Un número dual es a + bε con ε² = 0 (ε ≠ 0). Guarda a la vez un valor (a)
  y una derivada (b). Al evaluar una función f sobre a + ε se obtiene
  exactamente

      f(a + ε) = f(a) + f'(a)·ε

  porque, al desarrollar f en serie de Taylor, todos los términos con ε² o
  más se anulan. La parte dual del resultado es la derivada.

  Cada operación propaga valor y derivada a la vez. Por ejemplo:

      (a + bε)(c + dε) = ac + (ad + bc)ε      ; regla del producto
      sin(a + bε)      = sin a + (cos a · b)ε  ; regla de la cadena

  Por eso la regla de la cadena aparece sola al componer funciones: la
  derivada de la salida de cada paso entra como derivada de la entrada del
  siguiente.

  Uso
  ---
  El usuario escribe sus funciones con las operaciones de este namespace
  (d+, d*, dsin...) y no con las de clojure.core, que no saben operar con
  duales. Los números normales se aceptan mezclados con duales y se tratan
  como constantes (derivada 0).

      (diff (fn [x] (d+ (d* x x) (dsin x))) 1.0)   ; => 2·1 + cos 1

  Dominios
  --------
  No se lanzan excepciones por problemas de dominio; se sigue el
  comportamiento de los doubles: dlog de un negativo da NaN y la división por
  cero da ±Infinity o NaN.

  Costo
  -----
  `diff` hace una pasada. `grad` hace una pasada por cada variable de entrada
  (n pasadas para n variables). Esa es la limitación del modo forward y la
  razón de ser del modo reverse (fpnum.ad.reverse).")

;; ---------------------------------------------------------------------------
;; Tipo Dual
;; ---------------------------------------------------------------------------

(defrecord Dual [a b])

(defn dual
  "Construye el número dual a + bε, con a y b convertidos a double."
  [a b]
  (->Dual (double a) (double b)))

(defn dual?
  "Verdadero si `x` es un Dual."
  [x]
  (instance? Dual x))

(defn lift
  "Convierte un número en el dual constante x + 0ε. Si ya es dual, lo
  devuelve tal cual."
  [x]
  (if (dual? x)
    x
    (->Dual (double x) 0.0)))

;; ---------------------------------------------------------------------------
;; Aritmética (binaria, privada)
;; ---------------------------------------------------------------------------

(defn- fdiv
  "División de doubles con semántica IEEE (x/0 da ±Infinity o NaN). Los
  argumentos primitivos son necesarios: la división de Clojure con números en
  caja lanza ArithmeticException al dividir por cero, incluso con doubles."
  ^double [^double x ^double y]
  (/ x y))

(defn- add2 [x y]
  (let [{a1 :a b1 :b} (lift x)
        {a2 :a b2 :b} (lift y)]
    (->Dual (+ a1 a2) (+ b1 b2))))

(defn- sub2 [x y]
  (let [{a1 :a b1 :b} (lift x)
        {a2 :a b2 :b} (lift y)]
    (->Dual (- a1 a2) (- b1 b2))))

(defn- mul2 [x y]
  (let [{a1 :a b1 :b} (lift x)
        {a2 :a b2 :b} (lift y)]
    (->Dual (* a1 a2) (+ (* a1 b2) (* b1 a2)))))

(defn- div2 [x y]
  (let [{a1 :a b1 :b} (lift x)
        {a2 :a b2 :b} (lift y)]
    (->Dual (fdiv a1 a2)
            (fdiv (- (* b1 a2) (* a1 b2)) (* a2 a2)))))

;; ---------------------------------------------------------------------------
;; Operaciones públicas
;; ---------------------------------------------------------------------------

(defn d+
  "Suma de duales (y/o números). Sin argumentos devuelve 0."
  ([] (->Dual 0.0 0.0))
  ([x] (lift x))
  ([x y] (add2 x y))
  ([x y & more] (reduce add2 (add2 x y) more)))

(defn d-
  "Resta de duales (y/o números). Con un argumento devuelve su opuesto."
  ([x] (sub2 0.0 x))
  ([x y] (sub2 x y))
  ([x y & more] (reduce sub2 (sub2 x y) more)))

(defn d*
  "Producto de duales (y/o números), con la regla del producto:
  (a + bε)(c + dε) = ac + (ad + bc)ε. Sin argumentos devuelve 1."
  ([] (->Dual 1.0 0.0))
  ([x] (lift x))
  ([x y] (mul2 x y))
  ([x y & more] (reduce mul2 (mul2 x y) more)))

(defn d-div
  "Cociente de duales (y/o números), con la regla del cociente:
  (a + bε)/(c + dε) = a/c + ((bc − ad)/c²)ε. Con un argumento devuelve su
  inverso.

  Se llama `d-div` y no `d/` porque en Clojure un símbolo no puede terminar
  en `/`."
  ([x] (div2 1.0 x))
  ([x y] (div2 x y))
  ([x y & more] (reduce div2 (div2 x y) more)))

(defn dsin
  "sin(a + bε) = sin a + (b·cos a)ε."
  [x]
  (let [{:keys [a b]} (lift x)]
    (->Dual (Math/sin a) (* b (Math/cos a)))))

(defn dcos
  "cos(a + bε) = cos a − (b·sin a)ε."
  [x]
  (let [{:keys [a b]} (lift x)]
    (->Dual (Math/cos a) (- (* b (Math/sin a))))))

(defn dexp
  "exp(a + bε) = eᵃ + (b·eᵃ)ε."
  [x]
  (let [{:keys [a b]} (lift x)
        e (Math/exp a)]
    (->Dual e (* b e))))

(defn dlog
  "log(a + bε) = log a + (b/a)ε. Para a < 0 da NaN, para a = 0 da -Infinity."
  [x]
  (let [{:keys [a b]} (lift x)]
    (->Dual (Math/log a) (fdiv b a))))

(defn dpow
  "Potencia x^y con base y exponente duales (o números).

  * Exponente constante n: (a + bε)ⁿ = aⁿ + (n·aⁿ⁻¹·b)ε. Se trata aparte para
    que funcione en a = 0 o con a negativo y n entero (por ejemplo x² en 0).
  * Caso general: d(xʸ) = xʸ·(y'·ln x + y·x'/x). Los términos cuyo factor de
    derivada es cero se omiten, para no producir 0·∞ = NaN espurios."
  [x y]
  (let [{a1 :a b1 :b} (lift x)
        {a2 :a b2 :b} (lift y)
        v (Math/pow a1 a2)]
    (if (zero? b2)
      (->Dual v (if (zero? a2)
                  0.0
                  (* a2 (Math/pow a1 (dec a2)) b1)))
      (let [t-exp  (* b2 (Math/log a1))
            t-base (if (zero? b1) 0.0 (fdiv (* a2 b1) a1))]
        (->Dual v (* v (+ t-exp t-base)))))))

;; ---------------------------------------------------------------------------
;; Derivadas
;; ---------------------------------------------------------------------------

(defn diff
  "Derivada de `f` (función de una variable) en el punto `x`.
  Evalúa f en x + ε y devuelve la parte dual del resultado. Si f devuelve un
  número normal (función constante), la derivada es 0."
  [f x]
  (:b (lift (f (->Dual (double x) 1.0)))))

(defn grad
  "Gradiente de `f` en el punto `xs` (vector de números).

  `f` recibe un vector de duales y devuelve un dual (o número). Para obtener
  la derivada parcial i se siembra ε solo en la componente i, así que se hacen
  n evaluaciones de f para n variables."
  [f xs]
  (let [n (count xs)]
    (mapv (fn [i]
            (let [semilla (vec (map-indexed
                                (fn [j x] (->Dual (double x) (if (= i j) 1.0 0.0)))
                                xs))]
              (:b (lift (f semilla)))))
          (range n))))
