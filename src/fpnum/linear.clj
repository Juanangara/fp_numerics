(ns fpnum.linear
  "Vectores y matrices inmutables.

  Representación
  --------------
  * Vector: vector de Clojure de doubles, por ejemplo [1.0 2.0 3.0].
  * Matriz: vector de filas (vectores de doubles de igual longitud), por
    ejemplo [[1.0 2.0] [3.0 4.0]] es una matriz de 2 filas y 2 columnas.

  Al ser datos inmutables, toda operación devuelve una estructura nueva y
  nunca modifica sus argumentos. Todas las funciones son puras.

  Errores
  -------
  Si las dimensiones no son compatibles, las operaciones lanzan `ex-info` con
  las dimensiones en el mapa de datos. Sin esa validación, `map` cortaría en
  la secuencia más corta y daría resultados incorrectos sin avisar.

  Rendimiento
  -----------
  Los doubles van en caja y el producto de matrices cuesta O(n³), así que esto
  prioriza la claridad sobre la velocidad. Esta capa es siempre inmutable. La
  versión rápida de la SVD (fpnum.svd.core/svd) usa arrays primitivos de doubles
  internamente, encapsulados dentro de una llamada, y recibe y devuelve las
  estructuras de este namespace.")

;; ---------------------------------------------------------------------------
;; Dimensiones y validación
;; ---------------------------------------------------------------------------

(defn dims
  "Dimensiones [filas columnas] de la matriz `a`.
  Una matriz sin filas se considera de dimensiones [0 0]."
  [a]
  [(count a) (count (first a))])

(defn- check-same-length!
  [op x y]
  (when (not= (count x) (count y))
    (throw (ex-info (str "Vectores de distinta longitud en " op)
                    {:op op :n1 (count x) :n2 (count y)}))))

(defn- check-mat-vec!
  [a x]
  (let [[m n] (dims a)]
    (when (not= n (count x))
      (throw (ex-info "Dimensiones incompatibles en mat-vec"
                      {:op 'mat-vec :matriz [m n] :vector (count x)})))))

(defn- check-mat-mat!
  [a b]
  (let [[m n] (dims a)
        [p q] (dims b)]
    (when (not= n p)
      (throw (ex-info "Dimensiones incompatibles en mat-mat"
                      {:op 'mat-mat :a [m n] :b [p q]})))))

;; ---------------------------------------------------------------------------
;; Constructores
;; ---------------------------------------------------------------------------

(defn zeros
  "Vector de `n` ceros."
  [n]
  (vec (repeat n 0.0)))

(defn zero-matrix
  "Matriz de ceros de `m` filas y `n` columnas."
  [m n]
  (vec (repeat m (zeros n))))

(defn identity-matrix
  "Matriz identidad de `n` x `n`."
  [n]
  (vec (for [i (range n)]
         (vec (for [j (range n)]
                (if (= i j) 1.0 0.0))))))

;; ---------------------------------------------------------------------------
;; Vectores
;; ---------------------------------------------------------------------------

(defn v+
  "Suma de dos vectores, componente a componente."
  [x y]
  (check-same-length! 'v+ x y)
  (mapv + x y))

(defn v-
  "Resta de dos vectores, componente a componente."
  [x y]
  (check-same-length! 'v- x y)
  (mapv - x y))

(defn scale
  "Multiplica el vector `x` por el escalar `k`."
  [k x]
  (mapv #(* k %) x))

(defn dot
  "Producto escalar x·y = Σ xᵢ yᵢ."
  [x y]
  (check-same-length! 'dot x y)
  (reduce + 0.0 (map * x y)))

(defn norm2
  "Norma euclídea ‖x‖₂ = √(x·x).

  Nota: no se reescala para evitar desbordamiento con valores enormes; basta
  para el rango de valores con que trabaja este proyecto."
  [x]
  (Math/sqrt (dot x x)))

;; ---------------------------------------------------------------------------
;; Matrices
;; ---------------------------------------------------------------------------

(defn transpose
  "Transpuesta de `a`: la fila i de Aᵀ es la columna i de A.
  (Aᵀ)ᵢⱼ = Aⱼᵢ."
  [a]
  (if (empty? a)
    []
    (apply mapv vector a)))

(defn mat-vec
  "Producto matriz-vector A·x. La componente i es el producto escalar de la
  fila i de A con x. El resultado tiene tantas componentes como filas tiene A."
  [a x]
  (check-mat-vec! a x)
  (mapv #(dot % x) a))

(defn mat-mat
  "Producto matriz-matriz A·B. El elemento (i, j) es el producto escalar de la
  fila i de A con la columna j de B.

  Se transpone B una sola vez, de modo que cada columna de B pasa a ser una
  fila y el producto escalar recorre dos vectores contiguos en vez de buscar
  una columna en cada fila."
  [a b]
  (check-mat-mat! a b)
  (let [bt (transpose b)]
    (mapv (fn [fila] (mapv #(dot fila %) bt)) a)))

(defn norm-fro
  "Norma de Frobenius ‖A‖_F = √(Σᵢⱼ aᵢⱼ²).
  Equivale a la norma euclídea de la matriz vista como un solo vector largo."
  [a]
  (Math/sqrt (reduce + 0.0 (for [fila a x fila] (* x x)))))

;; ---------------------------------------------------------------------------
;; Sistemas lineales
;; ---------------------------------------------------------------------------

(defn solve
  "Resuelve el sistema cuadrado A·x = b y devuelve x.

  Método: eliminación gaussiana con pivoteo parcial.
  1. Eliminación hacia adelante: para cada columna k se elige como pivote la
     fila (desde k hacia abajo) con mayor valor absoluto en esa columna, se
     intercambia con la fila k y se resta un múltiplo de ella a las filas de
     abajo para anular la columna. Elegir el mayor pivote evita dividir por
     números pequeños, que amplificarían el error de redondeo.
  2. Sustitución hacia atrás sobre el sistema triangular resultante.

  Lanza `ex-info` si A no es cuadrada, si las dimensiones no coinciden con b,
  o si A es singular. Se considera singular cuando un pivote es, en valor
  absoluto, menor o igual que 1e-12 veces el mayor |aᵢⱼ| de A (tolerancia
  relativa a la escala de la matriz)."
  [a b]
  (let [n (count a)]
    (when (or (zero? n) (some #(not= n (count %)) a))
      (throw (ex-info "solve requiere una matriz cuadrada"
                      {:op 'solve :dims (dims a)})))
    (when (not= n (count b))
      (throw (ex-info "Dimensiones incompatibles en solve"
                      {:op 'solve :matriz [n n] :b (count b)})))
    (let [escala (reduce max 0.0 (for [fila a x fila] (Math/abs (double x))))
          tol    (* 1e-12 escala)
          aug    (mapv conj a b)
          tri    (reduce
                  (fn [m k]
                    (let [p      (apply max-key #(Math/abs (double (get-in m [% k])))
                                        (range k n))
                          pivote (get-in m [p k])]
                      (when (<= (Math/abs (double pivote)) tol)
                        (throw (ex-info "Matriz singular: no se puede resolver el sistema"
                                        {:op 'solve :columna k :pivote pivote})))
                      (let [m (-> m (assoc k (m p)) (assoc p (m k)))]
                        (reduce (fn [m i]
                                  (let [factor (/ (get-in m [i k]) pivote)]
                                    (assoc m i (v- (m i) (scale factor (m k))))))
                                m
                                (range (inc k) n)))))
                  aug
                  (range n))]
      (reduce (fn [x i]
                (let [fila (tri i)
                      resto (dot (subvec fila (inc i) n) (subvec x (inc i) n))]
                  (assoc x i (/ (- (peek fila) resto) (fila i)))))
              (zeros n)
              (range (dec n) -1 -1)))))
