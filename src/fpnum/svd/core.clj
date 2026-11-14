(ns fpnum.svd.core
  "Descomposición en valores singulares (SVD) por el método de Jacobi de un
  solo lado (Hestenes).

  Dos implementaciones, mismo algoritmo y misma API
  --------------------------------------------------
  * `svd-reference`: implementación de referencia. Estructuras inmutables
    (vectores de vectores) y funciones puras; es la más fácil de leer y de
    verificar, pero lenta (256×256 tarda decenas de segundos).
  * `svd`: implementación rápida. Hace exactamente las mismas operaciones en
    el mismo orden, pero sobre arrays `double[]` privados y sin recalcular los
    productos escalares que ya conoce (ver abajo). Los arrays nacen y mueren
    dentro de la llamada: la función sigue siendo pura de cara al exterior, y
    recibe y devuelve vectores inmutables.

  Las pruebas diferenciales exigen que ambas coincidan.

  Teoría
  ------
  Toda matriz A de m×n se descompone como A = U Σ Vᵀ, con U y V de columnas
  ortonormales y Σ diagonal con los valores singulares σ₁ ≥ σ₂ ≥ ... ≥ 0.
  Aquí se usa la forma económica: con k = min(m, n), U es m×k, V es n×k y hay
  k valores singulares.

  Idea del método (supongamos m ≥ n). Se busca una matriz ortogonal V tal que
  las columnas de B = A·V sean ortogonales entre sí. Entonces

      B = U Σ   con   σⱼ = ‖bⱼ‖  y  uⱼ = bⱼ/σⱼ,   así que   A = U Σ Vᵀ.

  V se construye con rotaciones planas. Para un par de columnas aₚ, a_q:

      α = aₚ·aₚ,   β = a_q·a_q,   γ = aₚ·a_q
      ζ = (β − α)/(2γ),   t = sgn(ζ)/(|ζ| + √(1 + ζ²)),   c = 1/√(1 + t²),   s = c·t

      aₚ' =  c·aₚ − s·a_q
      a_q' = s·aₚ + c·a_q

  El valor de t es la raíz de menor módulo de t² + 2ζt − 1 = 0, y con él las
  dos columnas nuevas quedan ortogonales (aₚ'·a_q' = 0). Las mismas rotaciones
  se aplican a las columnas de V, que parte de la identidad. Se barren todos
  los pares (p, q) repetidamente hasta que ninguno necesite rotación. La
  convergencia final es cuadrática, y en la práctica bastan unos 6 a 13
  barridos.

  Por qué no se forma AᵀA
  -----------------------
  Los valores singulares de A son las raíces de los autovalores de AᵀA, pero
  formar AᵀA eleva al cuadrado el número de condición: κ(AᵀA) = κ(A)². Con
  κ(A) = 1e10 serían 1e20, más allá de la precisión de un double, y los
  valores singulares pequeños se perderían por completo. Jacobi de un solo
  lado rota las columnas de A directamente, sin formar AᵀA, y conserva buena
  precisión relativa en los valores singulares pequeños.

  Qué hace más rápida a `svd`
  ---------------------------
  1. Arrays `double[]` primitivos: sin números en caja ni asignación de
     vectores nuevos en cada rotación.
  2. Normas de columna en caché. El método necesita α = ‖aₚ‖² y β = ‖a_q‖² en
     cada par, y un barrido tiene n(n−1)/2 pares: recalcularlas serían dos
     productos escalares de largo m por par, además de γ. En cambio se guardan
     las n normas², y al rotar un par se calculan las nuevas en el mismo bucle
     que escribe las columnas (suma de los cuadrados de las entradas recién
     escritas). Así cada par cuesta un producto escalar (γ) y una pasada de
     rotación, en vez de tres productos escalares y dos pasadas. Al sumarse en
     el mismo orden que `la/dot`, el valor en caché es idéntico al que
     recalcularía la versión de referencia, de modo que no hay deriva numérica
     (a diferencia de la fórmula α' = α − t·γ, que sufre cancelación).

  Casos especiales
  ----------------
  * Escala: antes de iterar se multiplica A por una potencia de 2 (exacto) que
    deja su mayor entrada en [1, 2), y los σ se reescalan al final. Sin eso,
    entradas muy pequeñas (≈1e-76) hacen que α·β se desborde por abajo a 0, y
    muy grandes que α se desborde a infinito; en ambos casos la medida de
    ortogonalidad se calcula mal y las columnas no se rotan.
  * Si m < n se calcula la SVD de Aᵀ y se intercambian U y V.
  * Si una columna queda numéricamente nula (σⱼ ≤ 4·max(m, n)·ε·‖A‖_F, es
    decir, del tamaño del ruido de redondeo; ocurre con matrices de rango
    deficiente), uⱼ = bⱼ/σⱼ no está determinada, y rotar ese ruido contra las
    demás columnas no converge. Esas columnas no se rotan, y su uⱼ se completa
    con un vector ortonormal a las demás para que U siga teniendo columnas
    ortonormales. El σⱼ reportado sigue siendo la norma real de la columna.
    Consecuencia: en matrices con valores singulares por debajo de ese umbral,
    los vectores singulares asociados no están determinados.

  Representación: igual que en fpnum.linear (vector de filas)."
  (:require [fpnum.linear :as la]))

(set! *warn-on-reflection* true)

(def ^:private eps 2.220446049250313E-16)

;; ---------------------------------------------------------------------------
;; Fórmulas compartidas por las dos implementaciones
;; ---------------------------------------------------------------------------

(defn- parametros-rotacion
  "Dado α = ‖aₚ‖², β = ‖a_q‖² y γ = aₚ·a_q (γ ≠ 0), devuelve [c s] de la
  rotación que ortogonaliza las dos columnas."
  [alpha beta gamma]
  (let [zeta (/ (- beta alpha) (* 2.0 gamma))
        signo (if (neg? zeta) -1.0 1.0)
        t (/ signo (+ (Math/abs zeta) (Math/sqrt (+ 1.0 (* zeta zeta)))))
        c (/ 1.0 (Math/sqrt (+ 1.0 (* t t))))]
    [c (* c t)]))

(defn- umbral-nulo
  "Norma por debajo de la cual una columna se considera numéricamente nula:
  4·max(m, n)·ε·‖A‖_F. Es del orden del error de redondeo absoluto con que
  cualquier algoritmo estable (LAPACK incluido) conoce un valor singular, así
  que lo que queda por debajo es indistinguible de cero. El factor 4 deja
  margen: en la matriz de Hilbert 10×10, κ ≈ 1.6e13, el valor singular menor
  es unas 8 veces mayor que este umbral y se conserva."
  [a-mat]
  (let [[m n] (la/dims a-mat)]
    (* 4.0 (max m n) eps (la/norm-fro a-mat))))

;; ---------------------------------------------------------------------------
;; Ensamblado de U, Σ y V (compartido)
;; ---------------------------------------------------------------------------

(defn- quitar-proyecciones
  "Resta a `x` su proyección sobre cada vector (ortonormal) de `base`."
  [base x]
  (reduce (fn [r u] (la/v- r (la/scale (la/dot u r) u))) x base))

(defn- completar-columna
  "Vector unitario ortogonal a todos los de `base` (ortonormales). Parte del
  vector canónico eₖ con mayor componente fuera del subespacio de `base`, y
  reortogonaliza dos veces para controlar el redondeo."
  [base m]
  (let [canonico (fn [k] (assoc (la/zeros m) k 1.0))
        mejor (apply max-key
                     (fn [k] (la/norm2 (quitar-proyecciones base (canonico k))))
                     (range m))
        r (->> (canonico mejor)
               (quitar-proyecciones base)
               (quitar-proyecciones base))]
    (la/scale (/ 1.0 (la/norm2 r)) r)))

(defn- columnas-u
  "Columnas de U: bⱼ/σⱼ, y para σⱼ numéricamente nulo (≤ `umbral`) una
  columna completada. `sigmas` y `cols` ya están ordenados."
  [cols sigmas m umbral]
  (reduce (fn [acc [col sigma]]
            (conj acc (if (> sigma umbral)
                        (la/scale (/ 1.0 sigma) col)
                        (completar-columna acc m))))
          []
          (map vector cols sigmas)))

(defn- ensamblar
  "A partir de las columnas ya ortogonalizadas `a` (B = A·V) y las de `v`
  (vectores de vectores columna), ordena y arma {:u :s :v :sweeps :converged?}."
  [a v sweeps converged? m n umbral]
  (let [sigmas (mapv la/norm2 a)
        orden (vec (sort-by #(- (sigmas %)) (range n)))
        s (mapv sigmas orden)]
    {:u (la/transpose (columnas-u (mapv a orden) s m umbral))
     :s s
     :v (la/transpose (mapv v orden))
     :sweeps sweeps
     :converged? converged?}))

;; ---------------------------------------------------------------------------
;; Implementación de referencia (inmutable)
;; ---------------------------------------------------------------------------

(defn- rotar
  "Devuelve [x' y'] con x' = c·x − s·y e y' = s·x + c·y."
  [c s x y]
  [(mapv (fn [xi yi] (- (* c xi) (* s yi))) x y)
   (mapv (fn [xi yi] (+ (* s xi) (* c yi))) x y)])

(defn- procesar-par
  "Si las columnas p y q de `a` no son ortogonales (con tolerancia relativa
  `tol`), las rota, y aplica la misma rotación a las de `v`.
  Devuelve [a' v' medida], con medida = |aₚ·a_q| / (‖aₚ‖‖a_q‖) antes de rotar.

  Un par con una columna numéricamente nula (norma² ≤ `umbral2`) se deja
  como está: esa columna es ruido de redondeo, y su ortogonalidad relativa
  con las demás no tiene sentido ni converge."
  [[a v off] [p q] tol umbral2]
  (let [ap (a p)
        aq (a q)
        alpha (la/dot ap ap)
        beta (la/dot aq aq)
        gamma (la/dot ap aq)
        denom (Math/sqrt (* alpha beta))
        rel (if (and (pos? denom) (> alpha umbral2) (> beta umbral2))
              (/ (Math/abs (double gamma)) denom)
              0.0)]
    (if (<= rel tol)
      [a v (max off rel)]
      (let [[c s] (parametros-rotacion alpha beta gamma)
            [ap' aq'] (rotar c s ap aq)
            [vp' vq'] (rotar c s (v p) (v q))]
        [(assoc a p ap' q aq')
         (assoc v p vp' q vq')
         (max off rel)]))))

(defn- barrido
  "Un barrido completo sobre todos los pares (p, q) con p < q.
  Devuelve [a' v' medida-máxima]."
  [a v tol umbral2]
  (let [n (count a)
        pares (for [p (range n) q (range (inc p) n)] [p q])]
    (reduce #(procesar-par %1 %2 tol umbral2) [a v 0.0] pares)))

(defn- jacobi
  "Itera barridos hasta que todas las columnas sean ortogonales (medida ≤ tol)
  o se alcance `max-sweeps`. `a` y `v` son vectores de columnas.
  Devuelve {:a :v :sweeps :converged?}."
  [a v tol umbral2 max-sweeps]
  (loop [a a v v sweeps 0]
    (if (>= sweeps max-sweeps)
      {:a a :v v :sweeps sweeps :converged? false}
      (let [[a' v' off] (barrido a v tol umbral2)]
        (if (<= off tol)
          {:a a' :v v' :sweeps (inc sweeps) :converged? true}
          (recur a' v' (inc sweeps)))))))

(defn- svd-alta-reference
  "SVD para m ≥ n, versión inmutable. Devuelve {:u :s :v :sweeps :converged?}."
  [a-mat tol max-sweeps]
  (let [[m n] (la/dims a-mat)
        umbral (umbral-nulo a-mat)
        {:keys [a v sweeps converged?]}
        (jacobi (la/transpose a-mat) (la/identity-matrix n) tol (* umbral umbral) max-sweeps)]
    (ensamblar a v sweeps converged? m n umbral)))

;; ---------------------------------------------------------------------------
;; Implementación rápida (arrays double[] encapsulados)
;; ---------------------------------------------------------------------------
;; Convención: un double[][] guarda una columna por posición (a[j] es la
;; columna j). Estos arrays son locales a una llamada de `svd` y se mutan
;; solo aquí; nunca salen de este namespace.

(defn- sumar-cuadrados
  "Σ xᵢ², sumando de izquierda a derecha como `la/dot` (mismo orden, mismo
  resultado exacto)."
  ^double [^doubles x]
  (let [m (alength x)]
    (loop [i 0 acc 0.0]
      (if (< i m)
        (let [xi (aget x i)]
          (recur (inc i) (+ acc (* xi xi))))
        acc))))

(defn- producto-escalar
  "x·y sumando de izquierda a derecha como `la/dot`. Es una función aparte, con
  retorno primitivo, porque un `loop` usado como expresión se compila como una
  función anónima que devuelve un objeto en caja y el compilador pierde el
  tipo double."
  ^double [^doubles x ^doubles y]
  (let [m (alength x)]
    (loop [i 0 acc 0.0]
      (if (< i m)
        (recur (inc i) (+ acc (* (aget x i) (aget y i))))
        acc))))

(defn- par!
  "Procesa el par de columnas (p, q): calcula γ, y si no son ortogonales las
  rota EN SITIO (columnas de `a` y de `v`) y actualiza sus normas² en `nrm2`.
  Devuelve la medida de ortogonalidad relativa antes de rotar.

  Solo se calcula un producto escalar (γ). α y β salen de `nrm2`; las normas
  nuevas se acumulan en el mismo bucle que escribe las columnas rotadas."
  [^"[[D" a ^"[[D" v ^doubles nrm2 p q tol umbral2]
  (let [p (long p)
        q (long q)
        tol (double tol)
        umbral2 (double umbral2)
        ^doubles ap (aget a p)
        ^doubles aq (aget a q)
        m (alength ap)
        alpha (aget nrm2 p)
        beta (aget nrm2 q)
        gamma (producto-escalar ap aq)
        denom (Math/sqrt (* alpha beta))
        rel (if (and (> denom 0.0) (> alpha umbral2) (> beta umbral2))
              (/ (Math/abs gamma) denom)
              0.0)]
    (when (> rel tol)
      (let [zeta (/ (- beta alpha) (* 2.0 gamma))
            signo (if (neg? zeta) -1.0 1.0)
            t (/ signo (+ (Math/abs zeta) (Math/sqrt (+ 1.0 (* zeta zeta)))))
            c (/ 1.0 (Math/sqrt (+ 1.0 (* t t))))
            s (* c t)]
        ;; Rotar las columnas de A y acumular sus normas² nuevas.
        (loop [i 0 sa 0.0 sb 0.0]
          (if (< i m)
            (let [x (aget ap i)
                  y (aget aq i)
                  nx (- (* c x) (* s y))
                  ny (+ (* s x) (* c y))]
              (aset ap i nx)
              (aset aq i ny)
              (recur (inc i) (+ sa (* nx nx)) (+ sb (* ny ny))))
            (do (aset nrm2 p sa)
                (aset nrm2 q sb))))
        ;; Misma rotación sobre las columnas de V.
        (let [^doubles vp (aget v p)
              ^doubles vq (aget v q)
              nv (alength vp)]
          (loop [i 0]
            (when (< i nv)
              (let [x (aget vp i)
                    y (aget vq i)]
                (aset vp i (- (* c x) (* s y)))
                (aset vq i (+ (* s x) (* c y))))
              (recur (inc i)))))))
    rel))

(defn- barrido!
  "Un barrido completo sobre los pares (p, q) con p < q, mutando `a`, `v` y
  `nrm2`. Devuelve la mayor medida de ortogonalidad vista."
  [^"[[D" a ^"[[D" v ^doubles nrm2 tol umbral2]
  (let [n (alength a)]
    (loop [p 0 off 0.0]
      (if (< p n)
        (recur (inc p)
               (double
                (loop [q (inc p) off off]
                  (if (< q n)
                    (recur (inc q)
                           (Math/max off (double (par! a v nrm2 p q tol umbral2))))
                    off))))
        off))))

(defn- columnas->arrays
  "Matriz (vector de filas) → double[][] con una columna por posición."
  ^"[[D" [a-mat]
  (into-array (map double-array (la/transpose a-mat))))

(defn- identidad-columnas
  "Identidad n×n como double[][] de columnas."
  ^"[[D" [n]
  (into-array (map (fn [j]
                     (let [^doubles col (double-array n)]
                       (aset col (int j) 1.0)
                       col))
                   (range n))))

(defn- svd-alta-rapida
  "SVD para m ≥ n, versión con arrays. Mismo resultado que `svd-alta-reference`."
  [a-mat tol max-sweeps]
  (let [[m n] (la/dims a-mat)
        umbral (umbral-nulo a-mat)
        umbral2 (* umbral umbral)
        a (columnas->arrays a-mat)
        v (identidad-columnas n)
        nrm2 (double-array (map #(sumar-cuadrados %) a))]
    (loop [sweeps 0]
      (if (>= sweeps max-sweeps)
        (ensamblar (mapv vec a) (mapv vec v) sweeps false m n umbral)
        (let [off (double (barrido! a v nrm2 tol umbral2))]
          (if (<= off tol)
            (ensamblar (mapv vec a) (mapv vec v) (inc sweeps) true m n umbral)
            (recur (inc sweeps))))))))

;; ---------------------------------------------------------------------------
;; API
;; ---------------------------------------------------------------------------

(defn- exponente-de-escala
  "Exponente binario e del mayor |aᵢⱼ| (0 si A es cero)."
  [a]
  (let [amax (reduce max 0.0 (for [fila a x fila] (Math/abs (double x))))]
    (if (pos? amax) (Math/getExponent (double amax)) 0)))

(defn- escalar-por-potencia-de-2
  "Multiplica cada entrada por 2^k. Es exacto (solo cambia el exponente), y
  evita desbordamientos y desbordamientos por abajo al elevar al cuadrado
  entradas muy grandes o muy pequeñas: tras escalar, el mayor |aᵢⱼ| queda en
  [1, 2)."
  [a k]
  (if (zero? k)
    a
    (mapv (fn [fila] (mapv #(Math/scalb (double %) (int k)) fila)) a)))

(defn- validar!
  [a]
  (let [[_ n] (la/dims a)]
    (when (or (empty? a) (zero? n) (some #(not= n (count %)) a))
      (throw (ex-info "svd requiere una matriz no vacía con todas las filas de igual longitud"
                      {:op 'svd :dims (la/dims a)})))))

(defn- svd*
  "Parte común de `svd` y `svd-reference`: validación, escalado, transposición
  si m < n y reescalado de los σ. `svd-alta` es la implementación para m ≥ n."
  [svd-alta a {:keys [tol max-sweeps] :or {max-sweeps 60}}]
  (validar! a)
  (let [[m n] (la/dims a)
        tol (or tol (* 2.0 (max m n) eps))
        e (exponente-de-escala a)
        a' (escalar-por-potencia-de-2 a (- e))
        r (if (>= m n)
            (svd-alta a' tol max-sweeps)
            (let [{:keys [u v] :as r} (svd-alta (la/transpose a') tol max-sweeps)]
              (assoc r :u v :v u)))]
    ;; U y V no dependen de la escala; solo los valores singulares.
    (update r :s (fn [s] (mapv #(Math/scalb (double %) (int e)) s)))))

(defn svd
  "SVD económica de la matriz `a` (m×n, vector de filas) por Jacobi de un solo
  lado. Con k = min(m, n) devuelve un mapa:

    :u           matriz m×k de columnas ortonormales
    :s           vector de k valores singulares, no negativos y en orden
                 decreciente
    :v           matriz n×k de columnas ortonormales
    :sweeps      barridos realizados
    :converged?  true si las columnas quedaron ortogonales dentro de la
                 tolerancia antes de agotar `:max-sweeps`

  de modo que A = U·diag(s)·Vᵀ (ver `reconstruct`).

  Es la versión rápida: trabaja internamente con arrays double[] y mantiene en
  caché las normas de las columnas (ver la documentación del namespace). De
  cara al exterior es pura: no modifica `a` y devuelve estructuras inmutables.
  Da los mismos resultados que `svd-reference`.

  Opciones:
    :tol         tolerancia de ortogonalidad relativa entre columnas
                 (por defecto 2·max(m, n)·ε, con ε ≈ 2.2e-16; no puede ser más
                 pequeña que el error de redondeo de un producto escalar)
    :max-sweeps  máximo de barridos (por defecto 60)

  Lanza `ex-info` si la matriz está vacía o sus filas no tienen igual longitud."
  ([a] (svd a {}))
  ([a opts] (svd* svd-alta-rapida a opts)))

(defn svd-reference
  "Implementación de referencia de `svd`: mismo algoritmo, misma API y mismo
  formato de salida, pero con estructuras inmutables (vectores de vectores) y
  recalculando los productos escalares en cada par. Es más lenta, y sirve
  para leer el algoritmo sin distracciones y para validar `svd` con pruebas
  diferenciales. Ver `svd` para el detalle de la salida y las opciones."
  ([a] (svd-reference a {}))
  ([a opts] (svd* svd-alta-reference a opts)))

(defn truncate
  "Aproximación de rango `k`: conserva las k primeras columnas de U y V y los k
  primeros valores singulares (los mayores, porque `s` viene ordenado).

  Por el teorema de Eckart–Young–Mirsky, U_k·Σ_k·V_kᵀ es la mejor aproximación
  de rango k a A, tanto en norma 2 (‖A − A_k‖₂ = σ_{k+1}) como en norma de
  Frobenius (‖A − A_k‖_F² = Σ_{i>k} σᵢ²).

  Lanza `ex-info` si k no está entre 1 y el número de valores singulares."
  [{:keys [u s v] :as r} k]
  (when-not (and (integer? k) (<= 1 k (count s)))
    (throw (ex-info "k debe estar entre 1 y el número de valores singulares"
                    {:op 'truncate :k k :max (count s)})))
  (assoc r
         :u (mapv #(subvec % 0 k) u)
         :s (subvec (vec s) 0 k)
         :v (mapv #(subvec % 0 k) v)))

(defn reconstruct
  "Reconstruye U·diag(s)·Vᵀ a partir del resultado de `svd`."
  [{:keys [u s v]}]
  (la/mat-mat (mapv (fn [fila] (mapv * fila s)) u)
              (la/transpose v)))

;; ---------------------------------------------------------------------------
;; Número de condición, rango numérico y pseudoinversa
;; ---------------------------------------------------------------------------
;; Las funciones siguientes aceptan una matriz o un resultado de `svd` (si ya
;; se calculó la SVD, se reutiliza en vez de repetirla).

(defn- ->svd
  [a-or-svd]
  (if (map? a-or-svd) a-or-svd (svd a-or-svd)))

(defn condition-number
  "Número de condición en norma 2: κ(A) = σmax/σmin. Es +Infinity si σmin = 0.

  Mide cuánto amplifica A los errores relativos: resolver A·w = y con datos
  perturbados en ε da una solución perturbada en torno a κ(A)·ε. Si κ ≈ 1e16
  la matriz es numéricamente singular en doble precisión."
  [a-or-svd]
  (let [s (:s (->svd a-or-svd))
        s-min (peek (vec s))]
    (if (zero? s-min)
      Double/POSITIVE_INFINITY
      (/ (first s) s-min))))

(defn- tolerancia-de-rango
  "σ por debajo de la cual un valor singular se trata como cero:
  rcond·σmax, con rcond = max(m, n)·ε por defecto (el mismo criterio que
  NumPy)."
  [{:keys [u v s]} rcond]
  (* (or rcond (* (max (count u) (count v)) eps)) (first s)))

(defn numerical-rank
  "Rango numérico: cantidad de valores singulares mayores que rcond·σmax
  (por defecto rcond = max(m, n)·ε).

  `rcond` es opcional (segundo argumento). Con aritmética de punto flotante el
  rango exacto no se puede decidir; lo que se mide es cuántas direcciones son
  distinguibles del ruido de redondeo."
  ([a-or-svd] (numerical-rank a-or-svd nil))
  ([a-or-svd rcond]
   (let [r (->svd a-or-svd)
         tol (tolerancia-de-rango r rcond)]
     (count (filter #(> % tol) (:s r))))))

(defn inverse-singular-values
  "Vector con 1/σᵢ para los σᵢ mayores que rcond·σmax y 0 para el resto: es la
  diagonal de Σ⁺. Ignorar los σ pequeños es lo que regulariza la pseudoinversa
  (ver `pinv`)."
  ([a-or-svd] (inverse-singular-values a-or-svd nil))
  ([a-or-svd rcond]
   (let [r (->svd a-or-svd)
         tol (tolerancia-de-rango r rcond)]
     (mapv #(if (> % tol) (/ 1.0 %) 0.0) (:s r)))))

(defn pinv
  "Pseudoinversa de Moore–Penrose A⁺ = V·Σ⁺·Uᵀ, de n×m si A es m×n.

  Σ⁺ invierte los valores singulares por encima de la tolerancia
  rcond·σmax (por defecto rcond = max(m, n)·ε) y deja en cero el resto. Cumple
  las cuatro condiciones de Moore–Penrose:

      A·A⁺·A = A,   A⁺·A·A⁺ = A⁺,   (A·A⁺)ᵀ = A·A⁺,   (A⁺·A)ᵀ = A⁺·A

  Si A tiene columnas independientes, A⁺ = (AᵀA)⁻¹Aᵀ y A⁺·y es la solución de
  mínimos cuadrados. Si A tiene rango deficiente, A⁺·y es, entre todas las
  soluciones de mínimos cuadrados, la de menor norma ‖w‖."
  ([a-or-svd] (pinv a-or-svd nil))
  ([a-or-svd rcond]
   (let [r (->svd a-or-svd)
         inv (inverse-singular-values r rcond)]
     (la/mat-mat (mapv (fn [fila] (mapv * fila inv)) (:v r))
                 (la/transpose (:u r))))))
