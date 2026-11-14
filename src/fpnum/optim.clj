(ns fpnum.optim
  "Descenso de gradiente y mínimos cuadrados, con el gradiente calculado por
  diferenciación automática.

  Contenido: `gradient-descent` (paso fijo o backtracking), el objetivo de
  Rosenbrock, los objetivos y soluciones de mínimos cuadrados (ecuaciones
  normales con `normal-equations`, y pseudoinversa vía SVD con
  `least-squares-svd`), la cota de paso seguro `safe-step` y datos de
  regresión reproducibles (`noisy-regression-data`).

  Teoría
  ------
  Para minimizar f se repite

      x ← x − α·∇f(x)

  El gradiente apunta en la dirección de máximo crecimiento de f, así que
  avanzar en sentido contrario la hace decrecer. El tamaño del paso α es
  delicado: si es muy pequeño la convergencia es lenta, y si es muy grande el
  método oscila o diverge.

  Dos estrategias para elegir α:
  * Paso fijo (`:alpha`). Simple; para mínimos cuadrados hay una cota segura
    (ver `safe-step`).
  * Backtracking con la condición de Armijo (`:backtracking? true`). En cada
    iteración se parte de un α grande y se reduce (α ← β·α) hasta que

        f(x − α·g) ≤ f(x) − c·α·‖g‖²,   con g = ∇f(x)

    es decir, hasta que f baje al menos una fracción c de lo que predice la
    aproximación lineal. Así cada paso aceptado hace decrecer f.

  Dos opciones más de `gradient-descent`:
  * `:grad-fn` elige cómo se calcula el gradiente: por defecto
    fpnum.ad.forward/grad (n pasadas), o fpnum.ad.reverse/grad (una pasada hacia
    adelante y una hacia atrás). Con el mismo objetivo dan el mismo resultado.
  * `:stop-when` permite parar con un criterio propio, por ejemplo la distancia
    a una solución conocida: con un problema mal condicionado, ‖∇f‖ pequeño no
    garantiza estar cerca del mínimo.

  Las funciones objetivo se escriben con las operaciones genéricas (d+, d*,
  etc., de fpnum.ad.forward o fpnum.ad.ops) y reciben un vector de duales o de
  nodos, según el modo.

  El núcleo es puro. Los datos ruidosos usan java.util.Random con semilla fija
  solo como generador determinista: la misma semilla da siempre los mismos
  datos."
  (:require [fpnum.ad.forward :as ad]
            [fpnum.linear :as la]
            [fpnum.svd.core :as svd]))

;; ---------------------------------------------------------------------------
;; Funciones objetivo
;; ---------------------------------------------------------------------------

(defn rosenbrock
  "Función de Rosenbrock f(x, y) = (1 − x)² + 100(y − x²)².
  Recibe un vector de duales [x y]. Su mínimo global es f(1, 1) = 0, en el
  fondo de un valle curvo y estrecho, lo que la hace difícil para el descenso
  de gradiente."
  [[x y]]
  (let [u (ad/d- 1 x)
        v (ad/d- y (ad/d* x x))]
    (ad/d+ (ad/d* u u) (ad/d* 100 v v))))

(defn rosenbrock-gradient
  "Gradiente analítico de Rosenbrock en [x y], para validar la AD:
  [−2(1 − x) − 400x(y − x²), 200(y − x²)]."
  [[x y]]
  (let [v (- y (* x x))]
    [(- (* -2.0 (- 1.0 x)) (* 400.0 x v))
     (* 200.0 v)]))

(defn least-squares-objective
  "Función objetivo de mínimos cuadrados ‖X·w − y‖², lista para `grad` o
  `gradient-descent`: devuelve una función de un vector de duales w.
  X es una matriz (vector de filas) e y un vector de números."
  [x-mat y]
  (fn [w]
    (apply ad/d+
           (map (fn [fila yi]
                  (let [r (ad/d- (apply ad/d+ (map ad/d* fila w)) yi)]
                    (ad/d* r r)))
                x-mat y))))

(defn least-squares-gradient
  "Gradiente analítico de ‖X·w − y‖²: 2·Xᵀ(X·w − y)."
  [x-mat y w]
  (la/scale 2.0 (la/mat-vec (la/transpose x-mat)
                            (la/v- (la/mat-vec x-mat w) y))))

(defn normal-equations
  "Solución de mínimos cuadrados por las ecuaciones normales:
  (XᵀX)·w = Xᵀy, resuelto con `fpnum.linear/solve`.

  Es la solución exacta cuando X tiene columnas linealmente independientes.
  Sirve de referencia para el descenso de gradiente. Lanza error si XᵀX es
  singular."
  [x-mat y]
  (let [xt (la/transpose x-mat)]
    (la/solve (la/mat-mat xt x-mat) (la/mat-vec xt y))))

(defn least-squares-svd
  "Solución de mínimos cuadrados por pseudoinversa vía SVD:

      w = X⁺·y = Σᵢ (uᵢᵀy / σᵢ)·vᵢ    (solo los σᵢ por encima de la tolerancia)

  sin formar X⁺ explícitamente. A diferencia de `normal-equations`:
  * no eleva al cuadrado el número de condición (el error crece como κ·ε y
    no como κ²·ε);
  * no falla si X tiene rango deficiente: devuelve la solución de mínima norma
    ‖w‖ entre todas las de mínimos cuadrados.

  `rcond` (opcional) fija la tolerancia rcond·σmax; por defecto max(m, n)·ε."
  ([x-mat y] (least-squares-svd x-mat y nil))
  ([x-mat y rcond]
   (let [r (svd/svd x-mat)
         inv (svd/inverse-singular-values r rcond)
         uty (la/mat-vec (la/transpose (:u r)) y)]
     (la/mat-vec (:v r) (mapv * inv uty)))))

(defn safe-step
  "Paso fijo seguro para mínimos cuadrados con matriz X.

  La función ‖Xw − y‖² tiene Hessiano 2·XᵀX, y el descenso con paso fijo α
  converge si α < 2/L, con L = 2·λmax(XᵀX) la constante de Lipschitz del
  gradiente. Como λmax(XᵀX) ≤ ‖X‖²_F (la traza de XᵀX es ‖X‖²_F y todos sus
  autovalores son no negativos), basta cualquier

      α < 1/‖X‖²_F

  Se devuelve α = 0.9/‖X‖²_F. Es una cota conservadora: el paso óptimo
  teórico es mayor, pero este se calcula sin conocer λmax."
  [x-mat]
  (let [f (la/norm-fro x-mat)]
    (/ 0.9 (* f f))))

;; ---------------------------------------------------------------------------
;; Datos ruidosos reproducibles
;; ---------------------------------------------------------------------------

(defn noisy-regression-data
  "Genera un problema de regresión lineal con ruido gaussiano, reproducible.

  Opciones:
    :n       número de observaciones (filas de X)
    :w-true  vector de coeficientes verdaderos (define las columnas de X)
    :sigma   desviación estándar del ruido
    :seed    semilla de java.util.Random

  Las entradas de X son gaussianas estándar y y = X·w-true + sigma·ruido.
  Devuelve {:x X :y y :w-true w-true}. La misma semilla produce siempre los
  mismos datos (el algoritmo de java.util.Random está especificado)."
  [{:keys [n w-true sigma seed]}]
  (let [rnd (java.util.Random. (long seed))
        p (count w-true)
        x-mat (vec (repeatedly n #(vec (repeatedly p (fn [] (.nextGaussian rnd))))))
        y (mapv (fn [fila]
                  (+ (la/dot fila w-true) (* sigma (.nextGaussian rnd))))
                x-mat)]
    {:x x-mat :y y :w-true (vec w-true)}))

;; ---------------------------------------------------------------------------
;; Descenso de gradiente
;; ---------------------------------------------------------------------------

(defn- f-value
  "Valor de f en `x` (se evalúa con duales constantes y se toma la parte a)."
  [f x]
  (:a (ad/lift (f (mapv ad/lift x)))))

(defn- armijo-step
  "Busca α por backtracking. Devuelve {:alpha α :x x'} con el primer α que
  cumple la condición de Armijo, o nil si α cae por debajo de `min-alpha`."
  [f x fx g g2 {:keys [alpha0 beta c min-alpha]}]
  (loop [alpha alpha0]
    (when (>= alpha min-alpha)
      (let [x' (la/v- x (la/scale alpha g))
            fx' (f-value f x')]
        (if (<= fx' (- fx (* c alpha g2)))
          {:alpha alpha :x x'}
          (recur (* beta alpha)))))))

(defn gradient-descent
  "Minimiza `f` por descenso de gradiente desde `x0`.

  `f` recibe un vector de duales y devuelve un dual (o número).

  Opciones (todas opcionales):
    :alpha         paso fijo (por defecto 0.01)
    :max-iter      máximo de iteraciones (por defecto 1000)
    :tol           se detiene si ‖∇f(x)‖ ≤ tol (por defecto 1e-8)
    :grad-fn       función (f x) → vector con el gradiente de f en x. Por
                   defecto fpnum.ad.forward/grad (n pasadas hacia adelante);
                   se puede pasar fpnum.ad.reverse/grad (una pasada hacia
                   adelante y una hacia atrás). `f` debe estar escrita con las
                   operaciones genéricas (fpnum.ad.ops o fpnum.ad.forward) para
                   que sirva con cualquiera de los dos.
    :stop-when     función de una fila {:iter :x :f :grad-norm :alpha}; si
                   devuelve verdadero, se detiene con :reason :stopped. Sirve
                   para parar con un criterio propio, por ejemplo la distancia
                   a una solución conocida: ‖∇f‖ pequeño no garantiza estar
                   cerca del mínimo si f está mal condicionada.
    :backtracking? si es true, elige α con la condición de Armijo en vez de
                   usar el paso fijo
    :alpha0        α inicial del backtracking (por defecto 1.0)
    :beta          factor de reducción de α (por defecto 0.5)
    :c             constante de Armijo (por defecto 1e-4)
    :min-alpha     α mínimo antes de rendirse (por defecto 1e-16)

  Devuelve un mapa:
    :x           punto final
    :f           f en el punto final
    :iterations  número de pasos dados
    :converged?  true si terminó por ‖∇f‖ ≤ tol (no si terminó por :stop-when)
    :reason      :converged, :stopped (lo pidió :stop-when), :max-iter,
                 :diverged (valores no finitos) o :line-search-failed (el
                 backtracking no encontró α)
    :trajectory  vector de filas {:iter :x :f :grad-norm :alpha}, una por
                 punto visitado. :alpha es el paso que se usó para salir de
                 ese punto; en la última fila es nil porque ya no se da
                 ningún paso."
  ([f x0] (gradient-descent f x0 {}))
  ([f x0 opts]
   (let [{:keys [alpha max-iter tol backtracking? stop-when grad-fn]
          :or {alpha 0.01 max-iter 1000 tol 1e-8 backtracking? false grad-fn ad/grad}} opts
         bt-opts (merge {:alpha0 1.0 :beta 0.5 :c 1e-4 :min-alpha 1e-16} opts)]
     (loop [x (mapv double x0)
            i 0
            traj []]
       (let [fx (f-value f x)
             g (grad-fn f x)
             g2 (la/dot g g)
             gn (Math/sqrt g2)
             fila (fn [a] {:iter i :x x :f fx :grad-norm gn :alpha a})
             terminar (fn [reason]
                        {:x x :f fx :iterations i
                         :converged? (= reason :converged)
                         :reason reason
                         :trajectory (conj traj (fila nil))})]
         (cond
           (not (and (Double/isFinite fx) (Double/isFinite gn))) (terminar :diverged)
           (<= gn tol) (terminar :converged)
           (and stop-when (stop-when (fila nil))) (terminar :stopped)
           (>= i max-iter) (terminar :max-iter)
           :else
           (let [paso (if backtracking?
                        (armijo-step f x fx g g2 bt-opts)
                        {:alpha alpha :x (la/v- x (la/scale alpha g))})]
             (if (nil? paso)
               (terminar :line-search-failed)
               (recur (:x paso) (inc i) (conj traj (fila (:alpha paso))))))))))))
