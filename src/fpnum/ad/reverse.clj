(ns fpnum.ad.reverse
  "Diferenciación automática en modo reverse.

  Teoría
  ------
  El modo forward (fpnum.ad.forward) obtiene una derivada parcial por pasada:
  un gradiente de n variables cuesta n evaluaciones de f. El modo reverse
  obtiene el gradiente completo con una pasada hacia adelante y una pasada
  hacia atrás, por unas 3 a 5 veces el costo de f, sin depender de n.

  1. Pasada hacia adelante: se evalúa f sobre `Node`s. Cada operación
     (d+, d*, dsin...) calcula su valor y devuelve un nodo nuevo que guarda la
     operación y sus argumentos. El resultado es un grafo acíclico de cálculo.
  2. Pasada hacia atrás: el adjunto de un nodo es ∂salida/∂nodo. El de la
     salida es 1. Por la regla de la cadena, cada nodo reparte su adjunto entre
     sus argumentos multiplicándolo por la derivada local de su operación, y
     el adjunto de un nodo usado en varios sitios es la SUMA de lo que recibe de
     todos. Los nodos se procesan en orden topológico inverso (padres antes que
     hijos), así que cuando le toca a un nodo ya ha recibido todo. Ese orden se
     obtiene gratis de los ids (ver abajo): basta procesar en orden decreciente
     de id.

  Identidad de los nodos: contador (enfoque 1b)
  ---------------------------------------------
  Para sumar bien los adjuntos hay que saber cuándo dos referencias son el mismo
  nodo. Cada nodo lleva un `id` entero, tomado de un contador global creciente
  (un `atom`) al crearlo. Es la única mutación del namespace, y es benigna: los
  ids solo SIRVEN PARA DISTINGUIR nodos, y el gradiente no depende de qué
  números se asignen. Igualdad y hash de un `long` son O(1), sin recursión.

  Por qué no una identidad estructural (id = [op id-arg ...], que sería pura y
  fusionaría subexpresiones iguales). Se implementó y se midió (docs/decisiones.md,
  docs/data/ad-reverse-sharing-1a.csv): la igualdad de vectores anidados es
  recursiva y sin memoria, y las tablas hash persistentes la invocan al buscar
  una clave aunque la clave encontrada sea otra distinta. Con ella, una suma de
  ~4000 términos o más desbordaba la pila (una suma de n términos es una cadena
  de n nodos, porque d+ pliega a la izquierda), y la misma expresión muy
  compartida construida dos veces era exponencial en la profundidad.

  INVARIANTE (del que depende la pasada hacia atrás): el id de un nodo es mayor
  que el de cada uno de sus argumentos, porque el contador crece y un nodo se
  crea después que sus argumentos. Se cumple también entre hilos (`swap!` es
  atómico). Una prueba lo verifica (fpnum.ad.reverse-sharing-test). Ver
  `hacia-atras` y docs/decisiones.md (D4).

  Consecuencia de 1b: dos subexpresiones iguales construidas por separado son
  nodos distintos (no se fusionan). Es más trabajo repetido, nunca un resultado
  distinto, y el costo sigue siendo lineal en el número de nodos creados.
  Uso
  ---
  Las funciones se escriben con las operaciones genéricas (fpnum.ad.ops o
  fpnum.ad.forward) y sirven igual en los dos modos:

      (require '[fpnum.ad.ops :as ad] '[fpnum.ad.reverse :as rev])
      (rev/grad (fn [[x y]] (ad/d+ (ad/d* x x y) (ad/dsin y))) [1.0 2.0])

  `f` recibe un vector de nodos y devuelve un nodo o un número. Si devuelve un
  número (o un dual constante), el gradiente es cero. Los nodos se comparan por
  identidad (no tienen igualdad estructural); para mostrarlos hay
  `print-method`, que imprime solo la operación y el valor y no recorre el
  grafo.

  Regla heurística: un Dual con derivada 0
  ----------------------------------------
  Si una operación mezcla un nodo con un `Dual` de derivada exactamente 0, el
  Dual se trata como una constante (ver fpnum.ad.forward/valor-de-dual-constante
  y docs/decisiones.md). Es una REGLA HEURÍSTICA: dentro de reverse no hay
  semillas de forward, así que un Dual solo puede venir de operar con números
  (por ejemplo (d* 2.0 3.0) o (d+)); pero un Dual con derivada cero también
  podría ser una variable de forward cuya derivada se anula en ese punto, y en
  ese caso la mezcla se aceptaría en silencio como constante. Mezclar modos
  (diferenciación anidada) no está soportado.

  Costo
  -----
  Memoria O(operaciones), porque se guarda todo el grafo. La pasada hacia atrás
  hace O(log n) por inserción en la cola de prioridad, pero en la práctica se
  mide casi lineal (docs/data/ad-forward-vs-reverse-notas.md)."
  (:require [fpnum.ad.protocol :as proto]))

(set! *warn-on-reflection* true)

(declare make-node const-node ->node)

(defn- fdiv
  "División de doubles con semántica IEEE (x/0 da ±Infinity o NaN). Ver la
  nota equivalente en fpnum.ad.forward: con números en caja Clojure lanza
  ArithmeticException al dividir por cero."
  ^double [^double x ^double y]
  (/ x y))

;; ---------------------------------------------------------------------------
;; Nodo del grafo de cálculo
;; ---------------------------------------------------------------------------

(def ^:private contador
  "Contador global de ids de nodos (ver la sección sobre identidad en el
  docstring del namespace)."
  (atom 0))

(defn- siguiente-id
  "Id nuevo, mayor que todos los anteriores. Un nodo se crea siempre después de
  sus argumentos, así que su id es mayor que el de cada uno de ellos."
  ^long []
  (long (swap! contador inc)))

(deftype Node [^long id op args ^double value]
  proto/ADValue
  (const-like [_ n] (const-node n))
  (add [x y] (let [^Node y (->node y)] (make-node :add [x y] (+ (.-value x) (.-value y)))))
  (sub [x y] (let [^Node y (->node y)] (make-node :sub [x y] (- (.-value x) (.-value y)))))
  (mul [x y] (let [^Node y (->node y)] (make-node :mul [x y] (* (.-value x) (.-value y)))))
  (div [x y] (let [^Node y (->node y)] (make-node :div [x y] (fdiv (.-value x) (.-value y)))))
  (pow [x y] (let [^Node y (->node y)] (make-node :pow [x y] (Math/pow (.-value x) (.-value y)))))
  (sin* [x] (make-node :sin [x] (Math/sin (.-value x))))
  (cos* [x] (make-node :cos [x] (Math/cos (.-value x))))
  (exp* [x] (make-node :exp [x] (Math/exp (.-value x))))
  (log* [x] (make-node :log [x] (Math/log (.-value x)))))

(defn node?
  "Verdadero si `x` es un nodo del grafo de reverse."
  [x]
  (instance? Node x))

(defn value
  "Valor (pasada hacia adelante) del nodo."
  ^double [^Node n]
  (.-value n))

(defmethod print-method Node
  [^Node n ^java.io.Writer w]
  (.write w (str "#Node[" (name (.-op n)) " " (.-value n) "]")))

(defn- make-node
  "Crea un nodo de operación con un id nuevo (ver `siguiente-id`)."
  ^Node [op args value]
  (Node. (siguiente-id) op args (double value)))

(defn- const-node
  "Nodo constante (no depende de ninguna entrada) con el valor del número `n`."
  ^Node [n]
  (Node. (siguiente-id) :const [] (double n)))

(defn- var-node
  "Nodo de una entrada, con el valor `x`. El índice `_` no se guarda: el orden
  de las entradas lo conserva el vector de hojas de `value-and-grad`."
  ^Node [_ x]
  (Node. (siguiente-id) :var [] (double x)))

(defn- ->node
  "Convierte un número en un nodo constante; los nodos pasan tal cual."
  ^Node [y]
  (cond
    (instance? Node y) y
    (number? y) (const-node y)
    :else (throw (ex-info "Operando no válido para AD en modo reverse: se esperaba un nodo o un número"
                          {:operando-clase (some-> y class .getName)}))))

;; ---------------------------------------------------------------------------
;; Pasada hacia atrás
;; ---------------------------------------------------------------------------

(defn- contribuciones
  "Cuánto aporta el adjunto `adj` de `nodo` al adjunto de cada uno de sus
  argumentos (vector en el mismo orden que :args): adj por la derivada local
  de la operación."
  [^Node nodo ^double adj]
  (let [args (.-args nodo)]
    (case (.-op nodo)
      :add [adj adj]
      :sub [adj (- adj)]
      :mul (let [[^Node a ^Node b] args]
             [(* adj (.-value b)) (* adj (.-value a))])
      :div (let [[^Node a ^Node b] args
                 av (.-value a)
                 bv (.-value b)]
             [(fdiv adj bv) (- (fdiv (* adj av) (* bv bv)))])
      :pow (let [[^Node a ^Node b] args
                 av (.-value a)
                 n (.-value b)]
             ;; d(aⁿ)/da = n·aⁿ⁻¹ (0 si n = 0, para que x⁰ no dé NaN en x = 0);
             ;; d(aⁿ)/dn = aⁿ·ln a, que se omite si el exponente es una constante.
             [(if (zero? n) 0.0 (* adj n (Math/pow av (- n 1.0))))
              (if (= :const (.-op b)) 0.0 (* adj (.-value nodo) (Math/log av)))])
      :sin [(* adj (Math/cos (.-value ^Node (first args))))]
      :cos [(- (* adj (Math/sin (.-value ^Node (first args)))))]
      :exp [(* adj (.-value nodo))]
      :log [(fdiv adj (.-value ^Node (first args)))]
      ;; hojas
      [])))

(defn- hacia-atras
  "Mapa id → adjunto de las HOJAS (entradas y constantes) del grafo que sale de
  `raiz`, cuyo adjunto es 1.

  INVARIANTE en el que se apoya: el id de un nodo es MAYOR que el de cada uno
  de sus argumentos (los ids salen de un contador creciente y un nodo se crea
  siempre después que sus argumentos; ver `siguiente-id`). Por tanto, procesar
  los nodos en orden DECRECIENTE de id es un orden topológico inverso: cuando
  le toca a un nodo, todos sus padres (que tienen un id mayor y están
  alcanzados desde la raíz) ya se procesaron, así que su adjunto está completo.

  La estructura de trabajo es una cola de prioridad (un `sorted-map` por id
  decreciente) de id → [nodo adjunto-acumulado]. Se saca siempre el id mayor, se
  reparte su adjunto entre sus argumentos y estos se insertan en la cola o, si
  ya estaban (nodo compartido), suman la contribución a su adjunto. No hace
  falta un recorrido previo, ni un conjunto de visitados, ni una pila explícita.

  Si se cambiara la estrategia de ids (por ejemplo a ids estructurales o
  asignados de forma perezosa), este invariante dejaría de cumplirse y habría
  que volver a un orden topológico explícito."
  [^Node raiz]
  (loop [cola (sorted-map-by > (.-id raiz) [raiz 1.0])
         hojas {}]
    (if (empty? cola)
      hojas
      (let [[id [^Node nodo adj]] (first cola)
            cola (dissoc cola id)
            args (.-args nodo)]
        (if (empty? args)
          (recur cola (assoc hojas id adj))
          (recur (reduce (fn [q [^Node a c]]
                           (update q (.-id a)
                                   (fn [e] (if e
                                             [(first e) (+ (double (second e)) (double c))]
                                             [a c]))))
                         cola
                         (map vector args (contribuciones nodo (double adj))))
                 hojas))))))
;; ---------------------------------------------------------------------------
;; API
;; ---------------------------------------------------------------------------

(defn- resultado-numerico
  "Valor numérico de lo que devuelve f: un nodo, un número o un dual constante."
  [salida]
  (cond
    (node? salida) (value salida)
    (number? salida) (double salida)
    (and (map? salida) (contains? salida :a)) (double (:a salida))
    :else (throw (ex-info "f debe devolver un nodo o un número"
                          {:clase (some-> salida class .getName)}))))

(defn value-and-grad
  "Valor y gradiente de `f` en el punto `xs` (vector de números), en una sola
  pasada hacia adelante y una hacia atrás.

  `f` recibe un vector de nodos (uno por componente de xs) y devuelve un nodo o
  un número. Devuelve {:value v :grad [∂f/∂x₀ ... ∂f/∂xₙ₋₁]}. Si f no depende de
  los nodos (devuelve un número), el gradiente es cero."
  [f xs]
  (let [hojas (vec (map-indexed var-node xs))
        salida (f hojas)]
    {:value (resultado-numerico salida)
     :grad (if (node? salida)
             (let [adj (hacia-atras salida)]
               (mapv (fn [^Node h] (double (get adj (.-id h) 0.0))) hojas))
             (vec (repeat (count xs) 0.0)))}))

(defn grad
  "Gradiente de `f` en el punto `xs` (vector de números) por reverse mode: una
  pasada hacia adelante y una hacia atrás, sin importar cuántas variables
  haya. Misma convención que `fpnum.ad.forward/grad`: `f` recibe un vector."
  [f xs]
  (:grad (value-and-grad f xs)))

(defn diff
  "Derivada de `f` (función de una variable) en el punto `x`. Misma convención
  que `fpnum.ad.forward/diff`: `f` recibe un solo argumento."
  [f x]
  (first (grad (fn [[v]] (f v)) [x])))
