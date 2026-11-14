(ns fpnum.ad.syntax
  "Sintaxis natural para funciones diferenciables: la macro `ad-fn`.

  Las funciones que se pasan a `fpnum.ad.forward/grad` o `fpnum.ad.reverse/grad`
  deben usar las operaciones de AD (`d+`, `d*`, `dsin`...) y no las de
  clojure.core, que no saben operar con duales ni con nodos. `ad-fn` permite
  escribirlas como matemática normal y reescribe el cuerpo:

      (ad-fn [[x y]] (+ (* x x y) (Math/sin y)))
      ;; equivale a
      (fn [[x y]] (fpnum.ad.ops/d+ (fpnum.ad.ops/d* x x y) (fpnum.ad.ops/dsin y)))

  Como las operaciones de AD son genéricas, la función resultante sirve igual en
  forward que en reverse.

  Cómo funciona (incremento 1)
  ----------------------------
  1. Expande todas las macros del cuerpo con `clojure.walk/macroexpand-all`, de
     modo que `->`, `->>`, `let`, `for`, `when`... quedan en formas básicas.
  2. Recorre el resultado y sustituye símbolos de una tabla fija:
     * `+  -  *  /` (y sus formas `clojure.core/+` ...) por `d+  d-  d*  d-div`.
       Se sustituyen también cuando aparecen como VALOR y no solo como cabeza
       de una llamada: `(reduce + ...)`, `(apply * ...)`.
     * `Math/sin  Math/cos  Math/exp  Math/log  Math/pow  Math/sqrt` (y con
       `java.lang.Math/`) por `dsin  dcos  dexp  dlog  dpow` y una raíz cuadrada
       `x^0.5`. También la forma `(. Math sin x)`.
  3. Devuelve `(fn params cuerpo-reescrito)`. Los símbolos generados están
     calificados (fpnum.ad.ops/d+...), así que no dependen de los alias del
     namespace del usuario.

  LÍMITES (este incremento es deliberadamente básico)
  ---------------------------------------------------
  1. SOLO REESCRIBE LO QUE VE DENTRO DEL CUERPO de la macro (incluidas las
     funciones anónimas `#(...)` y `(fn ...)` escritas dentro). NO reescribe las
     funciones definidas fuera, aunque se llamen desde dentro: una función
     externa que use `+` o `*` de clojure.core con valores de AD fallará con
     ClassCastException. Hay que escribirla con `d+`/`d*` o con `ad-fn`.
  2. TODA la aritmética `+ - * /` del cuerpo se reescribe, también la de índices
     y contadores: `(nth xs (+ i 1))` falla, porque `(d+ i 1)` devuelve un Dual.
     Por eso `inc` y `dec` NO se reescriben (suelen ser contadores); para sumar
     1 a un valor diferenciable, escribir `(+ x 1)`.
  3. NO hay manejo de sombreado: se sustituyen todos los símbolos de la tabla,
     aunque sean locales del usuario. Si el cuerpo define un local llamado `+`
     (por ejemplo `(let [+ ...] ...)`), el código resultante no compila.
     Tampoco se comprueba que el `+` sea el de clojure.core: si el namespace lo
     excluyó con `:refer-clojure :exclude` y define el suyo, se reescribe igual.
  4. Las formas citadas (`quote`) también se reescriben, y sus macros se
     expanden: `'(+ 1 2)` dentro del cuerpo se convierte en `'(fpnum.ad.ops/d+ 1 2)`.
  5. Solo las funciones de la tabla. Las comparaciones (`<`, `=`), los
     condicionales sobre valores diferenciables y otras funciones numéricas
     (`min`, `max`, `abs`, `Math/abs`...) no se reescriben y fallan si reciben
     valores de AD. (Con números normales funcionan igual que siempre.)
  6. La función devuelta es anónima (`fn` sin nombre) y la macro no admite la
     forma `(ad-fn nombre [params] ...)`.
  7. LLAMADA CON NÚMEROS NORMALES: la función devuelta NO se comporta como una
     función numérica corriente. Las operaciones de AD con solo números
     devuelven un Dual constante ((d+ 1 2) es Dual(3, 0), como en el forward original), así
     que `((ad-fn [x] (* x x)) 3.0)` devuelve `Dual(9, 0)` y no `9.0`. Para el
     valor numérico, `(:a resultado)`. Con `fwd/grad`, `rev/grad` y
     `gradient-descent` no hay que hacer nada: ya lo manejan."
  (:require [clojure.walk :as walk]
            [fpnum.ad.ops :as ad]))

(defn ad-sqrt
  "Raíz cuadrada como x^0.5 (a la que se reescribe `Math/sqrt`). En x = 0 la
  derivada es infinita, como la de la raíz."
  [x]
  (ad/dpow x 0.5))

(def ^:private operadores
  {'+ `ad/d+, 'clojure.core/+ `ad/d+
   '- `ad/d-, 'clojure.core/- `ad/d-
   '* `ad/d*, 'clojure.core/* `ad/d*
   '/ `ad/d-div, 'clojure.core// `ad/d-div})

(def ^:private clases-math #{"Math" "java.lang.Math"})

(def ^:private funciones-math
  {"sin" `ad/dsin, "cos" `ad/dcos, "exp" `ad/dexp
   "log" `ad/dlog, "pow" `ad/dpow, "sqrt" `ad-sqrt})

(defn- reescribir-forma
  "Reescribe UNA forma (el recorrido en postorden ya reescribió las interiores)."
  [f]
  (cond
    ;; (. Math sin x)
    (and (seq? f) (= '. (first f))
         (clases-math (str (second f)))
         (funciones-math (str (nth f 2 nil))))
    (cons (funciones-math (str (nth f 2))) (drop 3 f))

    (symbol? f)
    (or (operadores f)
        ;; Math/sin como cabeza de una llamada: (Math/sin x)
        (when-let [clase (namespace f)]
          (when (clases-math clase) (funciones-math (name f))))
        f)

    :else f))

(defn ad-rewrite
  "Expande las macros de `form` y sustituye los símbolos de la tabla de la macro
  `ad-fn`. Devuelve la forma reescrita (sin evaluarla). Es la parte de la macro
  que se puede probar por separado; ver los límites en el docstring del
  namespace."
  [form]
  (walk/postwalk reescribir-forma (walk/macroexpand-all form)))

(defmacro ad-fn
  "Como `fn`, pero reescribe el cuerpo para que `+ - * /` y `Math/sin`, `cos`,
  `exp`, `log`, `pow` y `sqrt` usen las operaciones de AD. La función resultante
  sirve en forward y en reverse. Para varias variables, recibir un vector:

      (rev/grad (ad-fn [[x y]] (+ (* x x y) (Math/sin y))) [1.0 2.0])

  LÍMITES: solo reescribe lo que ve dentro de su cuerpo (no las funciones
  definidas fuera), reescribe también la aritmética de índices, no maneja el
  sombreado de nombres, no reescribe `inc`, `dec` ni las comparaciones, y
  llamada con números normales devuelve un Dual constante, no un número. Ver el
  docstring de fpnum.ad.syntax para la lista completa."
  [params & body]
  `(fn ~params ~@(map ad-rewrite body)))
