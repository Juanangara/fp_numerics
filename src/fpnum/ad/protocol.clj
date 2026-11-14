(ns fpnum.ad.protocol
  "Protocolo para que las operaciones de AD (d+, d*, dsin...) sirvan con más
  de un tipo de valor.

  Hay dos modos de diferenciación automática:

  * forward (fpnum.ad.forward): el valor es un `Dual` (a + bε).
  * reverse (fpnum.ad.reverse): el valor es un `Node` del grafo de cálculo.

  Para que una misma función de usuario, escrita con `d+`, `d*`..., sirva en
  los dos modos, las operaciones públicas (que viven en fpnum.ad.forward y se
  reexportan en fpnum.ad.ops) miran el tipo de sus argumentos:

  * `Dual` o número: camino directo de siempre (el más rápido).
  * cualquier otro valor que implemente este protocolo (un `Node`): delegan en
    el protocolo.

  Los métodos binarios (`add`, `sub`, `mul`, `div`, `pow`) reciben como primer
  argumento un valor del protocolo y como segundo, ese mismo tipo o un número.
  Si el primer operando de la operación pública es un número y el segundo un
  valor del protocolo, se convierte antes con `const-like`.

  `Dual` NO implementa este protocolo a propósito: así el camino de forward
  sigue siendo una comprobación de tipo y la aritmética de siempre, sin
  despacho.")

(defprotocol ADValue
  "Operaciones que debe ofrecer un tipo de valor de AD (un nodo de reverse, por
  ejemplo) para que las operaciones genéricas de fpnum.ad.forward (`d+`, `d*`,
  `dsin`...) deleguen en él. Ver el docstring del namespace."
  (const-like [this n]
    "Constante (derivada 0) del mismo tipo que `this` con el valor del número `n`.")
  (add [x y] "x + y")
  (sub [x y] "x − y")
  (mul [x y] "x · y")
  (div [x y] "x / y")
  (pow [x y] "x elevado a y")
  (sin* [x] "sin x")
  (cos* [x] "cos x")
  (exp* [x] "eˣ")
  (log* [x] "ln x"))
