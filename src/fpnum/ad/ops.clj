(ns fpnum.ad.ops
  "Operaciones de diferenciación automática independientes del modo.

  Es una fachada fina: reexporta las operaciones de fpnum.ad.forward con un
  nombre neutro, para que el código que sirve igual en forward y en reverse no
  tenga que decir `forward/d*`. Las funciones escritas con estas operaciones
  se pueden pasar tanto a `fpnum.ad.forward/grad` como a
  `fpnum.ad.reverse/grad`:

      (require '[fpnum.ad.ops :as ad]
               '[fpnum.ad.forward :as fwd]
               '[fpnum.ad.reverse :as rev])

      (def f (fn [[x y]] (ad/d+ (ad/d* x x y) (ad/dsin y))))
      (fwd/grad f [1.0 2.0])    ; n pasadas hacia adelante
      (rev/grad f [1.0 2.0])    ; una pasada hacia adelante y una hacia atrás

  Cada operación acepta números, duales (forward) y nodos (reverse); ver
  fpnum.ad.protocol. Mezclar duales y nodos en una misma operación es un
  error."
  (:require [fpnum.ad.forward :as fwd]))

(def ^{:doc "Suma. Ver fpnum.ad.forward/d+."} d+ fwd/d+)
(def ^{:doc "Resta (o negación con un argumento). Ver fpnum.ad.forward/d-."} d- fwd/d-)
(def ^{:doc "Producto. Ver fpnum.ad.forward/d*."} d* fwd/d*)
(def ^{:doc "Cociente (o inverso con un argumento). Ver fpnum.ad.forward/d-div."} d-div fwd/d-div)
(def ^{:doc "Seno."} dsin fwd/dsin)
(def ^{:doc "Coseno."} dcos fwd/dcos)
(def ^{:doc "Exponencial."} dexp fwd/dexp)
(def ^{:doc "Logaritmo natural."} dlog fwd/dlog)
(def ^{:doc "Potencia x^y (base y exponente variables o constantes)."} dpow fwd/dpow)
