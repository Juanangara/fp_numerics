# Decisiones de diseño

Registro de las decisiones de diseño no obvias, con la evidencia medida y las
alternativas que se descartaron. Los datos están en `docs/data/` y los scripts que
los generan en `scripts/`.

| # | Decisión | Estado |
|---|---|---|
| D1 | Identidad de los nodos en reverse mode: contador (1b), tras el resultado negativo de la identidad estructural (1a) | Aplicada |
| D2 | Un `Dual` con derivada exactamente 0 se trata como constante al mezclarlo con un nodo (regla heurística) | Aplicada, con límite documentado |
| D3 | Despacho de las operaciones de forward por clase y no por interfaz | Aplicada |
| D4 | Pasada hacia atrás en orden decreciente de id (`sorted-map`) | Aplicada |
| D5 | Macro `ad-fn`: se descarta el segundo incremento (manejo de sombreado) | **Descartada** |

---

## D1. Identidad de los nodos en reverse mode

### Contexto
La pasada hacia atrás de reverse mode acumula adjuntos por nodo, y el adjunto de un
nodo usado en varios sitios es la **suma** de lo que recibe de todos. Para sumar bien
hay que saber cuándo dos referencias son *el mismo* nodo, es decir, hace falta una
**identidad de nodo** que sirva de clave en un mapa.

### Alternativas
- **1a — identidad estructural.** `id = [op id-arg1 id-arg2 ...]` (vectores anidados).
  Pura, sin estado, y de regalo fusiona subexpresiones iguales (eliminación de
  subexpresiones comunes).
- **1b — contador.** `id` = entero tomado de un `atom` creciente al crear el nodo.
- **1a′ — huella de 128 bits.** `id` = hash de 128 bits de `[op huellas-de-los-padres]`.
  Pura, con igualdad en O(1).
- **1c — identidad del objeto.** Sin ids: un mapa mutable local (`IdentityHashMap`) en
  la pasada hacia atrás. La mutación queda encerrada en esa función, como un
  `transient`.

### Resultado negativo de 1a (medido)
Se implementó 1a y se midió con `scripts/measure_reverse_sharing.clj`
(`docs/data/ad-reverse-sharing-1a.csv`). Tres casos:

- **L1**: cadena lineal `acc ← 0.999·acc + x` de n términos, construida una vez.
- **S1**: cadena de cuadrados `((x²)²)²…` de profundidad d, construida **una** vez.
- **S2**: la misma cadena de cuadrados construida **dos** veces y sumada.

| Caso | Tamaño | Con 1a |
|---|---|---|
| L1 | 500 / 1 000 / 2 000 / 3 000 | 0.041 / 0.049 / 0.121 / 0.298 s (costo superlineal) |
| L1 | **5 000 / 10 000 / 20 000 / 50 000** | **`StackOverflowError`** |
| S1 | d = 1 000 | 0.022 s, correcto (lineal) |
| S2 | d = 16 / 18 / 20 / 22 | 0.0050 / 0.0195 / 0.0786 / 0.3188 s |
| S2 | d = 24 / 26 / **28** | 1.26 / 5.07 / **20.12 s** |
| S2 | d = 30 en adelante | no medido: cada +2 de profundidad cuadruplica el tiempo, es decir, ×2 por nivel. Extrapolando, d = 30 ≈ 80 s y d = 40 ≈ 23 horas |

Una consecuencia práctica que no se había anticipado: como `d+` pliega de izquierda a
derecha, `(apply d+ xs)` sobre n términos es una **cadena de n sumas anidadas**. Por
tanto, con 1a **cualquier suma de ~4 000 términos o más fallaba**, por ejemplo el
objetivo de mínimos cuadrados con 5 000 observaciones.

### Por qué ocurre
La igualdad de vectores de Clojure es recursiva y no tiene memoria. Las tablas hash
persistentes la invocan al buscar una clave **aunque la clave encontrada sea otra
distinta**, porque solo comprueban los bits de la ruta del árbol y no el hash
completo. En una cadena larga, dos ids distintos comparten casi todo el prefijo
anidado, así que la comparación recorre toda la profundidad: de ahí el
`StackOverflowError` (L1) y el costo superlineal previo.

Para dos ids **iguales pero distintos como objetos** (S2), `[:mul c c]` frente a
`[:mul c' c']` compara `c` con `c'` **dos veces** por nivel, lo que da 2^d.
Dentro de un mismo grafo no ocurre, porque cada id se compara consigo mismo
(`identical?`, instantáneo): por eso S1 es lineal.

Forzar el hash al crear cada nodo (que se probó) evita el desborde en el *hash*,
pero no el de la *igualdad*.

### Comparación de las alternativas

| | **1b contador** | **1a′ huella 128 bits** | **1c identidad de objeto** |
|---|---|---|---|
| Igualdad y hash | O(1) | O(1) | no hay ids |
| Pureza | una mutación benigna: un `atom` global (como `gensym`) | pura | mutación local en la pasada hacia atrás |
| Subexpresiones iguales construidas dos veces | nodos distintos (más trabajo, nunca un resultado distinto) | **se fusionan** (CSE) | nodos distintos |
| Riesgo de corrección | ninguno: los ids solo distinguen nodos | **colisión**: ~1e-27 con un millón de nodos, pero si ocurre el gradiente sale mal **sin avisar** | ninguno |
| Depuración | ids legibles y secuenciales | huellas opacas | sin ids que inspeccionar |
| Código | ~5 líneas | ~25 líneas de mezcla de bits, que habría que validar | ~15 líneas |
| Efecto lateral útil | los ids crecen en orden de creación, lo que da un **orden topológico gratis** (ver D4) | — | — |

### Decisión: 1b
- **No tiene modos de fallo sutiles.** El único riesgo de 1a′ es un error numérico
  silencioso, el peor tipo de error en un sistema cuyo único trabajo es dar derivadas
  correctas.
- **La mutación no afecta al resultado.** El gradiente no depende de qué números de id
  se asignen; se verifica con una prueba (50 llamadas dan exactamente el mismo
  gradiente) y con otra de concurrencia.
- **Es la más simple** y la más fácil de explicar y depurar.
- **Lo que se pierde** (pureza estricta y fusión de subexpresiones iguales) no se
  necesita: una subexpresión repetida solo cuesta más nodos, y el costo sigue siendo
  lineal.
- 1c habría sido igualmente válida; se prefirió 1b porque permite D4.

### Resultado con 1b
`docs/data/ad-reverse-sharing-1b.csv`: **lineal en los tres casos y sin desbordes**.

| Caso | Con 1b |
|---|---|
| L1, 20 000 términos | 0.106 s (con 1a: `StackOverflowError`) |
| L1, 100 000 términos | 0.649 s |
| S1, d = 1 000 | 0.003 s |
| S2, d = 28 | 0.0002 s (con 1a: 20.12 s) |
| S2, d = 30 | 0.0002 s (con 1a: ~80 s extrapolado) |
| S2, d = 100 / 300 / 1 000 | 0.6 / 1.8 / 5.7 ms (lineal en la profundidad) |

Las pruebas que fallaban con 1a (9 aserciones en 6 pruebas, entre ellas las cadenas
de 20 000 y 50 000 términos y S2 con d = 30) pasan con 1b y están en
`test/fpnum/ad/reverse_sharing_test.clj`.

### Consecuencias
- Hay **un estado global**, el contador de ids. Es seguro entre hilos (`swap!` es
  atómico, y un nodo siempre se crea después que sus argumentos, así que su id es
  mayor que el de ellos).
- Los nodos de dos llamadas a `grad` distintas no deben mezclarse (nunca fue
  soportado).

---

## D2. Un `Dual` con derivada 0 se trata como constante (regla heurística)

### Problema
Las operaciones genéricas devuelven un `Dual` cuando todos sus argumentos son
números: `(d* 2.0 3.0)` da `Dual(6, 0)`, y `(d+)` sin argumentos da `Dual(0, 0)`. Eso
es necesario para no romper el comportamiento original de forward. Pero dentro de
una función que se diferencia con **reverse**, esa aritmética de constantes se
mezclaba después con un nodo (`(d-div x (d+ 1 (d* 1.0 1.0)))`) y fallaba con «no se
pueden mezclar duales y nodos». Lo encontró una propiedad con `test.check`, cuyo
`shrinking` redujo el caso a `x / (1 + 1·1)`.

### La regla
En `fpnum.ad.forward/valor-de-dual-constante`: si un `Dual` se combina con un valor de
otro modo de AD (un nodo de reverse) y su derivada es **exactamente 0**
(`(zero? b)`), se trata como la constante `a`. Si su derivada es distinta de 0, es
una mezcla real de modos y lanza `ex-info`.

### Por qué es una regla heurística, y no un hecho
El argumento a favor es: *dentro de reverse no hay semillas de forward, así que un
`Dual` solo puede venir de operar con números, y por tanto es una constante*. Pero la
regla comprueba una **consecuencia** («derivada 0») y no la **causa** («vino de
operar con números»). Son distintas:

- Un `Dual` con derivada 0 también puede ser **una variable de forward cuya derivada
  se anula en ese punto**. Ejemplo: forward `diff` de `x²` en `x = 0` produce
  `Dual(0, 0)`; si ese valor se mezclara con un nodo de reverse, la regla lo
  aceptaría **en silencio como constante**, y el resultado sería una derivada
  equivocada sin ningún aviso.
- Un `Dual` con derivada distinta de 0 sí se detecta como mezcla real. Detecta, pues,
  solo una parte de los usos incorrectos.

### Por qué se aceptó igualmente
- Mezclar modos (diferenciación anidada, por ejemplo forward sobre reverse para
  calcular Hessianos) **no está soportado**. La situación en la que la regla
  podría equivocarse solo aparece si alguien hace eso a propósito.
- Una alternativa más rigurosa es un tipo `Constant` distinto de `Dual` para los
  resultados de operar solo con números. Eso es más correcto, pero cambia el tipo de
  retorno de `(d+ 1 2)` y rompería las pruebas de forward y su API. No compensa
  para un caso fuera del alcance del proyecto.
- Se dejó documentado en `fpnum.ad.forward`, en `fpnum.ad.reverse` y aquí, y hay una
  prueba que fija el límite: `Dual(5, 0)` se acepta como constante, y `Dual(5, 1)`
  lanza error.

---

## D3. Despacho de las operaciones de forward por clase, no por interfaz

### Contexto
Para que una misma función de usuario sirva en forward y en reverse, `d+`, `d*`,
`dsin`... delegan en el protocolo `ADValue` cuando un argumento no es `Dual` ni
número. El requisito era que forward no se volviera más lento.

### Medición
`scripts/bench_forward_dispatch.clj`: comparación A/B con rondas intercaladas en la
misma JVM, frente a una copia congelada del forward original (`scripts/baseline/`). Criterio:
solo es regresión si los rangos no se solapan **y** la mediana empeora más de un 10 %.

- **Primera versión** (`(instance? ADValue x)`, comprobación contra una *interfaz*):
  `diff` de una función compuesta **1.113× con rangos disjuntos** (57.0–60.8 ms frente
  a 62.3–69.6 ms) y `grad` de Rosenbrock 1.078×. Es una regresión según el criterio.
  *Estos números vienen del registro de la sesión de trabajo: el CSV conserva solo la
  corrida final.*
- **Causa probable**: `instance?` contra una interfaz obliga a recorrer los supertipos
  secundarios de `Dual`, que implementa unas diez interfaces. Contra una clase es una
  comparación de punteros.
- **Versión final** (primero `(or (instance? Dual x) (instance? Number x))`, es decir,
  por clase, y delegación solo en el caso contrario): `docs/data/ad-forward-dispatch.csv`,
  cociente 1.006 a 1.033 en los cuatro casos, con todos los rangos solapados y resultados
  idénticos.
- **Repetición posterior** (mismo script, otro momento): cocientes **1.028, 1.103, 1.032 y
  1.092**, también con los rangos solapados en los cuatro casos (por tanto, sin regresión según
  el criterio), pero con dos casos cerca del 10 %. *Estos números vienen del registro de la
  sesión: el CSV guardado es el de la primera corrida y no se sobrescribió.*

### Qué se puede afirmar
El costo del despacho genérico en forward es **pequeño y variable: entre ~1 % y ~10 % según la
corrida y el caso, siempre con los rangos solapados**. Es mayor en las funciones diminutas (pocas
operaciones por llamada), donde la comprobación de tipo pesa más frente al trabajo útil. No se
puede afirmar una cifra más precisa con este equipo: el rendimiento del portátil varía entre
corridas (ver `docs/data/ad-forward-vs-reverse-notas.md`).

---

## D4. Pasada hacia atrás en orden decreciente de id

### Idea
Con ids de contador, cada nodo tiene un id **mayor** que el de todos sus argumentos,
porque se crea después. Por tanto, procesar los nodos en orden **decreciente de id** es
ya un orden topológico inverso válido (padres antes que hijos). Una cola de prioridad
(`sorted-map` por id) sustituye al recorrido en profundidad con pila explícita, al
conjunto de visitados y al vector de orden: se saca siempre el id mayor, se reparte su
adjunto entre sus argumentos, y estos se insertan en la cola o, si ya estaban (nodo
compartido), suman la contribución.

### Invariante (del que depende la corrección)
> **El id de un nodo es mayor que el de cada uno de sus argumentos.**

- **Por qué se cumple:** el contador (`siguiente-id`) crece, y un nodo se crea siempre
  después que sus argumentos, porque los necesita para calcular su valor. También se
  cumple entre hilos, porque `swap!` es atómico: el id de un argumento se asigna antes de
  que empiece a crearse el nodo que lo usa.
- **Dónde está documentado:** en el docstring del namespace `fpnum.ad.reverse` (sección
  sobre identidad) y en el de `hacia-atras`, y aquí.
- **Cómo se verifica:** `el-id-de-un-nodo-es-mayor-que-el-de-sus-argumentos` y
  `llamadas-concurrentes-no-se-interfieren`, en `test/fpnum/ad/reverse_sharing_test.clj`.
- **Cuándo deja de cumplirse:** si se cambiara la estrategia de ids (por ejemplo a ids
  estructurales o a huellas, 1a/1a′ de D1, o a ids asignados de forma perezosa). Entonces
  habría que volver a un orden topológico explícito. **Es una dependencia entre D1 y D4:**
  cambiar D1 obliga a revisar D4.

### Evaluación previa a aplicarlo
Se midió con una implementación de prueba fuera de `src/`: misma pasada hacia atrás sobre el
mismo grafo, mediana de 9 repeticiones tras calentar.

| Grafo | Anterior (DFS + visitados) | `sorted-map` por ids | Cociente | Diferencia relativa del gradiente |
|---|---|---|---|---|
| Mínimos cuadrados 100×100 (~20 k nodos) | 51.5 ms | 25.8 ms | 0.50 | 2.7e-15 |
| Mínimos cuadrados 100×200 (~40 k nodos) | 104.5 ms | 56.7 ms | 0.54 | 4.5e-15 |
| Cadena lineal de 20 000 | 100.4 ms | 33.7 ms | 0.34 | 0 |
| Suma de 20 000 productos | 106.5 ms | 42.7 ms | 0.40 | 0 |
| S2, profundidad 1 000 | 3.2 ms | 1.1 ms | 0.35 | 0 |

Es de 2 a 3 veces más rápida y da el mismo gradiente (las diferencias de ~1e-15 son del
distinto orden de suma de adjuntos).

### Explicación de la mejora de velocidad: HIPÓTESIS NO VERIFICADA
Se esperaba lo contrario (un `sorted-map` persistente cuesta O(log n) por operación, y el
recorrido en profundidad es O(n)), y salió más rápido. La explicación que se propone es que
**el recorrido anterior asignaba más memoria**: un par `[nodo expandido?]` por nodo en la
pila, un conjunto persistente de visitados, un vector con el orden y una segunda pasada con
un mapa de adjuntos, es decir, tres pasadas frente a una.

**No está verificada.** No se hizo un perfil de la JVM ni se midió la asignación de memoria
ni el tiempo del recolector de basura de cada variante. Lo único medido es el efecto: el
tiempo de la pasada hacia atrás y la igualdad de los gradientes. Podría haber otros factores
(por ejemplo, localidad de caché o el costo de hashear claves `long` en el conjunto de
visitados).

### Lo que NO resuelve
El riesgo de desborde no venía de la pasada hacia atrás (que ya era iterativa, con pila
explícita) sino de la igualdad de ids anidados, y eso ya lo resolvió D1. El cambio es una
**simplificación** (unas 15 líneas en lugar de unas 35) y una mejora de velocidad, no una
corrección de seguridad.

### Resultado después de aplicarlo
`docs/data/ad-forward-vs-reverse*.csv` y `ad-forward-vs-reverse-notas.md`, con muestras
intercaladas y razones pareadas:

- Con n = 200 en mínimos cuadrados, reverse pasó de ~102 ms a ~53–60 ms (unas 1.9 veces más
  rápido).
- **Punto de cruce forward/reverse: de n ≈ 25 a n ≈ 10–11** (razón 0.89 en n = 10 y 1.09 en
  n = 12 en mínimos cuadrados; interpolando, 11.1 y 10.4 en las dos funciones medidas, ver
  `docs/data/ajustes-log-log.csv`). El «≈ 25» anterior era una interpolación entre n = 10 y
  n = 50, no una medida directa.
- Todas las pruebas de reverse (incluidas las 500 expresiones aleatorias contra forward y las
  cadenas de 20 000 y 100 000 términos) pasan sin cambios.

---

## D5. Macro `ad-fn`: se descarta el incremento 2 (manejo de sombreado)

### Qué era
El incremento 1 de `ad-fn` (hecho) reescribe `+ - * /` y `Math/sin cos exp log pow sqrt`
del cuerpo a las operaciones de AD. El incremento 2 habría añadido **manejo de sombreado**:
no reescribir un símbolo si el usuario lo redefine como local (`(let [+ ...] ...)`,
`(fn [*] ...)`), y comprobar con `resolve` que el `+` sea el de `clojure.core` (y no uno
propio, si el namespace lo excluyó con `:refer-clojure :exclude`).

### Decisión
**Se descarta.** El incremento 1 se queda como la versión final de la macro, con sus límites
documentados en el docstring de `fpnum.ad.syntax` y fijados con pruebas
(`test/fpnum/ad/syntax_test.clj`, pruebas `limite-*`).

### Por qué
1. **El fallo es ruidoso, no silencioso.** Sin manejo de sombreado, un local llamado `+` se
   reescribe también en su *definición*, y queda un binder calificado
   (`fpnum.ad.ops/d+`) que el compilador rechaza («Can't let qualified name»). Hay una
   prueba que lo comprueba. El usuario ve el error de inmediato; nunca obtiene un resultado
   equivocado por este motivo.
2. **Es un caso raro.** Llamar a un local `+`, `-`, `*` o `/` no es idiomático en código
   numérico. En un proyecto donde las funciones diferenciables son expresiones
   matemáticas cortas, la probabilidad de encontrarlo es muy baja.
3. **El costo no es pequeño.** Hace falta seguir el alcance léxico a través de `let*`, `fn*`
   (con sus aridades), `loop*`, `letfn*` y `catch`, y decidir qué hacer con los
   destructuramientos ya expandidos. Es bastante más código y pruebas que todo el
   incremento 1, para un caso que falla de forma segura.
4. **No ataca los límites que sí pesan.** Los más relevantes son otros: que no se
   reescriben las funciones definidas fuera del cuerpo (límite 1), que la aritmética de
   índices también se reescribe (límite 2) y que llamada con números la función devuelve un
   `Dual` y no un número (límite 7). Ninguno depende del sombreado.
5. **La colisión con un `+` propio del namespace** (el caso de `resolve`) tampoco es un caso
   de uso del proyecto: quien excluye `+` de `clojure.core` ya está escribiendo con sus
   propias operaciones, y no necesita `ad-fn`.

### Lo que sigue sin cubrirse (y está documentado)
- **Casos silenciosos:** un símbolo usado como *dato* dentro del cuerpo (por ejemplo
  `'(+ 1 2)`, o la constante de un `case`) se reescribe sin aviso. Es el límite 4. Es el único
  modo de colisión de nombres que no produce un error de compilación, y también se deja tal
  cual por la misma razón (rareza), con su prueba (`limite-4-las-formas-citadas-...`).
- Si en el futuro hiciera falta, el incremento 2 se puede retomar sin rehacer nada: es una
  mejora del recorrido de `reescribir-forma`, aislada en `fpnum.ad.syntax`.