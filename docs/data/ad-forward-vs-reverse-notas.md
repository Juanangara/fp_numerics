# AD forward frente a reverse: tiempo de un gradiente contra n

Reproducir: `clojure -M -m fpnum.demo.reverse` (~25 s).

| Archivo | Contenido |
|---|---|
| `ad-forward-vs-reverse.csv` | los cuatro tamaños pedidos (n = 10, 50, 100, 200), corrida 2 |
| `ad-forward-vs-reverse-repeticion.csv` | los mismos tamaños, corrida 1 (repetición independiente) |
| `ad-forward-vs-reverse-cruce.csv` | tamaños intermedios (n = 4 a 30) para fijar el punto de cruce |
| `ajustes-log-log.csv` | exponentes y punto de cruce, calculados desde los CSV anteriores por `scripts/make_figures.py` (los que usan las figuras) |

Cada función está escrita **una sola vez** con las operaciones genéricas (`fpnum.ad.ops`) y se
diferencia tal cual en los dos modos. Los gradientes de los dos modos coinciden (diferencia
relativa máxima 5.7e-15, y exactamente 0 en la suma de cuadrados). Reverse usa la pasada hacia
atrás por ids decrecientes (D4 en `docs/decisiones.md`).

## Método y entorno
Intel Core i5-1235U (portátil), Windows 11, Java 25 (HotSpot), con otros programas abiertos.

- Cada muestra repite el gradiente hasta durar al menos ~20 ms; se calienta la JVM ~1.5 s
  antes de medir; se reporta el tiempo **por gradiente** (mediana y rango de 7 muestras).
- **Las muestras de forward y reverse se intercalan** (ronda par: forward y reverse; ronda
  impar: reverse y forward) y la **razón se calcula por rondas pareadas** (se reporta la
  mediana y el rango de las razones). Es necesario porque el rendimiento del equipo cambia
  entre corridas y dentro de una corrida: en una corrida previa, sin intercalar, forward salió
  **2.5 veces más lento** sin que su código hubiera cambiado (1 522 ms en lugar de ~600 ms con
  n = 200). **Solo las razones pareadas son comparables entre corridas; los tiempos absolutos
  no lo son.**

## Resultados (corrida 2, `ad-forward-vs-reverse.csv`)

Mínimos cuadrados ‖X·w − y‖² con m = 100 filas (costo de la función: O(m·n)):

| n | forward (ms) | reverse (ms) | forward / reverse |
|---|---|---|---|
| 10 | 1.75 | 2.22 | **0.78** [0.67–0.81] |
| 50 | 38.4 | 11.9 | 3.15 [2.91–3.40] |
| 100 | 148.1 | 25.5 | 5.77 [5.25–6.05] |
| 200 | 587.8 | 59.5 | **9.96** [9.33–10.35] |

Suma de cuadrados Σ wᵢ² (costo de la función: O(n)):

| n | forward (ms) | reverse (ms) | forward / reverse |
|---|---|---|---|
| 10 | 0.017 | 0.021 | 0.81 [0.52–1.09] |
| 50 | 0.328 | 0.115 | 2.86 [2.50–3.45] |
| 100 | 1.258 | 0.259 | 4.77 [4.21–5.09] |
| 200 | 4.897 | 0.559 | **8.78** [8.71–9.93] |

Corrida 1 (repetición independiente): razones 0.91 / 3.94 / 7.06 / **12.27** en mínimos
cuadrados y 0.94 / 3.56 / 6.06 / **11.00** en la suma de cuadrados, para n = 10 / 50 / 100 / 200.
Entre las dos corridas las razones varían alrededor de un 15–20 %.

## Qué muestra
- **Forward crece como n² y reverse como n.** Exponentes del ajuste log-log por mínimos
  cuadrados sobre los cuatro tamaños (n = 10, 50, 100, 200; `ajustes-log-log.csv`, dos
  corridas): forward 1.92–1.94 y reverse 1.05–1.09 en mínimos cuadrados (1.89–1.92 y
  1.10–1.11 en la suma de cuadrados). Es la ley esperada: forward necesita n pasadas,
  reverse una.
- **Punto de cruce: n ≈ 10–11.** La razón forward/reverse se mide directamente en n = 10 y
  n = 12 (`ad-forward-vs-reverse-cruce.csv`): en mínimos cuadrados pasa de 0.89 a 1.09, y en
  la suma de cuadrados de 0.97 a 1.11. Interpolando en escala log-log entre esos dos puntos,
  la razón vale 1 en n = 11.1 (mínimos cuadrados) y n = 10.4 (suma de cuadrados). Por debajo,
  reverse pierde: con n = 4 es unas 2.3 veces más lento (razón 0.43),
  porque cada nodo de reverse cuesta más que una operación dual (asigna un nodo, un vector de
  argumentos y entradas en una cola de prioridad persistente). Con n = 20 la razón es 1.7 y
  con n = 30, 2.4.
- **Efecto de la nueva pasada hacia atrás (D4):** con n = 200, reverse pasó de ~102 ms (recorrido
  en profundidad) a ~53–60 ms: unas 1.9 veces más rápido, en línea con la evaluación previa
  (2 a 3 veces sobre la pasada aislada). El cruce bajó de n ≈ 25 a n ≈ 10–11. Aquel «≈ 25» era una
  interpolación entre n = 10 y n = 50, no una medida directa.

## Cómo citarlo
> Un gradiente por reverse mode es más rápido que por forward a partir de unas 10 o 11 variables.
> Con 200 variables tarda ~55–60 ms frente a ~590–650 ms (unas 10 a 12 veces menos) en un
> problema de mínimos cuadrados con 100 observaciones. El tiempo de forward crece como n² y
> el de reverse como n. Con muy pocas variables (n = 4), reverse es ~2.3 veces más lento.

Tomar estas cifras como órdenes de magnitud de este equipo, no como constantes; las razones
son más fiables que los tiempos absolutos.
