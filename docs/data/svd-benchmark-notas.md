# Benchmark de la SVD: antes y después

Datos: matrices aleatorias n×n con entradas gaussianas (semilla fija, `fpnum.demo.svd/random-matrix`).
Algoritmo: Jacobi de un solo lado. Ambas implementaciones hacen **los mismos barridos** (10, 13 y 13 para
n = 64, 128 y 256) y dan **resultados idénticos bit a bit** (σ, U y V), con error relativo de reconstrucción
de 1.2e-14, 2.6e-14 y 5.6e-14.

- **Antes:** `fpnum.svd.core/svd-reference` (vectores inmutables) → `svd-benchmark-reference.csv`
- **Después:** `fpnum.svd.core/svd` (arrays `double[]` + normas de columna en caché) → `svd-benchmark-fast.csv`

Reproducir: `clojure -M -m fpnum.demo.svd`

## Entorno
Intel Core i5-1235U (portátil), Windows 11, Java 25 (HotSpot). Con otros programas abiertos durante las
mediciones (VS Code, Docker, VPN, navegador). Los tiempos de un portátil varían con la carga y la
temperatura, así que se dan con su rango.

## Tiempos (segundos)

| n | `svd-reference` (antes) | `svd` (después) |
|---|---|---|
| 64 | 0.8 – 1.7 | 0.01 – 0.06 |
| 128 | 5.8 – 15.7 | 0.04 – 0.24 |
| 256 | 44 – 113 | 0.57 – 1.0 |

Cómo se obtuvieron los rangos:

- **CSV de la demo** (una sola corrida, con un calentamiento de la JVM sobre una matriz de 32×32):
  reference 0.83 / 15.70 / 112.79 s, fast 0.055 / 0.237 / 1.005 s.
- **Medición anterior a la refactorización** (misma implementación de referencia, mismo equipo):
  0.78 / 5.76 / 45.0 s.
- **Repeticiones en una misma JVM ya caliente**: fast n = 256 → 0.61, 0.70, 0.78, 0.57, 0.71 s;
  fast n = 128 → 0.04 a 0.12 s; reference n = 128 → 15.4, 15.5, 14.9 s; reference n = 256 → 54.3 y 44.3 s.

## Cómo citarlo
La aceleración depende de qué corrida se compare, porque `svd-reference` varía mucho (genera mucha basura
y el recolector de basura y la carga del equipo pesan). Una afirmación defendible:

> Para 256×256, la versión de referencia tarda entre ~45 y ~113 s y la rápida entre ~0.6 y ~1 s: una
> mejora de **dos órdenes de magnitud** (unas 50 a 110 veces), con resultados idénticos.

Evitar citar una sola cifra de aceleración sin su rango.
