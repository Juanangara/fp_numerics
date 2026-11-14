# Cifras del README y su fuente

Cada cifra del README sale de un CSV de `docs/data/` (regenerable con el comando indicado) o de una prueba. Los tiempos se citan como razones y exponentes porque dependen de la máquina.

| Cifra | Valor | Fuente | Se regenera con |
|---|---|---|---|
| Error mínimo de diferencias finitas | 3.0·10⁻⁹ en h = 10⁻⁸ | `data/error-vs-h.csv` (`fd_error`) | `fpnum.demo.ad` |
| Error de AD forward | 0 | `data/error-vs-h.csv` (`ad_error`) | `fpnum.demo.ad` |
| Descenso en mínimos cuadrados | 37 iteraciones; distancia a ecuaciones normales 3.0·10⁻¹¹ | `data/descent-resumen.csv` | `fpnum.demo.descent` |
| Rosenbrock | 13 756 iteraciones | `data/descent-resumen.csv` | `fpnum.demo.descent` |
| Diferencia de σ frente a NumPy | ≤ 7.8·10⁻¹⁶·σ₁ (12 matrices) | `data/svd-accuracy.csv` | `scripts/generate_accuracy_csvs.clj` |
| Reconstrucción y ortonormalidad | ≤ 1.7·10⁻¹⁵ y ≤ 5.3·10⁻¹⁵; 1 a 9 barridos | `data/svd-accuracy.csv` | `scripts/generate_accuracy_csvs.clj` |
| Jacobi directo frente a vía AᵀA | tabla de la sección SVD: precisión | `data/svd-ata-vs-directo.csv` | `scripts/generate_accuracy_csvs.clj` |
| Determinante de Hilbert | 4.5·10⁻¹², 1.0·10⁻⁷, 2.5·10⁻⁵ | `data/svd-hilbert-determinante.csv` | `scripts/generate_accuracy_csvs.clj` |
| Razón referencia/rápida de la SVD | 15×, 66×, 112× (n = 64, 128, 256; una corrida) | `data/svd-benchmark-reference.csv`, `data/svd-benchmark-fast.csv`, `data/svd-benchmark-notas.md` | `fpnum.demo.svd benchmark` |
| Barridos de la SVD en el benchmark | 10, 13, 13 | `data/svd-benchmark-reference.csv` (`sweeps`) | ídem |
| Energía, error, PSNR y razón de compresión | tabla de compresión (k = 20: 99.82 %, 32.5 dB, 3.83×) | `data/svd-image-metrics.csv` | `fpnum.demo.svd` |
| Valores singulares de la imagen | Figura 2 | `data/svd-image-singular-values.csv` | `fpnum.demo.svd` |
| Iteraciones predichas y medidas según κ | 22/24, 113, 921, 7 312, 76 444; κ = 1000 no converge (error final 0.29, predicción ≈ 7.4 millones) | `data/conditioning-gd.csv` | `fpnum.demo.conditioning` |
| Error real con ‖∇f‖ ≤ 10⁻⁸ | 4.9·10⁻⁹ … 5.0·10⁻⁵ | `data/conditioning-gd.csv` (`dist_at_grad`) | `fpnum.demo.conditioning` |
| Error SVD frente a ecuaciones normales | tabla de κ | `data/conditioning-accuracy.csv` | `fpnum.demo.conditioning` |
| Rango deficiente: 81 iteraciones, 2.8·10⁻¹⁴, rango 3 | | `data/conditioning-resumen.csv` | `fpnum.demo.conditioning` |
| Razones forward/reverse | tabla de forward frente a reverse (dos corridas) | `data/ad-forward-vs-reverse.csv`, `data/ad-forward-vs-reverse-repeticion.csv`, `data/ad-forward-vs-reverse-notas.md` | `fpnum.demo.reverse` |
| Cruce forward/reverse | n ≈ 11.1 y 10.4 | `data/ajustes-log-log.csv`, `data/ad-forward-vs-reverse-cruce.csv` | `scripts/make_figures.py` |
| Exponentes forward y reverse | 1.89–1.94 y 1.05–1.11 | `data/ajustes-log-log.csv` | `scripts/make_figures.py` |
| Exponente de κ | 1.92 (κ ≥ 10) | `data/ajustes-log-log.csv` | `scripts/make_figures.py` |
| Coste del despacho genérico (D3) | ~1–10 % | `data/ad-forward-dispatch.csv` y `docs/decisiones.md` | `scripts/bench_forward_dispatch.clj` |
| Identidad de nodos: 1a frente a 1b (D1) | desbordamiento hacia 5 000 términos; ×4 por cada +2 de profundidad; 1b lineal | `data/ad-reverse-sharing-1a.csv`, `data/ad-reverse-sharing-1b.csv` | `scripts/measure_reverse_sharing.clj` |
| Cruce y pasada por ids (D4) | n ≈ 10–11; hipótesis no verificada | `docs/decisiones.md` | — |
| Pruebas | 230 pruebas, 1160 aserciones | `clojure -M:test` | — |

Todas las figuras salen de `scripts/make_figures.py` a partir de estos CSV.
