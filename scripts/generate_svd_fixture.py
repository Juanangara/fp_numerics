"""Genera el fixture de referencia para las pruebas de la SVD.

Usa NumPy (numpy.linalg.svd, que llama a LAPACK) para calcular los valores
singulares de varias matrices y los guarda en test/fpnum/svd/svd_reference.edn.

Es solo una herramienta para validar: el código de src/ NO usa NumPy ni
ninguna librería de álgebra lineal. El fixture ya generado está en el
repositorio, así que ejecutar este script solo hace falta para regenerarlo.

Dependencias: python -m pip install -r scripts/requirements.txt
Uso (desde la raíz del proyecto):
    python scripts/generate_svd_fixture.py
"""
import pathlib

import numpy as np

SALIDA = pathlib.Path("test/fpnum/svd/svd_reference.edn")
SEMILLA = 12345


def hilbert(n):
    i = np.arange(1, n + 1)
    return 1.0 / (i[:, None] + i[None, :] - 1.0)


def con_valores_singulares(rng, m, n, sigmas):
    """Matriz m x n con los valores singulares dados (U y V ortogonales
    aleatorias obtenidas con QR)."""
    k = min(m, n)
    u, _ = np.linalg.qr(rng.standard_normal((m, k)))
    v, _ = np.linalg.qr(rng.standard_normal((n, k)))
    return (u * np.asarray(sigmas)) @ v.T


def edn_vector(v):
    return "[" + " ".join(repr(float(x)) for x in v) + "]"


def edn_matriz(a):
    return "[" + "\n     ".join(edn_vector(fila) for fila in a) + "]"


def main():
    rng = np.random.default_rng(SEMILLA)
    casos = [
        ("diagonal-3x3", np.diag([3.0, -2.0, 1.0])),
        ("cuadrada-3x3", np.array([[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 4.0]])),
        ("alta-5x3", rng.standard_normal((5, 3))),
        ("ancha-3x5", rng.standard_normal((3, 5))),
        ("aleatoria-10x10", rng.standard_normal((10, 10))),
        ("aleatoria-20x8", rng.standard_normal((20, 8))),
        ("rango-1-4x4", np.outer([1.0, 2.0, 3.0, 4.0], [1.0, -1.0, 2.0, 0.5])),
        ("rango-2-6x5", rng.standard_normal((6, 2)) @ rng.standard_normal((2, 5))),
        ("decaimiento-geometrico-12x12",
         con_valores_singulares(rng, 12, 12, np.logspace(0, -10, 12))),
        ("hilbert-5", hilbert(5)),
        ("hilbert-8", hilbert(8)),
        ("hilbert-10", hilbert(10)),
    ]

    lineas = [
        ";; GENERADO por scripts/generate_svd_fixture.py. No editar a mano.",
        f";; numpy {np.__version__}, semilla {SEMILLA}.",
        ";; Valores singulares calculados con numpy.linalg.svd (LAPACK), en orden",
        ";; decreciente. LAPACK solo garantiza error absoluto ~ eps*sigma_max, así que",
        ";; los valores singulares muy pequeños de la referencia (p. ej. en Hilbert)",
        ";; tienen poca precisión relativa: comparar con tolerancia absoluta.",
        "[",
    ]
    for nombre, a in casos:
        s = np.linalg.svd(a, compute_uv=False)
        lineas.append(f' {{:name "{nombre}"')
        lineas.append(f"  :cond {repr(float(s[0] / s[-1])) if s[-1] > 0 else 'nil'}")
        lineas.append(f"  :matrix {edn_matriz(a)}")
        lineas.append(f"  :s {edn_vector(s)}}}")
    lineas.append("]")

    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    SALIDA.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    print(f"Escrito {SALIDA} con {len(casos)} casos")


if __name__ == "__main__":
    main()
