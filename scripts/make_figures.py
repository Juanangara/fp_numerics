"""Genera las figuras del README a partir de los CSV de docs/data/.

Escribe en docs/img/ (SVG y PNG de cada figura, salvo el montaje, que es solo PNG)
y escribe docs/data/ajustes-log-log.csv con los ajustes que usan las figuras
(exponentes de las leyes de potencia y punto de cruce forward/reverse).

Figuras
-------
  error-vs-h               error de la derivada por diferencias finitas frente a h, y por AD
  valores-singulares       valores singulares de cada canal de la imagen
  montaje-svd              la imagen reconstruida con distintos k      (solo PNG)
  forward-vs-reverse       tiempo de un gradiente frente al número de variables n
  iteraciones-vs-kappa     iteraciones del descenso de gradiente frente a κ(X)

Diseño
------
* Paleta Okabe-Ito (apta para daltonismo). Forward y reverse tienen SIEMPRE el
  mismo color (azul y bermellón) y el mismo marcador (círculo y cuadrado), y esos
  dos colores no se usan para nada más en ninguna figura.
* Tipografía grande, pensada para diapositivas. PNG a 200 dpi.
* Determinista: sal de hash fija y sin fecha en el SVG, así que regenerar las
  figuras con las mismas versiones no cambia los archivos.

Uso (desde la raíz del proyecto):
    python -m pip install -r scripts/requirements.txt
    python scripts/make_figures.py
"""
import csv
import math
import pathlib
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

RAIZ = pathlib.Path(__file__).resolve().parent.parent
DATOS = RAIZ / "docs" / "data"
IMG = RAIZ / "docs" / "img"
IMG_DEMO = RAIZ / "docs" / "demo-svd" / "img"

# ---------------------------------------------------------------------------
# Paleta Okabe-Ito
# ---------------------------------------------------------------------------
NEGRO = "#000000"
AZUL = "#0072B2"        # FORWARD
BERMELLON = "#D55E00"   # REVERSE
VERDE = "#009E73"
CELESTE = "#56B4E9"
PURPURA = "#CC79A7"
GRIS = "#5A5A5A"
FORWARD, REVERSE = AZUL, BERMELLON
MARCA_FORWARD, MARCA_REVERSE = "o", "s"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 16,
    "axes.labelsize": 17,
    "axes.titlesize": 17,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "legend.fontsize": 14.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.28,
    "grid.linewidth": 1.0,
    "lines.linewidth": 2.6,
    "lines.markersize": 9,
    "axes.linewidth": 1.4,
    "svg.fonttype": "path",        # el texto va como trazos: se ve igual en cualquier equipo
    "svg.hashsalt": "fp_numerics",  # ids estables: regenerar no cambia el SVG
    "figure.dpi": 100,
})


def guardar(fig, nombre, svg=True):
    """Guarda la figura como PNG (200 dpi) y, salvo que svg=False, como SVG."""
    IMG.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMG / f"{nombre}.png", dpi=200, bbox_inches="tight",
                metadata={"Software": None})
    if svg:
        fig.savefig(IMG / f"{nombre}.svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)
    print("escrito", f"docs/img/{nombre}.png" + (" y .svg" if svg else ""))


# ---------------------------------------------------------------------------
# Lectura de CSV y ajustes
# ---------------------------------------------------------------------------
def _convertir(v):
    if v is None or v == "":
        return None
    try:
        return float(v)
    except ValueError:
        return v


def leer(nombre):
    with open(DATOS / nombre, newline="", encoding="utf-8") as f:
        return [{k: _convertir(v) for k, v in fila.items()} for fila in csv.DictReader(f)]


def ajuste_potencia(x, y):
    """Ajuste y = c * x**p por mínimos cuadrados en escala log-log. Devuelve (p, c)."""
    p, b = np.polyfit(np.log10(x), np.log10(y), 1)
    return float(p), float(10 ** b)


def cruce_en_uno(n, razon):
    """n en que la razón (forward/reverse) cruza 1, por interpolación log-log entre
    los dos puntos medidos que la rodean. None si no cruza."""
    for (n1, r1), (n2, r2) in zip(zip(n, razon), zip(n[1:], razon[1:])):
        if (r1 - 1.0) * (r2 - 1.0) < 0:
            t = (0.0 - math.log(r1)) / (math.log(r2) - math.log(r1))
            return math.exp(math.log(n1) + t * (math.log(n2) - math.log(n1)))
    return None


def por_funcion(filas, funcion):
    return sorted((r for r in filas if r["function"] == funcion), key=lambda r: r["n"])


FUNCIONES = [("minimos-cuadrados", "Mínimos cuadrados (m = 100)"),
             ("suma-de-cuadrados", "Suma de cuadrados Σ wᵢ²")]


def calcular_ajustes():
    """Calcula los ajustes de las figuras y los escribe en docs/data/ajustes-log-log.csv."""
    filas = []  # (fit, function, source, metric, value)
    for archivo in ("ad-forward-vs-reverse.csv", "ad-forward-vs-reverse-repeticion.csv"):
        datos = leer(archivo)
        for funcion, _ in FUNCIONES:
            f = por_funcion(datos, funcion)
            n = [r["n"] for r in f]
            for modo, columna in (("forward", "forward_ms"), ("reverse", "reverse_ms")):
                p, _ = ajuste_potencia(n, [r[columna] for r in f])
                filas.append(("exponente-n", funcion, archivo, modo, p))
    cruce = leer("ad-forward-vs-reverse-cruce.csv")
    for funcion, _ in FUNCIONES:
        f = por_funcion(cruce, funcion)
        n_cruce = cruce_en_uno([r["n"] for r in f], [r["ratio_forward_over_reverse"] for r in f])
        filas.append(("cruce-forward-reverse", funcion, "ad-forward-vs-reverse-cruce.csv",
                      "n_donde_la_razon_vale_1", n_cruce))
    gd = leer("conditioning-gd.csv")
    medidas = [r for r in gd if r["iters_dist"] is not None and r["kappa"] >= 10]
    p, _ = ajuste_potencia([r["kappa"] for r in medidas], [r["iters_dist"] for r in medidas])
    filas.append(("exponente-kappa", "descenso-de-gradiente", "conditioning-gd.csv",
                  "iteraciones_vs_kappa_con_kappa_ge_10", p))
    destino = DATOS / "ajustes-log-log.csv"
    with open(destino, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["fit", "function", "source", "metric", "value"])
        for fila in filas:
            w.writerow([fila[0], fila[1], fila[2], fila[3], repr(fila[4])])
    print("escrito docs/data/ajustes-log-log.csv", f"({len(filas)} filas)")
    return {(f[0], f[1], f[2], f[3]): f[4] for f in filas}


def pot10(v):
    """1e-08 -> 10⁻⁸ (para potencias exactas de 10)."""
    e = round(math.log10(v))
    return "10" + str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))


def sci(v):
    """Formato 3.0e-9 -> '3.0·10⁻⁹' legible en una figura."""
    m, e = f"{v:.1e}".split("e")
    sup = str(int(e)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{m}·10{sup}"


# ---------------------------------------------------------------------------
# 1. error-vs-h
# ---------------------------------------------------------------------------
def fig_error_vs_h():
    filas = leer("error-vs-h.csv")
    h = np.array([r["h"] for r in filas])
    fd = np.array([r["fd_error"] for r in filas])
    mejor = int(np.argmin(fd))
    fig, ax = plt.subplots(figsize=(11.5, 6.6))
    ax.loglog(h, fd, "D-", color=GRIS, label="diferencias finitas hacia adelante")
    ax.axhline(1e-17, color=FORWARD, lw=3.0)
    ax.text(1.5e-16, 2.2e-17, "AD (forward): error = 0, exacto y sin depender de h",
            color=FORWARD, fontsize=15, va="bottom", ha="left", fontweight="bold")
    ax.plot([h[mejor]], [fd[mejor]], "o", ms=15, mfc="none", mec=NEGRO, mew=2.6, zorder=5)
    ax.annotate(f"mínimo: h = {pot10(h[mejor])}\nerror {sci(fd[mejor])}",
                xy=(h[mejor], fd[mejor]), xytext=(2.5e-12, 2e-11),
                fontsize=15, ha="left", va="center",
                arrowprops=dict(arrowstyle="->", color=NEGRO, lw=2.0))
    ax.text(2e-16, 6e-13, "cancelación\n(error ∝ ε/h)", fontsize=15, color=GRIS, va="center")
    ax.text(0.3, 2.5e-5, "truncamiento\n(error ∝ h)", fontsize=15, color=GRIS, ha="right", va="top")
    ax.set_xlabel("paso h")
    ax.set_ylabel("error absoluto de la derivada")
    ax.set_xlim(8e-17, 0.4)
    ax.set_ylim(3e-18, 3)
    ax.legend(loc="upper center", frameon=False)
    guardar(fig, "error-vs-h")


# ---------------------------------------------------------------------------
# 2. valores-singulares
# ---------------------------------------------------------------------------
def fig_valores_singulares():
    filas = leer("svd-image-singular-values.csv")
    i = np.array([r["i"] for r in filas])
    fig, ax = plt.subplots(figsize=(11.5, 6.6))
    canales = [("r", "canal R", PURPURA, "-"), ("g", "canal G", VERDE, "--"),
               ("b", "canal B", CELESTE, ":")]
    for col, etiqueta, color, estilo in canales:
        ax.semilogy(i, [r[col] for r in filas], estilo, color=color, label=etiqueta, lw=3.0)
    ks = [1, 5, 10, 20, 50]
    for k in ks:
        ax.axvline(k, color=GRIS, lw=1.2, ls=(0, (2, 3)), alpha=0.8)
    ax.set_xlabel("índice i")
    ax.set_ylabel("valor singular σᵢ (escala log)")
    ax.set_xlim(0, i.max() + 1)
    abajo = ax.get_ylim()[0]
    for k in ks:
        ax.text(k + 1.5, abajo * 1.12, f"k = {k}", fontsize=13.5, color=GRIS, rotation=90,
                va="bottom", ha="left")
    ax.legend(loc="upper right", frameon=False)
    guardar(fig, "valores-singulares")


# ---------------------------------------------------------------------------
# 3. montaje-svd (solo PNG)
# ---------------------------------------------------------------------------
def fig_montaje():
    metricas = {int(r["k"]): r for r in leer("svd-image-metrics.csv")}
    rango = max(metricas)
    paneles = [("original.png", "Original", f"{rango} valores singulares")]
    for k in (1, 5, 10, 20, 50):
        m = metricas[k]
        paneles.append((f"k-{k:03d}.png", f"k = {k}",
                        f"{m['compression_ratio']:.1f}× menos números · {m['psnr_db']:.1f} dB"))
    fig, ejes = plt.subplots(2, 3, figsize=(15.5, 7.6))
    for ax, (archivo, titulo, detalle) in zip(ejes.ravel(), paneles):
        ax.imshow(plt.imread(IMG_DEMO / archivo), interpolation="nearest")
        ax.set_title(titulo, fontsize=18, fontweight="bold", pad=22)
        ax.text(0.5, 1.02, detalle, transform=ax.transAxes, ha="center", va="bottom", fontsize=14)
        ax.axis("off")
    fig.subplots_adjust(wspace=0.06, hspace=0.34)
    fig.text(0.5, 0.045, "Imagen sintética reproducible: las cifras de compresión y de PSNR dependen de la imagen",
             ha="center", va="top", fontsize=13.5, color=GRIS)
    guardar(fig, "montaje-svd", svg=False)


# ---------------------------------------------------------------------------
# 4. forward-vs-reverse
# ---------------------------------------------------------------------------
def fig_forward_vs_reverse(ajustes):
    principal = leer("ad-forward-vs-reverse.csv")
    cruce = leer("ad-forward-vs-reverse-cruce.csv")
    fig, ejes = plt.subplots(1, 2, figsize=(16, 6.9))
    for ax, (funcion, titulo) in zip(ejes, FUNCIONES):
        f = por_funcion(principal, funcion)
        n = np.array([r["n"] for r in f])
        fw = np.array([r["forward_ms"] for r in f])
        rv = np.array([r["reverse_ms"] for r in f])
        ax.loglog(n, fw, MARCA_FORWARD + "-", color=FORWARD, label="forward")
        ax.loglog(n, rv, MARCA_REVERSE + "-", color=REVERSE, label="reverse")
        extra = por_funcion(cruce, funcion)
        ne = np.array([r["n"] for r in extra])
        ax.loglog(ne, [r["forward_ms"] for r in extra], MARCA_FORWARD, mfc="none", mec=FORWARD,
                  mew=2.0, ms=8)
        ax.loglog(ne, [r["reverse_ms"] for r in extra], MARCA_REVERSE, mfc="none", mec=REVERSE,
                  mew=2.0, ms=8)
        xs = np.logspace(math.log10(3.5), math.log10(260), 50)
        for y, color in ((fw, FORWARD), (rv, REVERSE)):
            p, c = ajuste_potencia(n, y)
            ax.loglog(xs, c * xs ** p, ":", color=color, lw=1.8, alpha=0.9)
        n_cruce = ajustes[("cruce-forward-reverse", funcion, "ad-forward-vs-reverse-cruce.csv",
                           "n_donde_la_razon_vale_1")]
        ax.axvline(n_cruce, color=GRIS, lw=1.6, ls="--")
        ax.text(n_cruce * 1.08, ax.get_ylim()[0] * 1.5, f"cruce\nn ≈ {n_cruce:.0f}", fontsize=14,
                color=GRIS, va="bottom", ha="left")
        ax.set_title(titulo, pad=12)
        ax.set_xlabel("número de variables n")
        ax.set_xticks([4, 10, 20, 50, 100, 200])
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
    ejes[0].set_ylabel("tiempo de un gradiente (ms)")
    # exponentes medidos, junto a las rectas
    for ax, (funcion, _) in zip(ejes, FUNCIONES):
        pf = ajustes[("exponente-n", funcion, "ad-forward-vs-reverse.csv", "forward")]
        pr = ajustes[("exponente-n", funcion, "ad-forward-vs-reverse.csv", "reverse")]
        ax.text(0.97, 0.20, f"forward ∝ n^{pf:.1f}", transform=ax.transAxes, color=FORWARD,
                fontsize=16, fontweight="bold", va="top", ha="right")
        ax.text(0.97, 0.11, f"reverse ∝ n^{pr:.1f}", transform=ax.transAxes, color=REVERSE,
                fontsize=16, fontweight="bold", va="top", ha="right")
    leyenda = [Line2D([], [], color=FORWARD, marker=MARCA_FORWARD, label="forward"),
               Line2D([], [], color=REVERSE, marker=MARCA_REVERSE, label="reverse"),
               Line2D([], [], color=GRIS, marker="o", mfc="none", ls="none", mew=2.0,
                      label="n = 4 a 30 (corrida aparte)")]
    fig.legend(handles=leyenda, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3,
               frameon=False)
    fig.text(0.5, -0.095, "Tiempos de esta implementación (Clojure, números en caja), medidos en un portátil: "
             "las razones y los exponentes son más fiables que los ms absolutos",
             ha="center", va="top", fontsize=13.5, color=GRIS)
    fig.subplots_adjust(wspace=0.18)
    guardar(fig, "forward-vs-reverse")


# ---------------------------------------------------------------------------
# 5. iteraciones-vs-kappa
# ---------------------------------------------------------------------------
def fig_iteraciones_vs_kappa(ajustes):
    gd = leer("conditioning-gd.csv")
    kappa = np.array([r["kappa"] for r in gd])
    pred = np.array([r["iters_pred"] for r in gd])
    alcanzaron = [r for r in gd if r["iters_dist"] is not None]
    fallaron = [r for r in gd if r["iters_dist"] is None]
    fig, ax = plt.subplots(figsize=(11.8, 7.0))
    ax.loglog(kappa, pred, "--", color=PURPURA, lw=2.8, zorder=2,
              label="predicción  ln(tol·√n) / ln|ρ|")
    ks = np.logspace(0, 3, 30)
    ax.loglog(ks, 7.7 * ks ** 2, ":", color=GRIS, lw=2.2, zorder=1, label="referencia ∝ κ(X)²")
    ax.loglog([r["kappa"] for r in alcanzaron], [r["iters_dist"] for r in alcanzaron], "o",
              color=NEGRO, ms=11, zorder=4, label="medidas (el descenso llegó a la tolerancia)")
    for r in fallaron:
        tope = r["iterations"]
        ax.axhline(tope, color=GRIS, lw=1.6, ls=(0, (6, 3)), zorder=1)
        ax.text(1.12, tope * 1.18, f"tope de {tope:,.0f} iteraciones".replace(",", " "),
                color=GRIS, fontsize=14, ha="left", va="bottom")
        ax.plot([r["kappa"]], [tope], "^", ms=17, mfc="white", mec=NEGRO, mew=3.0, zorder=6)
        ax.annotate("κ = 1000: NO CONVERGIÓ\nen 300 000 iteraciones\n(error relativo final "
                    f"{r['final_dist']:.2f}, falta mucho)",
                    xy=(r["kappa"], tope), xytext=(120, 2.4e3), fontsize=14.5,
                    ha="left", va="center", color=NEGRO, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.35", fc="#FFFFFF", ec=NEGRO, lw=1.8),
                    arrowprops=dict(arrowstyle="->", color=NEGRO, lw=2.2))
    for r in fallaron:
        ax.annotate(f"la predicción dice ≈ {r['iters_pred']/1e6:.1f} millones", xy=(r["kappa"], r["iters_pred"]),
                    xytext=(12, 1.6e6), fontsize=14.5, color=PURPURA, fontweight="bold",
                    ha="left", va="center",
                    arrowprops=dict(arrowstyle="->", color=PURPURA, lw=2.0))
    p = ajustes[("exponente-kappa", "descenso-de-gradiente", "conditioning-gd.csv",
                 "iteraciones_vs_kappa_con_kappa_ge_10")]
    ax.text(0.03, 0.97, f"medido: iteraciones ∝ κ(X)^{p:.2f}  (κ ≥ 10)", transform=ax.transAxes,
            fontsize=16, va="top", fontweight="bold")
    ax.text(0.03, 0.895, "esperado: iteraciones ∝ κ(XᵀX) = κ(X)², exponente 2",
            transform=ax.transAxes, fontsize=14.5, va="top", color=GRIS)
    ax.set_xlabel("número de condición κ(X) = σmax/σmin")
    ax.set_ylabel("iteraciones hasta error relativo ≤ 10⁻⁶")
    ax.set_xlim(0.8, 1500)
    ax.set_ylim(8, 5e7)
    ax.set_xticks([1, 3, 10, 30, 100, 1000])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
    ax.legend(loc="lower right", frameon=False)
    guardar(fig, "iteraciones-vs-kappa")


def main():
    ajustes = calcular_ajustes()
    fig_error_vs_h()
    fig_valores_singulares()
    fig_montaje()
    fig_forward_vs_reverse(ajustes)
    fig_iteraciones_vs_kappa(ajustes)
    print()
    for clave, valor in ajustes.items():
        if valor is not None:
            print(f"  {clave[0]:24s} {clave[1]:24s} {clave[2]:36s} {clave[3]:38s} {valor:.3f}")


if __name__ == "__main__":
    main()
