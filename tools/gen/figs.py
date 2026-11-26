"""Helper común: convierte una figura de matplotlib en un <svg> inline liviano para el campo `fig`.
Uso:
    from figs import nueva, svg
    fig, ax = nueva()            # o nueva(2, 1) para subplots
    ax.plot(...)
    s = svg(fig)                 # string "<svg ...>...</svg>"
"""
import io, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "svg.fonttype": "none", "font.family": "sans-serif", "font.size": 10,
    "axes.grid": True, "grid.alpha": 0.3, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 72, "lines.linewidth": 1.6,
})

def nueva(filas=1, cols=1, ancho=6.4, alto=None):
    alto = alto or (3.0 * filas if filas > 1 else 3.2)
    fig, ax = plt.subplots(filas, cols, figsize=(ancho, alto))
    return fig, ax

def svg(fig, max_kb=150):
    buf = io.StringIO()
    fig.tight_layout()
    fig.savefig(buf, format="svg", metadata={"Date": None, "Creator": None})
    plt.close(fig)
    s = buf.getvalue()
    s = s[s.index("<svg"):]
    s = re.sub(r"<metadata>.*?</metadata>", "", s, flags=re.S)
    s = re.sub(r"<!--.*?-->", "", s, flags=re.S)
    s = re.sub(r'\s(width|height)="[^"]*pt"', "", s, count=2)   # que escale con el contenedor
    s = re.sub(r">\s+<", "><", s)
    s = re.sub(r"(\d+\.\d{2})\d+", r"\1", s)                      # recorta decimales de coordenadas
    if len(s) > max_kb * 1000:
        raise ValueError(f"SVG de {len(s)//1000} KB: simplificá (menos puntos, rasterized=True en imshow)")
    return s.strip()
