#!/usr/bin/env python3
# Genera data/m3.js — U3: Aproximación polinómica trigonométrica, Fourier, FFT, Nyquist y filtros.
# Todos los números salen de verif_u3.py (se importa); las figuras, de figs.py.
import json, os, sys, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import verif_u3 as V
from figs import nueva, svg

R = V.R
PI = np.pi


# ------------------------------------------------------------------ helpers de formato
def f(x, d=3):
    """Número con coma decimal, sin ceros de más."""
    s = f"{float(x):.{d}f}".rstrip("0").rstrip(".")
    return s.replace(".", ",").replace("-", "−")


def SI(b):
    return "sí, coincide" if b else "NO coincide"
def PK(lst):
    return ", ".join(f"({a}, {b})" if isinstance(a, int) else f"{f(a, 3)} Hz → {f(b, 2)}" for a, b in lst)


def D(term, txt, src):
    return f'<div class="def"><b>{term}</b><p>{txt}</p><span class="src">{src}</span></div>'
def AGR(txt):
    return f'<div class="agregado"><b>Agregado</b><p>{txt}</p></div>'
def TRAMPA(txt):
    return f'<div class="trampa"><b>Así lo toma la cátedra</b><p>{txt}</p></div>'
def TRAMPA_P(txt):
    return f'<div class="trampa"><b>Así lo tomaría la cátedra</b><p>{txt}</p></div>'
def EJ(txt):
    return f'<div class="ej"><b>Ejemplo</b><p>{txt}</p></div>'
def TBL(head, rows):
    h = "".join(f"<th>{x}</th>" for x in head)
    b = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<div class="tbl"><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'
def UL(*items): return "<ul>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>"
def OL(*items): return "<ol>" + "".join(f"<li>{i}</li>" for i in items) + "</ol>"
P = lambda t: f"<p>{t}</p>"
H4 = lambda t: f"<h4>{t}</h4>"
CMD = lambda t: f'<pre class="cmd">{t}</pre>'
FIGHTML = lambda s: f'<figure class="fig">{s}</figure>'

# ------------------------------------------------------------------ figuras
FIG = {}
GRIS = dict(cmap="gray", interpolation="nearest", rasterized=True)


def fig_muestras():
    fig, ax = nueva(2, 1, alto=5.0)
    for nom, t_, y_, est, lab in (("m1", V.t1, V.y1, "-", "Muestra 1"), ("m2", V.t2, V.y2, "--", "Muestra 2")):
        fr, a = V.espectro(t_, y_)
        ax[0].plot(fr, a, est, label=lab, lw=1.2)
        m = fr <= 2.5
        ax[1].plot(fr[m], a[m], est, label=lab, lw=1.4)
    ax[0].set_xlim(0, 25); ax[1].set_xlim(0, 2.5)
    for a_ in ax:
        a_.set_xlabel("f [Hz]"); a_.set_ylabel("|FFT| normalizada"); a_.legend(loc="upper right", fontsize=8)
    ax[0].set_title("Espectros de muestra1.csv y muestra2.csv (arriba 0–25 Hz, abajo zoom 0–2,5 Hz)", fontsize=9)
    return svg(fig)


import re as _re
def seguro(fn):
    """Reintenta con otro ancho si el base64 del PNG rasterizado cae por azar en el patrón on\\w+= del validador."""
    for k in range(20):
        s_ = fn(5.6 + 0.02 * k)
        if not _re.search(r"<script|on\w+=|<image[^>]+href=\"http", s_, _re.I):
            return s_
    raise ValueError("no se pudo generar una figura limpia")


def fig_par_img(a, b, ta, tb):
    def hacer(ancho):
        fig, ax = nueva(1, 2, ancho=ancho, alto=2.9)
        ax[0].imshow(a, vmin=0, vmax=1, **GRIS); ax[0].set_title(ta, fontsize=9)
        ax[1].imshow(b, **GRIS); ax[1].set_title(tb, fontsize=9)
        for x in ax: x.axis("off")
        return svg(fig)
    return seguro(hacer)


def fig_tablero4():
    pb, spec = V.aplicar(V.TB, V.filtro_pasa_bajo(R["D0"][0], V.TB.shape))
    pa, _ = V.aplicar(V.TB, V.filtro_pasa_alto(R["D0"][1], V.TB.shape))
    fig, ax = nueva(1, 4, ancho=8.4, alto=2.6)
    ax[0].imshow(V.TB, vmin=0, vmax=1, **GRIS); ax[0].set_title("Original", fontsize=9)
    ax[1].imshow(spec, **GRIS); ax[1].set_title("log(1+|fftshift(fft2)|)", fontsize=8)
    ax[2].imshow(pb, **GRIS); ax[2].set_title(f"Pasa bajo D0={R['D0'][0]}", fontsize=9)
    ax[3].imshow(pa, **GRIS); ax[3].set_title(f"Pasa alto D0={R['D0'][1]}", fontsize=9)
    for x in ax: x.axis("off")
    return svg(fig)


def fig_sintetica():
    m = V.filtro_pasa_alto(10, V.FS.shape)
    res, spec = V.aplicar(V.FS, m)
    fig, ax = nueva(1, 5, ancho=9.0, alto=2.3)
    ims = [(V.FS, "Imagen"), (spec, "Espectro centrado"), (m, "Máscara (centrada)"),
           (np.fft.ifftshift(m), "Máscara sin centrar"), (res, "Resultado")]
    for a_, (im, tt) in zip(ax, ims):
        a_.imshow(im, **GRIS); a_.set_title(tt, fontsize=8); a_.axis("off")
    return svg(fig)


def fig_senal(t, a, b, ta, tb):
    fig, ax = nueva(2, 1, alto=4.2)
    ax[0].plot(t[::3], a[::3], lw=0.7, label=ta); ax[1].plot(t[::8], b[::8], lw=1.4, color="C1", label=tb)
    for x in ax: x.legend(loc="upper right", fontsize=8)
    ax[1].set_xlabel("tiempo [s]")
    return svg(fig)


def fig_muestreo_fs2():
    tc = np.linspace(0, 2, 300)
    fig, ax = nueva(2, 1, alto=4.0)
    for a_, fq in zip(ax, (1, 2)):
        ts = np.arange(0, 2 + 1 / 2, 1 / 2)
        a_.plot(tc, np.sin(2 * PI * fq * tc + PI / 2), color="0.7", lw=1, label=f"señal continua f = {fq} Hz")
        a_.plot(ts, np.sin(2 * PI * fq * ts + PI / 2), "o-", label="muestras con fs = 2 Hz")
        a_.legend(loc="lower right", fontsize=7); a_.set_ylim(-1.3, 1.3)
    ax[1].set_xlabel("tiempo [s]")
    return svg(fig)


def fig_muestreo_grid():
    tc = np.linspace(0, 2, 240)
    fsl = [10, 4, 2, 1]
    fig, ax = nueva(4, 2, ancho=7.0, alto=6.4)
    for i, fsm in enumerate(fsl):
        ts = np.arange(0, 2 + 1 / fsm, 1 / fsm)
        for j, fq in enumerate((1, 2)):
            a_ = ax[i][j]
            a_.plot(tc, np.sin(2 * PI * fq * tc + PI / 2), color="0.75", lw=0.8)
            a_.plot(ts, np.sin(2 * PI * fq * ts + PI / 2), "o-", ms=3, lw=1.1)
            a_.set_title(f"fs = {fsm} Hz, f = {fq} Hz", fontsize=8); a_.set_ylim(-1.3, 1.3)
            a_.tick_params(labelsize=7)
    return svg(fig)


def fig_ej_s():
    fig, ax = nueva(1, 2, ancho=7.0, alto=2.8)
    t = np.linspace(0, 13, 300)
    ax[0].plot(t, np.sin(t) + 1.2 * np.cos(2 * t)); ax[0].set_title("s(t) = sin(t) + 1,2·cos(2t)", fontsize=9)
    ax[0].set_xlabel("t [s]")
    for (fr, am), c in zip(R["ej_s"]["picos_fft"], ("C2", "C0")):
        ax[1].vlines(fr, 0, am, colors=c, lw=3)
    ax[1].set_xlim(0, 1); ax[1].set_ylim(0, 1.4); ax[1].set_xlabel("f [Hz]"); ax[1].set_ylabel("amplitud")
    ax[1].set_title("Espectro de frecuencia", fontsize=9)
    return svg(fig)


def fig_rayas():
    yy, xx = np.indices((128, 128))
    img = 0.5 + 0.5 * np.sin(2 * PI * xx / 16)
    sp = np.log(1 + np.abs(np.fft.fftshift(np.fft.fft2(img))))
    return fig_par_img(img, sp, "Imagen: franjas verticales (período 16 px)", "log(1+|fftshift(fft2)|)")


def fig_ev5():
    fig, ax = nueva(2, 1, alto=4.6)
    for nom, fs_, est in (("Muestra 1", 10.0, "--"), ("Muestra 2", 60.0, "-")):
        tt = np.arange(0, 20, 1 / fs_)
        fr, a = V.espectro(tt, V.senial_p5(tt))
        ax[0].plot(fr, a, est, label=nom, lw=1.2)
        m = fr <= 4
        ax[1].plot(fr[m], a[m], est, label=nom, lw=1.3)
    ax[0].set_xlim(0, 25); ax[1].set_xlim(0, 4)
    for a_ in ax:
        a_.set_xlabel("f [Hz]"); a_.legend(loc="upper right", fontsize=8)
    return svg(fig)


def fig_ev8():
    tc = np.linspace(0, 0.5, 500)
    k = np.arange(0, 11) / 20
    fig, ax = nueva(3, 1, alto=5.6)
    for a_, fq in zip(ax, (10, 20, 100)):
        a_.plot(tc, np.cos(2 * PI * fq * tc), color="0.7", lw=0.7, label=f"cos(2π·{fq}·t)")
        a_.plot(k, np.cos(2 * PI * fq * k), "o-", label="muestras fs = 20 Hz")
        a_.legend(loc="lower right", fontsize=7); a_.set_ylim(-1.4, 1.4); a_.tick_params(labelsize=7)
    ax[2].set_xlabel("tiempo [s]")
    return svg(fig)


def fig_tiempos():
    fig, ax = nueva(alto=2.8)
    labs = [f"{d['n']}" for d in R["tiempos"]]
    ax.bar(labs, [d["ms"] for d in R["tiempos"]], color=["C0", "C0", "C2", "C2", "C3"])
    ax.set_ylabel("ms por fft (esta máquina)"); ax.tick_params(axis="x", labelsize=8)
    return svg(fig)


def fig_escalon():
    x = np.linspace(-PI, PI, 400)
    fig, ax = nueva(alto=2.9)
    ax.plot(x, np.sign(np.sin(x)), color="k", lw=1, label="escalón")
    for K in (1, 3, 10):
        s = sum(4 / ((2 * m - 1) * PI) * np.sin((2 * m - 1) * x) for m in range(1, K + 1))
        ax.plot(x, s, lw=1.1, label=f"{K} término(s) impares")
    ax.legend(fontsize=7, loc="lower right")
    return svg(fig)


def fig_alias():
    fo = np.linspace(0, 3000, 1201)
    fig, ax = nueva(alto=2.8)
    ax.plot(fo, [V.alias(v, 1000) for v in fo])
    ax.axvline(500, ls="--", color="C3", lw=1); ax.text(520, 520, "fN = fs/2", color="C3", fontsize=8)
    ax.set_xlabel("frecuencia original [Hz]"); ax.set_ylabel("frecuencia reconstruida [Hz]")
    ax.set_title("fs = 1000 Hz: la frecuencia se refleja en [0, fs/2]", fontsize=9)
    return svg(fig)


def fig_resol():
    fig, ax = nueva(1, 2, ancho=7.2, alto=2.8)
    for a_, (b_, fsm) in zip(ax, ((2 * PI, 10), (6 * PI, 6))):
        tt = np.arange(0, b_, 1 / fsm); n = tt.size
        y = sum(np.sin(2 * PI * q * tt) for q in (0.5, 1, 2))
        Y = np.abs(np.fft.fft(y)); fr = np.fft.fftfreq(n, 1 / fsm); c = fr > 0
        a_.plot(fr[c], Y[c], ".-", ms=3)
        a_.set_title(f"[0, {'2π' if fsm == 10 else '6π'}), fs = {fsm}, n = {n}, Δf = {f(fsm / n)}", fontsize=8)
        a_.set_xlabel("f [Hz]")
    return svg(fig)


FIG["muestras"] = fig_muestras()
tb_pb, _ = V.aplicar(V.TB, V.filtro_pasa_bajo(R["D0"][0], V.TB.shape))
tb_pa, _ = V.aplicar(V.TB, V.filtro_pasa_alto(R["D0"][1], V.TB.shape))
FIG["tab_pb"] = fig_par_img(V.TB, tb_pb, "Imagen original", "Imagen procesada")
FIG["tab_pa"] = fig_par_img(V.TB, tb_pa, "Imagen original", "Imagen procesada")
FIG["tab4"] = fig_tablero4()
FIG["sint"] = fig_sintetica()
FIG["ej5"] = fig_senal(V.tS, V.ej5_orig, V.ej5_pb, "muestra original", "señal procesada")
FIG["ev6"] = fig_senal(V.tS, V.ev6_orig, V.ev6_pa, "muestra original", "señal procesada")
FIG["fs2"] = fig_muestreo_fs2()
FIG["mgrid"] = fig_muestreo_grid()
FIG["ej_s"] = fig_ej_s()
FIG["rayas"] = fig_rayas()
FIG["ev5"] = fig_ev5()
FIG["ev8"] = fig_ev8()
FIG["tiempos"] = fig_tiempos()
FIG["escalon"] = fig_escalon()
FIG["alias"] = fig_alias()
FIG["resol"] = fig_resol()


# ================================================================== números que se citan
ES = R["ej_s"]
BC = {n: b for (n, a, b, ex) in R["escalon_coef"]}
GIBBS = R["escalon_gibbs"][-1][2]
X4 = R["tdf4_X"]
def cplx(re, im):
    if abs(im) < 1e-9: return f(re, 4)
    return f"{f(re, 4)} {'−' if im < 0 else '+'} {f(abs(im), 4)}j"
X4s = [cplx(*z) for z in X4]
C = R["costos"]
MU = {d["fs"]: d for d in R["muestreo"]}
RES = R["resolucion"]
M1, M2 = R["muestras"]["m1"], R["muestras"]["m2"]
TI = {d["n"]: d["ms"] for d in R["tiempos"]}
MULT = {m[0]: m for m in R["mults"]}
def fact(d):
    sup = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    return "·".join(f"{p}{str(e).translate(sup) if e > 1 else ''}" for p, e in d.items())
FACT = {d["n"]: fact(d["fact"]) for d in R["tiempos"]}
# los textos sobre tiempos asumen este orden (se mide en cada corrida: si cambia, revisar el texto)
assert TI[65537] > 2 * TI[65536] and TI[131072] > TI[65536], TI
ANIOS = f(C[2]["anios_N2"], 1)
SEG_FFT = f(C[2]["seg_NlogN"], 1)

S = []
def sec(i, titulo, estrellas, pags, *html):
    S.append({"id": f"m3-s{i}", "titulo": titulo, "estrellas": estrellas, "pags": pags, "html": "".join(html)})

# ------------------------------------------------------------------ s1
sec(1, "Aproximación polinómica trigonométrica: de mínimos cuadrados a senos y cosenos", 2, "p2–5",
    P("Cuando los datos se repiten periódicamente, una recta o un polinomio no los representa bien: hace falta otra familia de funciones base. La idea de la unidad es usar el mismo esquema del ajuste por mínimos cuadrados, pero con senos y cosenos."),
    H4("El punto de partida (repaso de MNI, p2)"),
    P("Con una nube de datos se puede <b>interpolar</b> (la curva pasa por todos los puntos) o <b>ajustar</b> (la curva pasa cerca, minimizando el error). El ajuste polinomial por mínimos cuadrados busca"),
    CMD("F = α₁f₁(x) + α₂f₂(x) + … + αₙfₙ(x)"),
    P("Con F = α₁f₁(x) + α₂f₂(x) (y α₁ = 1 en la diapositiva) es el <b>ajuste lineal</b>. Las incógnitas son los coeficientes α; las fᵢ son la base elegida."),
    H4("Datos periódicos (p3–4)"),
    P("La p3 muestra datos que oscilan y pregunta <i>“En este caso, cual sería un mejor ajuste? Los valores parecen ser periodicos (al cabo de ciertos valores se repiten periodicamente)”</i>. La respuesta de la p4 es elegir funciones trigonométricas como base:"),
    CMD("F = sin(x) + 1.2·cos(2x)      α₁ = 1, f₁ = sin(x)      α₂ = 1.2, f₂ = cos(2x)"),
    P("“Se pueden aproximar los datos con una suma de funciones senoidales! (polinomio de funciones trigonométricas)”."),
    D("Teorema (p4–5)", "Cualquier función continua y periódica puede ser aproximada mediante la suma de senos y cosenos", "U3 p4"),
    H4("¿Cuántas sumas para el escalón? (p5)"),
    P("La p5 muestra sumas parciales cada vez más parecidas a una onda cuadrada y pregunta <i>“Cuantas sumas se requieren para aproximar de manera exacta la función escalon?”</i>. Respuesta: <b>infinitas</b>. Con una cantidad finita de términos siempre queda error, sobre todo en los saltos."),
    FIGHTML(FIG["escalon"]),
    AGR(f"Cuenta con verif_u3.py para el escalón sign(sin x) en [−π, π]: a₀ = {f(R['escalon_a0'], 4)}, aₙ = 0 y bₙ = 4/(nπ) para n impar y 0 para n par: b₁ = {f(BC[1], 4)}, b₃ = {f(BC[3], 4)}, b₅ = {f(BC[5], 4)}, b₇ = {f(BC[7], 4)}. El escalón es discontinuo, así que no cumple la hipótesis del teorema (continua): la serie converge en los puntos de continuidad, pero junto al salto la suma parcial se pasa del valor 1 en ≈{f(GIBBS, 3)} aunque se agreguen términos (fenómeno de Gibbs; con 1 término el máximo es {f(R['escalon_gibbs'][0][1], 4)}, con 10 es {f(R['escalon_gibbs'][4][1], 4)}, con 200 es {f(R['escalon_gibbs'][-1][1], 4)})."),
    H4("Polinómica lineal vs polinómica trigonométrica (WA25)"),
    TBL(["", "Aproximación polinómica (lineal o de grado n)", "Aproximación polinómica trigonométrica"], [
        ["Base fᵢ(x)", "1, x, x², …", "1, cos(nω₀x), sin(nω₀x)"],
        ["Forma", "F = α₁f₁(x) + … + αₙfₙ(x)", "La misma forma: suma de senos y cosenos (serie de Fourier)"],
        ["Datos adecuados", "Tendencias (crecen, decrecen, curvan)", "Datos que se repiten periódicamente (señales)"],
        ["Incógnitas", "Los coeficientes α (ecuaciones normales)", "a₀, aₙ y bₙ (integrales de la función, p6)"],
        ["Fuera del intervalo", "El polinomio diverge", "La aproximación se repite con período T"],
    ]),
    TRAMPA("WA25 pidió la <i>“Diferencia aproximación polinómica lineal vs polinómica trigonométrica”</i>. Lo que suma puntos: las dos son combinaciones lineales F = Σαᵢfᵢ(x) y lo que cambia es la base; la trigonométrica sirve <b>porque</b> los datos son periódicos y el teorema garantiza que una función continua y periódica se aproxima con senos y cosenos."),
)

# ------------------------------------------------------------------ s2
sec(2, "Serie de Fourier: coeficientes, período T y forma exponencial", 2, "p6–9",
    P("La serie de Fourier es la forma concreta de la aproximación trigonométrica: dice qué senos y cosenos usar y cómo calcular sus coeficientes."),
    D("Serie de Fourier (p6)", "Si la función a la que queremos aproximar es continua y periódica en R, […] La serie de Fourier de una función f(x) viene dada por: f(x) = ½a₀ + Σₙ₌₁^∞ aₙcos(nx) + Σₙ₌₁^∞ bₙsin(nx)", "U3 p6"),
    P("“Las funciones f₁ = cos(x) y f₂ = sin(x) son ortogonales en [−Π, Π]”. “Como en cualquier tipo de aproximación, el objetivo es encontrar los coeficiente a₀, aₙ y bₙ”:"),
    CMD("a₀ = (1/Π) ∫₋Π^Π f(x) dx\naₙ = (1/Π) ∫₋Π^Π f(x)·cos(nx) dx\nbₙ = (1/Π) ∫₋Π^Π f(x)·sin(nx) dx"),
    AGR(f"La p6 dice “El preriódo de estas funciones es T=Π”. El período de sin(x) y cos(x) es <b>2π</b> (= {f(R['periodo_sin'], 4)}); por eso se integra en [−π, π], que mide 2π. La ortogonalidad: ∫₋π^π sin(x)cos(x)dx = {f(R['orto'], 4)} (verif)."),
    H4("Preguntas de la p6"),
    UL("<b>¿Cuántos coeficientes para aproximar de manera exacta cualquier f?</b> Infinitos: a₀ y todos los aₙ, bₙ con n = 1, 2, 3, … Si se corta en N armónicos quedan 2N + 1 incógnitas (a₀, a₁…a_N, b₁…b_N) y la aproximación ya no es exacta (agregado: la cuenta 2N + 1).",
       "<b>¿Qué pasa si el período T ≠ Π?</b> Se generaliza replanteando los límites de integración (p7).",
       "<b>Obs (textual):</b> “Para calcular los coef necesito la función (para calcular las integrales)”. Con datos muestreados no hay función: ahí entra la TDF."),
    H4("Cualquier período T (p7–8)"),
    P("“Generalizar para cualquier periodo T replanteando los límites de integración. Definir ω₀ = 2Πf₀, f₀ = frecuencia fundamental y f₀ = 1/T”."),
    CMD("f(x) = a₀ + Σₙ₌₁^∞ aₙcos(nω₀x) + Σₙ₌₁^∞ bₙsin(nω₀x)      n = 1, 2, 3 …\na₀ = (1/T) ∫₀ᵀ f(x) dx\naₙ = (2/T) ∫₀ᵀ f(x)·cos(nω₀x) dx\nbₙ = (2/T) ∫₀ᵀ f(x)·sin(nω₀x) dx"),
    P("Notá que acá el término constante es a₀ (no ½a₀) porque a₀ ya lleva 1/T. “Se aplica en el procesamiento de señales”, a una <b>señal continua y periódica en el dominio del tiempo</b>: s(t) = s(t + T) ∀T ∈ R⁺."),
    D("Señal senoidal (p7)", "s(t) = A sin(ω₀t + θ); A = amplitud, θ = fase, ω₀ = frecuencia angular, f₀ = 1/T [Hertz][1/seg]", "U3 p7"),
    EJ(f"p8: s(t) = sin(t) → A = 1, ω₀ = 2πf₀ = 1, f₀ = 1/2π (= {f(ES['f1'], 4)} Hz), θ = 0. “Obs2: 1 señal de 1 Hz es mas lenta que 1 señal de 2 Hz (demora mas en completar el ciclo)”. “Si T&lt; entonces F&gt;”."),
    H4("Forma exponencial (p9)"),
    P("“Otra manera de representar una señal”: e^{jθ} = cosθ + j sinθ, y las fórmulas de Euler"),
    CMD("cos nω₀t = (e^{jnω₀t} + e^{−jnω₀t}) / 2\nsin nω₀t = (e^{jnω₀t} − e^{−jnω₀t}) / 2      (así en la p9)\n\nf(t) = a₀ + Σ aₙ·(e^{jnω₀t} + e^{−jnω₀t})/2 + Σ bₙ·(e^{jnω₀t} − e^{−jnω₀t})/2"),
    AGR("En la fórmula del seno falta la j: sin nω₀t = (e^{jnω₀t} − e^{−jnω₀t}) / (2j). Lo mismo dentro de bₙ y del segundo Σ de la serie de la p9. Esta forma exponencial es el puente hacia la transformada, que usa e^{±jωt}."),
    TRAMPA("IT25-T2 P2 b: <i>“¿Cuáles y cuántas son las incógnitas que se deben encontrar para aplicar la serie de Fourier? Explicar.”</i> Respuesta que suma: las incógnitas son a₀, aₙ y bₙ; para una representación exacta son infinitas (n = 1…∞), <b>porque</b> cada armónico nω₀ aporta un coseno y un seno; en la práctica se trunca y se calculan con las integrales, para lo cual hace falta conocer f."),
)

# ------------------------------------------------------------------ s3
sec(3, "Interpretación de la serie: dominio de la frecuencia y espectro", 3, "p10–13",
    P("Calcular los coeficientes no es el fin: lo que importa es leerlos como <b>qué frecuencias tiene la señal y con cuánta amplitud</b>. Eso es pasar del dominio del tiempo al dominio de la frecuencia."),
    P("La p10 muestra la misma señal vista de dos lados: <b>dominio del tiempo (amplitud vs tiempo)</b> → SF → <b>dominio de la frecuencia (amplitud vs frecuencia)</b>."),
    D("Espectro de frecuencia (p10)", "Esta representación permite discriminar las frecuencias presentes en dicha señal (o función), su distribución y la energía presente en cada frecuencia. Esta composición en frecuencia se denomina espectro de frecuencia. El espectro resultante al aplicar serie de Fourier es un espectro discreto (un valor de amplitud por cada frecuencia presente)", "U3 p10"),
    H4("Ejemplo de la p12–13: s(t) = sin(t) + 1.2·cos(2t)"),
    P("Se pide amplitud, fase y frecuencia fundamental de cada componente y cuál es más “rápida”. La resolución de la p13:"),
    TBL(["Componente", "Amplitud", "Fase", "ω", "f₀ = ω/2π", "Período"], [
        ["sin(t)", "A₁ = 1", "θ₁ = 0", "1", f"1/2Π = {f(ES['f1'], 4)} Hz", f"2π = {f(ES['T1'], 4)} s"],
        ["1.2·cos(2t)", "A₂ = 1.2", "θ₂ = 0", "2", f"1/Π = {f(ES['f2'], 4)} Hz", f"π = {f(ES['T2'], 4)} s"],
    ]),
    P("La más rápida es 1.2·cos(2t): tiene el doble de frecuencia (período más corto). El espectro tiene <b>dos barras</b>: altura 1 en 1/2π Hz y altura 1,2 en 1/π Hz."),
    FIGHTML(FIG["ej_s"]),
    P(f"(verif: la FFT de esa señal da picos en {f(ES['picos_fft'][0][0])} Hz con amplitud {f(ES['picos_fft'][0][1], 2)} y en {f(ES['picos_fft'][1][0])} Hz con amplitud {f(ES['picos_fft'][1][1], 2)}.)"),
    H4("Ventajas de pasar al dominio de la frecuencia"),
    UL("Se ven directamente <b>qué frecuencias</b> hay, cómo se distribuyen y cuánta energía lleva cada una (p10).",
       "Una señal que en el tiempo parece complicada se reduce a pocos números (amplitud por frecuencia).",
       "Permite <b>filtrar</b>: anular frecuencias (ruido de alta frecuencia, tendencia lenta) y volver al tiempo con la inversa (p23–24, p27).",
       "Aplicaciones de la p24–25: procesamiento de señales e imágenes, restauración, segmentación, reconocimiento, compresión."),
    TRAMPA("IT25-T3 P2 c: <i>“¿Cuáles son las ventajas de pasar del dominio del tiempo al dominio de la frecuencia?”</i> y b: <i>“¿Cuál es la interpretación de la Serie de Fourier? ¿y de la transformada rápida de Fourier?”</i>. La interpretación pedida es la de la p10 (discriminar frecuencias, distribución y energía), y hay que decir <b>de qué tipo es el espectro</b>: discreto en la serie, discreto (N valores, k = 0…N−1) en la FFT."),
)

# ------------------------------------------------------------------ s4
sec(4, "Transformadas: señal no periódica y transformada continua de Fourier", 3, "p14–15",
    P("La serie exige periodicidad. Para una señal continua que no se repite, la p14 propone: <b>“Suponemos que el período es infinito!”</b>. Eso lleva a la transformada continua."),
    D("Transformada continua de Fourier (p15)", "F(ω) = ∫₋∞^∞ v(t) e^{jωt} dt. Se aplica a señales (o funciones) contínuas pero no periódicas. Se consideta T como ∞. El resultado es un espectro contínuo. La transformada se realiza mediante el cálculo de una integral y no de una suma", "U3 p15"),
    P("D. Tpo v(t) → D. Frec F(ω). El ejemplo de la p15 es un pulso (forma de campana acotada) cuya transformada es una curva continua oscilante con máximo en 0. Sobre ese espectro la p15 anota: <b>“Espectro es simétrico respecto a x=0. Espectro es contínuo (T=∞). Cómo se interpreta este espectro?”</b>"),
    H4("Cómo se interpreta el espectro continuo"),
    UL("Como T → ∞, la frecuencia fundamental f₀ = 1/T → 0: las “barras” de la serie quedan infinitamente juntas y el espectro pasa de discreto a <b>continuo</b>: hay contenido en <b>todas</b> las frecuencias, con amplitud F(ω).",
       "Un pico alto en ω indica que esa frecuencia pesa mucho en la señal; donde F(ω) ≈ 0 esa frecuencia casi no está.",
       "Simetría respecto de 0: para una señal real, |F(−ω)| = |F(ω)|; la mitad negativa no agrega información (agregado; por eso en Python se grafican sólo las f &gt; 0, ver nyquist p10)."),
    AGR("La p15 escribe e^{+jωt}; la convención más usada es F(ω) = ∫ v(t)e^{−jωt}dt (con −j, igual que la TDF de la p16). Cambiar el signo sólo refleja el espectro (ω → −ω); el módulo |F(ω)| es el mismo. Si te piden escribirla, copiá la de la cátedra y aclarás la convención."),
    TRAMPA("IT25-T1 P2 a: <i>“¿Cuál es la principal diferencia entre la transformada discreta de Fourier y la transformada continua de Fourier?”</i>. La TCF trabaja sobre una función continua (integral de −∞ a ∞) y da un espectro continuo; la TDF trabaja sobre N muestras (una suma de N términos) y da un espectro discreto de N valores. Decí las dos cosas: <b>qué entra</b> (continua vs discreta) y <b>qué sale</b> (espectro continuo vs discreto)."),
)

# ------------------------------------------------------------------ s5
sec(5, "Transformada discreta de Fourier (TDF) y su costo", 3, "p16–19",
    P("En la computadora no hay funciones continuas: hay N muestras x[n]. La TDF es la versión que se calcula con una suma finita, y el problema que aparece es su costo."),
    D("Transformada discreta de Fourier (p16)", "X[k] = Σₙ₌₀^{N−1} x[n] e^{−j2πk n/N}, k = 0,1,…, N-1. Se aplica a señales discretas x[n]. No es necesaria la periodicidad. N es la cantidad total de puntos que tiene la señal (función). Cuanto &gt;&gt; N, &gt;&gt; resolución en frecuencia. El resultado es un espectro discreto", "U3 p16"),
    P("“Otra forma de escribir”: X[k] = Σₙ₌₀^{N−1} x[n]·ω_N^{nk}, con <b>ω_N = e^{−j2π/N}</b>. Las p16–17 muestran que vale tanto para una señal “Discrete, Non-Periodic” como “Discrete and Periodic”: el espectro X[k] es siempre <b>discreto</b> (N valores)."),
    AGR("La p16 termina la frase como “un espectro discreto y no de una suma” (arrastrada de la p15); lo que quiere decir es que el resultado es discreto y se calcula con una suma, no con una integral."),
    H4("Costo (p18–19), con los datos de la diapositiva"),
    P("N = 4, x[n] = 1.2, 1.4, 1.3, 1.0. Cada X[k] = Σₙ₌₀³ x[n]ω₄^{nk} es una incógnita. Con ω₄^m = e^{−j2πm/4} ∈ {1, −j, −1, j}:"),
    TBL(["k", "X[k] = Σ x[n]·ω₄^{nk}", "Resultado (verif)"], [
        ["0", "1.2 + 1.4 + 1.3 + 1.0", X4s[0]],
        ["1", "1.2 + 1.4(−j) + 1.3(−1) + 1.0(j)", X4s[1]],
        ["2", "1.2 − 1.4 + 1.3 − 1.0", X4s[2]],
        ["3", "1.2 + 1.4(j) + 1.3(−1) + 1.0(−j)", X4s[3]],
    ]),
    P(f"np.fft.fft da lo mismo ({'coincide' if R['tdf4_ok'] else 'NO coincide'}). “Supongamos que tenemos calculados y almacenados los coeficientes ω_N. Por cada X[k] necesitamos hacer N-1 sumas y N multiplicaciones complejas” → N valores de k × N multiplicaciones = <b>N²</b>: <b>Costo de O(N²)</b>. Para N = 4: {R['tdf4_ops']['mult']} multiplicaciones y {R['tdf4_ops']['sumas']} sumas."),
    TBL(["N", "1000", "10⁶", "10⁹"], [
        ["N²", "10⁶", "10¹²", "10¹⁸"],
        ["Tiempo a 1 ns/operación", f"{f(C[0]['seg_N2'], 3)} s", f"{f(C[1]['seg_N2'], 0)} s", f"{f(C[2]['seg_N2']/1e9, 0)}×10⁹ s ≈ {ANIOS} años"],
    ]),
    AGR(f"La p19 dice “10¹⁸ ns ~ 31.2 años”; la cuenta da 10⁹ s ≈ {ANIOS} años. El orden de magnitud (décadas) es lo que importa."),
    TRAMPA("IT25-T1 P2 c: <i>“¿Cuál es el tiempo de ejecución de la FFT? ¿y de la TDF tradicional?”</i> TDF directa: O(N²) <b>porque</b> cada uno de los N coeficientes X[k] necesita N multiplicaciones complejas y N−1 sumas. FFT: O(N log₂N). Dar la tabla con N = 10⁹ (décadas vs ~30 s) muestra que entendiste por qué importa."),
)

# ------------------------------------------------------------------ s6
cien = MULT[1024]
sec(6, "FFT: el algoritmo de Cooley–Tukey", 3, "p20–22",
    P("La TDF directa es inviable para N grande. La FFT calcula <b>exactamente la misma TDF</b> pero reutilizando cuentas, y baja el costo de N² a N log₂N."),
    D("FFT (p20)", "Algoritmo (FFT) para implementar la TDF. Desarrollado por Cooley y Tukey (1965). Baja el costo a N log₂(N) donde N es pot 2. Técnica D&amp;C", "U3 p20"),
    TBL(["N", "1000", "10⁶", "10⁹"], [
        ["N²", "10⁶", "10¹²", "10¹⁸ ns ≈ 31 años"],
        ["N log₂N (p20)", "10⁴", "20×10⁶", "30×10⁹ ns ~30 seg"],
        ["N log₂N (verif)", f(C[0]["NlogN"], 0), f(C[1]["NlogN"], 0), f"{f(C[2]['NlogN'], 0)} → {SEG_FFT} s"],
        ["Mejora N²/(N log₂N)", f(C[0]["cociente"], 0), f(C[1]["cociente"], 0), f(C[2]["cociente"], 0)],
    ]),
    H4("Qué aprovecha (p21)"),
    UL("<b>Periodicidad:</b> “ω_N^{nk} es periodico, con periodo fundamental N con respecto a n, y con respecto a k”: ω_N^{kn} = ω_N^{k(n+N)} = ω_N^{(k+N)n} (la diapositiva escribe m: debe ser múltiplo de N).",
       "<b>Simetría (complejos conjugados):</b> ω_N^{k(N−n)} = ω_N^{−kn} = (ω_N^{kn})*.",
       "<b>Identidad:</b> ω_N^{kN} = (e^{−j2π/N})^{kN} = e^{−j2πk} = 1."),
    H4("División en pares e impares (p22)"),
    P("Con “N (cant puntos) potencia de 2” se separa la suma según n par (2n) e impar (2n+1):"),
    CMD("X[k] = Σₙ₌₀^{N/2−1} [ x[2n]·ω_N^{k(2n)} + x[2n+1]·ω_N^{k(2n+1)} ]\n\nω_N^{k(2n)}   = e^{−j2π(2kn)/N} = e^{−j2πkn/(N/2)} = ω_{N/2}^{kn}\nω_N^{k(2n+1)} = ω_N^k · ω_{N/2}^{kn}\n\nX[k] = Σₙ₌₀^{N/2−1} x[2n]·ω_{N/2}^{kn} + ω_N^k · Σₙ₌₀^{N/2−1} x[2n+1]·ω_{N/2}^{kn}\n     = Σ pares + ω_N^k · Σ impares"),
    P("Cada Σ es una TDF de <b>N/2</b> puntos. Se repite la división (divide y vencerás) hasta TDF de 1 punto: el diagrama de la p22 lo muestra para 16 puntos: 1 señal de 16 → 2 de 8 → 4 de 4 → 8 de 2 → 16 de 1, es decir log₂16 = 4 niveles, con N operaciones por nivel → N log₂N."),
    AGR(f"La p22 escribe los exponentes como “k(2k)”; es k(2n). Verificado numéricamente con N = 8: X[k] = E[k] + ω₈^k·O[k] reproduce np.fft.fft ({SI(R['ct_ok'])}), y por periodicidad + simetría X[k + N/2] = E[k] − ω_N^k·O[k] ({SI(R['ct_mitad'])}): con una sola multiplicación ω_N^k·O[k] salen dos coeficientes. Contando multiplicaciones complejas en una FFT radix-2 recursiva con N = 1024: {cien[2]} (= N/2·log₂N = {cien[3]}) contra N² = {cien[1]} de la TDF directa."),
    H4("Ventajas de Cooley–Tukey (F18-P5)"),
    UL("Costo O(N log₂N) en lugar de O(N²): con N = 10⁹ pasa de décadas a segundos.",
       "Es <b>exacta</b>: da la misma TDF (no es una aproximación), sólo cambia el orden de las cuentas.",
       "Divide y vencerás: subproblemas iguales e independientes (recursiva, paralelizable, se puede hacer en el lugar).",
       "Con la FFT inversa (IFFT) se vuelve al tiempo con el mismo costo, lo que hace práctico el filtrado."),
    TRAMPA("IT25-T4 P2 c: <i>“¿En qué se basa la FFT para bajar el costo computacional?”</i> Las palabras que hay que escribir: <b>divide y vencerás</b> (pares/impares, dos TDF de N/2), <b>periodicidad</b> y <b>simetría</b> de ω_N^{nk}, y el resultado <b>O(N log₂N)</b>. F18-P5 pidió las ventajas del algoritmo de Cooley–Tukey."),
)

# ------------------------------------------------------------------ s7
sec(7, "N no potencia de 2: relleno con ceros, e IFFT", 3, "p23",
    P("La división pares/impares necesita que N se pueda partir a la mitad hasta llegar a 1. La p23 pregunta qué hacer si no es así, y si se puede volver al tiempo."),
    P("<b>“Que pasa si N (cant de puntos) no es potencia de 2?”</b> La respuesta de la diapositiva es el árbol: a la señal {x₀,…,x₇} se le agregan ceros hasta 16 puntos, {x₀,…,x₇,0,0,0,0,0,0,0,0}, y se sigue dividiendo en pares e impares: <b>relleno con ceros (zero padding)</b> hasta la siguiente potencia de 2."),
    TBL(["N original", "Siguiente potencia de 2"], [[str(a), str(b)] for a, b in R["pad"]]),
    H4("Qué cambia y qué no con el padding"),
    UL("Cambia N (y por lo tanto la grilla de frecuencias: el paso pasa a f_s/N′).",
       "<b>No cambia</b> f_s ni la frecuencia de Nyquist f_s/2: no genera ni corrige aliasing.",
       "No agrega información: los ceros no traen frecuencias nuevas; el espectro queda más “interpolado”."),
    P("<b>“Puedo volver desde el dominio de la frecuencia al dominio del tiempo?”</b> Sí: x[n] —FFT→ X(k) —IFFT→ x[n]. "
      f"Con los datos de la p18, ifft(fft(x)) devuelve {'; '.join(f(v, 1) for v in R['tdf4_ifft'])}."),
    H4("¿Y la fft de NumPy/SciPy? (TP3 ej. 2 y 6)"),
    AGR(f"NumPy y SciPy (pocketfft) <b>no exigen</b> potencia de 2 ni rellenan solas: aceptan cualquier N, descomponiendo N en factores pequeños (Cooley–Tukey de base mixta 2, 3, 4, 5, …) y, si queda un factor primo grande, usan otro algoritmo (Bluestein) que sigue siendo O(N log N) pero con más constante. Sólo rellenan si se lo pedís con <code>np.fft.fft(x, n)</code>. Tiempos medidos acá con NumPy {R['numpy']} (dependen de la máquina):"),
    TBL(["N", "Factorización", "ms por fft"], [[str(d["n"]), FACT[d["n"]], f(d["ms"], 2)] for d in R["tiempos"]]),
    FIGHTML(FIG["tiempos"]),
    P(f"Lectura: 2¹⁶ = 65536 es el más rápido; 65537 (un punto más, pero <b>primo</b>) tarda unas {f(TI[65537]/TI[65536], 1)} veces más. 166400 = {FACT[166400]} y 166464 = {FACT[166464]} son más grandes que 2¹⁷ pero sólo tienen factores chicos, así que su tiempo queda en el orden de 2¹⁷ (entre ellos dos la diferencia es chica y cambia de una corrida a otra). El tiempo depende de la factorización de N, no sólo de su tamaño."),
    TRAMPA("IT25-T1 P2 b: <i>“Dada la funcionalidad de fft provista por Scipy, ¿es necesario usar un número de puntos en potencia de dos para poder usarla? ¿Qué ocurre si no son potencia de dos?”</i> No es necesario; funciona con cualquier N, pero es más lenta cuanto peores son los factores de N (caso extremo: N primo). Mencioná además lo de la p23: el Cooley–Tukey radix-2 de la clase sí la necesita, y la salida es el relleno con ceros. WA25 V/F “FFT (n log n) sólo para n potencia de 2”: <b>Falso</b> por el “sólo”. ET25-T1 usó como distractor que el padding “generó aliasing, ya que modificó la frecuencia de Nyquist”: el padding no toca f_s."),
)

# ------------------------------------------------------------------ s8
E5, E6 = R["ej5"], R["EV_P6"]
sec(8, "Aplicaciones: procesamiento de señales y filtros", 3, "p24",
    P("El uso práctico de la FFT es modificar una señal en el dominio de la frecuencia (quitar ruido, quedarse con una banda) y volver al tiempo."),
    P("La p24 (captura de un programa de filtrado) lista: <b>Procesamiento digital de señales; Análisis de Fourier; Filtros (pasa baja, pasa alta, pasa banda)</b>. Se ve una señal ruidosa, su espectro de potencia con unos pocos picos, y la señal filtrada."),
    H4("Receta del filtrado 1D"),
    OL("Y = fft(y) y freqs = fftfreq(n, dt).",
       "Máscara: <b>pasa bajo</b> = 1 donde |f| &lt; f_c, 0 en el resto; <b>pasa alto</b> = 1 donde |f| ≥ f_c; <b>pasa banda</b> = 1 sólo entre f₁ y f₂.",
       "y_filtrada = ifft(Y·máscara).real."),
    TBL(["Original", "Procesada", "Filtro", "Por qué"], [
        ["Oscilación lenta + ruido fino", "Curva suave, misma forma lenta", "Pasa bajo", "El ruido es alta frecuencia; se conserva lo lento"],
        ["Oscilación lenta + oscilación rápida", "Sólo la oscilación rápida, centrada en 0", "Pasa alto", "Se quitó la componente lenta (baja frecuencia)"],
        ["Señal con muchos picos", "Una sola banda", "Pasa banda", "Se conserva un rango [f₁, f₂]"],
    ]),
    EJ(f"Recreado en verif_u3.py: señal lenta (0,15 y 0,4 Hz) + ruido de desvío {f(E5['std_ruido'], 3)}. Pasa bajo con f_c = {f(E5['fc'])} Hz: el error contra la señal lenta queda en {f(E5['err_pb'], 3)} (RMS). Con pasa alto el resultado se aleja ({f(E5['err_pa'], 3)}): queda casi sólo el ruido. Otra: 2·sin(2π·0,12t) + sin(2π·3t); pasa alto con f_c = 1 Hz deja la de 3 Hz con amplitud {f(E6['amp_pa'], 2)}; pasa bajo dejaría la lenta con amplitud {f(E6['amp_pb'], 2)}."),
    AGR("Un <b>filtro antialiasing</b> es un pasa bajo que se aplica a la señal <b>antes de muestrear</b> para quitar lo que esté por encima de f_s/2. Después de muestrear no corrige nada: el alias ya está mezclado con las frecuencias reales. Por eso “se eliminaron las señales de aliasing con un filtro antialiasing” no explica una señal procesada que sólo se suavizó."),
    TRAMPA("ET25-T3 ej. 5 (original lenta con ruido → procesada suave) es <b>pasa bajo</b>; EV-P6 (lenta + rápida → sólo la rápida) es <b>pasa alto</b>. Los distractores son la operación inversa, “filtro antialiasing” y “padding de ceros para modificar la frecuencia de Nyquist” (el padding no filtra ni cambia f_N)."),
)

# ------------------------------------------------------------------ s9
FB = R["frac_baja"]
sec(9, "Imágenes: imagen digital y espectro 2D", 3, "p25–26, N9",
    P("Una imagen también es una señal, en dos dimensiones. Su espectro 2D dice si tiene zonas suaves, bordes o patrones repetidos."),
    P("<b>Procesamiento digital de imágenes</b> (p25): Análisis de Fourier; Filtros (pasa baja, pasa alta, etc); Operaciones de restauración (interpolación, etc); Cambios en la morfología (erosión, dilatación, etc); Segmentación (detección de puntos, bordes, etc); Reconocimiento de objetos; Compresión de imágenes."),
    D("Imagen digital (p25)", "Imagen Digital: función f(x,y) de coordenadas espaciales", "U3 p25"),
    P("La p25 muestra imagen analógica → muestreo digital → cuantización de píxeles, y una matriz aleatoria con <code>np.random.random((100, 100))</code> y <code>imshow(X, cmap=\"gray\")</code>."),
    H4("Espectro 2D (p26)"),
    P("Dominio espacio 2D ⇄ FFT/IFFT ⇄ dominio frecuencia. El esquema de la p26 es el espectro <b>sin centrar</b>, como sale de fft2: (0,0) en la esquina superior izquierda, las <b>bajas</b> frecuencias en las cuatro esquinas y las <b>altas</b> en el centro. <code>np.fft.fftshift</code> reordena los cuadrantes para que la frecuencia cero quede en el centro (es como se grafica casi siempre); <code>np.fft.ifftshift</code> deshace ese reordenamiento."),
    P("“En el ejemplo: Las áreas más claras en el espectro muestran mayor energía (&gt; amplitud). Se observa mucha más energía concentrada en las frecuencias más altas” (se refiere a esa imagen hecha de trazos finos)."),
    D("Lectura del espectro de una imagen (p26)", "Zonas homogéneas en la imagen darán lugar a que la energía del espectro esté concentrada mayoritariamente en las bajas frecuencias. Zonas con muchos bordes, transiciones frecuentes de niveles de gris dará lugar a un espectro con componentes de alta frecuencia. Si en la imagen existe alguna regularidad (patrones que se repiten) dará lugar a picos de igual intensidad separados a una distancia regular", "U3 p26"),
    EJ(f"verif_u3.py: una mancha gaussiana suave tiene el {f(100*FB['suave'], 1)} % de la energía (sin la media) dentro de radio 8 del centro; el tablero sintético, sólo el {f(100*FB['tablero'], 1)} %: el resto son picos regulares lejos del centro. Unas franjas verticales de período 16 px (128 px de ancho) dan dos picos sobre el eje horizontal, en ±128/16 = ±8 (verif: picos en {PK(R['rayas_picos'])}): el espectro marca la dirección <b>perpendicular</b> a las franjas."),
    FIGHTML(FIG["rayas"]),
    H4("Preguntas para hacerse (nyquist N9)"),
    UL("<b>¿Qué hace fft2?</b> La TDF en dos dimensiones de la matriz f(x,y): da una matriz de coeficientes complejos indexada por dos frecuencias (horizontal y vertical).",
       f"<b>¿Podría replicar fft2 usando solamente fft?</b> Sí: fft a cada fila y después fft a cada columna del resultado (la TDF 2D es separable). Verificado: <code>np.fft.fft(np.fft.fft(A, axis=1), axis=0)</code> = <code>np.fft.fft2(A)</code> → {SI(R['fft2_filas_cols'])}."),
    TRAMPA_P("TP3 ej. 4 empareja imágenes con espectros. Regla: picos alineados en ciertas direcciones = bordes rectos perpendiculares a esas direcciones (un triángulo da una estrella de 3 rectas; un código QR, una cuadrícula horizontal/vertical); energía concentrada en el centro y casi nada afuera = imagen suave (una foto con niebla). ET25-T2 usó el distractor “la imagen de la derecha es el espacio de frecuencias”: un espectro centrado tiene simetría central y un máximo en el centro; una imagen filtrada no."),
)

# ------------------------------------------------------------------ s10
sec(10, "Filtrado de imágenes: pasa bajo y pasa alto", 3, "p27 + imagen_con_filtro.ipynb",
    P("Filtrar una imagen es multiplicar su espectro por una máscara y volver al espacio. Lo que se borra del espectro determina si la imagen queda suavizada o sólo con bordes."),
    D("Pasa baja (p27)", "Los filtros pasa-baja atenúan las componentes de medias-altas frecuencias y dejan intactas las bajas en función a la frecuencia de corte que se elija. Se usan para eliminar ruido de alta frecuencia, o eliminar todo lo que no sean variaciones suaves del nivel de gris.", "U3 p27"),
    D("Pasa alta (p27)", "Los filtros pasa-alta atenúan las componentes de baja frecuencia y dejan intactas las de medias-altas en función de la frecuencia de corte que se elija. Se usan para quedarnos con las propiedades de la imagen en los que los niveles de gris varían bruscamente, por los bordes de la imagen.", "U3 p27"),
    P("En el ejemplo de la p27 el “Filtro” es el espectro centrado con un <b>cuadrado negro en el centro</b> (se anulan las bajas) y la “Transformada inversa” muestra sólo contornos. A la pregunta “¿Que filtro se uso? Pasa alta o pasa baja?” la respuesta es <b>pasa alta</b>."),
    H4("El notebook imagen_con_filtro.ipynb (TP3 ej. 5)"),
    CMD("img_c2 = np.fft.fft2(img)                  # espectro de frecuencias\nimg_c3 = np.fft.fftshift(img_c2)           # espectro centrado (f = 0 al centro)\nimg_c4 = np.fft.ifftshift(img_c3 * filtro) # deshace el centrado\nimg_c5 = np.fft.ifft2(img_c4)              # vuelve al espacio (complejo)\n# se grafican np.log(1+np.abs(espectro)) y np.abs(img_c5)\n\ndef filtro_pasa_bajo(D0, imgShape):        # dado: 1 si distancia al centro &lt; D0\n    ...\ndef filtro_pasa_alto(D0, imgShape):\n    return 1 - filtro_pasa_bajo(D0, imgShape)   # o: 1 donde distancia &gt;= D0"),
    P("El círculo de radio D0 es la <b>frecuencia de corte</b>: un D0 chico en pasa bajo deja sólo lo muy suave (más borroso); un D0 grande en pasa alto deja sólo detalles finos."),
    H4("Tablero sintético (como el de ET25-T2/T3)"),
    P(f"128×128, casillas de 8 px en gris claro y blanco, con diagonales de casillas oscuras. Su espectro centrado tiene picos regulares; los más cercanos al centro están a radio {f(R['tablero_picos'][0][2], 2)} en (±4, ∓4) — las diagonales oscuras — y el damero aparece recién a radio {f(R['tablero_picos'][2][2], 2)}. Con <b>pasa bajo D0 = {R['D0'][0]}</b> sólo sobrevive la componente diagonal: la imagen queda como <b>franjas diagonales suaves</b> (correlación con un coseno diagonal puro = {f(R['tablero_pb_corr_diag'], 3)}, el signo indica que lo oscuro coincide con el valle). Con <b>pasa alto D0 = {R['D0'][1]}</b> se borra el nivel medio y lo suave: queda una imagen oscura con los <b>bordes</b> de las casillas marcados (brillo medio en los píxeles de borde {f(R['tablero_pa_borde_vs_int'][0], 3)} vs {f(R['tablero_pa_borde_vs_int'][1], 3)} en el interior)."),
    FIGHTML(FIG["tab4"]),
    AGR(f"Con padding de ceros (fft2 con s = 256×256) y después ifft2 se recupera la misma imagen con ceros alrededor (verif: {SI(R['pad_img_ok'])}): el padding no suaviza ni da franjas. El pasa alto quita la componente de frecuencia cero (el promedio), por eso el resultado queda oscuro: en una imagen sintética de media {f(R['sint_media'], 3)}, el pasa alto da media {f(R['sint_pa_media'], 3)}."),
    FIGHTML(FIG["sint"]),
    TRAMPA("ET25-T2 ej. 1 y ET25-T3 ej. 1: tablero → franjas diagonales suaves = <b>fft + pasa bajo + ifft</b>. ET25-T3 ej. 2: tablero → sólo los bordes de las casillas = <b>pasa alto</b>. ET25-T3 ej. 3 mostraba la secuencia imagen / espectro / máscara con un disco oscuro en el centro / máscara descentrada (cuartos de disco en las esquinas) / resultado oscuro con contornos = <b>pasa alto</b>. Distractores: la operación inversa, “es el espacio de frecuencias”, “padding de ceros”, “filtro antialiasing”."),
)

# ------------------------------------------------------------------ s11
LF = dict(R["limite_fase"])
sec(11, "Frecuencia de muestreo y frecuencia de Nyquist", 3, "N2–N3",
    P("Al muestrear una señal se pierde lo que pasa entre muestra y muestra. La frecuencia de Nyquist dice hasta qué frecuencia se puede recuperar información."),
    D("Muestreo (nyquist p2)", "f_s: Frecuencia de muestreo (o sampleo). T_s: Periodo de muestreo, T_s = 1/f_s. f_N: Frecuencia de Nyquist, f_N = f_s/2. No podemos obtener información de frecuencias mayores que f_N! Si tenemos una señal con f = 2 Hz, tenemos que tener una f_s de por lo menos 4 Hz (4 muestras en un segundo).", "Nyquist p2"),
    P("<b>De qué depende:</b> sólo de f_s (TP3 ej. 1 a). Para representar una frecuencia f hace falta f_s ≥ 2f, o sea f ≤ f_N. Con f_s = 10 Hz, f_N = 5 Hz (EV-P2)."),
    H4("El ejemplo de la p3"),
    P("Código con <code>tiempo = 6</code>, f = 0,5, 1 y 2 Hz (desfases 0, 0,1·2π y 0), <code>f_nyquist = f_muestreo / 2</code>. Con f_muestreo = 10 (f_N = 5) las tres se ven bien. Con f_muestreo = 2 (f_N = 1): 0,5 Hz aparece como un triángulo (4 muestras por ciclo); 1 Hz = f_N aparece como un zigzag de amplitud ≈0,7 (depende de la fase); 2 Hz = f_s queda plana (≈10⁻¹⁴): desaparece."),
    H4("El caso límite f = f_s/2"),
    P("La cátedra dice “por lo menos 4 Hz” para f = 2 Hz: admite f_s = 2f (f = f_N). Pero en f = f_N hay sólo 2 muestras por ciclo y lo que se ve depende de la fase:"),
    TBL(["Señal a 10 Hz, f_s = 20 Hz", "Máximo |muestra| (verif)"], [[f"sin(2π·10t + {fi}°)", f(v, 3)] for fi, v in R["limite_fase"]]),
    P("Con fase 0 (seno puro) todas las muestras dan 0: la frecuencia se pierde; con fase 90° (coseno) se ve ±1. Por eso, estrictamente, la condición segura es f &lt; f_N (f_s &gt; 2f)."),
    TRAMPA("EV-P8 / ET25-T3 ej. 4: f₁ = 10, f₂ = 20 y f₃ = 100 Hz con f_s = 20 Hz (f_N = 10). Con la convención de la cátedra (f ≤ f_N, “por lo menos”): <b>sólo f₁ = 10 Hz</b>. f₂ = f_s y f₃ = 5·f_s caen en 0 Hz (se ven constantes). Con la lectura estricta (f &lt; f_N) no se reconstruye ninguna, pero esa opción no estaba: marcá “Solo f1” y justificá el límite. ET25-T3 ej. 6: con 10, 15 y 100 Hz la única f_s de las opciones que cumple f_s ≥ 2·100 es <b>200 Hz</b>. EV-P2: “sólo puede representar correctamente frecuencias hasta” <b>f_s/2</b> (distractores 2f_s, f_s², f_s)."),
)

# ------------------------------------------------------------------ s12
sec(12, "Resolución de las frecuencias", 2, "N4–N6, N10",
    P("Además de hasta qué frecuencia se ve (f_N), importa qué tan juntas pueden estar dos frecuencias para distinguirlas: eso es la resolución Δf."),
    D("Resolución (nyquist p4)", "Δf = f_s/n", "Nyquist p4"),
    P("Con n = cantidad de muestras. Como n = f_s·(b − a), queda Δf = 1/(b − a): <b>sólo depende de la duración del registro</b>. Por eso la diapositiva concluye: <b>“Si mantenemos el rango (a, b), aumentar f_s no mejora la resolución!”</b>"),
    TBL(["Intervalo", "f_s", "n", "Δf", "Última f graficada", "f_N"], [
        ["[0, 2π)" if abs(d["b"] - 6.283) < 0.01 else "[0, 6π)", str(d["fs"]), str(d["n"]), f(d["df"], 3), f(d["fmax"], 3), f(d["fN"], 1)] for d in RES]),
    P("Los tres casos de la cátedra (verificados): con f_s = 10 en [0, 2π) n = 63 y Δf = 0,159; con f_s = 100, n = 629 pero Δf sigue en 0,159; recién al alargar a [0, 6π) con f_s = 6 baja a Δf = 0,053 y los picos de 0,5, 1 y 2 Hz quedan finos y en su lugar. La flecha “Que frecuencia es?” apunta al final del eje: es la máxima frecuencia graficable, ≈ f_N (5 Hz y 3 Hz)."),
    FIGHTML(FIG["resol"]),
    P("Código de la p4 y del “Extra” (p10, “solo ploteen los valores positivos de f”):"),
    CMD("y_hat = np.fft.fft(y)\nfreqs = np.fft.fftfreq(n, delta)\ncond =  0 &lt; freqs\nplt.plot(freqs[cond], np.abs(y_hat[cond]))"),
    AGR("fftfreq(n, d) devuelve las frecuencias de cada coeficiente: 0, Δf, 2Δf, … y después las negativas. Para una señal real, la mitad negativa es el espejo de la positiva, por eso alcanza con graficar f &gt; 0. Con la TDF de la p16, “Cuanto &gt;&gt; N, &gt;&gt; resolución” vale si los N puntos cubren más tiempo."),
)

# ------------------------------------------------------------------ s13
AT = R["alias_tabla"]
sec(13, "Aliasing y patrones de Moiré", 3, "N7–N8",
    P("Si una frecuencia supera f_N no desaparece: aparece disfrazada de otra más baja. Ese es el error más tomado en la evaluación práctica."),
    D("Aliasing (nyquist p7)", "La frecuencia se “refleja” en el rango [0, fs/2]", "Nyquist p7"),
    P("El gráfico de la p7 (f_s = 1 kHz) muestra la frecuencia reconstruida vs la original: sube hasta 500 Hz (la frecuencia de Nyquist), baja a 0 en 1 kHz, vuelve a subir hasta 500 en 1,5 kHz, etc. (diente de sierra)."),
    FIGHTML(FIG["alias"]),
    AGR("Fórmula del reflejo (no está escrita en la diapositiva, reproduce su gráfico): f_alias = |f − f_s·round(f/f_s)|. Ejemplos calculados: " + "; ".join(f"{f(a)} Hz con f_s = {f(b)} → {f(c)} Hz" for a, b, c in AT) + "."),
    H4("Cuándo hay aliasing (EV-P3)"),
    P("Cuando la señal tiene frecuencias mayores que f_N, es decir, cuando <b>la frecuencia de Nyquist es menor que la de la señal original</b>. No lo produce calcular la FFT con pocos puntos (eso afecta Δf), ni un filtro pasa bajo (que, aplicado antes de muestrear, lo evita)."),
    H4("Patrones de Moiré (p8)"),
    P("La p8 muestra una pared de ladrillos que, al bajar la resolución de la foto, se llena de ondas que no existen: <b>patrones de Moiré</b>. Es el aliasing en 2D: el patrón de ladrillos tiene una frecuencia espacial mayor que la mitad de la nueva frecuencia de muestreo (píxeles) y aparece reflejado como franjas de baja frecuencia."),
    H4("Las dos muestras del TP3 / ET25-T1"),
    P(f"muestra1.csv: n = {M1['n']}, dt = {f(M1['dt'])} s → f_s = {f(M1['fs'])} Hz, f_N = {f(M1['fN'])} Hz. muestra2.csv: n = {M2['n']}, dt = {f(M2['dt'])} s → f_s = {f(M2['fs'])} Hz, f_N = {f(M2['fN'])} Hz. La muestra 1 tiene un pico en ≈0,8 Hz que la muestra 2 no tiene: es el reflejo de una componente de 1,2 Hz (&gt; f_N = 1) → aliasing en la muestra 1. Detalle en la práctica de ET25-T1."),
    FIGHTML(FIG["muestras"]),
    TRAMPA("ET25-T1 ej. 1 y EV-P5 usan el mecanismo “fenómeno correcto, muestra equivocada”: el pico raro está siempre en la muestra de <b>menor</b> f_s (la de f_N más baja), y es el reflejo de una frecuencia real que la otra muestra sí ve en su lugar. El distractor del padding (“modificó la frecuencia de Nyquist”) es falso: f_N depende sólo de f_s. Y un pico no “se divide en dos” por aliasing."),
)

# ------------------------------------------------------------------ s14
sec(14, "Repaso para el integral: serie vs TCF vs TDF vs FFT", 3, "p6–23",
    P("Las cuatro preguntas de Fourier del integral 2025 piden, con distintas palabras, lo mismo: qué recibe cada herramienta, qué devuelve, cómo se lee su espectro y cuánto cuesta."),
    TBL(["", "Serie de Fourier", "Transformada continua (TCF)", "Transformada discreta (TDF)", "FFT"], [
        ["Entrada", "Función continua y periódica (período T)", "Función continua, no periódica (T → ∞)", "N muestras x[n]; no es necesaria la periodicidad", "Las mismas N muestras"],
        ["Cálculo", "Integrales para a₀, aₙ, bₙ (necesita la función)", "Integral ∫₋∞^∞ v(t)e^{jωt}dt", "Suma Σₙ₌₀^{N−1} x[n]e^{−j2πkn/N}", "La misma suma reorganizada (D&amp;C)"],
        ["Espectro", "Discreto: un valor por armónico nf₀", "Continuo, simétrico respecto de 0", "Discreto: N valores X[k], separados Δf = f_s/N", "Idéntico al de la TDF"],
        ["Costo", "—", "—", "O(N²)", "O(N log₂N)"],
        ["Cuándo", "Señal periódica conocida analíticamente", "Señal no periódica conocida analíticamente", "Señal muestreada (datos)", "Siempre que se calcule una TDF en la computadora"],
    ]),
    H4("Cómo se interpretan los espectros"),
    UL("Serie: barras en f₀, 2f₀, 3f₀…; la altura es la amplitud de cada armónico (ej. p13: barras de 1 y 1,2).",
       "TCF: curva continua; la altura en cada ω es cuánto aporta esa frecuencia; para señal real es simétrica respecto de 0.",
       "TDF/FFT: N valores complejos; |X[k]| a la frecuencia k·f_s/N; sólo tiene sentido hasta f_N = f_s/2 (la otra mitad es espejo); la resolución es Δf = f_s/N; picos anchos o corridos indican poca resolución; picos donde no los esperás pueden ser alias."),
    H4("¿Cualquier función? (IT25-T2 a y c)"),
    UL("Serie: según el teorema de la cátedra, funciones <b>continuas y periódicas</b>; además, para calcular los coeficientes hay que conocer la función (integrales). Una función no periódica va por la TCF.",
       "TDF: se aplica a <b>cualquier</b> conjunto de N muestras (no pide periodicidad), pero sólo describe bien frecuencias por debajo de f_N; si la señal tiene componentes más altas, aparecen alias.",
       "FFT: se aplica a cualquier señal muestreada, porque calcula la TDF. La versión de la clase pide N potencia de 2 (se rellena con ceros si no); las librerías aceptan cualquier N. Ventaja: O(N log₂N) en lugar de O(N²)."),
    TRAMPA("IT25-T4 P2: <i>“a. ¿Cuál es la diferencia entre la TCF, la TDF y la FFT? b. ¿Cómo son y cuál es la interpretación de los espectros obtenidos en cada caso? c. ¿En qué se basa la FFT para bajar el costo computacional?”</i> F20-P4: <i>“en qué casos aplicaría la Serie de Fourier, la TCF y la TDF”</i>. F18-P3: objetivo de la transformada, serie vs transformada, TDF vs FFT. El error que resta puntos: decir que la FFT es “otra transformada” o “una aproximación”. Es un <b>algoritmo</b> que calcula exactamente la TDF."),
)

# ------------------------------------------------------------------ s15
sec(15, "TP3: los notebooks y los ejercicios", 2, "TP3 ej. 1–6",
    P("La evaluación práctica de TP2/TP3 reutiliza los gráficos del TP3. Acá está qué da cada notebook, calculado."),
    H4("muestreo.ipynb (ej. 1)"),
    P("f₁ = 1 Hz y f₂ = 2 Hz con desfase π/2 (o sea cosenos), t ∈ [0, 2]. Al bajar f_s:"),
    TBL(["f_s", "f_N", "Muestras/ciclo de f₁", "Cómo se ve f₁ = 1 Hz", "Muestras/ciclo de f₂", "Cómo se ve f₂ = 2 Hz"], [
        ["60", f(MU[60]["fN"]), f(MU[60]["f1_por_periodo"]), "Perfecta", f(MU[60]["f2_por_periodo"]), "Perfecta"],
        ["10", f(MU[10]["fN"]), f(MU[10]["f1_por_periodo"]), "Bien (poligonal)", f(MU[10]["f2_por_periodo"]), "Reconocible, con picos recortados"],
        ["4", f(MU[4]["fN"]), f(MU[4]["f1_por_periodo"]), "Triángulo (0, ±1)", f(MU[4]["f2_por_periodo"]), "Zigzag ±1: f₂ = f_N, caso límite"],
        ["2", f(MU[2]["fN"]), f(MU[2]["f1_por_periodo"]), "Zigzag ±1: f₁ = f_N, límite", f(MU[2]["f2_por_periodo"]), f"Constante 1: alias a {f(MU[2]['f2_alias'])} Hz"],
        ["1", f(MU[1]["fN"]), f(MU[1]["f1_por_periodo"]), f"Constante 1: alias a {f(MU[1]['f1_alias'])} Hz", f(MU[1]["f2_por_periodo"]), f"Constante 1: alias a {f(MU[1]['f2_alias'])} Hz"],
    ]),
    P(f"Celda 2 (suma de 1, 2 y 5 Hz, f_s = 60, n = {R['muestreo_c2']['n']}): la FFT da picos en {', '.join(f(a, 2) for a, _ in R['muestreo_c2']['picos'])} Hz (Δf = 60/{R['muestreo_c2']['n']} = {f(R['muestreo_c2']['df'], 3)} Hz, por eso no caen exactos en 1, 2, 5)."),
    H4("tiempo_fft.ipynb (ej. 2)"),
    P(f"Mide <code>np.fft.fft(a, len(a))</code> para 166464 = {FACT[166464]}, 166400 = {FACT[166400]}, 2¹⁶, 2¹⁷ y 65537 (primo). Resultado (sección 7): las potencias de 2 son las más rápidas; el primo es el más lento aunque sea chico. Explicación: el algoritmo divide N en sus factores; cuanto más pequeños y repetidos (2), menos trabajo."),
    H4("muestra1 / muestra2 (ej. 3)"),
    P(f"Ambos CSV se reproducen exactamente (error máx. ≈{M2['err_modelo']:.0e}) con y(t) = 5cos(2π·1t) + 5cos(2π·1,2t) + 10sin(2π·3t) + cos(2π·20t) − 2,5sin(2π·0,125t) (encontrada con mínimos cuadrados en verif_u3.py; la cátedra no la da). Muestra 2 (f_s = 100) ve las cinco: 0,125, 1, 1,2, 3 y 20 Hz. Muestra 1 (f_s = 2, f_N = 1) ve 0,125 Hz, un pico en ≈0,8 Hz (alias de 1,2) y uno en el borde ≈1 Hz (1 Hz en f_N; el de 3 Hz cae en 1 Hz pero con fase seno da muestras nulas); el de 20 Hz cae en 0 Hz (suma una constante 1). Mejor muestreo: el 2."),
    H4("imagen_con_filtro.ipynb (ej. 5) y fft en librerías (ej. 6)"),
    P("Completar con fft2 → fftshift → (× filtro) → ifftshift → ifft2, graficar np.log(1+np.abs(·)) de los espectros y np.abs del resultado; pasa alto = 1 − pasa bajo (sección 10). Ej. 6: NumPy implementa la fft con pocketfft (Cooley–Tukey de base mixta + Bluestein para primos grandes, agregado en la sección 7); otras librerías: SciPy (scipy.fft, también pocketfft, con opciones de workers), FFTW (C, planifica el algoritmo según N y la máquina)."),
    AGR("Lo de pocketfft y FFTW no está en el material de la cátedra: es la respuesta al “Averiguar” del ej. 6."),
)


intro = (P("La unidad responde cómo representar una señal (o datos que se repiten) como suma de senos y cosenos, y cómo calcular eso rápido sobre datos muestreados: serie de Fourier → transformada continua → transformada discreta → FFT. Después, cómo muestrear sin perder frecuencias (Nyquist, aliasing) y cómo filtrar señales e imágenes en el dominio de la frecuencia.")
         + P("En el integral teórico siempre hay una pregunta de Fourier (serie vs TCF vs TDF vs FFT, espectros, costo N² vs N log N, potencia de 2). En la evaluación práctica de TP3: elegir entre aliasing / padding / pasa alto / pasa bajo mirando un gráfico de espectros, una señal o una imagen procesada, y elegir f_s."))

conexiones = UL(
    "<b>Con la unidad 2 (Ajuste de curvas):</b> la aproximación trigonométrica es el mismo esquema F = Σαᵢfᵢ(x) con otra base; el “ajuste trigonométrico” es una de las opciones que da ET25-T2 ej. 2.",
    "<b>Con la unidad 4 (EDO):</b> nada directo en el material; la idea de “período de muestreo” T_s se parece al paso h, y elegir h chico vs costo es el mismo compromiso que f_s.",
    "<b>Dentro de la unidad:</b> la serie necesita la función → con datos se usa la TDF → la TDF cuesta N² → la FFT la hace en N log₂N → la FFT sólo ve hasta f_s/2 → por encima, aliasing.")

fidelidad = {
    "fuentes": "trig.pdf (27 diapositivas, todas en imagen, renderizadas y leídas; la p11 es una captura de animación sin texto) y nyquist.pdf (10 diapositivas, renderizadas; p3–6 a 130 dpi para leer el código). TP3 (tp3.pdf, notebooks-TP3.md, muestra1.csv, muestra2.csv). Exámenes según 02-patron-evaluacion.md; gráficos de ET25-T1, EV-P5 y ET25-T3 ej. 3 vistos en parciales.pdf. Todo número sale de tools/gen/verif_u3.py.",
    "agregado": "Coeficientes y Gibbs del escalón; correcciones de erratas (período 2π en p6, falta de j en p9, signo de la TCF en p15, 31,7 años en p19, exponentes k(2n) en p22); fórmula explícita del alias; filtro antialiasing; NumPy/SciPy aceptan cualquier N (pocketfft, Bluestein) y tiempos medidos; fft2 = fft por filas y columnas; la fórmula exacta de la señal de muestra1/muestra2 (hallada por mínimos cuadrados); imágenes y señales sintéticas propias para recrear los ejercicios con gráfico.",
    "huecos": "No hay diapositiva de objetivos. La p11 (animación) no tiene contenido legible más allá del título. La figura de ET25-T1 ej. 1 no coincide en valores con los CSV del TP3 (el mecanismo sí); EV-P5 se resolvió inferiendo f_s de la muestra 1 desde el gráfico. Las fotos del TP3 ej. 4 y del póster de ET25-T3 ej. 3 no se reproducen (terceros).",
}

auditoria = [
    {"tema": "p2–5 Aproximación Polinómica Trigonométrica (ajuste por MC, datos periódicos, teorema, escalón)", "seccion": "m3-s1", "estado": "cubierto"},
    {"tema": "p6–9 Serie de Fourier (coeficientes, período T, señal senoidal, forma exponencial)", "seccion": "m3-s2", "estado": "cubierto"},
    {"tema": "p10–13 Serie de Fourier (interpretación): dominio de la frecuencia, espectro, ejemplo", "seccion": "m3-s3", "estado": "cubierto"},
    {"tema": "p14 Transformadas (función continua no periódica)", "seccion": "m3-s4", "estado": "cubierto"},
    {"tema": "p15 Transformada continua de Fourier", "seccion": "m3-s4", "estado": "cubierto"},
    {"tema": "p16–19 Transformada discreta de Fourier y su costo", "seccion": "m3-s5", "estado": "cubierto"},
    {"tema": "p20–22 Algoritmo: transformada rápida de Fourier (FFT)", "seccion": "m3-s6", "estado": "cubierto"},
    {"tema": "p23 FFT: N no potencia de 2, IFFT", "seccion": "m3-s7", "estado": "cubierto"},
    {"tema": "p24 APLICACIONES: señales y filtros", "seccion": "m3-s8", "estado": "cubierto"},
    {"tema": "p25–26 APLICACIONES: procesamiento digital de imágenes, espectro 2D", "seccion": "m3-s9", "estado": "cubierto"},
    {"tema": "p27 APLICACIONES: filtrado pasa baja / pasa alta", "seccion": "m3-s10", "estado": "cubierto"},
    {"tema": "Nyquist p2 Frecuencia de Muestreo y Frecuencia de Nyquist", "seccion": "m3-s11", "estado": "cubierto"},
    {"tema": "Nyquist p3 Ejemplo", "seccion": "m3-s11", "estado": "cubierto"},
    {"tema": "Nyquist p4–6 Resolución de la Frecuencias", "seccion": "m3-s12", "estado": "cubierto"},
    {"tema": "Nyquist p7 Aliasing", "seccion": "m3-s13", "estado": "cubierto"},
    {"tema": "Nyquist p8 Patrones de Moiré", "seccion": "m3-s13", "estado": "cubierto"},
    {"tema": "Nyquist p9 Preguntas para hacerse (fft2)", "seccion": "m3-s9", "estado": "cubierto"},
    {"tema": "Nyquist p10 Extra (graficar sólo f positivas)", "seccion": "m3-s12", "estado": "cubierto"},
    {"tema": "TP3 ej. 1–6", "seccion": "m3-s15", "estado": "cubierto"},
]

glosario = [
    {"termino": "Aproximación polinómica trigonométrica", "def": "F = Σαᵢfᵢ(x) con senos y cosenos como base; para datos que se repiten periódicamente.", "seccion": "m3-s1"},
    {"termino": "Serie de Fourier", "def": "f(x) = a₀ + Σaₙcos(nω₀x) + Σbₙsin(nω₀x) para una función continua y periódica; coeficientes por integrales.", "seccion": "m3-s2"},
    {"termino": "Frecuencia fundamental f₀", "def": "f₀ = 1/T; ω₀ = 2πf₀ es la frecuencia angular. Menor T, mayor f.", "seccion": "m3-s2"},
    {"termino": "Espectro de frecuencia", "def": "Composición en frecuencia de una señal: qué frecuencias tiene, su distribución y la energía en cada una (amplitud vs frecuencia).", "seccion": "m3-s3"},
    {"termino": "Transformada continua (TCF)", "def": "Integral de −∞ a ∞ para señales continuas no periódicas (T = ∞); espectro continuo y simétrico respecto de 0.", "seccion": "m3-s4"},
    {"termino": "Transformada discreta (TDF)", "def": "X[k] = Σx[n]e^{−j2πkn/N} sobre N muestras; no pide periodicidad; espectro discreto; costo O(N²).", "seccion": "m3-s5"},
    {"termino": "FFT", "def": "Algoritmo de Cooley y Tukey (1965) que calcula exactamente la TDF en O(N log₂N), por divide y vencerás (pares/impares).", "seccion": "m3-s6"},
    {"termino": "ω_N (factor de giro)", "def": "ω_N = e^{−j2π/N}; su periodicidad y simetría (conjugados) son lo que aprovecha la FFT.", "seccion": "m3-s6"},
    {"termino": "Zero padding (relleno con ceros)", "def": "Agregar ceros hasta una potencia de 2 (p23). Cambia N y la grilla Δf; no cambia f_s ni f_N, no crea ni corrige aliasing.", "seccion": "m3-s7"},
    {"termino": "IFFT", "def": "Transformada inversa: vuelve del dominio de la frecuencia al del tiempo (o del espacio).", "seccion": "m3-s7"},
    {"termino": "Frecuencia de muestreo f_s", "def": "Muestras por segundo; T_s = 1/f_s es el período de muestreo.", "seccion": "m3-s11"},
    {"termino": "Frecuencia de Nyquist f_N", "def": "f_N = f_s/2: máxima frecuencia de la que se obtiene información. Depende sólo de f_s.", "seccion": "m3-s11"},
    {"termino": "Resolución Δf", "def": "Δf = f_s/n = 1/duración: separación entre frecuencias del espectro. Aumentar f_s sin alargar el registro no la mejora.", "seccion": "m3-s12"},
    {"termino": "Aliasing", "def": "Una frecuencia mayor que f_N aparece reflejada como otra menor en [0, f_s/2].", "seccion": "m3-s13"},
    {"termino": "Patrones de Moiré", "def": "Aliasing en imágenes: al bajar la resolución aparecen franjas que no existen en el original.", "seccion": "m3-s13"},
    {"termino": "Filtro antialiasing", "def": "Pasa bajo aplicado antes de muestrear para quitar lo que supera f_s/2; después de muestrear ya no sirve.", "seccion": "m3-s8"},
    {"termino": "Filtro pasa bajo", "def": "Deja las bajas frecuencias y atenúa las medias-altas: quita ruido, suaviza (imagen borrosa, señal suave).", "seccion": "m3-s10"},
    {"termino": "Filtro pasa alto", "def": "Atenúa las bajas y deja las medias-altas: quedan los bordes (imagen oscura con contornos) o la oscilación rápida.", "seccion": "m3-s10"},
    {"termino": "fft2 / fftshift / ifftshift / ifft2", "def": "Espectro 2D / centrar la frecuencia cero / deshacer el centrado / volver a la imagen.", "seccion": "m3-s9"},
    {"termino": "Frecuencia de corte D0", "def": "Radio de la máscara circular del notebook: define qué frecuencias pasan.", "seccion": "m3-s10"},
]


# ================================================================== QUIZ
Q = []
def q(sec_, tema, pide, sobre, ops, expl, fig=None):
    it = {"id": f"m3-q{len(Q)+1:02d}", "seccion": f"m3-s{sec_}", "tema": tema, "pide": pide, "sobre": sobre,
          "opciones": [{"t": t, "v": v} for t, v in ops], "explicacion": expl}
    if fig: it["fig"] = fig
    Q.append(it)

# --- con figura
q(13, "Aliasing entre dos muestreos (ET25-T1)", "una",
  "La figura muestra los espectros de muestra1.csv (dt = 0,5 s) y muestra2.csv (dt = 0,01 s), dos muestreos de la misma señal en [0, 6π]. ¿A qué se debe que la muestra 1 tenga un pico cerca de 0,8 Hz que no está en la muestra 2?",
  [("Es aliasing: una frecuencia mayor que la de Nyquist de la muestra 1 aparece reflejada como una frecuencia menor.", True),
   ("En la muestra 1 se rellenó con ceros hasta potencia de 2 y eso generó aliasing al modificar su frecuencia de Nyquist.", False),
   ("Por la baja frecuencia de Nyquist de la muestra 1, el pico de 1 Hz se dividió en dos picos; eso se llama aliasing.", False),
   ("Por aliasing, la muestra 2 no detecta un pico que sí aparece en la muestra 1, debido a su frecuencia de Nyquist.", False)],
  f"La muestra 1 tiene f_s = {f(M1['fs'])} Hz y f_N = {f(M1['fN'])} Hz; la componente real de 1,2 Hz supera f_N y se refleja en |1,2 − 2| = 0,8 Hz. La muestra 2 (f_N = {f(M2['fN'])} Hz) ve 1,2 Hz en su lugar. Lo del padding es falso: rellenar con ceros cambia N, no f_s ni f_N. Un pico no se divide en dos por aliasing. Y decir que es la muestra 2 la que falla invierte el caso: la de mayor f_s es la confiable.",
  FIG["muestras"])
q(10, "Filtro sobre el tablero: franjas", "una",
  "A la izquierda, la imagen original; a la derecha, la misma imagen procesada. ¿Cómo se obtuvo la imagen procesada?",
  [("Se usó fft, se aplicó un filtro pasa bajo y luego ifft para obtener la imagen de la derecha.", True),
   ("Se usó fft, se aplicó un filtro pasa alto y luego ifft para obtener la imagen de la derecha.", False),
   ("La imagen de la derecha es el espacio de frecuencias obtenido al aplicar fft a la original.", False),
   ("Se usó fft haciendo un padding de ceros, y luego ifft para obtener la imagen procesada de la derecha.", False)],
  f"Las franjas diagonales suaves son la componente de frecuencia más baja del tablero (las diagonales oscuras, picos a radio {f(R['tablero_picos'][0][2], 2)} del centro); un pasa bajo de radio {R['D0'][0]} deja sólo esa y borra los bordes de las casillas. El pasa alto haría lo contrario: imagen oscura con bordes. Un espectro centrado tendría un máximo brillante en el centro y puntos aislados, no franjas. El padding seguido de ifft devuelve la misma imagen con ceros alrededor.",
  FIG["tab_pb"])
q(10, "Filtro sobre el tablero: bordes", "una",
  "A la izquierda, la imagen original; a la derecha, la misma imagen procesada. ¿Cómo se obtuvo la imagen procesada?",
  [("Se usó fft, se aplicó un filtro pasa alto y luego ifft para obtener la imagen de la derecha.", True),
   ("Se usó fft, se aplicó un filtro pasa bajo y luego ifft para obtener la imagen de la derecha.", False),
   ("La imagen de la derecha es el espacio de frecuencias obtenido al aplicar fft a la original.", False),
   ("Se aplicó un filtro antialiasing sobre la imagen y luego ifft para obtener la imagen de la derecha.", False)],
  f"La procesada es oscura (se quitó la frecuencia cero, el promedio) y sólo brillan los bordes de las casillas: lo que queda son las transiciones bruscas de gris, que la p27 asocia al pasa alta (brillo medio en bordes {f(R['tablero_pa_borde_vs_int'][0], 3)} vs {f(R['tablero_pa_borde_vs_int'][1], 3)} en el interior). Un pasa bajo daría una imagen suave sin bordes. No es un espectro: no tiene el máximo central ni la simetría de un espectro. El antialiasing se aplica antes de muestrear y no produce bordes.",
  FIG["tab_pa"])
q(10, "Qué filtro se utilizó (secuencia completa)", "una",
  "La figura muestra una imagen, su espectro centrado, la máscara usada (centrada y sin centrar; negro = 0) y el resultado. ¿Qué filtro se utilizó?",
  [("Un filtro pasa alto: la máscara anula el círculo central del espectro centrado.", True),
   ("Un filtro pasa bajo: la máscara conserva sólo el círculo central del espectro.", False),
   ("Un filtro antialiasing: la máscara elimina las frecuencias reflejadas del espectro.", False),
   ("Ningún filtro: se rellenó el espectro con ceros en el centro para volver con ifft.", False)],
  f"En la máscara centrada el disco del centro es negro (0) y el resto blanco (1): se eliminan las bajas frecuencias, que tras fftshift están en el centro (sin centrar quedan en las esquinas, por eso los cuartos de disco). El resultado es oscuro (media {f(R['sint_pa_media'], 3)} contra {f(R['sint_media'], 3)} de la original) con contornos: pasa alto. Si el disco fuera blanco y el resto negro sería pasa bajo. Antialiasing no es lo que muestra una máscara aplicada después de muestrear, y llamar “padding” a anular el centro confunde el relleno de N con un filtrado.",
  FIG["sint"])
q(8, "Señal con ruido procesada", "una",
  "A partir de un muestreo se obtuvo la muestra original (arriba) y se procesó para obtener la señal de abajo. ¿Cómo se obtuvo la señal procesada?",
  [("Se usó fft, se aplicó un filtro pasa bajo y después ifft.", True),
   ("Se usó fft, se aplicó un filtro pasa alto y después ifft.", False),
   ("Se eliminaron las señales de aliasing con un filtro antialiasing.", False),
   ("Se usó fft con padding de ceros para modificar la frecuencia de Nyquist.", False)],
  f"El ruido fino es alta frecuencia; la curva procesada conserva la forma lenta: se dejaron las bajas y se anularon las altas (pasa bajo, f_c = {f(E5['fc'])} Hz en la recreación, error RMS {f(E5['err_pb'], 3)} contra la señal lenta). Un pasa alto dejaría justamente el ruido. El antialiasing es previo al muestreo y no se puede aplicar a la muestra ya tomada para limpiar ruido. El padding no cambia la frecuencia de Nyquist ni filtra nada.",
  FIG["ej5"])
q(8, "Señal lenta + rápida procesada (EV-P6)", "una",
  "La muestra original (arriba) se procesó para obtener la señal de abajo. ¿Cómo se obtuvo?",
  [("Se aplicó un filtro pasa alto.", True),
   ("Se aplicó un filtro pasa bajo.", False),
   ("Se aplicó un filtro antialiasing.", False),
   ("Se hizo padding de ceros y se usó fft.", False)],
  f"La procesada es sólo la oscilación rápida, centrada en 0 y con amplitud ≈{f(E6['amp_pa'], 2)}: se quitó la oscilación lenta (baja frecuencia) y se dejó la alta, o sea pasa alto. El pasa bajo daría la onda lenta de amplitud ≈{f(E6['amp_pb'], 2)}. El filtro antialiasing actúa antes del muestreo y quitaría lo rápido, no lo lento. El padding de ceros no elimina ninguna componente de la señal.",
  FIG["ev6"])
q(13, "Tercer pico distinto en dos muestreos (EV-P5)", "una",
  "Una misma señal se muestreó con dos frecuencias de muestreo distintas. Los espectros coinciden en los dos primeros picos, pero el tercero es distinto en cada caso. ¿Por qué?",
  [("El tercer pico de la muestra 2 es real y el de la muestra 1 es resultado del aliasing.", True),
   ("El tercer pico de la muestra 1 es real y el de la muestra 2 es resultado del aliasing.", False),
   ("El tercer pico de la muestra 1 es real y el de la muestra 2 es resultado del padding.", False),
   ("El tercer pico de la muestra 2 es real y el de la muestra 1 es resultado del padding.", False)],
  f"La muestra 1 sólo llega hasta 5 Hz (f_s = 10, f_N = 5): no puede ver 23 Hz, y 23 Hz con f_s = 10 se refleja en |23 − 2·10| = {f(R['EV_P5_alias'][1][1])} Hz. La muestra 2 (f_s = 60, f_N = 30) ve 23 Hz en su lugar. Atribuirle el alias a la muestra 2 es el mecanismo de “muestra equivocada”: con f_N = 30 no hay reflejo de 23 Hz. El padding cambia N, no crea picos en otra frecuencia.",
  FIG["ev5"])
q(15, "muestreo.ipynb con f_s = 2 Hz", "una",
  "En muestreo.ipynb se muestrean f₁ = 1 Hz y f₂ = 2 Hz (desfase π/2) con f_s = 2 Hz. ¿Qué describe lo que muestra la figura?",
  [("f₁ está justo en f_N y se ve como zigzag ±1; f₂ = f_s se ve constante (alias en 0 Hz).", True),
   ("Ambas se reconstruyen bien porque f_s = 2 Hz es el doble de la frecuencia f₁ = 1 Hz.", False),
   ("f₂ se ve como zigzag de 2 Hz y f₁ queda constante porque es menor que la de Nyquist.", False),
   ("Ninguna cambia: con desfase π/2 las muestras caen siempre sobre los máximos de la onda.", False)],
  f"Con f_s = 2, f_N = 1: f₁ = 1 Hz está en el límite y sólo hay {f(MU[2]['f1_por_periodo'])} muestras por ciclo (valores {', '.join(f(v) for v in MU[2]['f1_muestras'])}), un zigzag que recién con fase coseno conserva la amplitud. f₂ = 2 Hz = f_s: todas las muestras dan {', '.join(f(v) for v in MU[2]['f2_muestras'])}, una constante (alias en {f(MU[2]['f2_alias'])} Hz). Decir que ambas salen bien ignora f₂; invertir cuál queda constante es el error de caso; el desfase no salva a f₂.",
  FIG["fs2"])
q(3, "Espectro de s(t) = sin t + 1,2 cos 2t", "una",
  "La figura muestra s(t) = sin(t) + 1,2·cos(2t) y su espectro. ¿Qué afirmación es correcta?",
  [("Hay dos barras: altura 1 en 1/2π ≈ 0,159 Hz y altura 1,2 en 1/π ≈ 0,318 Hz.", True),
   ("Hay dos barras: altura 1 en 1 Hz y altura 1,2 en 2 Hz, porque ω vale 1 y 2.", False),
   ("La barra más alta corresponde a sin(t), que es la componente más rápida de s(t).", False),
   ("El espectro es continuo porque s(t) es una función continua del tiempo t real.", False)],
  f"f₀ = ω/2π: para sin(t), ω = 1 → {f(ES['f1'], 4)} Hz con amplitud 1; para 1,2cos(2t), ω = 2 → {f(ES['f2'], 4)} Hz con amplitud 1,2 (p13). Poner 1 y 2 Hz confunde frecuencia angular con frecuencia en Hz. La más rápida es 1,2cos(2t), no sin(t). La señal es periódica y se descompone con la serie: el espectro es discreto, un valor por frecuencia presente.",
  FIG["ej_s"])
q(9, "Espectro de franjas verticales", "una",
  "Imagen de 128×128 con franjas verticales de período 16 píxeles y su espectro centrado. ¿Qué describe el espectro?",
  [("Dos picos sobre el eje horizontal, a ±8 del centro, además del valor central.", True),
   ("Dos picos sobre el eje vertical, a ±8 del centro, paralelos a las franjas.", False),
   ("Energía repartida en todas las frecuencias porque la imagen tiene muchos bordes.", False),
   ("Un único pico en el centro porque la imagen tiene un solo nivel de gris medio.", False)],
  f"Las franjas varían sólo en horizontal (x), con 128/16 = 8 ciclos: aparecen picos en (0, ±8) del centro, sobre el eje horizontal (verif: picos en {PK(R['rayas_picos'])}), perpendicular a las franjas. Ubicarlos en el eje vertical confunde la dirección. Es una senoidal suave, sin bordes bruscos, así que no hay energía en todas las frecuencias; y el pico central sólo es el promedio, no describe la oscilación (p26: patrón regular → picos separados a distancia regular).",
  FIG["rayas"])

# --- opción única sin figura
q(11, "Hasta qué frecuencia (EV-P2)", "una",
  "Una señal muestreada con frecuencia de muestreo f_s = 10 Hz sólo puede representar correctamente frecuencias hasta:",
  [("f_s/2 = 5 Hz", True), ("f_s = 10 Hz", False), ("2·f_s = 20 Hz", False), ("f_s² = 100 Hz", False)],
  "La frecuencia de Nyquist es f_N = f_s/2 = 5 Hz (nyquist p2: “No podemos obtener información de frecuencias mayores que f_N!”). f_s = 10 Hz es lo que se suele confundir: una senoidal de 10 Hz muestreada a 10 Hz da todas las muestras iguales. 2·f_s confunde la condición f_s ≥ 2f con su inversa, y f_s² es un múltiplo sin sentido físico (mecanismo de múltiplo equivocado).")
q(13, "Cuándo hay aliasing (EV-P3)", "una",
  "¿En cuál de estos casos se observará el fenómeno de aliasing al realizar el muestreo de una señal?",
  [("La frecuencia de Nyquist es menor que la de la señal original", True),
   ("La frecuencia de Nyquist es mayor que la de la señal original", False),
   ("Cuando la FFT se calcula con pocos puntos de la señal", False),
   ("Al aplicar un filtro pasa-bajo sobre la señal muestreada", False)],
  "Hay aliasing cuando la señal tiene frecuencias por encima de f_N = f_s/2: se reflejan en [0, f_s/2] (nyquist p7). Si f_N es mayor que la frecuencia de la señal, el muestreo es correcto. Pocos puntos empeoran la resolución Δf = f_s/n, pero no generan alias. Un pasa bajo (antes de muestrear) es justamente lo que lo evita; después no lo crea.")
q(11, "f₁ = 10, f₂ = 20, f₃ = 100 Hz con f_s = 20 (EV-P8)", "una",
  "Una señal tiene tres senoidales de f₁ = 10 Hz, f₂ = 20 Hz y f₃ = 100 Hz y se muestrea con f_s = 20 Hz. ¿Qué frecuencias se podrán reconstruir correctamente?",
  [("Solo f1 = 10 Hz", True), ("Solo f2 = 20 Hz", False), ("f1 = 10 Hz y f2 = 20 Hz", False), ("Solamente f3 = 100 Hz", False)],
  "f_N = 20/2 = 10 Hz. Con la convención de la cátedra (“f_s de por lo menos 4 Hz” para 2 Hz, es decir f ≤ f_N), sólo f₁ = 10 Hz, que está justo en el límite (y depende de la fase: con fase seno sus muestras dan 0). f₂ = 20 Hz = f_s y f₃ = 100 Hz = 5f_s se reflejan en 0 Hz: se ven como una constante. Elegir f₂ o f₁ y f₂ confunde f_s con f_N; f₃ es la más lejana del rango.")
q(11, "Qué f_s elegir (ET25-T3 ej. 6)", "una",
  "Una señal tiene tres senoidales de f₁ = 10 Hz, f₂ = 15 Hz y f₃ = 100 Hz. ¿Qué frecuencia de muestreo f_s debería usar para reproducirlas correctamente?",
  [("f_s = 200 Hz", True), ("f_s = 100 Hz", False), ("f_s = 20 Hz", False), ("f_s = 10 Hz", False)],
  "Manda la frecuencia más alta: f_s ≥ 2·100 = 200 Hz para que f_N = 100 Hz la alcance. Con 100 Hz, f_N = 50 y f₃ cae en 0 Hz; con 20 Hz, f_N = 10 y se pierden 15 y 100 Hz (15 aparece en 5 Hz); con 10 Hz, f_N = 5 y no se reconstruye ninguna. Las opciones mal son el mecanismo de múltiplo equivocado: f_s igual a la frecuencia, o el doble de la menor.")
q(13, "Calcular un alias", "una",
  "Una senoidal de 23 Hz se muestrea con f_s = 10 Hz. ¿En qué frecuencia aparece en el espectro (eje de 0 a f_s/2)?",
  [("3 Hz", True), ("13 Hz", False), ("5 Hz", False), ("7 Hz", False)],
  f"La frecuencia se refleja en [0, f_s/2] = [0, 5]: f_alias = |23 − 10·round(23/10)| = |23 − 20| = {f(R['EV_P5_alias'][1][1])} Hz. 13 Hz es restar f_s una sola vez y queda fuera del rango visible. 5 Hz es f_N, el borde del eje, no el alias. 7 Hz sería 10 − 3, reflejar una vez de más.")
q(5, "Costo de la TDF directa", "una",
  "Para la TDF directa con N = 10⁹ puntos, suponiendo 1 ns por operación, ¿qué afirmación es correcta?",
  [(f"Cuesta O(N²): unas 10¹⁸ operaciones, del orden de {ANIOS} años; la FFT tarda unos {SEG_FFT} s.", True),
   (f"Cuesta O(N log₂N): unos {SEG_FFT} s, igual que la FFT, porque las dos calculan la misma suma.", False),
   ("Cuesta O(N): cada X[k] es una sola suma de N términos que se calcula una única vez.", False),
   ("Cuesta O(N²) sólo si N no es potencia de 2; con 2³⁰ puntos ya se reduce a O(N log₂N).", False)],
  f"Por cada uno de los N coeficientes X[k] hay N multiplicaciones complejas y N − 1 sumas: N² = 10¹⁸ ns ≈ {ANIOS} años (la p19 dice 31,2). La FFT (N log₂N ≈ {f(C[2]['NlogN']/1e9, 1)}×10⁹) tarda ≈{SEG_FFT} s. Las dos calculan lo mismo, pero el costo es del algoritmo, no del resultado. O(N) olvida que hay N coeficientes. Y la potencia de 2 importa para la FFT radix-2, no para la TDF directa.")
q(7, "fft de SciPy sin potencia de 2 (IT25-T1)", "una",
  "Con la fft de SciPy/NumPy y un número de puntos que NO es potencia de 2, ¿qué ocurre?",
  [("La calcula igual sobre los N puntos; según los factores de N tarda más, y un N primo es el peor caso.", True),
   ("Da error: el algoritmo de Cooley–Tukey sólo está definido para N igual a una potencia de 2.", False),
   ("Rellena automáticamente con ceros hasta la siguiente potencia de 2 y devuelve ese espectro.", False),
   ("Descarta las últimas muestras hasta la potencia de 2 anterior y calcula la FFT con esas.", False)],
  f"Las librerías aceptan cualquier N (descomponen N en factores; para primos grandes usan otro algoritmo O(N log N)). Medido acá: 65536 = 2¹⁶ tarda {f(TI[65536], 2)} ms y 65537 (primo) {f(TI[65537], 2)} ms, la salida tiene 65537 valores. No da error. No rellena sola: sólo si pasás n con np.fft.fft(x, n). Tampoco descarta datos. El relleno con ceros es lo que muestra la p23 para el radix-2 de la clase.")
q(6, "En qué se basa la FFT (IT25-T4 c)", "una",
  "¿En qué se basa la FFT para bajar el costo computacional respecto de la TDF?",
  [("Divide la suma en índices pares e impares (dos TDF de N/2) y aprovecha periodicidad y simetría de ω_N.", True),
   ("Aproxima la integral de la transformada continua con trapecios, lo que reduce el número de términos.", False),
   ("Submuestrea la señal a la mitad en cada paso, descartando muestras hasta quedarse con potencias de 2.", False),
   ("Calcula sólo las frecuencias positivas, ya que las negativas son el espejo de las positivas.", False)],
  "La p20–22: “Técnica D&C”, X[k] = Σpares + ω_N^k·Σimpares, cada Σ es una TDF de N/2 y se repite log₂N veces; aprovecha la periodicidad y la simetría (conjugados) de ω_N^{nk}. Es exacta: no aproxima ninguna integral. No descarta muestras: usa todas, separadas en pares e impares. Graficar sólo f > 0 es una cuestión de visualización (nyquist p10), no lo que baja el costo.")
q(4, "TCF vs TDF (IT25-T1 a)", "una",
  "¿Cuál es la principal diferencia entre la transformada continua (TCF) y la discreta (TDF) de Fourier?",
  [("La TCF integra una función continua y da un espectro continuo; la TDF suma N muestras y da N valores.", True),
   ("La TCF se aplica a señales periódicas y la TDF a no periódicas; ambas devuelven un espectro discreto.", False),
   ("La TCF es exacta y la TDF es una aproximación que se calcula con un algoritmo de costo N log N.", False),
   ("La TCF trabaja en el dominio de la frecuencia y la TDF en el dominio del tiempo de la señal.", False)],
  "p15: TCF = integral, señal continua no periódica, espectro continuo. p16: TDF = suma sobre N puntos de una señal discreta, espectro discreto, no necesita periodicidad. Decir que la TCF es para periódicas confunde con la serie de Fourier, y el espectro de la TCF no es discreto. N log N es el costo de la FFT, no de la TDF (O(N²)). Ambas pasan del tiempo a la frecuencia: confundir dominio es el mecanismo 4.")
q(12, "Mejorar la resolución en frecuencia", "una",
  f"Con [0, 2π) y f_s = 10 Hz se obtiene n = {RES[0]['n']} y Δf = {f(RES[0]['df'])} Hz. ¿Qué cambio mejora la resolución en frecuencia?",
  [(f"Registrar más tiempo: con [0, 6π) y f_s = 6 Hz, Δf baja a {f(RES[2]['df'])} Hz.", True),
   (f"Subir f_s a 100 Hz en el mismo intervalo: n sube a {RES[1]['n']} y Δf baja.", False),
   ("Bajar f_s a 2 Hz en el mismo intervalo, para que f_N quede más cerca de los picos.", False),
   ("Graficar sólo las frecuencias positivas, como recomienda la diapositiva extra.", False)],
  f"Δf = f_s/n = 1/duración. La cátedra lo muestra: con f_s = 100, n = {RES[1]['n']} pero Δf sigue en {f(RES[1]['df'])} (“aumentar f_s no mejora la resolución”); con [0, 6π) Δf = {f(RES[2]['df'])}. Bajar f_s tampoco cambia Δf (baja n en la misma proporción) y además acerca el aliasing. Graficar f > 0 sólo oculta la mitad espejo del espectro.")
q(7, "Efecto del padding de ceros", "una",
  "Una señal de 1000 muestras (f_s = 100 Hz) se rellena con ceros hasta 1024 antes de la FFT. ¿Qué efecto tiene?",
  [("Cambia N y la grilla de frecuencias (f_s/1024); no cambia f_N ni agrega información.", True),
   ("Genera aliasing porque al cambiar N se modifica la frecuencia de Nyquist de la señal.", False),
   ("Duplica la frecuencia de Nyquist porque hay más puntos por segundo en la señal.", False),
   ("Elimina el ruido de alta frecuencia porque los ceros no tienen ninguna oscilación.", False)],
  "El relleno sólo agrega ceros al final: N pasa a 1024 (potencia de 2, p23), el paso de la grilla a f_s/1024, y el espectro queda interpolado. f_s sigue en 100 Hz, así que f_N = 50 Hz no cambia: no genera aliasing (es el distractor de ET25-T1) ni lo duplica (las muestras siguen separadas 0,01 s). Tampoco filtra: las componentes de la señal siguen todas.")

# --- verdaderas
q(14, "Serie, TCF, TDF y FFT", "verdaderas", "los espectros de Fourier",
  [("El espectro de la serie de Fourier es discreto: un valor de amplitud por cada frecuencia presente.", True),
   ("La FFT da el mismo resultado que la TDF; lo que cambia es el costo del cálculo.", True),
   ("La TCF da un espectro discreto porque se calcula sobre una función continua.", False),
   ("La TDF sólo puede aplicarse a señales periódicas, igual que la serie.", False)],
  "La p10 dice que la serie da un espectro discreto (un valor por frecuencia presente). La FFT es un algoritmo para implementar la TDF (p20): mismo resultado, O(N log₂N) en vez de O(N²). La TCF da un espectro continuo (p15: “El resultado es un espectro contínuo”). La TDF no requiere periodicidad (p16: “No es necesaria la periodicidad”), a diferencia de la serie.")
q(6, "Costo y potencia de 2 (WA25)", "verdaderas", "la FFT",
  [("Divide la TDF en las sumas de índices pares e impares y repite la división log₂N veces.", True),
   ("Con N = 10⁹ y 1 ns por operación tarda unos 30 s, contra décadas de la TDF directa.", True),
   ("Tiene costo N log N sólo si N es potencia de 2; para cualquier otro N vuelve a costar N².", False),
   ("Aproxima la TDF: ahorra cálculos a cambio de un error de truncamiento en X[k].", False)],
  f"Pares/impares y log₂N niveles es la técnica D&C de la p22. La p20 da 30×10⁹ ns ≈ 30 s contra ≈31 años (verif: {SEG_FFT} s y {ANIOS} años). “Sólo si N es potencia de 2” es falso por el absoluto: con relleno de ceros (p23) N′ es menor que 2N y sigue siendo O(N log N), y las librerías usan factores mixtos. La FFT es exacta: da la misma TDF, no una aproximación.")
q(11, "Nyquist y aliasing", "verdaderas", "el muestreo de una señal",
  [("La frecuencia de Nyquist es f_s/2 y depende sólo de la frecuencia de muestreo.", True),
   ("Una frecuencia mayor que f_N se refleja en el rango [0, f_s/2].", True),
   ("Aumentar f_s manteniendo el intervalo [a, b] mejora la resolución Δf.", False),
   ("El aliasing se corrige después de muestrear calculando la FFT con muchos más puntos.", False)],
  "f_N = f_s/2 (nyquist p2) y no depende de la señal ni de N. La p7 dice que la frecuencia se “refleja” en [0, f_s/2]. Δf = f_s/n con n = f_s·(b−a): f_s se cancela; la p5 lo dice textual. Una vez muestreado, el alias está mezclado con frecuencias reales: más puntos (padding) no lo quitan; se evita antes, con f_s mayor o un filtro antialiasing.")
q(10, "Filtros", "verdaderas", "el filtrado en el dominio de la frecuencia",
  [("Un pasa bajo deja las bajas frecuencias y atenúa las medias-altas; sirve para quitar ruido de alta frecuencia.", True),
   ("En el notebook, el pasa alto se obtiene como 1 − pasa bajo, o invirtiendo la condición de distancia D0.", True),
   ("Un pasa alto aplicado a una imagen conserva las zonas homogéneas y borra los bordes de los objetos.", False),
   ("Para filtrar se multiplica la imagen original por la máscara circular y después se aplica fft2.", False)],
  "Las dos primeras salen de la p27 y del notebook (filtro_pasa_bajo pone 1 si la distancia al centro es menor que D0; el pasa alto es el complemento). El pasa alto hace lo opuesto a lo dicho: se queda con los bordes, donde el gris varía bruscamente, y borra lo homogéneo. La máscara se multiplica por el espectro (fft2 + fftshift), no por la imagen: el orden es fft2, máscara, ifftshift, ifft2.")
q(9, "Espectro de una imagen (p26)", "verdaderas", "el espectro de una imagen",
  [("Las zonas homogéneas concentran la energía en las bajas frecuencias.", True),
   ("Los bordes y transiciones frecuentes de gris dan componentes de alta frecuencia.", True),
   ("Un patrón que se repite da picos de igual intensidad a distancia regular.", True),
   ("Sin fftshift, la frecuencia cero de fft2 queda en el centro de la matriz.", False)],
  "Las tres primeras son textuales de la p26. La cuarta está al revés: fft2 deja la frecuencia (0,0) en la esquina superior izquierda (las bajas en las cuatro esquinas, como el esquema de la p26); fftshift es el que la lleva al centro, e ifftshift deshace ese movimiento antes de ifft2.")
q(2, "Coeficientes de la serie", "verdaderas", "las incógnitas de la serie de Fourier",
  [("Son a₀, aₙ y bₙ; para una representación exacta hacen falta infinitos coeficientes.", True),
   ("Para calcularlos con las integrales hace falta conocer la función f.", True),
   ("Si el período T no es 2π, la serie de Fourier ya no se puede aplicar a la función f.", False),
   ("Se obtienen resolviendo las ecuaciones normales, como en mínimos cuadrados.", False)],
  "p6: “el objetivo es encontrar los coeficiente a₀, aₙ y bₙ” con n = 1…∞, y “Obs: Para calcular los coef necesito la función”. Para T cualquiera la p7 generaliza con ω₀ = 2π/T y límites 0…T, así que la serie se aplica igual. Los coeficientes salen de integrales (gracias a la ortogonalidad de senos y cosenos), no de un sistema de ecuaciones normales.")
q(1, "Polinómica vs trigonométrica (WA25)", "verdaderas", "la aproximación polinómica y la trigonométrica",
  [("Ambas tienen la forma F = α₁f₁(x) + … + αₙfₙ(x); cambia la familia de funciones base.", True),
   ("La trigonométrica es la indicada cuando los datos se repiten periódicamente.", True),
   ("Una señal de 1 Hz es más rápida que una de 2 Hz porque completa antes el ciclo.", False),
   ("La aproximación lineal ajusta bien datos periódicos si se suman suficientes términos de grado 1.", False)],
  "p2–4: la misma forma F = Σαᵢfᵢ(x), con base polinómica o con senos y cosenos; la p3 elige la trigonométrica porque los valores “parecen ser periódicos”. La p8 dice lo contrario de la tercera: 1 Hz es más lenta que 2 Hz (tarda más en completar el ciclo). Sumar términos de grado 1 sigue dando una recta: nunca oscila.")
q(9, "Pasos del notebook de imágenes", "verdaderas", "fft2 y las funciones del notebook imagen_con_filtro",
  [("fft2 puede replicarse aplicando fft a cada fila y luego a cada columna.", True),
   ("ifftshift deshace el centrado del espectro antes de aplicar ifft2.", True),
   ("ifft2 devuelve una matriz compleja; por eso se grafica np.abs del resultado.", True),
   ("fftshift modifica los valores del espectro, suavizándolo para graficarlo.", False)],
  f"La TDF 2D es separable (verif: fft por filas y columnas = fft2 → {SI(R['fft2_filas_cols'])}), que responde la pregunta de nyquist p9. En el notebook img_c4 = ifftshift (deshace el centrado) e img_c5 = ifft2, que devuelve complejos y se grafica np.abs. fftshift sólo reordena los cuadrantes (mueve la frecuencia cero al centro); no cambia ningún valor.")
q(7, "Padding e IFFT", "verdaderas", "el relleno con ceros y la IFFT",
  [("Con la IFFT se vuelve del dominio de la frecuencia al dominio del tiempo.", True),
   ("Relleno con fft2(s mayor) y luego ifft2 devuelve la imagen original con ceros alrededor.", True),
   ("El relleno con ceros modifica la frecuencia de Nyquist de la señal muestreada.", False),
   ("El relleno con ceros agrega a la señal componentes de frecuencia que antes no estaban en ella.", False)],
  f"p23: “Puedo volver desde el dominio de la frecuencia al dominio del tiempo?” → IFFT. Verificado: ifft2(fft2(tablero, s = 256×256)) da el tablero en la esquina y ceros en el resto ({SI(R['pad_img_ok'])}): el padding no suaviza ni genera franjas. f_N = f_s/2 no depende de N. Los ceros no traen frecuencias nuevas: el espectro sólo se evalúa en una grilla más fina.")


# ================================================================== DESARROLLO
DEV = []
def d(sec_, real, pregunta, modelo, fig=None, practica=False):
    it = {"id": f"m3-d{len(DEV)+1:02d}", "seccion": f"m3-s{sec_}", "real": real, "pregunta": pregunta, "modelo": modelo}
    if fig: it["fig"] = fig
    if practica: it["practica"] = True
    DEV.append(it)

# ---------------- teoría real
d(14, "IT25-T1 P2",
  "a. ¿Cuál es la principal diferencia entre la transformada discreta de Fourier y la transformada continua de Fourier? b. Dada la funcionalidad de fft provista por Scipy, ¿es necesario usar un número de puntos en potencia de dos para poder usarla? ¿Qué ocurre si no son potencia de dos? c. ¿Cuál es el tiempo de ejecución de la FFT? ¿y de la TDF tradicional?",
  P("<b>a)</b> La <b>TCF</b> se aplica a señales continuas y no periódicas (T → ∞): F(ω) = ∫₋∞^∞ v(t)e^{jωt}dt es una <b>integral</b> y el resultado es un <b>espectro continuo</b> (hay amplitud para todo ω). La <b>TDF</b> se aplica a una señal discreta de N muestras x[n]: X[k] = Σₙ₌₀^{N−1} x[n]e^{−j2πkn/N}, k = 0…N−1, es una <b>suma finita</b> y el resultado es un <b>espectro discreto</b> de N valores, separados Δf = f_s/N. Además la TDF no necesita periodicidad. La diferencia de fondo es qué entra (función continua vs muestras) y qué sale (espectro continuo vs discreto); la TDF es la que se puede calcular en una computadora.")
  + P(f"<b>b)</b> No es necesario. La fft de SciPy (y la de NumPy) acepta cualquier N: factoriza N y aplica Cooley–Tukey con esos factores (y otro algoritmo si queda un primo grande). Si N no es potencia de 2 el resultado es el mismo (la TDF de los N puntos), pero <b>tarda más</b> cuanto peores son los factores de N: medido acá, 65536 = 2¹⁶ tarda {f(TI[65536], 2)} ms y 65537 (primo) {f(TI[65537], 2)} ms. El Cooley–Tukey radix-2 de la clase (p22) sí necesita N potencia de 2; la diapositiva p23 muestra la solución: <b>rellenar con ceros</b> hasta la siguiente potencia de 2 (en SciPy, pasando n). Eso cambia N y la grilla de frecuencias, pero no f_s ni la frecuencia de Nyquist.")
  + P(f"<b>c)</b> FFT: <b>O(N log₂N)</b>. TDF tradicional: <b>O(N²)</b>, porque cada uno de los N coeficientes X[k] requiere N multiplicaciones complejas y N−1 sumas (p18). Con N = 10⁹ y 1 ns por operación: N² = 10¹⁸ ns ≈ {ANIOS} años (la p19 dice 31,2) contra N log₂N ≈ 3×10¹⁰ ns ≈ {SEG_FFT} s."))
d(14, "IT25-T2 P2",
  "a. ¿Cualquier función puede ser aproximada con la serie de Fourier y la transformada discreta de Fourier? Explicar. b. ¿Cuáles y cuántas son las incógnitas que se deben encontrar para aplicar la serie de Fourier? Explicar. c. ¿A cualquier función se le puede aplicar la FFT? ¿por qué? ¿Qué ventajas tiene?",
  P("<b>a)</b> No cualquiera. La <b>serie de Fourier</b> se aplica a funciones <b>continuas y periódicas</b> (teorema de la p4: “Cualquier función continua y periódica puede ser aproximada mediante la suma de senos y cosenos”), y además hay que conocer la función para calcular los coeficientes con integrales. Si la función no es periódica corresponde la transformada continua (T = ∞). Con discontinuidades (escalón) la serie converge, pero se necesitan infinitos términos y cerca del salto siempre queda un sobrepico. La <b>TDF</b> se aplica a cualquier señal <b>muestreada</b> (N valores; no necesita periodicidad), pero sólo representa bien frecuencias menores que la de Nyquist f_s/2: si la señal tiene componentes mayores, aparecen reflejadas (aliasing).")
  + P("<b>b)</b> Las incógnitas son los coeficientes <b>a₀, aₙ y bₙ</b> de f(x) = a₀ + Σaₙcos(nω₀x) + Σbₙsin(nω₀x): a₀ (valor medio), y un aₙ y un bₙ por cada armónico nω₀. Para una representación exacta son <b>infinitas</b> (n = 1…∞) porque cada frecuencia múltiplo de la fundamental aporta un coseno y un seno; si se trunca en N armónicos son 2N + 1. Se calculan con a₀ = (1/T)∫₀ᵀf, aₙ = (2/T)∫₀ᵀf·cos(nω₀x), bₙ = (2/T)∫₀ᵀf·sin(nω₀x), por eso hace falta conocer f.")
  + P("<b>c)</b> La FFT se aplica a cualquier señal <b>discreta</b> (N muestras), porque es sólo un algoritmo que calcula la TDF; no se aplica a una función continua sin muestrear. La versión de Cooley–Tukey de la clase pide N potencia de 2 (si no, se rellena con ceros); las librerías aceptan cualquier N. Ventajas: costo <b>O(N log₂N)</b> en lugar de O(N²) (con N = 10⁹: ~30 s contra décadas), resultado <b>exacto</b> (no es una aproximación de la TDF), y con la IFFT se vuelve al tiempo con el mismo costo, lo que hace práctico filtrar."))
d(14, "IT25-T3 P2",
  "a. Dada una señal de audio continua y periódica, ¿recomendaría aplicar la transformada rápida de Fourier o la serie de Fourier? Explicar por qué. b. ¿Cuál es la interpretación de la Serie de Fourier? ¿y de la transformada rápida de Fourier? c. ¿Cuáles son las ventajas de pasar del dominio del tiempo al dominio de la frecuencia?",
  P("<b>a)</b> Hay dos lecturas y conviene escribir las dos, empezando por la primera. <b>(1) Si la señal se conoce como función</b> (expresión analítica s(t) continua y periódica), la herramienta teórica es la <b>serie de Fourier</b>: es exactamente para señales continuas y periódicas (p7) y da los coeficientes por integrales. <b>(2) En la práctica un audio se graba</b>: es una señal muestreada, no una fórmula, y para calcular los coeficientes de la serie “necesito la función” (p6). Entonces se aplica la <b>FFT</b> sobre las muestras, que calcula la TDF en O(N log₂N), eligiendo f_s mayor que el doble de la máxima frecuencia del audio para no tener aliasing. Recomendación: FFT sobre la señal muestreada, porque es lo único computable; la serie es el modelo que justifica el resultado.")
  + P("<b>b)</b> <b>Serie:</b> la señal se descompone en una suma de senos y cosenos de frecuencias f₀, 2f₀, 3f₀…; cada coeficiente es la amplitud de esa frecuencia. Su espectro es <b>discreto</b>: un valor de amplitud por cada frecuencia presente (p10), lo que permite “discriminar las frecuencias presentes, su distribución y la energía en cada frecuencia”. <b>FFT:</b> da los N coeficientes X[k] de la TDF; |X[k]| es la amplitud de la frecuencia k·f_s/N. Es un espectro discreto con resolución Δf = f_s/N y sólo tiene sentido hasta f_s/2 (la otra mitad es el espejo); un pico indica una frecuencia presente en la señal (o un alias, si la señal tenía componentes por encima de f_s/2).")
  + P("<b>c)</b> Ventajas: (1) se ven directamente qué frecuencias tiene la señal y con cuánta energía, cosa que en el tiempo está mezclada; (2) se puede <b>filtrar</b> anulando frecuencias (ruido de alta frecuencia con pasa bajo, tendencia lenta con pasa alto, una banda) y volver al tiempo con la IFFT; (3) la señal se resume en pocos coeficientes (compresión, reconocimiento); (4) con la FFT todo eso cuesta O(N log N)."))
d(14, "IT25-T4 P2",
  "a. ¿Cuál es la diferencia entre la transformada continua de Fourier, la transformada discreta de Fourier y la Transformada Rápida de Fourier? b. ¿Cómo son y cuál es la interpretación de los espectros obtenidos en cada caso? c. ¿En qué se basa la FFT para bajar el costo computacional?",
  P("<b>a)</b> <b>TCF</b>: para señales continuas no periódicas; es una integral F(ω) = ∫₋∞^∞ v(t)e^{jωt}dt. <b>TDF</b>: para N muestras de una señal discreta; es una suma X[k] = Σx[n]e^{−j2πkn/N}, no requiere periodicidad y cuesta O(N²). <b>FFT</b>: no es otra transformada sino un <b>algoritmo</b> (Cooley y Tukey, 1965) que calcula exactamente la TDF con costo O(N log₂N).")
  + P("<b>b)</b> <b>TCF</b>: espectro <b>continuo</b> y simétrico respecto de 0; la altura en cada ω indica cuánto aporta esa frecuencia (como T = ∞ las frecuencias quedan infinitamente juntas). <b>TDF</b>: espectro <b>discreto</b> de N valores en las frecuencias k·f_s/N; se interpreta hasta f_N = f_s/2 y su resolución es Δf = f_s/N (más N con más duración = más resolución). <b>FFT</b>: el mismo espectro que la TDF, porque calcula lo mismo; un pico es una frecuencia presente, uno inesperado puede ser un alias.")
  + P("<b>c)</b> En la técnica <b>divide y vencerás</b>: separa la suma en índices pares e impares, X[k] = Σx[2n]ω_{N/2}^{kn} + ω_N^k·Σx[2n+1]ω_{N/2}^{kn}, es decir dos TDF de N/2 puntos, y repite hasta llegar a TDF de 1 punto (log₂N niveles con N operaciones cada uno). Funciona porque ω_N^{nk} es <b>periódico</b> (período N en n y en k), <b>simétrico</b> (ω_N^{k(N−n)} = (ω_N^{kn})*) y ω_N^{kN} = 1: con esto X[k + N/2] = E[k] − ω_N^k·O[k] reutiliza las mismas cuentas. Resultado O(N log₂N), con N potencia de 2."))
d(14, "F20-P4",
  "Explique en qué casos aplicaría la Serie de Fourier, la Transformada Continua de Fourier y la Transformada Discreta de Fourier para señales determinadas.",
  TBL(["Herramienta", "Cuándo", "Por qué"], [
      ["Serie de Fourier", "Señal continua y periódica conocida como función (s(t) = s(t+T))", "Se descompone en armónicos de f₀ = 1/T; los coeficientes salen de integrales sobre un período; espectro discreto"],
      ["TCF", "Señal continua <b>no</b> periódica conocida como función", "Se toma T → ∞; la integral de −∞ a ∞ da un espectro continuo"],
      ["TDF", "Señal <b>discreta</b> (muestreada, datos) de N puntos, periódica o no", "Es una suma finita: se calcula en la computadora (con la FFT); espectro discreto de N valores, válido hasta f_s/2"]])
  + P("Ejemplos: s(t) = sin(t) + 1,2cos(2t) → serie (dos coeficientes no nulos, en 1/2π y 1/π Hz). Un pulso aislado → TCF (espectro continuo). Una grabación de audio o los CSV del TP3 → TDF/FFT, con f_s suficiente para no tener aliasing. La justificación en cada caso es el tipo de señal: continua/discreta y periódica/no periódica."))
d(14, "F18-P3",
  "¿Cuál es el objetivo de la transformada de Fourier de una función? Diferencia entre la Serie y la transformada de Fourier. Diferencia entre la TDF y la FFT.",
  P("<b>Objetivo:</b> pasar la función del dominio del tiempo al de la frecuencia, para conocer qué frecuencias contiene, cómo se distribuyen y cuánta energía tiene cada una (su espectro), y poder operar ahí (filtrar, comprimir, analizar) y volver con la inversa.")
  + P("<b>Serie vs transformada:</b> la serie se aplica a funciones <b>periódicas</b>: se expresan como suma de senos y cosenos en las frecuencias nf₀ y el espectro es discreto (un valor por armónico). La transformada se aplica a funciones <b>no periódicas</b>: se considera T = ∞, la suma se vuelve una integral y el espectro es continuo.")
  + P("<b>TDF vs FFT:</b> la TDF es la definición X[k] = Σₙ x[n]e^{−j2πkn/N} para una señal discreta; calcularla directamente cuesta O(N²). La FFT es un algoritmo (Cooley–Tukey) que obtiene el mismo resultado en O(N log₂N), dividiendo la suma en pares e impares y aprovechando periodicidad y simetría de ω_N; en su versión básica requiere N potencia de 2."))
d(6, "F18-P5", "¿Cuáles son las ventajas del algoritmo de Cooley–Tukey?",
  UL(f"<b>Costo O(N log₂N)</b> en lugar de O(N²): con N = 1024, {MULT[1024][2]} multiplicaciones complejas contra {MULT[1024][1]}; con N = 10⁹ y 1 ns/operación, ≈{SEG_FFT} s contra ≈{ANIOS} años. Porque cada nivel de división reutiliza las cuentas: X[k] y X[k+N/2] salen de los mismos E[k] y ω_N^k·O[k].",
     "<b>Exacto:</b> da la misma TDF; no hay error de aproximación (sólo redondeo).",
     "<b>Divide y vencerás:</b> el problema de N se parte en dos idénticos de N/2 (pares e impares), recursivo y fácil de programar; los subproblemas son independientes (paralelizables).",
     "<b>Inversa con el mismo costo</b> (IFFT): hace práctico ir y volver del dominio de la frecuencia para filtrar.",
     "Hizo posible el procesamiento digital de señales e imágenes en tiempo razonable.")
  + P("Limitación a mencionar: la versión radix-2 pide N potencia de 2; si no, se rellena con ceros (p23), o se usa una variante con otros factores."))
d(1, "WA25", "¿Cuál es la diferencia entre la aproximación polinómica lineal y la aproximación polinómica trigonométrica?",
  P("Las dos son ajustes de la forma <b>F = α₁f₁(x) + α₂f₂(x) + … + αₙfₙ(x)</b>: se busca una combinación lineal de funciones base que represente los datos. La diferencia es la <b>base</b>: en la polinómica lineal las fᵢ son 1 y x (o potencias de x en un polinomio), y en la trigonométrica son senos y cosenos de frecuencias nω₀ (la serie de Fourier).")
  + P("<b>Cuándo cada una:</b> la lineal/polinómica describe tendencias (crecer, decrecer); la trigonométrica se usa cuando los datos <b>se repiten periódicamente</b> (p3), porque el teorema garantiza que una función continua y periódica se aproxima con senos y cosenos. Una recta nunca puede oscilar. <b>Incógnitas:</b> en la lineal, los α por ecuaciones normales de mínimos cuadrados; en la trigonométrica, a₀, aₙ y bₙ, que salen de integrales gracias a la ortogonalidad de senos y cosenos. <b>Interpretación:</b> los coeficientes trigonométricos son amplitudes de cada frecuencia, o sea el espectro de la señal."))
d(7, "WA25 (V/F)", "Verdadero o falso, justificando: “La FFT tiene costo N log N sólo para N potencia de 2”.",
  P("<b>Falso</b>, por el “sólo”. Hay que mostrar las dos razones:")
  + OL("Con el material de la cátedra: el Cooley–Tukey radix-2 (p20–22) se formula para N potencia de 2, pero si N no lo es se <b>rellena con ceros</b> hasta la siguiente potencia (p23). Esa potencia N′ es menor que 2N, así que el costo sigue siendo O(N′ log N′) = O(N log N).",
       f"Además, las implementaciones reales (NumPy/SciPy) factorizan N en primos pequeños (Cooley–Tukey de base mixta) y usan otro algoritmo O(N log N) para primos grandes: funcionan con cualquier N. Lo que cambia es la constante: 65537 (primo) tardó {f(TI[65537], 2)} ms y 65536 = 2¹⁶ {f(TI[65536], 2)} ms en esta máquina.")
  + P("Lo verdadero de la frase: el caso <b>más rápido</b> y el de la deducción de la clase es N potencia de 2."))

# ---------------- variantes rotadas
d(5, None, "a. ¿Por qué la serie de Fourier exige periodicidad y la TDF no? b. ¿Qué significa que “cuanto mayor N, mayor resolución en frecuencia”? c. Si para la TDF directa “por cada X[k] necesitamos N−1 sumas y N multiplicaciones”, deduzca su costo.",
  P("<b>a)</b> La serie representa la función como suma de armónicos de f₀ = 1/T: si la función no se repite no hay T ni f₀ (el caso T → ∞ lleva a la transformada continua). La TDF trabaja con N valores dados: la suma Σx[n]e^{−j2πkn/N} se puede calcular para cualquier secuencia finita, periódica o no (p16–17).")
  + P("<b>b)</b> Los N coeficientes X[k] están en las frecuencias k·f_s/N: el paso es Δf = f_s/N. Con más puntos (en más tiempo) el paso es más fino y se distinguen frecuencias cercanas. Si N crece sólo por subir f_s en el mismo intervalo, Δf no cambia (nyquist p5).")
  + P(f"<b>c)</b> Hay N coeficientes, cada uno con N multiplicaciones y N−1 sumas: N·N multiplicaciones y N(N−1) sumas → <b>O(N²)</b>. Para N = 4: {R['tdf4_ops']['mult']} multiplicaciones y {R['tdf4_ops']['sumas']} sumas."))
d(6, None, "Escriba la descomposición en pares e impares de la TDF que usa la FFT y explique qué propiedades de ω_N se aprovechan y por qué el costo resulta N log₂N.",
  CMD("X[k] = Σₙ₌₀^{N/2−1} x[2n]·ω_N^{k(2n)} + Σₙ₌₀^{N/2−1} x[2n+1]·ω_N^{k(2n+1)}\nω_N^{k(2n)} = ω_{N/2}^{kn},   ω_N^{k(2n+1)} = ω_N^k·ω_{N/2}^{kn}\nX[k] = E[k] + ω_N^k·O[k]     (E: TDF de los pares, O: TDF de los impares, de N/2 puntos)")
  + P("<b>Propiedades:</b> periodicidad (ω_N^{kn} = ω_N^{k(n+N)} = ω_N^{(k+N)n}), por la que E y O, de período N/2, sirven también para k ≥ N/2; simetría (ω_N^{k(N−n)} = (ω_N^{kn})*, o sea ω_N^{k+N/2} = −ω_N^k), de modo que X[k+N/2] = E[k] − ω_N^k·O[k]; y la identidad ω_N^{kN} = 1.")
  + P(f"<b>Costo:</b> cada nivel combina N/2 pares de valores con una multiplicación cada uno; se divide hasta TDF de 1 punto, en log₂N niveles: (N/2)·log₂N multiplicaciones → O(N log₂N). Verificado con una FFT recursiva: N = 8 → {MULT[8][2]}, N = 16 → {MULT[16][2]}, N = 1024 → {MULT[1024][2]} multiplicaciones (N² = {MULT[1024][1]}). Requiere poder partir a la mitad en cada nivel: N potencia de 2."))
d(13, None, "a. Defina frecuencia de muestreo y frecuencia de Nyquist. ¿De qué depende la frecuencia de Nyquist? b. ¿Qué ocurre con una componente de frecuencia mayor que f_N? c. ¿Qué es un filtro antialiasing y cuándo se aplica?",
  P("<b>a)</b> f_s es la cantidad de muestras por segundo (T_s = 1/f_s es el período de muestreo). La frecuencia de Nyquist es <b>f_N = f_s/2</b>: la máxima frecuencia de la que se puede obtener información. Depende <b>sólo de f_s</b>, porque hacen falta al menos dos muestras por ciclo para ver una oscilación (para f = 2 Hz, f_s de por lo menos 4 Hz).")
  + P("<b>b)</b> No desaparece: se <b>refleja</b> en [0, f_s/2] y aparece como una frecuencia menor que no existe en la señal (aliasing). Ej.: 23 Hz con f_s = 10 aparece en 3 Hz; 1,2 Hz con f_s = 2 aparece en 0,8 Hz. En imágenes da patrones de Moiré.")
  + P("<b>c)</b> Es un filtro <b>pasa bajo</b> que se aplica a la señal <b>antes de muestrear</b> para eliminar lo que esté por encima de f_s/2, de modo que no haya nada que se refleje. Después del muestreo no sirve: el alias ya se mezcló con frecuencias reales. (Agregado: la cátedra no lo define; ET25-T3 lo usa como opción y dejó incompleta la pregunta “El antialiasing ocurre cuando…”.)"))
d(9, None, "Explique cómo es el espectro 2D de una imagen con zonas homogéneas, con muchos bordes y con un patrón repetido. ¿Dónde están las bajas frecuencias en la salida de fft2 y para qué se usa fftshift?",
  P("Según la p26: <b>zonas homogéneas</b> → la energía del espectro se concentra en las <b>bajas</b> frecuencias; <b>bordes y transiciones frecuentes de gris</b> → componentes de <b>alta</b> frecuencia; <b>patrones que se repiten</b> → <b>picos</b> de igual intensidad separados a distancia regular, en la dirección perpendicular al patrón.")
  + P(f"Verificado con imágenes sintéticas: una mancha suave tiene el {f(100*FB['suave'], 1)} % de la energía dentro de radio 8 del centro, un tablero sólo el {f(100*FB['tablero'], 1)} %; franjas verticales de período 16 px dan picos en (0, ±8).")
  + P("La salida de fft2 deja la frecuencia (0,0) en la esquina superior izquierda: las bajas están en las <b>cuatro esquinas</b> y las altas en el centro (esquema p26). <b>fftshift</b> reordena los cuadrantes para llevar la frecuencia cero al centro, que es como se grafica y como se arma la máscara circular del filtro; <b>ifftshift</b> lo deshace antes de ifft2. Se grafica log(1 + |espectro|) porque el valor central es muchísimo mayor que el resto."))
d(7, None, "Compare aliasing y relleno con ceros (padding): qué los produce, qué efecto tienen en el espectro y si cambian la frecuencia de Nyquist.",
  TBL(["", "Aliasing", "Padding de ceros"], [
      ["Qué lo produce", "Muestrear con f_s &lt; 2·f_max (f_N menor que alguna frecuencia de la señal)", "Agregar ceros a las N muestras para llegar a una potencia de 2 (p23) u otra longitud"],
      ["Efecto en el espectro", "Aparecen picos en frecuencias falsas (reflejadas en [0, f_s/2]) y la componente real no se ve", "Grilla más fina (paso f_s/N′), espectro interpolado; no aparecen ni desaparecen componentes"],
      ["¿Cambia f_N?", "No: f_N = f_s/2 está fija; el problema es que la señal la supera", "No: f_s no cambia, así que f_N tampoco"],
      ["¿Se corrige después?", "No; se evita antes (más f_s o filtro antialiasing)", "No hay nada que corregir"]])
  + P("Por eso, en ET25-T1, “se rellenaron ceros … y esto generó aliasing, ya que modificó la frecuencia de Nyquist” es falso dos veces: el padding no modifica f_N y no genera aliasing."))

# ---------------- práctica real (evaluación TP3 y Moodle)
fb = M1["picos_bins"]
d(13, "ET25-T1 ej. 1 (= ET25-T3 ej. 7, TP3 ej. 3)",
  "Se realizaron 2 muestreos de una misma señal en [0, 6π]: muestra1.csv (dt = 0,5 s) y muestra2.csv (dt = 0,01 s). Con los espectros de la figura, responda a qué se debe que en el espacio de frecuencias la muestra 1 tiene un pico que no está presente en la muestra 2: a. En la muestra 1 se rellenaron ceros (padding) para obtener un valor que sea potencia de 2 y esto generó aliasing, ya que modificó la frecuencia de Nyquist original. b. El pico se debe al fenómeno de aliasing, resultado de una frecuencia mayor a la frecuencia de Nyquist en la muestra 1, que aparece reflejada como una frecuencia menor. c. Debido a la baja frecuencia de Nyquist de la muestra 1, el pico de 1.25 Hz se dividió en dos. Este fenómeno se conoce como Aliasing. d. Debido al fenómeno de aliasing, la muestra 2 no es capaz de detectar el pico que sí aparece claramente en la muestra 1. Esto se debe a la frecuencia de Nyquist.",
  P("<b>Respuesta: b.</b>")
  + OL(f"Muestra 1: n = {M1['n']}, dt = {f(M1['dt'])} s → f_s = {f(M1['fs'])} Hz, <b>f_N = {f(M1['fN'])} Hz</b>, Δf = {f(M1['df'], 4)} Hz. Muestra 2: n = {M2['n']}, dt = {f(M2['dt'])} s → f_s = {f(M2['fs'])} Hz, <b>f_N = {f(M2['fN'])} Hz</b>, Δf = {f(M2['df'], 4)} Hz.",
       f"<code>np.fft.fft</code> + <code>fftfreq(n, dt)</code>, f &gt; 0, normalizado. Muestra 2: picos en ≈0,125, 1, 1,2, 3 y 20 Hz (con zero padding sólo para ubicarlos: {', '.join(f(a, 3) for a, _ in M2['picos_finos'] if _ > 0.4)} Hz y uno chico en 20). Muestra 1: picos en ≈{f(fb[0][0], 2)} y ≈{f(fb[1][0], 2)} Hz y uno en el borde ≈1 Hz.",
       "La muestra 2 tiene f_N = 50 Hz: ve todas las componentes en su lugar. La muestra 1 tiene f_N = 1 Hz: la componente de 1,2 Hz la supera y se refleja en |1,2 − 2| = 0,8 Hz. Ese pico de 0,8 Hz es el que la muestra 2 no tiene: <b>aliasing en la muestra 1</b> (opción b).",
       f"Además (verif): 3 Hz cae en 1 Hz con muestras nulas (máx |10 sin(2π·3·0,5k)| = {R['m1_3Hz_max']:.0e}) y 20 Hz cae en 0 Hz (todas las muestras de cos(2π·20·0,5k) valen {', '.join(f(v) for v in R['m1_20Hz_vals'])}, una constante). La muestra 1 identifica 3 frecuencias (una falsa) y la 2, cinco.")
  + P("<b>Por qué no las otras:</b> a) el padding cambia N, no f_s; la frecuencia de Nyquist no se modifica y no hay aliasing por rellenar. c) el aliasing no “divide” un pico: lo traslada a otra frecuencia. d) invierte las muestras: la de mayor f_s (muestra 2) es la que ve bien; la que tiene el pico espurio es la 1.")
  + AGR("La señal que reproduce ambos CSV (error ~1e−14) es 5cos(2π·1t) + 5cos(2π·1,2t) + 10sin(2π·3t) + cos(2π·20t) − 2,5sin(2π·0,125t). La figura del examen tiene otros valores en el eje (picos en ≈0,5, 1,25, 1,5 y ≈18,5 Hz), así que no salió exactamente de estos CSV o de esta normalización; el razonamiento y la opción correcta son los mismos. En la foto el alumno marcó c, que es incorrecta."),
  fig=FIG["muestras"], practica=True)
d(10, "ET25-T2 ej. 1 (= ET25-T3 ej. 1)",
  "En la figura se ve una imagen original (tablero) a la izquierda y la misma imagen procesada a la derecha (franjas diagonales suaves). Responda cómo se obtuvo la imagen procesada: a. Se usó fft y se aplicó un filtro pasa alto, luego se usó ifft. b. Se usó fft y se aplicó un filtro pasa bajo, luego se usó ifft. c. La imagen de la derecha es el espacio de frecuencias obtenido al aplicar fft. d. Se usó fft haciendo un padding de ceros, y luego ifft.",
  P("<b>Respuesta: b (pasa bajo).</b> Recreado con un tablero sintético propio de 128×128 (casillas de 8 px, damero gris claro/blanco con diagonales oscuras):")
  + OL(f"fft2 + fftshift: los picos más cercanos al centro están a radio {f(R['tablero_picos'][0][2], 2)}, en (±4, ∓4): son las diagonales oscuras (período 32 px en diagonal). El damero aparece a radio {f(R['tablero_picos'][2][2], 2)} o más.",
       f"Máscara pasa bajo de radio D0 = {R['D0'][0]} (1 adentro, 0 afuera), ifftshift, ifft2, abs: sólo sobrevive la componente diagonal → <b>franjas diagonales suaves</b> (correlación con un coseno diagonal puro: {f(R['tablero_pb_corr_diag'], 3)}).",
       f"Con pasa alto (D0 = {R['D0'][1]}) sale lo contrario: imagen oscura con los bordes de las casillas.")
  + P("<b>Por qué no las otras:</b> a) el pasa alto deja bordes, no franjas suaves; c) un espectro centrado tiene un máximo brillante en el centro y puntos aislados con simetría central, no franjas; d) padding + ifft devuelve la misma imagen con ceros alrededor (verificado), no la suaviza. Justificación con la p27: el pasa bajo elimina “todo lo que no sean variaciones suaves del nivel de gris”."),
  fig=FIG["tab_pb"], practica=True)
d(10, "ET25-T3 ej. 2",
  "Imagen original: tablero. Imagen procesada: gris oscuro donde sólo quedan marcados los bordes de las casillas. ¿Cómo se obtuvo? (mismas opciones: fft + pasa alto + ifft / fft + pasa bajo + ifft / es el espacio de frecuencias / fft con padding de ceros + ifft).",
  P("<b>Respuesta: fft + filtro pasa alto + ifft.</b>")
  + P(f"Con el tablero sintético y pasa alto D0 = {R['D0'][1]} (1 − pasa bajo): se anula la frecuencia cero (el promedio, por eso el resultado es oscuro: media {f(R['tablero_pa_media'], 3)}) y las componentes suaves; quedan las transiciones bruscas de gris. Brillo medio en píxeles de borde de casilla: {f(R['tablero_pa_borde_vs_int'][0], 3)}; a 2 px o más del borde: {f(R['tablero_pa_borde_vs_int'][1], 3)}. Es lo que dice la p27: el pasa alta se usa para “quedarnos con las propiedades de la imagen en los que los niveles de gris varían bruscamente, por los bordes”.")
  + P("Pasa bajo daría franjas suaves (el ejercicio anterior); un espectro no se ve como la imagen con bordes; el padding no altera la imagen."),
  fig=FIG["tab_pa"], practica=True)
d(10, "ET25-T3 ej. 3",
  "Dada la siguiente figura (imagen, su espectro centrado, la máscara usada centrada y sin centrar, y el resultado), responda qué filtro se utilizó: se aplicó un filtro pasa alto / pasa bajo / antialiasing / de bordes. (La figura es un equivalente sintético propio de la del examen.)",
  P("<b>Respuesta: pasa alto.</b> Cómo leer la secuencia:")
  + OL("Espectro centrado (fftshift): la frecuencia cero está en el centro.",
       "Máscara centrada: disco <b>negro (0) en el centro</b> y blanco (1) alrededor → elimina las bajas frecuencias.",
       "Máscara sin centrar (ifftshift): el disco aparece partido en cuatro cuartos en las esquinas, que es donde están las bajas en la salida de fft2.",
       f"Resultado: oscuro (media {f(R['sint_pa_media'], 3)} contra {f(R['sint_media'], 3)} de la original) con sólo los contornos: típico del pasa alto.")
  + P("Si el disco central fuera blanco y el resto negro, sería pasa bajo y el resultado una imagen borrosa. El antialiasing se aplica antes de muestrear, no sobre el espectro de una imagen ya digital. “De bordes” describe el efecto, pero el filtro que la cátedra nombra (p24, p27) es <b>pasa alto</b>; si el enunciado aceptara ambos, la justificación es que el pasa alto es un detector de bordes."),
  fig=FIG["sint"], practica=True)
ev8 = {e["f"]: e for e in R["EV_P8"]}
d(11, "ET25-T3 ej. 4 / EV-P8",
  "Dada una señal compuesta de tres ondas senoidales con frecuencias f1 = 10 Hz, f2 = 20 Hz y f3 = 100 Hz. Supongamos que se realiza un muestreo con fs = 20 Hz. Responda qué frecuencias se podrán reconstruir correctamente: f1 = 10 Hz y f2 = 20 Hz / Solo f2 = 20 Hz / Solo f1 = 10 Hz / Solamente f3 = 100 Hz.",
  P("f_N = f_s/2 = 10 Hz. Alias (verif): 10 Hz → " + f(ev8[10]["alias"]) + " Hz (= f_N), 20 Hz → " + f(ev8[20]["alias"]) + " Hz, 100 Hz → " + f(ev8[100]["alias"]) + " Hz. f₂ y f₃ son múltiplos enteros de f_s: todas sus muestras valen lo mismo (constante, “0 Hz”). Quedan descartadas.")
  + P("<b>Dos lecturas para f₁ = 10 Hz = f_N</b> (escribí primero la de la cátedra):")
  + OL("<b>Convención de la cátedra (f ≤ f_N):</b> la diapositiva dice que para f = 2 Hz hace falta “una f_s de <b>por lo menos</b> 4 Hz”, o sea f_s = 2f alcanza. Entonces f₁ se reconstruye: <b>Solo f1 = 10 Hz</b>.",
       f"<b>Lectura estricta (f &lt; f_N):</b> en f = f_N hay 2 muestras por ciclo y lo que se ve depende de la fase: con sin(2π·10t) todas las muestras dan 0 (máx {f(ev8[10]['sin_max'], 3)}); con fase 30° la amplitud aparente es {f(LF[30], 3)}, con 45° {f(LF[45], 3)}, con 90° (coseno) {f(LF[90], 3)}. Estrictamente no se garantiza ninguna; como esa opción no está, la respuesta sigue siendo “Solo f1”, aclarando el caso límite.")
  + P(f"Para reconstruir las tres hace falta f_s &gt; 200 Hz. Con f_s = 21 Hz, 10 Hz ya queda bajo f_N, 20 Hz aparece en {f(R['fs21'][1][1])} Hz y 100 Hz en {f(R['fs21'][2][1])} Hz.")
  + TRAMPA("EV-P8 / ET25-T3 ej. 4: elegir “f1 y f2” confunde f_s con f_N (f₂ = 20 Hz es igual a f_s, no a f_N). La figura muestra las muestras con fase coseno: 10 Hz alterna ±1 y 20 y 100 Hz dan siempre 1."),
  fig=FIG["ev8"], practica=True)
d(8, "ET25-T3 ej. 5",
  "Muestra original: oscilación lenta con ruido de alta frecuencia superpuesto. Señal procesada: curva suave con la misma forma lenta. ¿Cómo se obtuvo la señal procesada? a. fft + filtro pasa alto; b. fft + filtro pasa bajo; c. se eliminaron las señales de aliasing con un filtro antialiasing; d. se usó fft haciendo un padding de ceros para modificar la frecuencia de Nyquist.",
  P("<b>Respuesta: b (pasa bajo).</b> Recreación: y = 2sin(2π·0,15t) + sin(2π·0,4t) + ruido gaussiano (desvío " + f(E5["std_ruido"], 3) + "), t ∈ [0, 20), dt = 0,01. Y = fft(y); máscara |f| &lt; " + f(E5["fc"]) + " Hz; ifft → curva suave. Error RMS contra la parte lenta: <b>" + f(E5["err_pb"], 3) + "</b>. Con pasa alto se obtiene casi sólo el ruido (error " + f(E5["err_pa"], 3) + ").")
  + P("Por qué: el ruido fino es alta frecuencia y la forma lenta es baja frecuencia; la p27 dice que el pasa baja se usa “para eliminar ruido de alta frecuencia”. c) el antialiasing actúa antes del muestreo y no se trata de alias; d) el padding no cambia la frecuencia de Nyquist ni elimina componentes."),
  fig=FIG["ej5"], practica=True)
d(8, "EV-P6",
  "Muestra original: oscilación lenta con una oscilación rápida superpuesta. Señal procesada: sólo una oscilación rápida de amplitud ≈1 centrada en 0. ¿Cómo se obtuvo? Se aplicó un filtro antialiasing / un filtro pasa-alto / un filtro pasa-bajo / se realizó un padding de ceros y se usó fft.",
  P("<b>Respuesta: filtro pasa-alto.</b> Recreación: 2sin(2π·0,12t) + sin(2π·3t). Pasa alto con f_c = " + f(E6["fc"]) + " Hz: queda la componente de 3 Hz con amplitud " + f(E6["amp_pa"], 2) + " (error RMS contra sin(2π·3t): " + f(E6["err_pa"], 4) + "). El pasa bajo dejaría la onda lenta de amplitud " + f(E6["amp_pb"], 2) + ".")
  + P("La procesada no tiene la tendencia lenta y está centrada en 0: se quitaron las bajas frecuencias. Un antialiasing (pasa bajo previo al muestreo) quitaría lo rápido, no lo lento; el padding no elimina componentes. En la captura del Moodle la opción marcada era pasa-alto, y coincide con este análisis."),
  fig=FIG["ev6"], practica=True)
e6 = {a: b for a, _, b in R["ET_ej6"]}
d(11, "ET25-T3 ej. 6",
  "Dada una señal compuesta de tres ondas senoidales con frecuencias f1 = 10 Hz, f2 = 15 Hz y f3 = 100 Hz, ¿qué frecuencia de muestreo fs debería usar para poder reproducirlas correctamente? fs = 10 Hz / 100 Hz / 20 Hz / 200 Hz.",
  P("Manda la máxima: f₃ = 100 Hz necesita f_N ≥ 100, o sea <b>f_s ≥ 200 Hz</b> (convención de la cátedra, “por lo menos” el doble). Respuesta: <b>200 Hz</b>.")
  + TBL(["f_s", "f_N", "10 Hz se ve en", "15 Hz se ve en", "100 Hz se ve en"],
        [[str(fs_), f(fs_ / 2), *[f(a) + " Hz" for _, a in e6[fs_]]] for fs_ in (10, 20, 100, 200)])
  + P("Con 200 Hz, 100 Hz está justo en f_N (caso límite, depende de la fase); en la práctica se elige algo mayor que 200. Las otras opciones aplican el doble a la frecuencia equivocada (20 = 2·10) o usan f_s igual a la máxima (100), que la manda a 0 Hz."),
  practica=True)
d(11, "EV-P2",
  "Una señal muestreada con frecuencia de muestreo fs = 10 Hz sólo puede representar correctamente frecuencias hasta: 2 fs / fs^2 / fs/2 / fs.",
  P("<b>f_s/2 = 5 Hz</b>, la frecuencia de Nyquist (nyquist p2: f_N = f_s/2, “No podemos obtener información de frecuencias mayores que f_N!”). Porque se necesitan al menos 2 muestras por ciclo: con 10 muestras por segundo, a lo sumo 5 ciclos por segundo.")
  + P("2f_s invierte la condición (f_s ≥ 2f); f_s² no tiene sentido físico; con f = f_s todas las muestras caen en la misma fase y la señal se ve constante (alias en 0 Hz)."),
  practica=True)
d(13, "EV-P3",
  "Responda en cuál de estos casos se observará el fenómeno de aliasing al realizar el muestreo de una señal: Cuando la FFT se calcula con pocos puntos / La frecuencia de Nyquist es menor que la de la señal original / Al aplicar un filtro pasa-bajo / La frecuencia de Nyquist es mayor que la de la señal original.",
  P("<b>Cuando la frecuencia de Nyquist es menor que la de la señal original.</b> Si la señal tiene componentes por encima de f_N = f_s/2, se reflejan en [0, f_s/2] (nyquist p7). Ej.: 1,2 Hz con f_s = 2 (f_N = 1) aparece en 0,8 Hz.")
  + P("Pocos puntos empeoran la resolución Δf = f_s/n (picos anchos), no crean frecuencias falsas. El pasa bajo, aplicado antes de muestrear, es lo que <b>evita</b> el aliasing. Si f_N es mayor que la frecuencia de la señal, el muestreo es correcto."),
  practica=True)
d(13, "EV-P5",
  "A partir de una misma señal se obtienen dos muestras con diferentes frecuencias de muestreo. Luego de aplicar FFT se obtiene la figura: ambas coinciden en los dos primeros picos (≈0,7 y ≈1,25 Hz), pero el tercero es distinto: la muestra 1 tiene un pico en ≈3 Hz y la muestra 2 en ≈23 Hz. ¿Por qué? 3er pico de la muestra 1 real y el de la muestra 2 por padding / 3er pico de la muestra 1 real y el de la muestra 2 por aliasing / 3er pico de la muestra 2 real y el de la muestra 1 por aliasing / 3er pico de la muestra 2 real y el de la muestra 1 por padding.",
  P("<b>Respuesta: el tercer pico de la muestra 2 es real y el de la muestra 1 es resultado del aliasing</b> (inferida del gráfico, con verificación numérica).")
  + OL("En el gráfico original la curva de la muestra 1 termina cerca de 5 Hz y la de la muestra 2 llega a 25 Hz: f_N₁ ≈ 5 (f_s₁ ≈ 10) y f_N₂ ≥ 25.",
       f"Una muestra con f_N = 5 no puede ver 3 Hz como reflejo de nada… salvo que haya algo por encima de 5: 23 Hz con f_s = 10 se refleja en |23 − 20| = {f(R['EV_P5_alias'][1][1])} Hz. Coincide exactamente con el pico de la muestra 1.",
       f"Recreado (figura): sin(2π·0,7t) + sin(2π·1,25t) + sin(2π·23t) con f_s = 10 da picos en {', '.join(f(a, 2) for a, _ in R['EV_P5']['m1']['picos'])} Hz; con f_s = 60 da {', '.join(f(a, 2) for a, _ in R['EV_P5']['m2']['picos'])} Hz.",
       "La muestra 2 no puede tener un alias en 23 Hz si su f_N es ≥ 25; el padding no traslada picos a otra frecuencia.")
  + P(f"Qué queda abierto: f_s₁ no está escrita; si fuera 20 Hz, 23 también caería en {f(R['EV_P5_alias'][2][1])} Hz y la conclusión es la misma. En la captura el alumno marcó “3er pico de la muestra 1 real y el de la 2 por aliasing”, que invierte el caso."),
  fig=FIG["ev5"], practica=True)

# ---------------- TP3
d(15, None, "TP3 ej. 1: muestreo.ipynb. a. Investigue la frecuencia de Nyquist: ¿de qué depende? Calcúlela en la notebook. b. Cambie la frecuencia de muestreo de 60 a 10, 4, 2 y 1 y explique qué ocurre con las señales (f1 = 1 Hz y f2 = 2 Hz, desfase π/2, t en [0, 2]). c. Aplique fft y fftfreq a la señal suma (1, 2 y 5 Hz con fs = 60).",
  P("<b>a)</b> f_N = f_s/2; depende sólo de f_s. Con f_s = 60, f_N = 30 Hz.")
  + P("<b>b)</b>")
  + TBL(["f_s", "f_N", "f₁ = 1 Hz", "f₂ = 2 Hz"], [
      ["60", "30", "60 muestras/ciclo: perfecta", "30 muestras/ciclo: perfecta"],
      ["10", "5", "10 muestras/ciclo: bien", "5 muestras/ciclo: reconocible, picos recortados"],
      ["4", "2", "4 muestras/ciclo: triángulo 1, 0, −1, 0", "2 muestras/ciclo (f = f_N): zigzag ±1"],
      ["2", "1", "f = f_N: zigzag ±1", f"f = f_s: constante 1 (alias {f(MU[2]['f2_alias'])} Hz)"],
      ["1", "0,5", f"f = f_s: constante 1 (alias {f(MU[1]['f1_alias'])} Hz)", "f = 2f_s: constante 1"]])
  + P("Cuanto más baja f_s, menos muestras por ciclo: primero la forma se vuelve poligonal, en f = f_N queda un zigzag que sólo conserva la amplitud por la fase coseno, y por encima de f_N la señal aparece con otra frecuencia (aquí 0 Hz: una constante). Las dos señales llegan a verse iguales: aliasing.")
  + P(f"<b>c)</b> <code>y_hat = np.fft.fft(y); freqs = np.fft.fftfreq(n, 1/60); cond = freqs &gt; 0; plt.plot(freqs[cond], np.abs(y_hat[cond]))</code>. Con n = {R['muestreo_c2']['n']}, Δf = {f(R['muestreo_c2']['df'], 3)} Hz; picos en {', '.join(f(a, 2) for a, _ in R['muestreo_c2']['picos'])} Hz (las frecuencias 1, 2 y 5 redondeadas a la grilla)."),
  fig=FIG["mgrid"], practica=True)
d(7, None, "TP3 ej. 2: ejecute tiempo_fft.ipynb (np.fft.fft(a, len(a)) para n = 166464, 166400, 2¹⁶, 2¹⁷ y 2¹⁶ + 1), analice el algoritmo y explique los resultados.",
  TBL(["n", "Factorización", "ms (esta máquina)"], [[str(dd["n"]), FACT[dd["n"]], f(dd["ms"], 2)] for dd in R["tiempos"]])
  + P(f"Los tiempos absolutos dependen de la máquina (NumPy {R['numpy']}); el orden es lo que importa. <b>2¹⁶ y 2¹⁷</b> son los más rápidos: el algoritmo divide a la mitad en cada nivel (Cooley–Tukey radix-2, p22), y 2¹⁷ tarda ≈{f(TI[131072]/TI[65536], 1)} veces 2¹⁶, como N log₂N (algo más de 2). <b>166400 = {FACT[166400]}</b> y <b>166464 = {FACT[166464]}</b> no son potencias de 2 pero se factorizan en primos chicos (2, 3, 5, 13, 17): el algoritmo de base mixta los divide igual y tardan lo esperable para su tamaño (≈{f((TI[166400] + TI[166464]) / 2 / TI[131072], 1)} veces 2¹⁷ en esta corrida); entre ellos la diferencia es chica y cambia entre corridas. <b>65537</b> es primo: no se puede dividir y la librería usa otro algoritmo (Bluestein) que pasa por FFT de mayor tamaño; tarda ≈{f(TI[65537]/TI[65536], 1)} veces lo de 65536 con un solo punto más.")
  + P("Conclusión: la fft acepta cualquier N, pero el costo depende de la factorización; por eso conviene N potencia de 2 o, si no, rellenar con ceros hasta una (p23)."),
  fig=FIG["tiempos"], practica=True)
d(15, None, "TP3 ej. 3: muestra1.csv y muestra2.csv, misma señal en [0, 6π]. Grafique el espectro de cada muestra. a. ¿Cuántas frecuencias logra identificar de cada muestreo? b. ¿Qué diferencias respecto a lo esperado encuentra? ¿Por qué? c. ¿Cuál de las dos opciones de muestreo es mejor? Concluya.",
  P(f"Datos: muestra 1 n = {M1['n']}, dt = {f(M1['dt'])} (f_s = {f(M1['fs'])}, f_N = {f(M1['fN'])}); muestra 2 n = {M2['n']}, dt = {f(M2['dt'])} (f_s = {f(M2['fs'])}, f_N = {f(M2['fN'])}).")
  + P("<b>a)</b> Muestra 2: <b>cinco</b>, en ≈0,125, 1, 1,2, 3 y 20 Hz. Muestra 1: <b>tres</b> picos, en ≈0,13, ≈0,8 y en el borde ≈1 Hz.")
  + P("<b>b)</b> La muestra 1 no puede ver nada por encima de 1 Hz. Reflejos (verif): " + "; ".join(f"{f(a)} Hz → {f(b)} Hz" for a, b in R["tp3_alias_m1"]) + ". El pico de 0,8 Hz no existe en la señal: es el alias de 1,2 Hz. El de 3 Hz cae en 1 Hz pero con fase seno sus muestras son nulas; el de 20 Hz cae en 0 Hz (suma una constante). Es <b>aliasing</b>, porque f_N = 1 Hz es menor que varias frecuencias de la señal.")
  + P("<b>c)</b> El muestreo 2: f_N = 50 Hz cubre la componente más alta (20 Hz). Las dos tienen resolución parecida (Δf ≈ 0,05 Hz) porque cubren el mismo intervalo: lo que las distingue es f_s.")
  + AGR("La señal exacta (hallada por mínimos cuadrados; la cátedra no la da): 5cos(2π·1t) + 5cos(2π·1,2t) + 10sin(2π·3t) + cos(2π·20t) − 2,5sin(2π·0,125t)."),
  fig=FIG["muestras"], practica=True)
d(9, None, "TP3 ej. 4: dadas cuatro imágenes (I triángulo de Sierpinski de trazos, II foto de un parque con niebla, III espiral fractal, IV código QR) y cuatro espectros (A: cruz horizontal/vertical fina con rayas; B: cuadrícula de bandas horizontales y verticales; C: estrella de tres rectas a 60°; D: mancha central difusa), determine qué espectro corresponde a cada imagen.",
  P("Las imágenes del TP3 son de terceros y no se reproducen; el emparejamiento se razona con la p26 (respuesta <b>inferida</b>, sin clave oficial):")
  + TBL(["Imagen", "Espectro", "Por qué"], [
      ["I (triángulo)", "C", "Tres familias de bordes rectos a 60° → tres rectas de energía perpendiculares a esos bordes (estrella)"],
      ["II (niebla)", "D (sin rótulo)", "Zonas homogéneas, pocos bordes → energía concentrada en las bajas (centro)"],
      ["III (espiral)", "A", "Muchas direcciones mezcladas, sin bordes rectos dominantes; la cruz suele venir de los bordes de la imagen"],
      ["IV (QR)", "B", "Bloques cuadrados alineados con los ejes → energía en bandas horizontales y verticales regulares"]])
  + P("Regla para el examen: bordes rectos → rectas de energía perpendiculares; patrón periódico → picos regulares; imagen suave → centro. El emparejamiento III–A es el menos seguro: se deja por descarte."),
  practica=True)
d(10, None, "TP3 ej. 5: complete imagen_con_filtro.ipynb (espectro, espectro centrado, reconstrucción del espectro y de la imagen), explique cada función, use el pasa bajo, implemente el pasa alto y explique las imágenes obtenidas.",
  CMD("img_c2 = np.fft.fft2(img)          # TDF 2D (FFT): espectro de la imagen\nimg_c3 = np.fft.fftshift(img_c2)   # lleva la frecuencia (0,0) al centro\nimg_c4 = np.fft.ifftshift(img_c3 * filtro)   # deshace el centrado (ya filtrado)\nimg_c5 = np.fft.ifft2(img_c4)      # vuelve al espacio; se grafica np.abs(img_c5)\n\ndef filtro_pasa_alto(D0, imgShape):\n    return 1 - filtro_pasa_bajo(D0, imgShape)\n# equivalente: base[y, x] = 1 si distancia((y, x), centro) &gt;= D0")
  + P("<b>c) Pasa bajo:</b> la imagen obtenida es una versión <b>suavizada / borrosa</b> de la original: se conservan las variaciones suaves del gris y se pierden bordes y detalles (lo que está fuera del círculo D0). Con D0 más chico, más borrosa.")
  + P(f"<b>e) Pasa alto:</b> se obtiene una imagen <b>oscura con los bordes</b> (contornos, texturas finas): se eliminó la frecuencia cero (brillo medio) y las variaciones suaves. Ejemplo sintético: media original {f(R['sint_media'], 3)}, tras pasa alto {f(R['sint_pa_media'], 3)}. Los espectros se grafican con np.log(1 + np.abs(·)) porque el valor central domina."),
  fig=FIG["tab4"], practica=True)
d(7, None, "TP3 ej. 6: a. Averiguar cómo es implementada la función fft en la librería NumPy. b. ¿Cómo es definida en otras librerías? Explicar.",
  AGR("No está en el material: es lo que se averigua.")
  + P(f"<b>a)</b> NumPy (desde 1.17; aquí {R['numpy']}) usa <b>pocketfft</b>: Cooley–Tukey de <b>base mixta</b> (descompone N en factores 2, 3, 4, 5, 7, … y aplica pasos de cada base) y, si queda un factor primo grande, el algoritmo de <b>Bluestein</b>, que convierte la TDF en una convolución que se calcula con FFT de tamaño potencia de 2. Así acepta cualquier N con costo O(N log N); por eso 65537 (primo) es más lento que 65536 pero funciona. No rellena con ceros salvo que se pase <code>n</code>. Define X[k] = Σx[n]e^{{−j2πkn/N}} (signo −, sin normalizar; ifft divide por N).")
  + P("<b>b)</b> <b>SciPy</b> (<code>scipy.fft</code>): también pocketfft (versión C++), con opciones extra (<code>workers</code> para paralelizar, <code>next_fast_len</code> para elegir un N rápido, transformadas reales y de coseno). <b>FFTW</b> (C, usado por MATLAB y otras): arma un “plan” que elige el algoritmo según N y la máquina. Todas calculan la misma TDF; difieren en velocidad, normalización y opciones."),
  practica=True)

# ---------------- práctica de predicción
d(5, None, "Calcule a mano la TDF de x[n] = 1.2, 1.4, 1.3, 1.0 (N = 4, datos de la diapositiva) con X[k] = Σₙ₌₀³ x[n]·ω₄^{nk}, ω₄ = e^{−j2π/4}. ¿Cuántas multiplicaciones complejas hizo?",
  P("ω₄^m para m = 0, 1, 2, 3: " + "; ".join(f"{f(a)} {'−' if b < 0 else '+'} {f(abs(b))}j" for _, a, b in R["tdf4_w"]) + " (o sea 1, −j, −1, j; después se repiten).")
  + TBL(["k", "Cuenta", "X[k]"], [
      ["0", "1,2 + 1,4 + 1,3 + 1,0", X4s[0]],
      ["1", "1,2 + 1,4(−j) + 1,3(−1) + 1,0(j)", X4s[1]],
      ["2", "1,2 + 1,4(−1) + 1,3(1) + 1,0(−1)", X4s[2]],
      ["3", "1,2 + 1,4(j) + 1,3(−1) + 1,0(−j)", X4s[3]]])
  + P(f"Comprobado con np.fft.fft. X[0] es la suma (componente de frecuencia 0); X[3] es el conjugado de X[1] porque la señal es real (simetría). Operaciones: 4 coeficientes × 4 multiplicaciones = {R['tdf4_ops']['mult']}, y {R['tdf4_ops']['sumas']} sumas → O(N²). La IFFT devuelve {'; '.join(f(v, 1) for v in R['tdf4_ifft'])}.")
  + TRAMPA_P("Que pidan una TDF chica a mano para mostrar el O(N²), o que pidan qué X[k] son conjugados entre sí."),
  practica=True)
d(3, None, "Para s(t) = sin(t) + 1,2·cos(2t) calcule la amplitud, la fase y la frecuencia fundamental de cada componente, diga cuál es más rápida y dibuje el espectro.",
  TBL(["Componente", "A", "θ", "ω [rad/s]", "f₀ = ω/2π [Hz]", "T [s]"], [
      ["sin(t)", "1", "0", "1", f"1/2π = {f(ES['f1'], 4)}", f(ES["T1"], 4)],
      ["1,2·cos(2t)", "1,2", "0 (coseno)", "2", f"1/π = {f(ES['f2'], 4)}", f(ES["T2"], 4)]])
  + P(f"Más rápida: 1,2cos(2t), por tener mayor frecuencia (menor período: “Si T &lt; entonces F &gt;”). Espectro: dos barras, altura 1 en {f(ES['f1'], 3)} Hz y 1,2 en {f(ES['f2'], 3)} Hz (verif con FFT: {PK(ES['picos_fft'])}).")
  + P("Sobre la fase: la diapositiva pone θ₂ = 0 tomando el coseno como referencia; si se escribe como seno, 1,2cos(2t) = 1,2sin(2t + π/2), fase π/2.")
  + TRAMPA_P("Confundir ω con f: decir que las barras están en 1 y 2 Hz."),
  fig=FIG["ej_s"], practica=True)
d(1, None, "Calcule los coeficientes de Fourier del escalón f(x) = −1 en (−π, 0), 1 en (0, π), y explique cuántos términos hacen falta para aproximarlo de manera exacta.",
  P("Con las fórmulas de la p6 (período 2π): a₀ = (1/π)∫f = 0 (simétrico); aₙ = (1/π)∫f·cos(nx) = 0 (f impar por cos par); bₙ = (1/π)∫f·sin(nx) = (2/π)∫₀^π sin(nx)dx = 2(1 − cos nπ)/(nπ): 4/(nπ) si n es impar, 0 si es par.")
  + TBL(["n", "bₙ (integración numérica)", "4/(nπ)"], [[str(n), f(b, 4), f(e, 4)] for n, a, b, e in R["escalon_coef"]])
  + P("f(x) ≈ (4/π)[sin x + sin 3x/3 + sin 5x/5 + …]. Hacen falta <b>infinitos</b> términos: con cualquier cantidad finita queda error. Máximo de la suma parcial cerca del salto: " + "; ".join(f"{K} término(s): {f(m, 4)}" for K, m, _ in R["escalon_gibbs"]) + ": el sobrepico tiende a ≈" + f(GIBBS, 3) + " y no desaparece (fenómeno de Gibbs, agregado).")
  + TRAMPA_P("Responder la pregunta de la p5 (“¿Cuántas sumas para aproximar de manera exacta la función escalón?”) con un número finito."),
  fig=FIG["escalon"], practica=True)
d(12, None, "Se suman senoidales de 0,5, 1 y 2 Hz. Calcule n y Δf = f_s/n para: (i) [0, 2π) con f_s = 10; (ii) [0, 2π) con f_s = 100; (iii) [0, 6π) con f_s = 6. ¿Qué cambio mejora la resolución? ¿Hasta qué frecuencia llega el eje en cada caso?",
  TBL(["Caso", "n", "Δf [Hz]", "Última f &gt; 0 [Hz]", "f_N", "Picos (verif) [Hz]"],
      [[c_, str(dd["n"]), f(dd["df"], 3), f(dd["fmax"], 3), f(dd["fN"], 1), ", ".join(f(a, 2) for a, _ in dd["picos"])] for c_, dd in zip(("(i)", "(ii)", "(iii)"), RES)])
  + P("(i) y (ii) tienen el mismo Δf: multiplicar f_s por 10 multiplica n por 10 y Δf = f_s/n no cambia (“Si mantenemos el rango (a, b), aumentar f_s no mejora la resolución!”). (iii) triplica la duración → Δf ≈ 1/3: los picos caen en su lugar y son finos. El eje llega hasta ≈f_N (5, 50 y 3 Hz): es la “frecuencia” que señala la diapositiva."),
  fig=FIG["resol"], practica=True)
d(13, None, "Calcule en qué frecuencia aparece cada senoidal al muestrear: 3 Hz con f_s = 4; 6 Hz con f_s = 10; 9 Hz con f_s = 10; 23 Hz con f_s = 10 y 20; 1,2, 3 y 20 Hz con f_s = 2; 100 Hz con f_s = 20; 60 Hz con f_s = 50.",
  P("Regla (reproduce el gráfico de la p7): f_alias = |f − f_s·round(f/f_s)|, que siempre cae en [0, f_s/2].")
  + TBL(["f [Hz]", "f_s [Hz]", "f_N", "Aparece en [Hz]"], [[f(a), f(b), f(b / 2), f(c)] for a, b, c in AT])
  + P("Las que ya están por debajo de f_N no cambian (ninguna en esta lista); las que son múltiplo de f_s caen en 0 Hz (se ven constantes).")
  + TRAMPA_P("Dar f − f_s una sola vez (23 − 10 = 13), que queda fuera del eje [0, f_s/2]."),
  fig=FIG["alias"], practica=True)
d(6, None, "Complete la tabla de costos de la TDF directa y de la FFT para N = 1000, 10⁶ y 10⁹, con 1 ns por operación, y el cociente entre ambos. ¿Cuántas multiplicaciones complejas hace una FFT radix-2 con N = 1024?",
  TBL(["N", "N²", "Tiempo TDF", "N log₂N", "Tiempo FFT", "N²/(N log₂N)"], [
      [f(c_["N"], 0), f"{c_['N2']:.0e}".replace("e+", "×10^"), f"{f(c_['seg_N2'], 3)} s" + (f" ≈ {f(c_['anios_N2'], 1)} años" if c_["N"] == 10**9 else ""),
       f(c_["NlogN"], 0), f"{f(c_['seg_NlogN'], 6)} s", f(c_["cociente"], 0)] for c_ in C])
  + P(f"La diapositiva redondea N log₂N a 10⁴, 20×10⁶ y 30×10⁹ y da “~31.2 años” (la cuenta da {ANIOS}). FFT radix-2 con N = 1024: (N/2)·log₂N = 512·10 = <b>{MULT[1024][2]}</b> multiplicaciones (contadas en una implementación recursiva), contra N² = {MULT[1024][1]}."),
  practica=True)


unit = {
    "id": 3, "titulo": "Aproximación polinómica trigonométrica: Fourier, FFT, Nyquist y filtros", "corto": "Fourier y FFT",
    "fuente": "trig.pdf (27 diapositivas) + nyquist.pdf (10 diapositivas) + TP3", "peso": 2,
    "intro": intro, "secciones": S, "conexiones": conexiones, "fidelidad": fidelidad,
    "auditoria": auditoria, "glosario": glosario, "quiz": Q, "desarrollo": DEV,
}
if __name__ == "__main__":
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "m3.js")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("ASO.addUnit(" + json.dumps(unit, ensure_ascii=False, indent=1) + ");\n")
    figs_items = [it for it in Q + DEV if it.get("fig")]
    print("ok", len(Q), "preguntas,", len(S), "secciones,", len(DEV), "desarrollo;",
          len(figs_items), "ítems con fig;", os.path.getsize(out) // 1000, "KB el archivo")
