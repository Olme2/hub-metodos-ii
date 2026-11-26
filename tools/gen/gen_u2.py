#!/usr/bin/env python3
# Genera data/m2.js — U2: Ajuste de curvas (Métodos Numéricos II 2026).
# Todo número sale de verif_u2.py (se importa R y los datasets de las figuras).
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verif_u2 import R, D
from figs import nueva, svg

def n(x, d=3):
    s = f"{x:.{d}f}"
    if s.startswith("-"): s = "−" + s[1:]
    return s

def D_(term, txt, src):
    return f'<div class="def"><b>{term}</b><p>{txt}</p><span class="src">{src}</span></div>'
def AGR(txt):
    return f'<div class="agregado"><b>Agregado</b><p>{txt}</p></div>'
def TRAMPA(txt, real=True):
    rot = "Así lo toma la cátedra" if real else "Así lo tomaría la cátedra"
    return f'<div class="trampa"><b>{rot}</b><p>{txt}</p></div>'
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
FIG = lambda s: f'<figure class="fig">{s}</figure>'

ej, cv, et, F = R["ej"], R["covid"], R["et1"], R["figs"]

# ------------------------------------------------------------------ figuras
def recta(ax, x, a0, a1, lab="Ajuste lineal (MC)"):
    xx = np.linspace(min(x), max(x), 2); ax.plot(xx, a0 + a1 * xx, color="#d62728", label=lab)

def fig_scatter(key, conrecta=True, titulo=None, xl="x", yl="y"):
    x, y = D[key]; fig, ax = nueva()
    ax.scatter(x, y, s=14, label="Datos")
    if conrecta: recta(ax, x, F[key]["a0"], F[key]["a1"])
    if titulo: ax.set_title(titulo)
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.legend(loc="best", fontsize=8)
    return svg(fig)

def fig_it3():
    return fig_scatter("it3", titulo="Nube de puntos y ajuste lineal por mínimos cuadrados")

def fig_et1():
    x = np.array([3, 6, 7, 8, 9.]); y = np.array([47, 493, 1215, 3132, 8292.])
    fig, ax = nueva(); ax.scatter(x, y, color="k", s=18, label="Datos")
    xx = np.linspace(3, 9, 2); ax.plot(xx, et["int_examen"] + et["coef_examen"] * xx, color="#333", label="Ajuste lineal")
    for a, b in zip(x, y): ax.annotate(f"({int(a)}, {int(b)})", (a, b), textcoords="offset points", xytext=(4, 4), fontsize=8)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.legend(loc="upper left", fontsize=8)
    return svg(fig)

def fig_et2():
    fig, axs = nueva(1, 2, alto=3.0)
    x, y = D["expo"]; axs[0].scatter(x, y, s=12); axs[0].set_title("Gráfico 1")
    x, y = D["parab"]; axs[1].scatter(x, y, s=12); axs[1].set_title("Gráfico 2")
    for a in axs: a.set_xlabel("x"); a.set_ylabel("y")
    return svg(fig)

def fig_tp4():
    t = R["tp2e4"]; fig, ax = nueva(); ax.scatter(t["x"], t["y"], s=8, label="datos.csv")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.legend(fontsize=8); ax.set_title("TP2 ej. 4: datos")
    return svg(fig)

def fig_tp4_ajuste():
    t = R["tp2e4"]; fig, ax = nueva(); ax.scatter(t["x"], t["y"], s=8, label="datos.csv")
    xx = np.linspace(0, 10, 100)
    ax.plot(xx, np.polyval(t["tab"][0]["c"], xx), color="#7f7f7f", label="grado 1")
    ax.plot(xx, np.polyval(t["tab"][1]["c"], xx), color="#d62728", label="grado 2")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.legend(fontsize=8)
    return svg(fig)

def fig_p11():
    fig, axs = nueva(2, 2, alto=5.2)
    tit = {"lin": "Tendencia lineal: recta OK", "ruido": "Sin tendencia: r ≈ 0",
           "expo": "Exponencial: recta mala", "seno": "Periódica: recta mala"}
    for ax, k in zip(axs.flat, ["lin", "ruido", "expo", "seno"]):
        x, y = D[k]; ax.scatter(x, y, s=8); recta(ax, x, F[k]["a0"], F[k]["a1"], "recta MC")
        ax.set_title(f"{tit[k]} (r = {n(F[k]['r'], 2)})", fontsize=9)
    return svg(fig)

def fig_lineal():
    x, y = D["expo"]; fig, axs = nueva(1, 2, alto=3.0)
    axs[0].scatter(x, y, s=10); axs[0].set_title("y vs x"); axs[0].set_xlabel("x")
    axs[1].scatter(x, np.log(y), s=10); axs[1].set_title("ln y vs x (recta)"); axs[1].set_xlabel("x")
    return svg(fig)

def fig_residuos_sec():
    x = np.arange(1, 7); y = np.array([0, 0.6, 1.77, 1.92, 3.31, 3.52])
    fig, ax = nueva(); ax.scatter(x, y, s=18, label="Datos (ejemplo p9)")
    recta(ax, x, ej["a0"], ej["a1"], f"F(x) = {n(ej['a0'])} + {n(ej['a1'])}x")
    for a, b in zip(x, y): ax.plot([a, a], [b, ej["a0"] + ej["a1"] * a], color="#1f77b4", lw=1)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.legend(fontsize=8)
    return svg(fig)

def fig_resid_q():
    x, y = D["curva"]; r = y - (F["curva"]["a0"] + F["curva"]["a1"] * x)
    fig, ax = nueva(); ax.axhline(0, color="k", lw=1); ax.scatter(x, r, s=16)
    ax.set_title("Residuos Rᵢ = yᵢ − (a₀ + a₁xᵢ) de un ajuste lineal"); ax.set_xlabel("x"); ax.set_ylabel("residuo")
    return svg(fig)

def fig_loglog():
    x, y = D["pot"]; fig, axs = nueva(1, 2, alto=3.0)
    axs[0].scatter(x, y, s=10); axs[0].set_title("y vs x"); axs[0].set_xlabel("x")
    axs[1].scatter(np.log(x), np.log(y), s=10); axs[1].set_title("ln y vs ln x"); axs[1].set_xlabel("ln x")
    return svg(fig)

def fig_recip():
    x, y = D["recip"]; fig, axs = nueva(1, 2, alto=3.0)
    axs[0].scatter(x, y, s=10); axs[0].set_title("y vs x"); axs[0].set_xlabel("x")
    axs[1].scatter(x, 1 / y, s=10); axs[1].set_title("1/y vs x"); axs[1].set_xlabel("x")
    return svg(fig)

def fig_expo_q():
    return fig_scatter("expo", titulo=f"Ajuste lineal por MC (r = {n(F['expo']['r'], 2)})")

def fig_ruido_q():
    return fig_scatter("ruido", titulo="Ajuste lineal por MC")

def fig_parab_q():
    return fig_scatter("parab", conrecta=False, titulo="Datos")

# ------------------------------------------------------------------ secciones
S = []
def sec(i, titulo, estrellas, pags, *html):
    S.append({"id": f"m2-s{i}", "titulo": titulo, "estrellas": estrellas, "pags": pags, "html": "".join(html)})

sec(1, "Aproximación de funciones: objetivo, interpolación vs ajuste", 2, "p2–4",
    P("Tenés una nube de puntos medidos y querés una función que la describa para entender la tendencia y predecir. La pregunta es si la función tiene que pasar por todos los puntos (interpolar) o sólo seguir la tendencia con cierto error (ajustar)."),
    D_("Objetivo", "Encontrar la función que describe mejor los datos.", "U2 p2"),
    H4("Exploración visual de los datos (p2)"),
    UL("¿Hay alguna tendencia?", "¿Cuáles son los valores mínimos y máximos presentes?", "¿Existen valores alejados (outliers)?",
       "¿Qué información me da esta nube de puntos?", "¿Puedo hacer alguna predicción?"),
    P("Siempre se grafica primero: el ejemplo de la p9 arranca con \"1- Graficamos los datos\" y recién después elige el modelo."),
    H4("Dos maneras de aproximar"),
    TBL(["", "Interpolación (p3)", "Ajuste de curvas (p4)"], [
        ["Qué hace", "\"Aseguro que la aproximación para las abcisas de los datos coinciden con la ordenadas\": pasa por todos los puntos.", "\"Ajuste de la curva con otra función\" (en el ejemplo, lineal) que no tiene por qué pasar por los puntos."],
        ["Tendencia", "La diapositiva pregunta \"¿Describe la tendencia de los datos?\": con datos ruidosos, no; la polilínea sigue cada salto.", "\"¿Puedo ver la tendencia de los datos?\": sí, la recta la muestra."],
        ["Error en los datos", "Error cero en los puntos dados (copia también el ruido).", "Acepta cierto error en cada punto a cambio de describir el comportamiento."],
        ["Pregunta clave", "—", "\"¿Cuánto me interesa que pase exactamente por los puntos dados?\""],
    ]),
    D_("Idea central de la p4", "Muchas veces es más útil ajustar con una función lineal (recta) que describa el comportamiento de los datos (con cierto error) y no por ejemplo interpolar con un polinomio de orden alto.", "U2 p4"),
    AGR("Por qué: los datos experimentales traen error de medición. Interpolar obliga al polinomio a reproducir ese error, y un polinomio de orden alto oscila entre los nodos (fenómeno de Runge). El ajuste por mínimos cuadrados usa muchos más datos que parámetros y así promedia el ruido."),
    TRAMPA("IT25-T1 1a preguntó la diferencia entre \"ajuste lineal\" y \"ajuste lineal por mínimos cuadrados\". Un ajuste lineal es cualquier recta propuesta para los datos (hay infinitas, p5 muestra tres). \"Por mínimos cuadrados\" fija el criterio que elige una sola: la que minimiza E = Σ(yₙ − (a₀ + a₁xₙ))²."),
)

sec(2, "Regresión lineal: el modelo y el error a minimizar", 3, "p5–6",
    P("Sobre una nube de puntos se pueden trazar muchas rectas (la p5 dibuja tres). Hace falta un criterio numérico para decidir cuál es \"la mejor\": ese criterio es el error."),
    CMD("F(x, aᵢ) = a₀ + a₁x        a₀: ordenada al origen    a₁: coef lineal   (p5)"),
    D_("Objetivo de la regresión lineal", "Encontrar una función \"lineal\" F que aproxime a los datos con el menor error posible. Encontrar los coeficientes de tal manera que el error (distancia de los datos a la recta) sea mínima.", "U2 p6"),
    P("Los datos son los puntos (xᵢ, yᵢ). El error que define la cátedra:"),
    CMD("E(α₁, α₂, …, αₙ) = [ Σₙ₌₁ᴺ ( yₙ − F(x, α₁, …, αₙ) )² ]^(1/2)       (forma general, p6)\n\nPara que sea mínimo:   E = Σₙ₌₁ᴺ ( yₙ − (a₀ + a₁xₙ) )²       (dato − recta)"),
    UL("La \"distancia\" es <b>vertical</b>: dato yₙ menos el valor de la recta en xₙ (el gráfico de la p6 marca segmentos verticales).",
       "Se elevan al cuadrado para que los errores positivos y negativos no se cancelen y para penalizar más los grandes.",
       "Minimizar la forma con raíz o sin raíz da los mismos coeficientes (la raíz es creciente); por eso se trabaja con la suma de cuadrados."),
    AGR("La forma general de la p6 vale para <b>cualquier</b> F(x, α): el criterio de mínimos cuadrados se puede plantear siempre. Lo que no siempre se puede es resolverlo con un sistema lineal: eso sólo pasa cuando F es lineal en los parámetros (recta, polinomio, plano). Para F no lineal en los parámetros (y = a₁e^{b₁x}) la cátedra linealiza (p12); la alternativa es mínimos cuadrados no lineales, iterativo (por ejemplo <code>scipy.optimize.curve_fit</code>), que no está en las diapositivas."),
    TRAMPA("F18-P4 pide el objetivo de aproximar por mínimos cuadrados en regresión lineal: es esta diapositiva. Escribí F, el error E con los cuadrados y la frase \"encontrar los coeficientes que hacen mínimo E\"."),
)

sec(3, "Regresión lineal por mínimos cuadrados: condición de mínimo y ecuaciones normales", 3, "p7–8",
    P("Ya está definido el error E(a₀, a₁). Falta encontrar los a₀, a₁ que lo hacen mínimo: se deriva respecto de cada coeficiente y se iguala a cero."),
    D_("Mínimos cuadrados", "Encontrar aᵢ, tal que E sea mínimo (mínimos cuadrados).", "U2 p7"),
    CMD("E = Σₙ₌₁ᴺ ( yₙ − (a₀ + a₁xₙ) )²\n\nCondición de mínimo:   ∂E/∂a₀ = 0   y   ∂E/∂a₁ = 0\n\n∂E/∂a₀ = −2 Σ ( yₙ − a₀ − a₁xₙ ) = 0\n∂E/∂a₁ = −2 Σ ( yₙ − a₀ − a₁xₙ ) xₙ = 0"),
    P("\"Al resolver este sistema de ecuaciones encuentro los coeficientes que buscaba\" (p7). Ordenado (p8):"),
    CMD("N a₀ + (Σxₙ) a₁ = Σyₙ\n(Σxₙ) a₀ + (Σxₙ²) a₁ = Σxₙyₙ\n\nNotación:  Sx = Σx,  Sy = Σy,  Sxy = Σxy,  Sxx = Σx²\n\na₀ = (Sxx·Sy − Sxy·Sx) / (N·Sxx − (Sx)²)\na₁ = (N·Sxy − Sx·Sy) / (N·Sxx − (Sx)²)"),
    UL("<b>2 incógnitas</b> (a₀, a₁) → <b>2 ecuaciones normales</b>. Regla: una ecuación por cada coeficiente.",
       "El sistema es <b>lineal</b> en a₀, a₁ aunque E sea cuadrática: por eso se resuelve directo, sin iterar."),
    AGR("¿Siempre existe la recta? El denominador N·Sxx − (Sx)² = N·Σ(xₙ − x̄)² es cero sólo si todos los xₙ son iguales (recta vertical, no es función de x). Con al menos dos abscisas distintas el sistema tiene solución única y es un mínimo (E es una parábola convexa en a₀, a₁). O sea: <b>calcularla</b> siempre se puede (salvo ese caso); que la recta sea un <b>buen modelo</b> es otra cosa (p11)."),
    TRAMPA("IT25-T1 1b-c: \"Escriba la ecuación que representa la condición de mínimos cuadrados en un ajuste lineal\" → ∂E/∂a₀ = 0 y ∂E/∂a₁ = 0 con E = Σ(yₙ − (a₀ + a₁xₙ))², y las dos derivadas desarrolladas. \"¿Siempre es posible encontrar un ajuste lineal usando MC?\" tiene dos lecturas: matemáticamente sí (salvo todos los x iguales); como modelo adecuado, no siempre (exponencial, periódico). Respondé las dos."),
)

sec(4, "Ejemplo de la cátedra y cálculo desde sumatorias", 2, "p9–10",
    P("Cómo se aplica la fórmula de la p8 con números, y cómo se resuelven los ejercicios de Moodle que dan las sumatorias ya hechas."),
    TBL(["x", "1", "2", "3", "4", "5", "6"], [["y", "0", "0.6", "1.77", "1.92", "3.31", "3.52"]]),
    OL("Graficamos los datos.",
       "\"Como los puntos parecen seguir una tendencia lineal, buscamos ajustar con una recta aplicando el método de regresión lineal por mínimos cuadrados.\"",
       f"N = {ej['N']}, Sx = {n(ej['Sx'],0)}, Sy = {n(ej['Sy'],2)}, Sx² = {n(ej['Sxx'],0)}, Sxy = {n(ej['Sxy'],2)}.",
       "a₀ = −0.734, a₁ = 0.739 → <b>F(x) = −0.734 + 0.739x</b> (p10)."),
    P(f"Verificado en Python: a₀ = {n(ej['a0'],5)}, a₁ = {n(ej['a1'],5)}; la diapositiva trunca a 3 decimales. Con esta recta, r = {n(ej['r'],4)} (fórmula de la p19)."),
    H4("Recta desde sumatorias y varianzas (formato Moodle, EV-P7)"),
    AGR("Los evaluativos dan Σxᵢ, Σyᵢ, Σxᵢyᵢ y las varianzas en lugar de los datos. Es la misma recta de MC escrita con medias: x̄ = Σx/n, ȳ = Σy/n, covarianza s_xy = Σxy/n − x̄·ȳ, pendiente a₁ = s_xy / s_x², ordenada a₀ = ȳ − a₁·x̄ (la recta pasa por (x̄, ȳ)). Convención: varianza y covarianza dividiendo por n, que es la que reproduce las opciones de EV-P7 (con n−1 las opciones no coinciden). El signo de s_xy da si la nube es creciente o decreciente."),
    EJ(f"EV-P7 (n = 8, Σx = 24, Σxy = 64, Σy = 40, s_x² = 6): x̄ = {n(R['evp7']['xm'],0)}, ȳ = {n(R['evp7']['ym'],0)}, s_xy = 64/8 − 3·5 = {n(R['evp7']['sxy'],0)}, a₁ = −7/6 = {n(R['evp7']['b'],4)}, a₀ = 5 + (7/6)·3 = {n(R['evp7']['a'],1)} → <b>Y = 8.5 − 1.167x</b>; en x = 4: {n(R['evp7']['p4'],3)}. Nube <b>decreciente</b>."),
    TRAMPA("EV-P7 pone como distractores la recta con los coeficientes intercambiados (Y = −1.167 + 8.5x), la predicción corrida en 1 (4.83 en lugar de 3.83) y \"nube creciente\" con pendiente negativa."),
)

sec(5, "¿Siempre es bueno ajustar con una recta?", 3, "p11",
    P("La recta de MC se puede calcular para cualquier nube, pero eso no quiere decir que sirva. La p11 muestra cuatro nubes: en una la recta va bien (✓), en una es dudosa y en dos está mal (✗)."),
    FIG(fig_p11()),
    TBL(["Nube (p11)", "Marca", "Por qué"], [
        ["Tendencia lineal con dispersión", "✓", "Los puntos se reparten a ambos lados de la recta sin patrón: la recta describe la tendencia."],
        ["Puntos alrededor de 0 sin tendencia", "duda (emoji pensativo)", f"La recta sale casi horizontal y r ≈ 0 (en esta recreación r = {n(F['ruido']['r'],2)}): no hay relación lineal que modelar."],
        ["Crecimiento casi nulo y después explosivo", "✗", "Es exponencial: la recta queda por encima en el medio y por debajo en los extremos. Se linealiza (p12)."],
        ["Oscilación", "✗", "Datos periódicos: ninguna recta ni polinomio bajo los describe. Ajuste trigonométrico (U3)."],
    ]),
    AGR(f"Ojo con r: en la nube exponencial recreada la recta da r = {n(F['expo']['r'],2)}, que \"parece bueno\", y sin embargo el modelo es incorrecto. Lo que delata el error es el <b>patrón de los residuos</b> (todos del mismo signo por tramos), no el valor de r. También indica problemas la <b>dispersión que crece con x</b> (heterocedasticidad): la recta pesa igual todos los puntos y los de mayor dispersión dominan E."),
    TRAMPA("IT25-T3 1 dio una nube que crece y se aplana con dispersión creciente, con la recta de MC encima: a) problemas, b) ¿r cercano a 1?, c) qué modelo. ET25-T2 ej. 2 pidió elegir entre lineal, cuadrática, exponencial, lineal con linealización y trigonométrico, explicando por qué sí y por qué no cada una."),
)

sec(6, "Linealización", 3, "p12–13",
    P("Si los datos siguen una curva conocida pero no lineal en los parámetros (exponencial, potencia, saturación), no se puede aplicar directo la recta de MC. Se transforman las variables para que la relación quede como recta y se usa la regresión lineal sobre los datos transformados."),
    D_("Linealización de la exponencial", "Parece tener la forma de: y = a₁·e^{b₁x}. ¿Cómo encontrar a₁ y b₁? ln y = ln(a₁·e^{b₁x}) → ln y = ln a₁ + b₁x, que tiene la forma F(x, aᵢ) = a₀ + a₁x (ln a₁: ordenada al origen). Se grafica ln y vs x.", "U2 p12"),
    FIG(fig_lineal()),
    H4("Los tres patrones de la p13"),
    TBL(["Si los datos siguen", "Se puede linealizar así", "Se visualiza", "Pendiente / ordenada"], [
        ["y = a₂·x^{b₂}", "ln y = ln a₂ + b₂ ln x", "ln y vs ln x", "pendiente b₂, ordenada ln a₂"],
        ["y = a₃·(x / (b₃ + x))", "b₃ + x = a₃x/y → (b₃ + x)/(a₃x) = 1/y → (b₃/a₃)(1/x) + 1/a₃ = 1/y", "1/y vs 1/x", "pendiente b₃/a₃, ordenada 1/a₃"],
        ["y = 1/(a₄x + b₄)", "1/y = a₄x + b₄", "1/y vs x", "pendiente a₄, ordenada b₄"],
        ["y = a₁·e^{b₁x} (p12)", "ln y = ln a₁ + b₁x", "ln y vs x", "pendiente b₁, ordenada ln a₁"],
    ]),
    P("Procedimiento: (1) transformar los datos (Yᵢ = ln yᵢ, 1/yᵢ, etc.); (2) recta de MC con la fórmula de la p8 sobre (Xᵢ, Yᵢ); (3) deshacer la transformación para recuperar los parámetros (a₁ = e^{ordenada}, a₃ = 1/ordenada, b₃ = pendiente·a₃)."),
    AGR("Dos cuidados que no están en la diapositiva: (1) ln y y 1/y exigen y > 0 (o y ≠ 0); con datos negativos o nulos no se puede linealizar así. (2) La recta minimiza el error en las variables transformadas (Σ(ln yᵢ − …)²), no Σ(yᵢ − a₁e^{b₁xᵢ})²: los coeficientes son buenos pero no son exactamente los de mínimos cuadrados sobre y. El r o R² que da la recta es el del espacio transformado."),
    TRAMPA("ET25-T2 ej. 2 d) pedía \"Regresión lineal con linealización (indicar cuál)\": para el gráfico exponencial es ln y = ln a₁ + b₁x, graficando ln y vs x. En ET25-T1 ej. 2 la corrección marginal pregunta \"¿cuánto vale R²?\" sobre una recta forzada a datos exponenciales: se espera R² bajo y proponer la linealización (inferido)."),
)

sec(7, "Regresión polinomial simple", 2, "p14–15",
    P("Cuando la tendencia es curva pero no corresponde a ninguno de los patrones linealizables, se ajusta un polinomio de grado mayor en una sola variable independiente."),
    D_("Regresión Polinomial Simple", "Se realiza el ajuste con un polinomio de una única variable independiente (simple). La regresión lineal es una regresión polinómica de grado 1. Del mismo modo se puede ajustar mediante regresión polinómica de grado >1.", "U2 p14"),
    CMD("Regresión polinomial de segundo orden (ajustamos con una parábola):\nF(x, aᵢ) = a₀ + a₁x + a₂x²\n\naᵢ se obtienen minimizando el error:\nE = Σₙ₌₁ᴺ ( yₙ − (a₀ + a₁xₙ + a₂xₙ²) )²\n\n∂E/∂a₀ = 0    ∂E/∂a₁ = 0    ∂E/∂a₂ = 0\n\nEcuaciones normales:\nN a₀    + (Σxₙ) a₁  + (Σxₙ²) a₂ = Σyₙ\n(Σxₙ) a₀  + (Σxₙ²) a₁ + (Σxₙ³) a₂ = Σxₙyₙ\n(Σxₙ²) a₀ + (Σxₙ³) a₁ + (Σxₙ⁴) a₂ = Σxₙ²yₙ"),
    UL("Grado 2: <b>3 incógnitas</b> (a₀, a₁, a₂) → <b>3 ecuaciones normales</b>. En general, grado m → m + 1 incógnitas y m + 1 ecuaciones.",
       "Sigue siendo un sistema <b>lineal</b> porque F es lineal en los aᵢ (no importa que sea cuadrática en x)."),
    D_("Observación", "Si se aumenta el grado del polinomio, comienza a aparecer inestabilidad (sistemas mal condicionados).", "U2 p15"),
    AGR("La diapositiva escribe el término cuadrático como \"a₂x_{n²}\" (subíndice); es a₂xₙ². Además de la inestabilidad numérica, subir el grado lleva a sobreajuste: el polinomio empieza a seguir el ruido y R² sube muy poco. En el TP2 ej. 4 (datos.csv) R² pasa de " + n(R['tp2e4']['tab'][0]['R2'],4) + " (grado 1) a " + n(R['tp2e4']['tab'][1]['R2'],4) + " (grado 2) y después casi no cambia hasta grado 6 (" + n(R['tp2e4']['tab'][5]['R2'],4) + "): se elige grado 2."),
    TRAMPA("\"¿Cuántas ecuaciones normales e incógnitas tiene una regresión polinomial de grado 3?\" → 4 y 4. Error típico: contestar N (el número de datos).", real=False),
)

sec(8, "Regresión polinomial múltiple", 2, "p16",
    P("Cuando y depende de más de una variable independiente (dureza del acero según cobre y temperatura), se ajusta un plano en lugar de una curva."),
    D_("Regresión Polinomial Múltiple", "Se realiza el ajuste mediante un plano (más de una variable independiente).", "U2 p16"),
    CMD("F(x₁, x₂, aᵢ) = a₀ + a₁x₁ + a₂x₂\n\nE = Σₙ₌₁ᴺ ( yₙ − (a₀ + a₁x₁ₙ + a₂x₂ₙ) )²\n\n∂E/∂a₀ = 0    ∂E/∂a₁ = 0    ∂E/∂a₂ = 0\n\n⎡ N      Σx₁ₙ      Σx₂ₙ    ⎤ ⎡a₀⎤   ⎡ Σyₙ    ⎤\n⎢ Σx₁ₙ   Σx₁ₙ²     Σx₁ₙx₂ₙ ⎥ ⎢a₁⎥ = ⎢ Σx₁ₙyₙ ⎥\n⎣ Σx₂ₙ   Σx₁ₙx₂ₙ   Σx₂ₙ²   ⎦ ⎣a₂⎦   ⎣ Σx₂ₙyₙ ⎦"),
    AGR("Erratas de la p16 (corregidas arriba): la diapositiva escribe E = Σ(yₙ − (a₀ + a₁x₁ₙ + <b>a₁</b>x₂ₙ))² (debe ser a₂x₂ₙ); en la matriz pone Σx₂ₙ en el lugar (2,2) (debe ser Σx₁ₙ²) y Σx^{2_{2n}} en el (3,3) (es Σx₂ₙ²)."),
    TBL(["", "Polinomial simple de 2.º orden (p14–15)", "Polinomial múltiple (p16)"], [
        ["Variables independientes", "Una: x", "Dos o más: x₁, x₂"],
        ["Modelo", "a₀ + a₁x + a₂x² (parábola)", "a₀ + a₁x₁ + a₂x₂ (plano)"],
        ["Incógnitas / ecuaciones normales", "3 / 3", "3 / 3 (con dos variables)"],
        ["Error a minimizar", "Σ(yₙ − (a₀ + a₁xₙ + a₂xₙ²))²", "Σ(yₙ − (a₀ + a₁x₁ₙ + a₂x₂ₙ))²"],
        ["Elementos del sistema", "Σxₙ hasta Σxₙ⁴", "Σx₁ₙ, Σx₂ₙ, Σx₁ₙ², Σx₂ₙ², Σx₁ₙx₂ₙ"],
    ]),
    P("Las dos son mínimos cuadrados lineales (F lineal en los aᵢ): se minimiza E y se resuelve un sistema lineal. Coinciden en número de incógnitas sólo por casualidad (grado 2 vs dos variables)."),
    EJ(f"EV-P1: y = 161.3365 + 32.9688u − 0.0855v (u = % de cobre, v = temperatura). Con u = 0.05 y v = 1125: 161.3365 + {n(R['evp1']['t_u'],5)} − {n(R['evp1']['t_v'],4)} = <b>{n(R['evp1']['dur'],3)}</b>. Verificado: ajustando el plano a las 12 durezas de la tabla se reobtienen los coeficientes ({n(R['evp1']['coef_tabla'][0],4)}, {n(R['evp1']['coef_tabla'][1],4)}, {n(R['evp1']['coef_tabla'][2],4)})."),
    TRAMPA("IT25-T2 1 pidió la diferencia entre polinomial de 2.º orden y polinomial múltiple, cuántas ecuaciones normales e incógnitas en cada caso, y la ecuación del error en cada caso. La tabla de arriba es la respuesta."),
)

sec(9, "Dados los datos, ¿qué método usamos?", 3, "p17",
    P("Con todas las herramientas en la mano (recta, linealización, polinomio simple, múltiple), la decisión sale de mirar los datos. La p17 lo plantea como escalera: regresión lineal → linealización → regresión polinomial simple → regresión polinomial múltiple."),
    TBL(["Lo que ves en el gráfico", "Modelo", "Incógnitas"], [
        ["Tendencia recta, dispersión pareja", "Regresión lineal F = a₀ + a₁x", "a₀, a₁"],
        ["Crecimiento casi nulo y después explosivo (o decaimiento rápido a 0)", "Exponencial y = a₁e^{b₁x}, linealizada: ln y vs x", "a₁, b₁"],
        ["Recta en papel log-log", "Potencia y = a₂x^{b₂}: ln y vs ln x", "a₂, b₂"],
        ["Crece y se satura hacia un valor máximo (y > 0)", "y = a₃x/(b₃ + x): 1/y vs 1/x", "a₃, b₃"],
        ["Parábola (un máximo o un mínimo)", "Polinomial de 2.º orden", "a₀, a₁, a₂"],
        ["y depende de dos variables", "Polinomial múltiple (plano)", "a₀, a₁, a₂"],
        ["Oscilación periódica", "Ajuste trigonométrico (U3)", "a₀, aₖ, bₖ"],
    ]),
    P("Criterio en orden: (1) graficar; (2) mirar forma, outliers y si la dispersión cambia; (3) proponer el modelo más simple que siga la forma; (4) evaluar con residuos, E_RMS y r (p18–19); (5) si un modelo más complejo apenas mejora r, quedarse con el simple."),
    AGR("Un modelo logarítmico y = a₀ + a₁ ln x (crece rápido y se aplana, admite valores negativos cerca de x = 0) no está en las diapositivas, pero se ajusta con la misma recta de MC usando ln x como variable. Es la opción natural para nubes como la de IT25-T3."),
    TRAMPA("ET25-T2 ej. 2: para la nube exponencial, la respuesta es exponencial (c) implementada como lineal con linealización ln y (d); para la parábola decreciente, cuadrática (b). Pedía además \"por qué podría ser elegido y por qué no\" cada caso: hay que justificar también los descartados."),
)

sec(10, "Evaluación de la calidad del ajuste: residuo, error y error cuadrático medio", 3, "p18",
    P("Después de ajustar hay que medir qué tan bien quedó. Lo primero es el residuo de cada punto; con los residuos se arma el error total y el error cuadrático medio."),
    FIG(fig_residuos_sec()),
    CMD("Residuo (regresión lineal simple):   Rᵢ = yᵢ − (a₀ + a₁x),   i = 1, …, n\n\nError:                     E = Σ ( yᵢ − F(x, aᵢ) )²\n\nError cuadrático medio:    E_RMS = ( Σ ( yᵢ − F(x, aᵢ) )² / (n + 1) )^(1/2)"),
    D_("Error cuadrático medio", "E_RMS = (Σ(yᵢ − F(x,aᵢ))² / (n+1))^{1/2}", "U2 p18"),
    UL("E depende de cuántos datos hay (sumar más puntos lo hace crecer). E_RMS lo divide por la cantidad y saca la raíz: queda en las <b>mismas unidades que y</b> y se lee como el tamaño típico del error.",
       "Sirve para comparar modelos sobre los mismos datos: menor E_RMS, mejor ajuste.",
       "Los residuos se miran también como gráfico: deben repartirse al azar alrededor de 0. Un patrón (curva en U, signos por tramos, abanico que se abre) indica que el modelo no capta la forma o que la dispersión no es pareja."),
    AGR(f"La cátedra divide por n + 1. Lo usual es dividir por n (RMSE) o por n − p (p = número de parámetros, estimador insesgado). Usá la de la cátedra y aclaralo. En el ejemplo de la p9: E = {n(ej['E'],4)}, E_RMS(cátedra) = √(E/7) = {n(ej['Erms_cat'],4)}; con n: {n(ej['Erms_n'],4)}."),
    TRAMPA("IT25-T4 1a: diferencia conceptual entre ECM y r. El ECM mide el tamaño absoluto del error, en unidades de y; r mide qué fracción de la dispersión de los datos explica el ajuste, sin unidades. Un ajuste puede tener r alto y ECM grande (datos de valores grandes) o al revés.", real=True),
)

sec(11, "Coeficiente de correlación", 3, "p19",
    P("E_RMS dice cuánto se equivoca el ajuste, pero no si eso es mucho o poco respecto de cuánto varían los datos. El coeficiente de correlación compara las dos dispersiones."),
    D_("Coeficiente de correlación", "Una manera de medir la relación entre las variables es el coeficiente de correlación. Este valor que varía entre 0 y 1 nos indica la dispersión de los datos respecto al ajuste.", "U2 p19"),
    CMD("S_T = Σ ( yᵢ − ȳ )²            dispersión de los datos experimentales\nS_R = Σ ( yᵢ − F(x, aᵢ) )²      dispersión de los datos respecto a F\n\nr = √( (S_T − S_r) / S_T )       Coef de correlación"),
    P("\"r cercano a 1 indica que existe una buena correlación entre el polinomio de ajuste, con los datos experimentales\" (p19). El código de la diapositiva:"),
    CMD("def St(yreal):\n    promedio=np.mean(yreal)\n    sum=np.sum((yreal-promedio)**2)\n    return(sum)\n\ndef Sr(yreal,yaprox):\n    sum=np.sum((yreal-yaprox)**2)\n    return (sum)\n\ndef r(Sr,St):\n    r_cuad=1-(Sr/St)\n    r=np.sqrt(r_cuad)\n    return(r)"),
    UL("S_T es el error de usar la media ȳ como modelo; S_R, el error del ajuste. (S_T − S_R)/S_T = r² es la fracción de la dispersión que el ajuste explica (<b>coeficiente de determinación R²</b>).",
       "Si F pasa por todos los puntos, S_R = 0 y r = 1. Si el ajuste no mejora a la media, S_R = S_T y r = 0.",
       "La definición usa sólo yᵢ y F(xᵢ): vale para <b>cualquier</b> modelo (recta, parábola, exponencial), no sólo lineal. La p19 dice \"polinomio de ajuste\" y en la p23 se compara r del lineal y del cuadrático."),
    AGR(f"El coeficiente de Pearson de estadística tiene signo (−1 ≤ r ≤ 1): negativo si la nube es decreciente (EV-P7 da r = {n(R['evp7']['r'],3)}). Para la recta de MC, |Pearson| coincide con el r de la p19 (ejemplo p9: ambos {n(ej['r'],4)}). La definición de la cátedra, con la raíz, siempre es ≥ 0. Si el modelo se ajustó linealizado, r se calcula en el espacio transformado salvo que se recalculen S_T y S_R con y original."),
    H4("Otras métricas (TP2 ej. 1)"),
    AGR(f"No están en las diapositivas: <b>R²</b> = 1 − S_R/S_T; <b>R² ajustado</b> = 1 − (1 − R²)(n − 1)/(n − p − 1), penaliza agregar parámetros; <b>MAE</b> = Σ|yᵢ − F(xᵢ)|/n, menos sensible a outliers que el ECM; <b>MSE</b> y <b>RMSE</b> (= E_RMS con n); <b>MAPE</b> = (100/n)Σ|(yᵢ − F(xᵢ))/yᵢ|, en %; análisis gráfico de residuos; AIC/BIC para comparar modelos con distinto número de parámetros. En el ejemplo de la p9: R² = {n(ej['R2'],4)}, R² ajustado = {n(ej['R2adj'],4)}, MAE = {n(ej['MAE'],4)}."),
    TRAMPA("IT25-T4 1b: \"Tanto el ECM como r, ¿sólo pueden ser usados en ajustes lineales?\" → No: los dos se calculan con yᵢ y F(xᵢ, aᵢ) de cualquier modelo. 1c: escribir r = √((S_T − S_r)/S_T) con S_T y S_R definidos."),
)

sec(12, "Ejemplo actual: casos de COVID-19 en Tucumán", 2, "p20–23",
    P("La cátedra cierra con datos reales para mostrar que un r \"aceptable\" no alcanza: hay que mirar más datos y comparar modelos antes de sacar conclusiones."),
    OL("<b>Datos</b>: primera semana de septiembre de 2020, capital (fuente lagaceta.com.ar): x = 1…7, y = 159, 156, 187, 150, 219, 198, 218.",
       "<b>Visualización</b>: ¿qué comportamiento esperábamos? ¿qué tendencia se observa? ¿qué modelo sería el más adecuado?",
       f"<b>Ajuste lineal</b>: r = 0.77. ¿Es un buen ajuste? ¿Es satisfactorio? ¿Qué podría pasar en el futuro? Verificado: F(x) = {n(cv['a0'],1)} + {n(cv['a1'],3)}x, r = {n(cv['r'],4)}.",
       "<b>Más datos</b> (1 al 19 de septiembre): r = 0.86. \"Mejoró mucho más el ajuste\"; \"si miramos estas semanas los casos crecen de manera lineal\". Obs: la serie inicia en x = 1 (primer día de medición, no la fecha). \"¿Pero qué pasa si miramos la serie con más datos hacia atrás?\"",
       "<b>Serie 1/jul–19/sep</b>: lineal r = 0.83; exponencial \"R = ?\" (la curva explota, escala 1e59); cuadrático r = 0.96.",
       "<b>Sacar conclusiones</b>."),
    P("Lecturas que se esperan: (1) con 7 datos ruidosos, r = 0.77 es pobre y la recta no permite extrapolar; (2) un r que sube con más datos no prueba que el fenómeno sea lineal, sólo describe esa ventana; (3) mirando la serie completa, la forma es curva (casi nula y después creciente) y el cuadrático (r = 0.96) describe mejor que la recta (0.83); (4) el exponencial ajustado mal (sin linealizar bien o con parámetros que explotan) puede ser peor que todos: elegir modelo exige graficar el ajuste, no sólo mirar un número."),
    AGR("La diapositiva usa R y r indistintamente para el coeficiente de correlación. Los datos de las p22–23 no están en el material: sólo se verificó r = 0.77 de la p21."),
)

intro = P("Ajustar es encontrar la función que describe mejor una nube de datos con error, sin obligarla a pasar por cada punto. La unidad arma todo sobre <b>mínimos cuadrados</b>: definir el error E = Σ(yₙ − F(xₙ))², derivar respecto de cada coeficiente, igualar a cero y resolver las ecuaciones normales; después extiende la idea a modelos linealizables, polinomios y planos, y la cierra con cómo medir la calidad (residuos, E_RMS, r).") + \
        P("En el integral teórico siempre hay una pregunta de ajuste (IT25 temas 1–4): condición de MC y su ecuación, ecuaciones normales e incógnitas, ECM vs r, y leer un gráfico de nube + recta para decir qué anda mal y qué modelo usar. En la evaluación práctica y Moodle: elegir el ajuste para un gráfico, R² de una recta forzada sobre datos exponenciales, predicciones con un ajuste dado y la recta desde sumatorias.")

conexiones = UL("Con la unidad 3: el ajuste trigonométrico (serie de Fourier discreta) es otro ajuste por mínimos cuadrados, con senos y cosenos como funciones base; es la respuesta para datos periódicos que la recta no puede seguir (p11).",
                "Con la unidad 3, otra vez: la diferencia \"aproximación polinómica vs trigonométrica\" (WA25) es la elección de base para el mismo criterio de MC.",
                "Con interpolación (Métodos Numéricos I): interpolar fuerza error cero en los nodos; ajustar minimiza el error global. La p4 recomienda ajustar antes que interpolar con polinomios de orden alto.",
                "Con álgebra lineal: las ecuaciones normales son un sistema lineal; la inestabilidad de los polinomios de grado alto (p15) es mal condicionamiento de esa matriz.")

fidelidad = {
    "fuentes": "ajuste.pdf (23 diapositivas; las 23 renderizadas y leídas, p6, p15, p16, p18 y p19 también a 120 dpi). TP2 (tp2.pdf) y datos.csv. Patrón 02-patron-evaluacion.md (IT25-T1..T4, ET25-T1/T2, EV-P1, P4, P7, F18-P4, WA25). Números: tools/gen/verif_u2.py.",
    "agregado": "Existencia de la recta (N·Sxx − Sx² ≠ 0); recta desde medias y varianzas (convención n) para EV-P7; cuidados de la linealización (y > 0, minimiza en el espacio transformado); modelo logarítmico; heterocedasticidad y patrón de residuos; E_RMS con n vs n+1; Pearson con signo; R², R² ajustado, MAE, MAPE, AIC/BIC (TP2 ej. 1); mínimos cuadrados no lineales (curve_fit); sobreajuste; fenómeno de Runge; erratas de p14–16 corregidas.",
    "huecos": "Ninguna diapositiva ilegible. Los datos de la serie COVID de las p22–23 no están en el material (r = 0.86, 0.83 y 0.96 no se pueden verificar). La convención de varianza (n o n−1) no está explicitada por la cátedra; se usa n porque reproduce las opciones de EV-P7.",
}

auditoria = [
    {"tema": "Objetivo: encontrar la función que describe mejor los datos (p2)", "seccion": "m2-s1", "estado": "cubierto"},
    {"tema": "Aproximación de funciones (p2–4)", "seccion": "m2-s1", "estado": "cubierto"},
    {"tema": "Regresión Lineal (p5–6)", "seccion": "m2-s2", "estado": "cubierto"},
    {"tema": "Regresión Lineal por mínimos cuadrados (p7–8)", "seccion": "m2-s3", "estado": "cubierto"},
    {"tema": "Regresión Lineal por mínimos cuadrados — EJEMPLO (p9–10)", "seccion": "m2-s4", "estado": "cubierto"},
    {"tema": "¿Siempre es bueno ajustar con una recta? (p11)", "seccion": "m2-s5", "estado": "cubierto"},
    {"tema": "Linealización (p12–13)", "seccion": "m2-s6", "estado": "cubierto"},
    {"tema": "Regresión Polinomial Simple (p14–15)", "seccion": "m2-s7", "estado": "cubierto"},
    {"tema": "Regresión Polinomial Múltiple (p16)", "seccion": "m2-s8", "estado": "cubierto"},
    {"tema": "Dados los datos, ¿qué método usamos? (p17)", "seccion": "m2-s9", "estado": "cubierto"},
    {"tema": "Evaluación de la calidad del ajuste: residuo, E, E_RMS (p18)", "seccion": "m2-s10", "estado": "cubierto"},
    {"tema": "Evaluación de la calidad del ajuste: coeficiente de correlación (p19)", "seccion": "m2-s11", "estado": "cubierto"},
    {"tema": "Veamos un ejemplo actual: COVID-19 Tucumán (p20–23)", "seccion": "m2-s12", "estado": "cubierto"},
]

glosario = [
    {"termino": "Interpolación", "def": "Aproximación que pasa exactamente por todos los puntos dados (error cero en los nodos); con datos ruidosos no describe la tendencia.", "seccion": "m2-s1"},
    {"termino": "Ajuste de curvas", "def": "Aproximar los datos con otra función que describe su comportamiento con cierto error, sin pasar necesariamente por los puntos.", "seccion": "m2-s1"},
    {"termino": "Ajuste lineal", "def": "Cualquier recta F = a₀ + a₁x propuesta para los datos; sin criterio hay infinitas.", "seccion": "m2-s2"},
    {"termino": "Mínimos cuadrados", "def": "Criterio que elige los coeficientes que hacen mínimo E = Σ(yₙ − F(xₙ, aᵢ))².", "seccion": "m2-s3"},
    {"termino": "Condición de mínimo", "def": "∂E/∂aᵢ = 0 para cada coeficiente; en la recta, ∂E/∂a₀ = 0 y ∂E/∂a₁ = 0.", "seccion": "m2-s3"},
    {"termino": "Ecuaciones normales", "def": "Sistema lineal que sale de la condición de mínimo; una ecuación por incógnita (recta: 2; parábola o plano con dos variables: 3).", "seccion": "m2-s3"},
    {"termino": "Residuo", "def": "Rᵢ = yᵢ − F(xᵢ): diferencia vertical entre el dato y el ajuste.", "seccion": "m2-s10"},
    {"termino": "Error E", "def": "Suma de los residuos al cuadrado, Σ(yᵢ − F(x, aᵢ))²; crece con la cantidad de datos.", "seccion": "m2-s10"},
    {"termino": "Error cuadrático medio (E_RMS)", "def": "(Σ(yᵢ − F)²/(n+1))^{1/2} en la cátedra; error típico en unidades de y.", "seccion": "m2-s10"},
    {"termino": "S_T", "def": "Σ(yᵢ − ȳ)²: dispersión de los datos experimentales respecto de su media.", "seccion": "m2-s11"},
    {"termino": "S_R", "def": "Σ(yᵢ − F(x, aᵢ))²: dispersión de los datos respecto al ajuste.", "seccion": "m2-s11"},
    {"termino": "Coeficiente de correlación r", "def": "r = √((S_T − S_r)/S_T), entre 0 y 1; cercano a 1 indica buena correlación entre el ajuste y los datos.", "seccion": "m2-s11"},
    {"termino": "Coeficiente de determinación R²", "def": "r² = 1 − S_R/S_T: fracción de la dispersión de y que explica el ajuste (agregado).", "seccion": "m2-s11"},
    {"termino": "Linealización", "def": "Transformar variables (ln y, ln x, 1/y, 1/x) para que un modelo no lineal en los parámetros quede como recta y usar la regresión lineal.", "seccion": "m2-s6"},
    {"termino": "Regresión polinomial simple", "def": "Ajuste con un polinomio de una única variable independiente; la regresión lineal es la de grado 1.", "seccion": "m2-s7"},
    {"termino": "Regresión polinomial múltiple", "def": "Ajuste mediante un plano: más de una variable independiente, F = a₀ + a₁x₁ + a₂x₂.", "seccion": "m2-s8"},
    {"termino": "Sistema mal condicionado", "def": "Sistema cuya solución cambia mucho ante pequeños cambios en los datos; aparece al subir el grado del polinomio.", "seccion": "m2-s7"},
    {"termino": "Heterocedasticidad", "def": "Dispersión de los datos que cambia a lo largo de x (abanico); la recta de MC pesa igual todos los puntos (agregado).", "seccion": "m2-s5"},
    {"termino": "Covarianza s_xy", "def": "Σxy/n − x̄ȳ (convención n); su signo indica nube creciente o decreciente; a₁ = s_xy/s_x².", "seccion": "m2-s4"},
]

# ------------------------------------------------------------------ quiz
Q = []
def q(sec_n, tema, pide, sobre, ops, expl, fig=None):
    d = {"id": f"m2-q{len(Q)+1:02d}", "seccion": f"m2-s{sec_n}", "tema": tema, "pide": pide, "sobre": sobre,
         "opciones": [{"t": t, "v": v} for t, v in ops], "explicacion": expl}
    if fig: d["fig"] = fig
    Q.append(d)

q(1, "Interpolación vs ajuste", "verdaderas", "la interpolación y el ajuste de curvas",
  [("La interpolación garantiza que la aproximación coincide con las ordenadas en las abscisas de los datos.", True),
   ("Con datos experimentales ruidosos, un ajuste suele describir mejor la tendencia que una interpolación.", True),
   ("El ajuste de curvas exige que la función pase exactamente por todos los puntos medidos.", False),
   ("Interpolar con un polinomio de orden alto es siempre preferible porque deja error nulo.", False)],
  "La p3 define la interpolación como la aproximación que coincide con las ordenadas en las abscisas de los datos: verdadera. La p4 dice que muchas veces es más útil ajustar una recta que describa el comportamiento con cierto error que interpolar con un polinomio de orden alto: por eso el ajuste describe mejor la tendencia con datos ruidosos. Exigir que el ajuste pase por todos los puntos es confundirlo con la interpolación. Y el error nulo en los nodos no hace preferible al polinomio de orden alto: copia el ruido y oscila.")

q(2, "Qué minimiza MC", "una", "En la regresión lineal por mínimos cuadrados de la cátedra, ¿qué cantidad se hace mínima?",
  [("La suma de los cuadrados de las diferencias verticales yₙ − (a₀ + a₁xₙ).", True),
   ("La suma de las distancias perpendiculares de cada punto a la recta.", False),
   ("La suma de los valores absolutos de yₙ − (a₀ + a₁xₙ), sin elevar.", False),
   ("La cantidad de puntos que quedan por fuera de la recta ajustada.", False)],
  "La p6 define el error como la distancia de los datos a la recta, E = Σ(yₙ − (a₀ + a₁xₙ))²: diferencias verticales (dato menos recta en el mismo x) al cuadrado. Las distancias perpendiculares corresponden a otro método (mínimos cuadrados totales), no al de la cátedra. Los valores absolutos dan otro criterio (mínima desviación absoluta) que no se resuelve con ecuaciones normales lineales. Contar puntos fuera de la recta no es un error diferenciable ni mide distancias.")

q(3, "Condición de mínimo", "verdaderas", "la condición de mínimos cuadrados para la recta F = a₀ + a₁x",
  [("Se plantea ∂E/∂a₀ = 0 y ∂E/∂a₁ = 0 con E = Σ(yₙ − (a₀ + a₁xₙ))².", True),
   ("∂E/∂a₁ = −2Σ(yₙ − a₀ − a₁xₙ)xₙ = 0 es una de las dos ecuaciones.", True),
   ("Las ecuaciones resultantes son no lineales en a₀ y a₁ y se resuelven iterando.", False),
   ("La condición se plantea derivando E respecto de cada dato xₙ e igualando a cero.", False)],
  "La p7 escribe la condición de mínimo como ∂E/∂a₀ = 0 y ∂E/∂a₁ = 0, y desarrolla ∂E/∂a₁ = −2Σ(yₙ − a₀ − a₁xₙ)xₙ = 0: las dos son verdaderas. El sistema que resulta (p8) es lineal en a₀ y a₁, por eso tiene solución cerrada; decir que es no lineal y hay que iterar es falso. Se deriva respecto de los coeficientes incógnita, no respecto de los datos xₙ, que son números conocidos.")

q(3, "Ecuaciones normales de la recta", "una", "¿Cuántas ecuaciones normales e incógnitas tiene la regresión lineal por mínimos cuadrados con N = 20 datos?",
  [("2 ecuaciones normales y 2 incógnitas (a₀ y a₁).", True),
   ("20 ecuaciones normales, una por cada dato, y 2 incógnitas.", False),
   ("1 ecuación normal (∂E/∂a₁ = 0) y 1 incógnita, la pendiente.", False),
   ("3 ecuaciones normales: dos derivadas y la de la ordenada al origen.", False)],
  "Hay una ecuación normal por cada coeficiente a determinar: la recta tiene a₀ y a₁, así que son 2 ecuaciones (∂E/∂a₀ = 0 y ∂E/∂a₁ = 0) y 2 incógnitas, sin importar cuántos datos haya. Los 20 datos entran dentro de las sumatorias, no como ecuaciones aparte. Una sola ecuación sólo alcanzaría si la ordenada estuviera fija. La ordenada al origen es a₀, que ya tiene su derivada: no agrega una tercera.")

q(3, "Existencia de la recta", "una", "¿En qué caso NO se puede calcular la recta de mínimos cuadrados con las fórmulas de la p8?",
  [("Cuando todos los datos tienen la misma abscisa x, porque N·Sxx − (Sx)² = 0.", True),
   ("Cuando los datos no siguen una tendencia lineal, porque el sistema no tiene solución.", False),
   ("Cuando hay valores de y negativos, porque las fórmulas usan logaritmos.", False),
   ("Cuando el coeficiente de correlación da menor que 0.5, porque el sistema diverge.", False)],
  "El denominador de a₀ y a₁ es N·Sxx − (Sx)² = N·Σ(xₙ − x̄)², que vale cero sólo si todas las abscisas son iguales: es el único caso en que el sistema no tiene solución única. Que los datos no sean lineales no impide calcular la recta, sólo la hace un mal modelo (p11). Los logaritmos aparecen en la linealización, no en la recta. r se calcula después del ajuste y no afecta la resolución del sistema.")

q(2, "MC y cualquier función", "verdaderas", "el método de mínimos cuadrados",
  [("El criterio de minimizar Σ(yₙ − F(xₙ, aᵢ))² se puede plantear para cualquier modelo F.", True),
   ("Con F lineal en los parámetros, la condición de mínimo da un sistema lineal de ecuaciones normales.", True),
   ("Para y = a₁e^{b₁x} las ecuaciones normales de la p8 se aplican directamente sobre y.", False),
   ("Aplicar MC garantiza que el modelo elegido describe bien cualquier conjunto de datos.", False)],
  "La forma general del error de la p6 vale para cualquier F: el criterio siempre se puede plantear. Cuando F es lineal en los aᵢ (recta, polinomio, plano) la condición ∂E/∂aᵢ = 0 da un sistema lineal. Para la exponencial no: hay que linealizar (ln y = ln a₁ + b₁x, p12) y recién ahí aplicar la recta. MC encuentra el mejor juego de coeficientes para el modelo elegido, pero si el modelo no tiene la forma de los datos (p11) el ajuste es malo igual: es el absoluto falso de WA25.")

v7 = R["evp7q"]
q(4, "Recta desde sumatorias", "una", f"De n = {v7['n']} observaciones: Σxᵢ = {v7['Sx']}, Σyᵢ = {v7['Sy']}, Σxᵢyᵢ = {v7['Sxy']}, s_x² = {v7['sx2']} (varianzas y covarianza dividiendo por n). La recta de regresión de Y sobre X es:",
  [(f"Y = {n(v7['a'],1)} − {n(-v7['b'],1)}x", True),
   (f"Y = −{n(-v7['b'],1)} + {n(v7['a'],1)}x", False),
   (f"Y = {n(v7['a'],1)} + {n(-v7['b'],1)}x", False),
   (f"Y = {n(v7['ym'],1)} − {n(-v7['b'],1)}x", False)],
  f"x̄ = {v7['Sx']}/{v7['n']} = {n(v7['xm'],0)}, ȳ = {n(v7['ym'],0)}, s_xy = {v7['Sxy']}/{v7['n']} − {n(v7['xm'],0)}·{n(v7['ym'],0)} = {n(v7['sxy'],1)}. Pendiente a₁ = s_xy/s_x² = {n(v7['b'],1)}; ordenada a₀ = ȳ − a₁x̄ = {n(v7['a'],1)}. La recta con los números intercambiados pone la pendiente como ordenada (distractor de EV-P7). La de pendiente positiva ignora el signo de la covarianza. La que usa ȳ como ordenada olvida corregir por a₁x̄: la recta pasa por (x̄, ȳ), no por (0, ȳ).")

q(4, "Predicción con el ejemplo", "una", "Con el ajuste del ejemplo de la cátedra, F(x) = −0.734 + 0.739x, ¿cuál es la predicción para x = 8?",
  [(f"{n(R['q_ej8']['p'],3)}", True), (f"{n(R['q_ej8']['swap'],3)}", False), ("5.912", False), ("0.005", False)],
  f"Se reemplaza en la recta: −0.734 + 0.739·8 = {n(R['q_ej8']['p'],3)}. El valor {n(R['q_ej8']['swap'],3)} sale de intercambiar ordenada y pendiente (0.739 − 0.734·8), el mismo mecanismo que EV-P7. 5.912 es sólo 0.739·8, olvidando la ordenada al origen. 0.005 es la suma −0.734 + 0.739, es decir evaluar en x = 1 en lugar de x = 8.")

e7 = R["evp7"]
q(4, "EV-P7", "verdaderas", "una muestra de 8 observaciones con Σxᵢ = 24, Σxᵢyᵢ = 64, Σyᵢ = 40, s_y² = 12, s_x² = 6 (convención: dividir por n)",
  [("La recta de regresión es Y = 8.5 − 1.167x.", True),
   ("La predicción para x = 4 es 3.83.", True),
   ("La predicción para x = 4 es 4.83.", False),
   ("Las observaciones forman una nube de puntos creciente.", False)],
  f"x̄ = 3, ȳ = 5, s_xy = 64/8 − 3·5 = {n(e7['sxy'],0)}, a₁ = −7/6 = {n(e7['b'],3)}, a₀ = 5 + (7/6)·3 = {n(e7['a'],1)}: la recta Y = 8.5 − 1.167x es correcta. En x = 4 da 8.5 − 4.667 = {n(e7['p4'],2)}, así que 4.83 es la predicción corrida en 1 (distractor del Moodle). La covarianza negativa (r = {n(e7['r'],3)}) indica nube decreciente, no creciente. Con n−1 en la covarianza la recta sería Y = 9 − 1.333x y no aparece en las opciones: por eso la convención es n.")

q(5, "Nube sin tendencia", "una", "Se ajustó una recta por MC a estos datos. ¿Qué conclusión es correcta?",
  [("No hay relación lineal apreciable: la recta sale casi horizontal y r ≈ 0.", True),
   ("La recta es un buen modelo porque pasa justo por el centro de la nube de puntos medidos.", False),
   ("Hay que interpolar con un polinomio de alto grado para captar todos los puntos.", False),
   ("El ajuste falló porque los datos negativos no admiten mínimos cuadrados.", False)],
  f"Es la nube dudosa de la p11: puntos repartidos alrededor de 0 sin tendencia. La recta de MC sale casi horizontal y r = {n(F['ruido']['r'],2)}: el ajuste no explica casi nada de la dispersión. Pasar por el centro de la nube no la hace un buen modelo, cualquier recta de MC pasa por (x̄, ȳ). Interpolar con grado alto copiaría el ruido (p4). Los valores negativos no son problema para la recta: sólo lo son al linealizar con logaritmos.", fig=fig_ruido_q())

q(5, "Recta sobre exponencial", "una", "La recta de MC sobre estos datos da el r indicado. ¿Qué afirmación es correcta?",
  [("La recta es inadecuada: los residuos tienen patrón y la forma es exponencial.", True),
   ("Como r está cerca de 1, la recta es el modelo correcto para estos datos.", False),
   ("La recta es adecuada pero conviene agregar más datos para que r llegue a 1.", False),
   ("El problema es la dispersión de los datos; con un filtro la recta sería válida.", False)],
  f"La nube crece casi nada y después explota: es el caso exponencial que la p11 marca con ✗. Aunque r = {n(F['expo']['r'],2)} parezca aceptable, la recta queda por encima de los datos en el medio y por debajo en los extremos (residuos con patrón) y predice valores negativos al principio. Un r relativamente alto no prueba que el modelo sea el correcto. Más datos con la misma forma no arreglan una forma equivocada. No es un problema de ruido sino de modelo: hay que linealizar.", fig=fig_expo_q())

q(6, "Linealización exponencial", "una", "Para estos datos, ¿qué transformación permite usar la regresión lineal?",
  [("Graficar ln y vs x: queda ln y = ln a₁ + b₁x.", True),
   ("Graficar ln y vs ln x: queda ln y = ln a₂ + b₂ ln x.", False),
   ("Graficar 1/y vs 1/x: queda 1/y = (b₃/a₃)(1/x) + 1/a₃.", False),
   ("Graficar y vs x²: queda y = a₀ + a₂x², sin término lineal.", False)],
  "Los datos crecen de forma exponencial, y = a₁e^{b₁x}. Tomando logaritmo, ln y = ln a₁ + b₁x (p12): ln y vs x es una recta con ordenada ln a₁ y pendiente b₁. La transformación log-log corresponde a la potencia y = a₂x^{b₂}; 1/y vs 1/x, a la saturación y = a₃x/(b₃ + x); ambas son de la p13 pero para otras formas. Graficar contra x² es un caso particular de parábola, que no sigue el crecimiento explosivo.", fig=fig_expo_q())

q(6, "Patrón potencia", "una", "Los datos (izquierda) quedan sobre una recta al graficar ln y vs ln x (derecha). ¿Qué modelo corresponde?",
  [("y = a₂·x^{b₂}, con pendiente b₂ y ordenada ln a₂.", True),
   ("y = a₁·e^{b₁x}, con pendiente b₁ y ordenada ln a₁.", False),
   ("y = 1/(a₄x + b₄), con pendiente a₄ y ordenada b₄.", False),
   ("y = a₃·x/(b₃ + x), con pendiente b₃/a₃ y ordenada 1/a₃.", False)],
  "La p13 asocia cada visualización a un patrón: ln y vs ln x recta ⇔ y = a₂x^{b₂}, porque ln y = ln a₂ + b₂ ln x. La exponencial se ve recta en ln y vs x (no en log-log). La recíproca 1/(a₄x + b₄) se ve recta en 1/y vs x. La saturación a₃x/(b₃ + x) se ve recta en 1/y vs 1/x. Los cuatro tienen pendiente y ordenada con nombre propio: el distractor cambia sólo la escala del eje.", fig=fig_loglog())

q(6, "Patrón recíproco", "una", "Los datos (izquierda) quedan sobre una recta al graficar 1/y vs x (derecha). ¿Qué modelo se ajusta?",
  [("y = 1/(a₄x + b₄), con 1/y = a₄x + b₄.", True),
   ("y = a₃·x/(b₃ + x), con 1/y = (b₃/a₃)(1/x) + 1/a₃.", False),
   ("y = a₁·e^{b₁x} con b₁ negativo, con ln y = ln a₁ + b₁x.", False),
   ("y = a₂·x^{b₂} con b₂ negativo, con ln y = ln a₂ + b₂ ln x.", False)],
  "La p13: si los datos siguen y = 1/(a₄x + b₄), entonces 1/y = a₄x + b₄ y se visualiza como recta en 1/y vs x. La saturación también usa 1/y pero contra 1/x. Una exponencial decreciente o una potencia con exponente negativo también bajan como estos datos, pero se ven rectas en ln y vs x y en ln y vs ln x respectivamente, no en 1/y vs x. Lo que decide es en qué gráfico transformado aparece la recta.", fig=fig_recip())

qe = R["q_exp"]
q(6, "Recuperar parámetros", "una", "Al ajustar ln y vs x por MC se obtuvo ln y = 0.7 + 0.3x. ¿Cuál es el modelo exponencial y = a₁e^{b₁x}?",
  [(f"a₁ = {n(qe['a'],3)}, b₁ = 0.3", True), ("a₁ = 0.7, b₁ = 0.3", False), (f"a₁ = {n(qe['a'],3)}, b₁ = {n(np.exp(0.3),3)}", False), ("a₁ = 0.3, b₁ = 0.7", False)],
  f"La ordenada de la recta es ln a₁, así que a₁ = e^0.7 = {n(qe['a'],3)}; la pendiente es directamente b₁ = 0.3 (p12). Dejar a₁ = 0.7 olvida deshacer el logaritmo. Exponenciar también la pendiente (e^0.3 = {n(np.exp(0.3),3)}) es deshacer de más: b₁ ya está en el exponente. Intercambiar 0.3 y 0.7 confunde ordenada con pendiente. Con el modelo correcto, y(5) = {n(qe['y5'],3)}.")

q(6, "Linealizaciones de la p13", "verdaderas", "las linealizaciones de la cátedra",
  [("y = a₃x/(b₃ + x) se linealiza graficando 1/y vs 1/x.", True),
   ("En ln y = ln a₁ + b₁x, la ordenada al origen de la recta es ln a₁.", True),
   ("y = a₂x^{b₂} se linealiza graficando ln y vs x.", False),
   ("1/y = a₄x + b₄ requiere además tomar logaritmo en ambos miembros.", False)],
  "La p13 lleva y = a₃x/(b₃ + x) a (b₃/a₃)(1/x) + 1/a₃ = 1/y: recta en 1/y vs 1/x. La p12 marca ln a₁ como la ordenada al origen de ln y = ln a₁ + b₁x. La potencia y = a₂x^{b₂} necesita logaritmo en los dos ejes (ln y vs ln x); con ln y vs x se linealiza la exponencial, no la potencia. 1/y = a₄x + b₄ ya es una recta en 1/y vs x: no hace falta ningún logaritmo.")

q(7, "Ecuaciones normales polinomio", "una", "Para una regresión polinomial simple de segundo orden con N datos, ¿cuántas ecuaciones normales e incógnitas hay?",
  [("3 ecuaciones normales y 3 incógnitas: a₀, a₁ y a₂.", True),
   ("2 ecuaciones normales y 3 incógnitas: falta un dato.", False),
   ("N ecuaciones normales, una por punto, y 3 incógnitas.", False),
   ("6 ecuaciones normales: una por cada potencia Σxₙ hasta Σxₙ⁴.", False)],
  "F = a₀ + a₁x + a₂x² tiene 3 coeficientes; se plantea ∂E/∂a₀ = ∂E/∂a₁ = ∂E/∂a₂ = 0 (p15), que da 3 ecuaciones normales con 3 incógnitas. No faltan ecuaciones: hay una por coeficiente. Los N datos entran en las sumatorias. Las sumatorias Σxₙ…Σxₙ⁴ son los coeficientes de la matriz del sistema, no ecuaciones separadas.")

q(7, "Grado del polinomio", "verdaderas", "subir el grado en la regresión polinomial simple",
  [("Según la p15, aparece inestabilidad (sistemas mal condicionados).", True),
   ("La regresión lineal es la regresión polinómica de grado 1.", True),
   ("Siempre mejora el modelo, porque el error E disminuye o en el peor caso queda igual.", False),
   ("Las ecuaciones normales dejan de ser lineales a partir de grado 3.", False)],
  "La OBS de la p15 dice que al aumentar el grado aparece inestabilidad (sistemas mal condicionados), y la p14 que la lineal es la polinómica de grado 1: ambas verdaderas. Que E baje con más parámetros no significa mejor modelo: el polinomio empieza a seguir el ruido (sobreajuste) y extrapola mal. Las ecuaciones normales siguen siendo lineales en los aᵢ para cualquier grado; lo que empeora es el condicionamiento de la matriz.")

t4 = R["tp2e4"]["tab"]
q(7, "Elegir grado (TP2 ej. 4)", "una", f"Para datos.csv del TP2 se obtuvo R² = {n(t4[0]['R2'],4)} (grado 1), {n(t4[1]['R2'],4)} (grado 2), {n(t4[2]['R2'],4)} (grado 3) y {n(t4[5]['R2'],4)} (grado 6). ¿Qué grado conviene?",
  [("Grado 2: el salto grande está de 1 a 2 y después R² casi no cambia.", True),
   ("Grado 6: es el de mayor R² y por lo tanto el mejor modelo.", False),
   ("Grado 1: todos superan 0.9, así que conviene el más simple.", False),
   ("Grado 3: es el primero que supera a todos los anteriores.", False)],
  f"De grado 1 a 2, R² sube de {n(t4[0]['R2'],4)} a {n(t4[1]['R2'],4)} y E_RMS (cátedra) baja de {n(t4[0]['Erms_cat'],2)} a {n(t4[1]['Erms_cat'],2)}; de 2 a 6 la mejora es en el cuarto decimal. Grado 6 tiene el R² más alto por definición (más parámetros) pero suma inestabilidad (p15) sin ganancia real. Grado 1 deja curvatura sin explicar: el salto a grado 2 muestra que la recta no alcanza. Grado 3 mejora a grado 2 en una diferencia despreciable.")

q(8, "Polinomial vs múltiple", "una", "¿Cuál es la diferencia entre una regresión polinomial de segundo orden y una regresión polinomial múltiple con dos variables?",
  [("La primera usa una variable (a₀ + a₁x + a₂x²); la segunda, dos (a₀ + a₁x₁ + a₂x₂).", True),
   ("La primera se resuelve por mínimos cuadrados y la segunda necesita interpolación.", False),
   ("La primera tiene 3 ecuaciones normales y la segunda tiene 2, una por variable.", False),
   ("La primera ajusta un plano y la segunda ajusta una parábola en el espacio.", False)],
  "La p14 define la polinomial simple con una única variable independiente (parábola a₀ + a₁x + a₂x²) y la p16 la múltiple con más de una variable (plano a₀ + a₁x₁ + a₂x₂). Las dos se resuelven por mínimos cuadrados con ecuaciones normales. Con dos variables el plano tiene 3 coeficientes, así que también son 3 ecuaciones, no 2. Invertir parábola y plano es el error de intercambio.")

qa = R["q_acero"]
q(8, "Predicción con ajuste múltiple", "una", "Para la dureza del acero se ajustó y = 161.3365 + 32.9688u − 0.0855v (u: % de cobre, v: temperatura en °F). ¿Qué dureza predice para u = 0.10 y v = 1200?",
  [(n(qa['dur'], 2), True), (n(qa['sin_u'], 2), False), (n(qa['signo'], 2), False), (n(R['evp1']['dur'], 2), False)],
  f"161.3365 + 32.9688·0.10 − 0.0855·1200 = 161.3365 + 3.29688 − 102.6 = {n(qa['dur'],2)}. {n(qa['sin_u'],2)} omite el término del cobre. {n(qa['signo'],2)} suma el término de la temperatura en lugar de restarlo. {n(R['evp1']['dur'],2)} es la respuesta de EV-P1 (u = 0.05, v = 1125): valores de otro enunciado.")

q(8, "Sistema de la regresión múltiple", "verdaderas", "la regresión polinomial múltiple F = a₀ + a₁x₁ + a₂x₂",
  [("Se obtienen 3 ecuaciones normales con ∂E/∂a₀ = ∂E/∂a₁ = ∂E/∂a₂ = 0.", True),
   ("El elemento (2,2) de la matriz de ecuaciones normales es Σx₁ₙ².", True),
   ("El error a minimizar es Σ(yₙ − (a₀ + a₁x₁ₙ + a₁x₂ₙ))², tal cual la p16.", False),
   ("El término Σx₁ₙx₂ₙ aparece sólo en el vector de la derecha del sistema.", False)],
  "Con tres coeficientes se plantean tres derivadas nulas y salen 3 ecuaciones normales (p16). Al derivar respecto de a₁ el coeficiente de a₁ es Σx₁ₙ², aunque la diapositiva escribe Σx₂ₙ (errata). La p16 también escribe a₁x₂ₙ en el error, que es una errata: el segundo coeficiente es a₂. Σx₁ₙx₂ₙ está en la matriz (posiciones (2,3) y (3,2)); el vector de la derecha tiene Σyₙ, Σx₁ₙyₙ y Σx₂ₙyₙ.")

q(9, "Parábola decreciente", "una", "¿Qué tipo de ajuste conviene para estos datos?",
  [("Regresión cuadrática: la caída se acelera, como una parábola con a₂ < 0.", True),
   ("Regresión lineal: el r de la recta es alto y la tendencia general es claramente decreciente.", False),
   ("Regresión exponencial decreciente, linealizando con ln y vs x.", False),
   ("Ajuste trigonométrico, porque la curva es parte de un coseno.", False)],
  f"La nube cae cada vez más rápido (concavidad negativa): es una parábola con a₂ < 0; la cuadrática da R² = {n(F['parab']['R2cuad'],3)}. La recta tiene r = {n(F['parab']['r'],2)}, alto, pero deja residuos con patrón (la recta queda por debajo de los datos en el centro y por encima en los extremos). Una exponencial decreciente es convexa: baja rápido al principio y se aplana, lo contrario de estos datos. Un tramo sin ninguna oscilación completa no justifica un ajuste trigonométrico (ET25-T2 ej. 2).", fig=fig_parab_q())

q(5, "Nube que se satura (IT25-T3)", "una", "¿Qué problemas muestra este ajuste lineal por MC?",
  [("La nube es curva y se aplana, y además la dispersión crece con x: la recta no sigue la forma.", True),
   ("Ninguno: la recta pasa por el medio de la nube y el r resulta cercano a 1.", False),
   ("El único problema es que hay pocos datos para que el ajuste sea confiable.", False),
   ("Hay que interpolar los puntos, porque la recta no pasa por ninguno de ellos.", False)],
  f"Es la nube de IT25-T3: sube rápido, se aplana, y la dispersión aumenta hacia la derecha. La recta queda por encima de los datos al principio, por debajo en la zona media y otra vez por encima donde se aplanan; los residuos tienen patrón. Acá r = {n(F['it3']['r'],2)}: no es cercano a 1 y, aunque lo fuera, no garantizaría el modelo. Hay muchos datos: el problema no es la cantidad. La recta no tiene por qué pasar por los puntos (eso sería interpolar, p3).", fig=fig_it3())

qr = R["q_rms"]
q(10, "E_RMS de la cátedra", "una", "Un ajuste sobre n = 4 datos dejó residuos 1, −2, 2 y −1. Con la fórmula de la cátedra, E_RMS = (Σ(yᵢ − F)²/(n + 1))^{1/2}, ¿cuánto vale?",
  [(n(qr['cat'], 3), True), (n(qr['usual'], 3), False), (n(qr['sin_raiz'], 3), False), (n(qr['E'], 3), False)],
  f"E = 1 + 4 + 4 + 1 = {n(qr['E'],0)}; la cátedra divide por n + 1 = 5 y saca raíz: √2 = {n(qr['cat'],3)}. {n(qr['usual'],3)} es dividir por n = 4 (RMSE usual, no la fórmula de la p18). {n(qr['sin_raiz'],3)} es E/(n+1) sin la raíz. {n(qr['E'],3)} es el error total E, que no se promedia.")

q(10, "Gráfico de residuos", "una", "Este es el gráfico de residuos de un ajuste lineal. ¿Qué indica?",
  [("La recta no capta una curvatura: residuos en U; falta un término cuadrático.", True),
   ("El ajuste es bueno, porque hay residuos positivos y negativos que se compensan entre sí.", False),
   ("Los residuos negativos del centro son outliers y hay que eliminarlos.", False),
   ("El error es aleatorio y sólo se reduce midiendo con más precisión.", False)],
  f"Residuos positivos en los extremos y negativos en el centro forman un patrón sistemático (U): la forma de los datos es curva y la recta no la sigue. Con una cuadrática R² pasa de {n(F['curva']['R2'],3)} a {n(F['curva']['R2cuad'],3)}. Que la suma se compense es propio de toda recta de MC (∂E/∂a₀ = 0 obliga a Σresiduos = 0), no prueba nada. No son outliers sueltos sino un tramo entero. Un error aleatorio no dibuja una U.", fig=fig_resid_q())

qr1 = R["q_r"]
q(11, "Calcular r", "una", "Para un ajuste se obtuvo S_T = 200 y S_R = 50. Con la definición de la cátedra, ¿cuánto vale r?",
  [(n(qr1['r'], 3), True), (n(qr1['R2'], 3), False), ("0.500", False), ("0.250", False)],
  f"r = √((S_T − S_R)/S_T) = √(150/200) = √0.75 = {n(qr1['r'],3)} (p19). 0.750 es (S_T − S_R)/S_T sin la raíz, o sea r² (el coeficiente de determinación). 0.500 es √(S_R/S_T), que mide la parte no explicada. 0.250 es S_R/S_T sin raíz. La definición de la diapositiva es la raíz de la fracción explicada.")

q(11, "ECM vs r (IT25-T4)", "verdaderas", "el error cuadrático medio y el coeficiente de correlación",
  [("El E_RMS tiene las unidades de y y mide el tamaño típico del error.", True),
   ("r compara la dispersión respecto del ajuste con la de los datos respecto de su media.", True),
   ("Los dos sólo pueden usarse para evaluar ajustes lineales, no polinomiales ni linealizados.", False),
   ("Un r cercano a 1 garantiza que el modelo elegido es el correcto.", False)],
  "El E_RMS (p18) promedia los residuos al cuadrado y saca raíz: queda en unidades de y. r = √((S_T − S_R)/S_T) (p19) compara S_R con S_T: es adimensional y relativo. Ambos se calculan con yᵢ y F(xᵢ) de cualquier modelo (la p23 da r de un cuadrático), así que no son exclusivos de la recta. Y un r alto puede convivir con un modelo equivocado (recta sobre exponencial, p11): hay que mirar los residuos.")

q(11, "Interpretar un ajuste (EV-P4)", "una", "En una población se ajustó, con r = 0.9, la relación peso (Y, kg) = −50 + 0.7·altura (X, cm). Entonces:",
  [("Hay tendencia a que, a mayor altura, mayor peso, y la nube es creciente.", True),
   ("Todas las personas que miden 170 cm pesan exactamente 69 kg.", False),
   ("El 90 % de las personas de la población está sobre la recta ajustada.", False),
   ("Hay más personas altas que personas con peso bajo en la población.", False)],
  "Pendiente positiva y r = 0.9 indican asociación creciente: en promedio, más altura va con más peso. Reemplazar 170 en la recta da 69 kg como valor esperado, no el peso de cada persona: es la sobreinterpretación determinista de EV-P4. r no es un porcentaje de puntos sobre la recta. La recta no dice nada sobre cuántas personas hay en cada rango de valores, sólo describe la relación entre las variables.")

q(12, "Caso COVID", "verdaderas", "el ejemplo de casos de COVID-19 en Tucumán de la cátedra",
  [("Con los 7 datos de la primera semana de septiembre, la recta dio r = 0.77.", True),
   ("Con la serie desde julio, el ajuste cuadrático dio r = 0.96 frente a 0.83 del lineal.", True),
   ("Que r subiera a 0.86 con más datos prueba que el contagio crece linealmente.", False),
   ("Un r alto permite extrapolar con seguridad los casos de las semanas siguientes.", False)],
  "La p21 da r = 0.77 para la primera semana (verificado: 0.7715) y la p23 compara lineal r = 0.83 con cuadrático r = 0.96. La p22 dice que mirando esas semanas los casos crecen de forma lineal, pero enseguida pregunta qué pasa con más datos hacia atrás: la serie completa es curva, así que un r mejor en una ventana no prueba linealidad. Extrapolar fuera del rango de los datos es riesgoso con cualquier r, y más si el modelo no tiene la forma correcta.")

# ------------------------------------------------------------------ desarrollo
dev = []
def d(sec_n, real, preg, modelo, practica=False, fig=None):
    it = {"id": f"m2-d{len(dev)+1:02d}", "seccion": f"m2-s{sec_n}", "real": real, "pregunta": preg, "modelo": modelo}
    if practica: it["practica"] = True
    if fig: it["fig"] = fig
    dev.append(it)

d(3, "IT25-T1 P1",
  "a. ¿Cuál es la diferencia entre hacer un ajuste lineal y un ajuste lineal por mínimos cuadrados? b. Escriba la ecuación que representa la condición de mínimos cuadrados en un ajuste lineal. c. ¿Siempre es posible encontrar un ajuste lineal usando mínimos cuadrados? Explicar.",
  P("<b>a)</b> Un <b>ajuste lineal</b> es aproximar los datos con una recta F(x) = a₀ + a₁x que describa su comportamiento con cierto error; sin criterio hay infinitas rectas posibles (p5). El ajuste lineal <b>por mínimos cuadrados</b> agrega el criterio que elige una sola: la recta cuyos coeficientes hacen mínimo el error E = Σₙ₌₁ᴺ (yₙ − (a₀ + a₁xₙ))², la suma de los cuadrados de las distancias verticales de los datos a la recta (p6). La diferencia es, entonces, que el segundo define qué significa \"mejor recta\" y da una forma de calcularla.")
  + P("<b>b)</b> Con E = Σ(yₙ − (a₀ + a₁xₙ))², la condición de mínimo es")
  + CMD("∂E/∂a₀ = 0   y   ∂E/∂a₁ = 0\n∂E/∂a₀ = −2 Σ (yₙ − a₀ − a₁xₙ) = 0\n∂E/∂a₁ = −2 Σ (yₙ − a₀ − a₁xₙ) xₙ = 0")
  + P("que ordenadas dan las ecuaciones normales N a₀ + (Σxₙ)a₁ = Σyₙ; (Σxₙ)a₀ + (Σxₙ²)a₁ = Σxₙyₙ, con solución a₀ = (Sxx·Sy − Sxy·Sx)/(N·Sxx − Sx²), a₁ = (N·Sxy − Sx·Sy)/(N·Sxx − Sx²) (p7–8).")
  + P("<b>c)</b> Hay dos lecturas; escribí las dos, primero la matemática. (1) <b>Calcularla: sí, casi siempre</b>, porque el sistema es lineal y su determinante N·Sxx − (Sx)² = N·Σ(xₙ − x̄)² sólo se anula si todos los xₙ son iguales; con al menos dos abscisas distintas hay solución única y es un mínimo (E es cuadrática convexa en a₀, a₁). (2) <b>Que sea un buen ajuste: no siempre</b>, porque si los datos tienen forma exponencial, periódica o curva la recta de MC existe pero describe mal la tendencia (p11 marca esos casos con ✗); ahí hay que linealizar (p12–13), usar un polinomio (p14) o un ajuste trigonométrico. Y para modelos no lineales en los parámetros MC no se aplica directo: se linealiza primero.")
  + AGR("La condición de existencia (determinante no nulo) y la convexidad de E no están en las diapositivas."))

d(8, "IT25-T2 P1",
  "a. ¿Cuál es la diferencia entre hacer un ajuste polinomial de segundo orden y una regresión polinomial múltiple? b. ¿Cuántas ecuaciones normales son necesarias en cada caso? ¿y cuántas incógnitas? c. ¿En ambos casos se puede implementar el método de mínimos cuadrados? Escribir la ecuación del error a minimizar en cada caso.",
  P("<b>a)</b> El ajuste polinomial de segundo orden es una <b>regresión polinomial simple</b>: \"se realiza el ajuste con un polinomio de una única variable independiente\" (p14), una parábola F(x, aᵢ) = a₀ + a₁x + a₂x². La <b>regresión polinomial múltiple</b> \"se realiza mediante un plano (más de una variable independiente)\" (p16): F(x₁, x₂, aᵢ) = a₀ + a₁x₁ + a₂x₂. La diferencia es el número de variables independientes: una (con potencias de x) frente a varias (cada una con su coeficiente). Por eso la primera modela curvatura en x y la segunda la influencia conjunta de dos factores (ej.: dureza del acero según cobre y temperatura, EV-P1).")
  + P("<b>b)</b> Una ecuación normal por cada coeficiente, porque se pide ∂E/∂aᵢ = 0 para cada uno. Polinomial de 2.º orden: <b>3 incógnitas</b> (a₀, a₁, a₂) y <b>3 ecuaciones normales</b>. Múltiple con dos variables: <b>3 incógnitas</b> (a₀, a₁, a₂) y <b>3 ecuaciones normales</b>; en general, con k variables, k + 1. Los sistemas son distintos:")
  + CMD("2.º orden:  N a₀ + (Σxₙ)a₁ + (Σxₙ²)a₂ = Σyₙ\n           (Σxₙ)a₀ + (Σxₙ²)a₁ + (Σxₙ³)a₂ = Σxₙyₙ\n           (Σxₙ²)a₀ + (Σxₙ³)a₁ + (Σxₙ⁴)a₂ = Σxₙ²yₙ\n\nMúltiple:   N a₀ + (Σx₁ₙ)a₁ + (Σx₂ₙ)a₂ = Σyₙ\n           (Σx₁ₙ)a₀ + (Σx₁ₙ²)a₁ + (Σx₁ₙx₂ₙ)a₂ = Σx₁ₙyₙ\n           (Σx₂ₙ)a₀ + (Σx₁ₙx₂ₙ)a₁ + (Σx₂ₙ²)a₂ = Σx₂ₙyₙ")
  + P("<b>c)</b> Sí, en los dos, porque ambos modelos son <b>lineales en los coeficientes</b>: al derivar E respecto de cada aᵢ e igualar a cero queda un sistema lineal que se resuelve directo. Errores a minimizar:")
  + CMD("2.º orden:  E = Σₙ₌₁ᴺ ( yₙ − (a₀ + a₁xₙ + a₂xₙ²) )²\nMúltiple:   E = Σₙ₌₁ᴺ ( yₙ − (a₀ + a₁x₁ₙ + a₂x₂ₙ) )²")
  + AGR("La p16 escribe el error múltiple con a₁x₂ₙ y la matriz con Σx₂ₙ en (2,2) y Σx^{2_{2n}} en (3,3): son erratas; lo correcto es lo de arriba (a₂x₂ₙ, Σx₁ₙ², Σx₂ₙ²). La OBS de la p15 vale para el polinomial: subir el grado lleva a sistemas mal condicionados."))

it3 = F["it3"]
d(5, "IT25-T3 P1",
  "Analice el gráfico siguiente, donde se observa una nube de puntos de datos y un ajuste lineal por mínimos cuadrados. a. ¿Qué problemas se detectan mirando el gráfico? Explicar en cualquier caso (aun si no encuentra ningún problema). b. De manera intuitiva, el coeficiente de correlación, ¿le daría cercano a 1? ¿por qué? c. ¿Qué modelo sugiere usar en este caso? Escribir la ecuación del modelo que considera mejor para ajustar y cuáles serían las incógnitas.",
  P("<i>(Gráfico recreado con datos propios de la misma forma que el del examen.)</i>")
  + P("<b>a)</b> Tres problemas. (1) <b>Forma</b>: los datos crecen rápido al principio y después se aplanan (concavidad negativa, saturación); la recta no puede seguir esa curvatura: queda por encima de los datos al principio (donde son muy negativos), por debajo en el tramo medio y otra vez por encima a la derecha donde se aplanan, o sea los residuos tienen un patrón sistemático en lugar de repartirse al azar. (2) <b>Dispersión creciente</b>: la nube se abre hacia la derecha; la recta de MC pesa igual todos los puntos y los de mayor dispersión dominan E. (3) Como consecuencia, la recta da malas predicciones en los extremos (al principio sobreestima los valores más negativos y, extrapolando, seguiría creciendo cuando los datos tienden a estabilizarse). Es el caso de la p11 \"¿siempre es bueno ajustar con una recta?\": no.")
  + P(f"<b>b)</b> No, no daría cercano a 1, porque r = √((S_T − S_R)/S_T) (p19) compara la dispersión de los datos respecto de la recta (S_R) con la dispersión respecto de su media (S_T); acá S_R es grande por dos razones: la recta no sigue la curva y la dispersión propia de los datos de la derecha es grande. Daría un valor intermedio (positivo, porque la tendencia general sube). En la recreación, r = {n(it3['r'],2)}. Aunque saliera más alto, r no prueba que el modelo sea el correcto: lo que delata el problema es el patrón de los residuos.")
  + P(f"<b>c)</b> Un modelo que crezca y se aplane. Opciones, en este orden: (1) <b>logarítmico</b> y = a₀ + a₁ ln x, incógnitas a₀ y a₁; se ajusta con la recta de MC usando X = ln x como variable (linealización). Es el que mejor representa una curva que sube rápido, se aplana y toma valores negativos cerca de x = 0. (2) <b>Polinomial de 2.º orden</b> F = a₀ + a₁x + a₂x² con a₂ < 0 (p14), incógnitas a₀, a₁, a₂ y 3 ecuaciones normales: sigue la curvatura dentro del rango pero termina bajando, así que no extrapola bien. (3) Saturación y = a₃x/(b₃ + x) (p13, 1/y vs 1/x, incógnitas a₃, b₃): tiene la forma de saturación, pero exige y > 0 y estos datos tienen valores negativos, así que no se puede linealizar así. En la recreación: R² lineal {n(it3['R2'],3)}, cuadrático {n(it3['R2cuad'],3)}, logarítmico {n(it3['R2log'],3)}. Para la dispersión creciente, agregado: mínimos cuadrados ponderados o transformar y.")
  + AGR("El modelo logarítmico, la heterocedasticidad y los mínimos cuadrados ponderados no están en las diapositivas."),
  fig=fig_it3())

d(11, "IT25-T4 P1",
  "a. Al momento de evaluar la calidad de un ajuste lineal por mínimos cuadrados, ¿cuál es la diferencia conceptual entre el error cuadrático medio y el coeficiente de correlación? b. Tanto el error cuadrático medio como el coeficiente de correlación, ¿sólo pueden ser usados en ajustes lineales? Explicar. c. Escriba la ecuación del coeficiente de correlación.",
  P("<b>a)</b> El <b>error cuadrático medio</b> E_RMS = (Σ(yᵢ − F(x, aᵢ))²/(n + 1))^{1/2} (p18) mide el <b>tamaño absoluto</b> del error típico del ajuste: está en las mismas unidades que y y depende de la escala de los datos (un E_RMS de 10 es enorme si y vale 1 y despreciable si y vale 10⁶). El <b>coeficiente de correlación</b> r (p19) mide la calidad <b>relativa</b>: qué parte de la dispersión total de los datos (S_T, respecto de su media) explica el ajuste, comparándola con la dispersión que queda respecto de F (S_R); es adimensional y va de 0 a 1, \"r cercano a 1 indica que existe una buena correlación entre el polinomio de ajuste con los datos experimentales\". Uno responde \"¿cuánto me equivoco?\" y el otro \"¿cuánto mejor que la media es este modelo?\".")
  + P("<b>b)</b> No. Los dos se calculan sólo con los datos yᵢ y los valores del ajuste F(xᵢ, aᵢ), sin importar qué forma tenga F, porque ni S_R ni S_T dependen de que F sea una recta. La p19 habla de \"polinomio de ajuste\" y en la p23 la cátedra compara r del ajuste lineal (0.83) con el del cuadrático (0.96). Sirven para comparar modelos sobre los mismos datos. Cuidado: si el modelo se ajustó linealizado (ln y), r calculado sobre la recta transformada mide la calidad en ln y, no en y.")
  + P("<b>c)</b>")
  + CMD("S_T = Σ ( yᵢ − ȳ )²          dispersión de los datos experimentales\nS_R = Σ ( yᵢ − F(x, aᵢ) )²    dispersión de los datos respecto a F\n\nr = √( (S_T − S_r) / S_T )")
  + AGR("r² = 1 − S_R/S_T es el coeficiente de determinación R². El coeficiente de Pearson de estadística tiene signo (−1..1); para la recta de MC su valor absoluto coincide con este r. La E_RMS usual divide por n (la cátedra por n + 1)."))

d(2, "F18-P4",
  "Regresión lineal: ¿cuál es el objetivo de aproximar por mínimos cuadrados?",
  P("El objetivo es encontrar la función lineal F(x, aᵢ) = a₀ + a₁x que aproxime a los datos (xᵢ, yᵢ) <b>con el menor error posible</b>, describiendo su tendencia aunque no pase por los puntos (p4–6). \"Mínimos cuadrados\" define qué error se minimiza: la suma de los cuadrados de las distancias verticales de los datos a la recta,")
  + CMD("E = Σₙ₌₁ᴺ ( yₙ − (a₀ + a₁xₙ) )²")
  + P("Se eleva al cuadrado porque así los errores positivos y negativos no se cancelan, los errores grandes pesan más y E queda derivable. Los coeficientes salen de la condición de mínimo ∂E/∂a₀ = 0 y ∂E/∂a₁ = 0, que da un sistema lineal de 2 ecuaciones normales (p7–8) con solución a₀ = (Sxx·Sy − Sxy·Sx)/(N·Sxx − Sx²), a₁ = (N·Sxy − Sx·Sy)/(N·Sxx − Sx²).")
  + P("Para qué sirve: describir la tendencia de datos con error de medición (en vez de interpolarlos con un polinomio de orden alto, p4), hacer predicciones y cuantificar la calidad del ajuste con el residuo, E_RMS y r (p18–19)."))

d(2, "WA25 (V/F)",
  "Verdadero o falso, justificar: \"El ajuste por mínimos cuadrados es aplicable a cualquier función.\"",
  P("<b>Falso</b> tal como está escrito (es un absoluto). Escribí primero el falso y después la salvedad.")
  + P("Por qué es falso: el procedimiento de la cátedra (plantear ∂E/∂aᵢ = 0 y resolver las ecuaciones normales) da un sistema lineal sólo cuando el modelo es <b>lineal en los parámetros</b>: recta (p8), polinomio (p15), plano (p16). Para un modelo como y = a₁e^{b₁x} las derivadas dan ecuaciones no lineales; la cátedra no aplica MC directo sino que primero <b>linealiza</b> (ln y = ln a₁ + b₁x, p12–13). Y aunque los coeficientes se puedan calcular, MC no garantiza un buen ajuste si el modelo no tiene la forma de los datos: la p11 muestra que la recta de MC sobre datos exponenciales o periódicos es un mal ajuste.")
  + P("La salvedad (segunda lectura): el <b>criterio</b> de MC, minimizar E = Σ(yₙ − F(xₙ, α))² (forma general de la p6), se puede plantear para cualquier F. Para F no lineal se resuelve con métodos iterativos (mínimos cuadrados no lineales) o linealizando. Si la cátedra entiende \"aplicable\" como \"se puede plantear\", la justificación anterior muestra el matiz.")
  + AGR("Los mínimos cuadrados no lineales (iterativos, curve_fit) no están en las diapositivas."))

d(1, None,
  "Variante: a) ¿Qué diferencia hay entre interpolar un conjunto de datos y ajustarlos con una curva? b) ¿Por qué la cátedra recomienda muchas veces ajustar con una recta en lugar de interpolar con un polinomio de orden alto? c) ¿Qué conviene hacer antes de elegir el modelo?",
  P("<b>a)</b> La <b>interpolación</b> asegura que la aproximación coincide con las ordenadas en las abscisas de los datos: pasa por todos los puntos (p3). El <b>ajuste de curvas</b> aproxima con otra función que describe el comportamiento de los datos con cierto error, sin pasar necesariamente por ellos (p4). Interpolar prioriza reproducir cada dato; ajustar prioriza describir la tendencia.")
  + P("<b>b)</b> Porque con datos experimentales (con error de medición) lo que interesa es la tendencia: \"muchas veces es más útil ajustar con una función lineal (recta) que describa el comportamiento de los datos (con cierto error) y no por ejemplo interpolar con un polinomio de orden alto\" (p4). El interpolante copia el ruido, oscila entre los puntos y predice mal fuera de ellos; además los polinomios de grado alto llevan a sistemas mal condicionados (p15). La pregunta clave de la p4: \"¿cuánto me interesa que pase exactamente por los puntos dados?\".")
  + P("<b>c)</b> Explorar visualmente los datos (p2): ver si hay tendencia, cuáles son los mínimos y máximos, si hay outliers, qué información da la nube y si se puede predecir. El ejemplo de la p9 arranca con \"1- Graficamos los datos\" y recién después elige la recta \"como los puntos parecen seguir una tendencia lineal\".")
  + AGR("La oscilación entre nodos de los polinomios de interpolación de grado alto (fenómeno de Runge) es de Métodos Numéricos I, no de esta unidad."))

d(6, None,
  "Variante: Los datos de un experimento parecen seguir alguno de estos modelos: y = a₂x^{b₂}, y = a₃x/(b₃ + x), y = 1/(a₄x + b₄). a) Linealice cada uno e indique qué se grafica en cada eje. b) ¿Cómo se recuperan los parámetros a partir de la recta? c) ¿La recta obtenida minimiza el mismo error que si se ajustara el modelo original?",
  P("<b>a)</b> (p13)")
  + TBL(["Modelo", "Linealización", "Ejes"], [
      ["y = a₂x^{b₂}", "ln y = ln a₂ + b₂ ln x", "ln y vs ln x"],
      ["y = a₃x/(b₃ + x)", "1/y = (b₃/a₃)(1/x) + 1/a₃", "1/y vs 1/x"],
      ["y = 1/(a₄x + b₄)", "1/y = a₄x + b₄", "1/y vs x"]])
  + P("En cada caso se transforman los datos (Xᵢ, Yᵢ) y se aplica la regresión lineal por MC de la p8 a los transformados, porque en esas variables el modelo tiene la forma F = a₀ + a₁X.")
  + P("<b>b)</b> Potencia: b₂ = pendiente, a₂ = e^{ordenada}. Saturación: ordenada = 1/a₃ → a₃ = 1/ordenada; pendiente = b₃/a₃ → b₃ = pendiente·a₃. Recíproca: a₄ = pendiente, b₄ = ordenada. Exponencial (p12): b₁ = pendiente, a₁ = e^{ordenada}.")
  + P("<b>c)</b> No exactamente: la recta minimiza Σ(Yᵢ − (a₀ + a₁Xᵢ))² en las variables transformadas (por ejemplo en ln y), no Σ(yᵢ − F(xᵢ))² en las originales. Los parámetros son buenos y suelen ser suficientes, pero el logaritmo comprime los valores grandes y da más peso relativo a los chicos. Por eso r o R² de la recta describe la calidad en el espacio transformado. Además ln y y 1/y exigen y > 0 o y ≠ 0.")
  + AGR("El inciso c y la condición y > 0 no están en las diapositivas."))

d(7, None,
  "Variante: a) Deduzca las ecuaciones normales de la regresión polinomial de segundo orden. b) Generalice: ¿cuántas ecuaciones e incógnitas tiene un polinomio de grado m? c) ¿Qué pasa al aumentar el grado?",
  P("<b>a)</b> Con F = a₀ + a₁x + a₂x², el error es E = Σ(yₙ − (a₀ + a₁xₙ + a₂xₙ²))². Derivando:")
  + CMD("∂E/∂a₀ = −2 Σ (yₙ − a₀ − a₁xₙ − a₂xₙ²) = 0\n∂E/∂a₁ = −2 Σ (yₙ − a₀ − a₁xₙ − a₂xₙ²) xₙ = 0\n∂E/∂a₂ = −2 Σ (yₙ − a₀ − a₁xₙ − a₂xₙ²) xₙ² = 0\n\nN a₀ + (Σxₙ)a₁ + (Σxₙ²)a₂ = Σyₙ\n(Σxₙ)a₀ + (Σxₙ²)a₁ + (Σxₙ³)a₂ = Σxₙyₙ\n(Σxₙ²)a₀ + (Σxₙ³)a₁ + (Σxₙ⁴)a₂ = Σxₙ²yₙ")
  + P("Es un sistema lineal de 3×3 porque F es lineal en a₀, a₁, a₂ (p15).")
  + P("<b>b)</b> Grado m: m + 1 coeficientes (a₀…aₘ), m + 1 derivadas nulas, m + 1 ecuaciones normales con sumatorias Σxₙ hasta Σxₙ^{2m}. La regresión lineal es el caso m = 1 (p14).")
  + P("<b>c)</b> \"Si se aumenta el grado del polinomio, comienza a aparecer inestabilidad (sistemas mal condicionados)\" (p15): las sumatorias Σxₙᵏ crecen muchísimo y las columnas de la matriz se parecen, así que pequeños errores de redondeo cambian mucho los coeficientes. Además, aunque E siempre baja con más grado, el polinomio empieza a seguir el ruido (sobreajuste) y extrapola mal.")
  + AGR("El sobreajuste y la explicación del mal condicionamiento por columnas casi paralelas no están en la diapositiva."))

d(10, None,
  "Variante: a) Defina residuo, error E y error cuadrático medio. b) ¿Por qué se usan los cuadrados de los residuos y no los residuos? c) ¿Qué información da el gráfico de residuos que no da E_RMS?",
  P("<b>a)</b> (p18) Residuo de cada dato: Rᵢ = yᵢ − (a₀ + a₁xᵢ), la diferencia vertical entre el dato y el ajuste. Error: E = Σ(yᵢ − F(x, aᵢ))². Error cuadrático medio: E_RMS = (Σ(yᵢ − F(x, aᵢ))²/(n + 1))^{1/2}, que queda en unidades de y y no crece con la cantidad de datos (la cátedra divide por n + 1; lo usual es n).")
  + P("<b>b)</b> Porque los residuos tienen signo: su suma se cancela (para la recta de MC, ∂E/∂a₀ = 0 obliga a ΣRᵢ = 0 siempre), así que no mide nada. Elevando al cuadrado todos suman positivo, los errores grandes pesan más y E es derivable, lo que permite plantear ∂E/∂aᵢ = 0 y obtener un sistema lineal.")
  + P("<b>c)</b> E_RMS es un solo número; el gráfico de residuos muestra <b>dónde y cómo</b> se equivoca el modelo. Si los residuos se reparten al azar alrededor de 0, el modelo captó la forma. Si forman una U o signos por tramos, falta curvatura (otro modelo); si se abren en abanico, la dispersión crece con x; un punto aislado muy lejos es un outlier.")
  + AGR("La propiedad ΣRᵢ = 0 de la recta de MC y la lectura del gráfico de residuos no están explícitas en las diapositivas."))

d(12, None,
  "Variante (ejemplo de la cátedra): con los casos de COVID de la capital de Tucumán, un ajuste lineal de la primera semana dio r = 0.77; con 19 días, r = 0.86; con la serie desde julio, el lineal dio r = 0.83 y el cuadrático r = 0.96. a) ¿r = 0.77 es satisfactorio? b) ¿El paso a 0.86 prueba que los casos crecen linealmente? c) ¿Qué conclusión se saca de la serie completa?",
  P(f"<b>a)</b> No: r = 0.77 (verificado {n(cv['r'],4)}, recta F = {n(cv['a0'],1)} + {n(cv['a1'],3)}x con los 7 datos de la p20) significa que la recta explica sólo r² ≈ {n(cv['R2'],2)} de la dispersión. Con 7 puntos tan dispersos no alcanza para describir ni para predecir: los datos suben y bajan día a día.")
  + P("<b>b)</b> No. La p22 dice que \"si miramos estas semanas los casos crecen de manera lineal\", pero enseguida pregunta \"¿qué pasa si miramos la serie con más datos hacia atrás?\": un r mejor en una ventana corta sólo describe esa ventana. r mide qué tan bien la recta sigue esos datos, no si el fenómeno es lineal.")
  + P("<b>c)</b> Con la serie completa (p23) los casos son casi nulos durante semanas y después crecen cada vez más rápido: la forma es curva. El cuadrático (r = 0.96) describe mucho mejor que la recta (0.83). El exponencial que muestra la diapositiva (\"R = ?\") explota fuera de escala: un modelo mal ajustado puede ser peor que todos. Conclusión: graficar siempre, mirar la serie completa, comparar modelos con r y con el gráfico, y no extrapolar.")
  + AGR("Los datos de las p22–23 no están en el material; sólo se verificó r = 0.77 de la p21."))

d(5, None,
  "Variante: ¿Siempre es bueno ajustar con una recta? Analice cuatro nubes: (1) tendencia lineal con dispersión; (2) puntos alrededor de 0 sin tendencia; (3) crecimiento casi nulo y luego explosivo; (4) oscilación periódica. Para cada una diga si la recta de MC sirve, por qué, y qué haría.",
  P("La recta de MC <b>siempre se puede calcular</b>, pero no siempre sirve (p11).")
  + P(f"<b>(1)</b> Sí sirve (✓): los puntos se reparten a ambos lados sin patrón; es el caso del ejemplo de la p9. r alto (en la recreación de la sección {n(F['lin']['r'],2)}).")
  + P(f"<b>(2)</b> Dudoso: la recta sale casi horizontal y r ≈ 0 (recreación: {n(F['ruido']['r'],2)}). No hay relación lineal: la recta sólo dice que y no depende de x. Hay que preguntarse si existe alguna relación o si es ruido.")
  + P(f"<b>(3)</b> No sirve (✗): es exponencial. La recta queda por encima en el medio y por debajo en los extremos, predice negativos al inicio. Aunque r parezca aceptable (recreación: {n(F['expo']['r'],2)}) el modelo es incorrecto. Se linealiza: ln y = ln a₁ + b₁x (p12).")
  + P(f"<b>(4)</b> No sirve (✗): ningún polinomio de grado bajo sigue una oscilación; la recta da r bajo (recreación: {n(F['seno']['r'],2)}). Se usa ajuste trigonométrico (U3).")
  + AGR("La lectura del caso (2) y los valores de r son de datos propios que recrean la diapositiva (que no da números)."))

# ----------------------- práctica
etx = et["exp"]
d(4, "ET25-T1 ej. 2",
  "Un programa en Python para realizar una regresión lineal por mínimos cuadrados arroja los siguientes resultados: coef_ [1146.59578596] intercept_ −4931.29430665. El gráfico muestra los puntos (3, 47), (6, 493), (7, 1215), (8, 3132), (9, 8292) y la recta. Analice el resultado: ¿es un buen modelo? ¿cuánto vale R²? ¿qué modelo propone? Calcúlelo.",
  P("<b>1. ¿Los puntos reproducen el ajuste?</b> Con los 5 puntos etiquetados (N = 5, Sx = 33, Sy = 13179, Sxx = 239, Sxy = 111288) las fórmulas de la p8 dan "
    f"a₁ = {n(et['a1'],5)} y a₀ = {n(et['a0'],5)}, frente a coef_ = 1146.59579 e intercept_ = −4931.29431 del enunciado: diferencias de {n(et['d_coef'],3)} y {n(et['d_int'],3)} (≈ 0.005 %). "
    "Los puntos reproducen el ajuste salvo el redondeo de las etiquetas (las ordenadas reales tienen decimales). En la foto no se ven otros puntos; un punto más sin etiquetar cambiaría los coeficientes mucho más que eso.")
  + P("<b>2. Diagnóstico de la recta.</b>")
  + TBL(["x", "y", "recta", "residuo"], [[str(int(a)), str(int(b)), n(c, 1), n(r_, 1)] for a, b, c, r_ in zip([3, 6, 7, 8, 9], [47, 493, 1215, 3132, 8292], et['yap'], et['res'])])
  + P(f"Los residuos siguen el patrón + − − − +: la recta queda por debajo en los extremos y por encima en el medio, que es la firma de datos curvos. Además predice {n(et['yap'][0],1)} en x = 3 para un dato positivo. "
      f"R² = 1 − S_R/S_T = 1 − {n(et['Sr'],0)}/{n(et['St'],0)} = <b>{n(et['R2'],3)}</b> (r = {n(et['r'],3)}); con los coeficientes del enunciado da lo mismo a 3 decimales ({n(et['R2_examen'],3)}). E_RMS (cátedra) = {n(et['Erms_cat'],0)}. Mal modelo.")
  + P("<b>3. Modelo propuesto: exponencial linealizado</b> y = a₁e^{b₁x} → ln y = ln a₁ + b₁x (p12). Recta de MC sobre (x, ln y):")
  + CMD(f"ln a₁ = {n(etx['lna'],4)}   →   a₁ = e^{n(etx['lna'],4)} = {n(etx['a'],4)}\nb₁    = {n(etx['b'],4)}            (cada unidad de x multiplica y por e^b₁ = {n(etx['factor'],3)})\n\ny ≈ {n(etx['a'],3)} · e^({n(etx['b'],4)} x)")
  + TBL(["x", "y", "exponencial"], [[str(int(a)), str(int(b)), n(c, 1)] for a, b, c in zip([3, 6, 7, 8, 9], [47, 493, 1215, 3132, 8292], etx['yap'])])
  + P(f"R² en ln y = {n(etx['R2_log'],4)}; recalculado en y original = {n(etx['R2_orig'],3)}; E_RMS (cátedra) = {n(etx['Erms_cat'],0)} frente a {n(et['Erms_cat'],0)} de la recta. "
      f"Una cuadrática daría R² = {n(et['cuad']['R2'],3)} pero predice {n(et['cuad']['yap'][1],1)} en x = 6 (negativo): tampoco sirve. "
      f"Extrapolando a x = 10: recta {n(et['p10_lin'],0)}, cuadrática {n(et['p10_cuad'],0)}, exponencial {n(et['p10_exp'],0)}: la elección del modelo cambia todo.")
  + AGR("R² y el análisis de residuos no están con ese nombre en las diapositivas (r² de la p19). La verificación de los puntos es una corrida propia."),
  practica=True, fig=fig_et1())

d(9, "ET25-T2 ej. 2",
  "¿Qué tipo de ajuste sugiere hacer en el caso de los datos graficados? Explicar cada caso por qué podría ser elegido y por qué no. a) Regresión lineal b) Regresión cuadrática c) Regresión exponencial d) Regresión lineal con linealización (indicar cuál) e) Ajuste trigonométrico. (Gráfico 1: crecimiento casi nulo y luego explosivo. Gráfico 2: decreciente cada vez más rápido.)",
  P("<i>(Gráficos recreados con datos propios de la misma forma que los del examen.)</i>")
  + P("<b>Gráfico 1 (exponencial).</b>")
  + UL(f"<b>a) Lineal</b> — podría: es la más simple y la tendencia general es creciente. No: la forma es curva, los residuos tienen patrón y predice valores negativos al principio (recreación: R² = {n(F['expo']['R2'],3)}).",
       "<b>b) Cuadrática</b> — podría: tiene curvatura y dentro del rango sigue bastante. No: una parábola tiene un mínimo y vuelve a subir hacia la izquierda; no reproduce el tramo plano ni el crecimiento multiplicativo.",
       "<b>c) Exponencial</b> — <b>sí</b>: y = a₁e^{b₁x} tiene exactamente esta forma (casi nula y después explosiva, p11–12). Por qué no directa: no es lineal en los parámetros, MC no se aplica directo.",
       f"<b>d) Lineal con linealización</b> — <b>la elegida</b>: ln y = ln a₁ + b₁x, graficando ln y vs x (p12), que implementa c) con la recta de MC (recreación: R² en y = {n(F['expo']['R2exp'],3)}, b₁ = {n(F['expo']['b_exp'],3)}). Condición: y > 0.",
       "<b>e) Trigonométrico</b> — no: no hay oscilación ni periodicidad.")
  + P("<b>Gráfico 2 (parábola decreciente).</b>")
  + UL(f"<b>a) Lineal</b> — podría: la tendencia es claramente decreciente y r es alto (recreación: r = {n(F['parab']['r'],2)}). No: la caída se acelera, los residuos forman un patrón.",
       f"<b>b) Cuadrática</b> — <b>la elegida</b>: F = a₀ + a₁x + a₂x² con a₂ < 0 sigue la concavidad negativa (recreación: R² = {n(F['parab']['R2cuad'],3)}, a₂ = {n(F['parab']['a2'],3)}); 3 ecuaciones normales.",
       "<b>c) Exponencial</b> — no: una exponencial decreciente (b₁ < 0) es convexa, cae rápido y se aplana; acá pasa lo contrario. Una exponencial con a₁ < 0 y b₁ > 0 sí es cóncava y decreciente, pero siempre negativa: para estos datos positivos haría falta una constante extra (y = c + a₁e^{b₁x}), que ya no se linealiza con ln y.",
       "<b>d) Linealización</b> — no hay patrón de la p13 que dé recta con esta forma; ln y no sirve si y llega a valores ≤ 0.",
       "<b>e) Trigonométrico</b> — podría parecer un tramo de coseno, pero sin una oscilación completa no hay base para un modelo periódico.")
  + AGR("La descripción de los gráficos sale del patrón; los R² son de los datos recreados, no del examen."),
  practica=True, fig=fig_et2())

ev1 = R["evp1"]
d(8, "EV-P1",
  "Veinte tipos de hojas de acero procesadas en frío tienen diferentes composiciones de cobre a diferentes temperaturas de templado. Dureza: Cobre 0.02 % → 78.9, 65.1, 55.2, 56.4; 0.10 % → 80.9, 69.7, 57.4, 55.4; 0.18 % → 85.3, 71.8, 60.7, 58.9 (para T = 1000, 1100, 1200, 1300 °F). Se realiza un ajuste lineal por mínimos cuadrados y se obtiene y = 161.3365 + 32.9688u − 0.0855v, donde u es el contenido de cobre y v la temperatura. Calcular la dureza cuando el cobre es de 0.05 % a una temperatura de 1125 °F.",
  P("Es una <b>regresión polinomial múltiple</b> (p16): un plano en dos variables independientes, y = a₀ + a₁u + a₂v. Se reemplaza:")
  + CMD(f"y = 161.3365 + 32.9688·0.05 − 0.0855·1125\n  = 161.3365 + {n(ev1['t_u'],5)} − {n(ev1['t_v'],4)}\n  = {n(ev1['dur'],3)}")
  + P(f"<b>Dureza ≈ {n(ev1['dur'],2)}</b>. Es una interpolación dentro del rango de los datos (u entre 0.02 y 0.18, v entre 1000 y 1300), así que la predicción es razonable.")
  + P(f"Verificación: armando las 3 ecuaciones normales con las 12 durezas de la tabla y resolviendo, se obtiene a₀ = {n(ev1['coef_tabla'][0],4)}, a₁ = {n(ev1['coef_tabla'][1],4)}, a₂ = {n(ev1['coef_tabla'][2],4)}: reproduce los coeficientes del enunciado con la celda leída como 57.4 (con 57.9 cambiarían a₀ y a₂). El enunciado dice \"veinte tipos\", pero la tabla tiene 12 valores.")
  + TRAMPA("Errores típicos en las opciones: olvidar el término del cobre, sumar el término de la temperatura en lugar de restarlo, o usar 0.05 como 5.", real=False),
  practica=True)

d(11, "EV-P4",
  "En una población se obtiene, con una bondad de ajuste de 0.9, que la relación entre nivel de glucemia (Y) y nivel de colesterol (X) es de Y = 20 + X/4. Entonces: a. Todos los individuos con un valor de colesterol 100, presentan glucemia 45. b. a y c son correctas. c. Las observaciones se muestran como una nube de puntos creciente. d. Hay más individuos con colesterol alto que con glucemia baja. e. Existe tendencia a que, a mayor nivel de glucemia, mayor nivel de colesterol.",
  P("La pregunta es <b>ambigua</b>: c y e describen la misma asociación positiva. Marcá <b>c</b> y, si hay lugar para justificar, explicá e.")
  + P(f"<b>a) Falsa.</b> Reemplazando X = 100 da Y = 20 + 25 = {n(R['evp4']['y100'],0)}, pero la recta da el valor <b>esperado</b> (promedio) de la glucemia para ese colesterol, no el de cada individuo: con bondad 0.9 hay dispersión alrededor de la recta. Es la sobreinterpretación determinista.")
  + P("<b>b) Falsa</b>, porque a es falsa.")
  + P("<b>c) Verdadera.</b> Pendiente 1/4 > 0 y bondad de ajuste alta: la nube es creciente y está bastante concentrada alrededor de la recta.")
  + P("<b>d) Falsa.</b> La recta no dice nada sobre cuántos individuos hay en cada rango de valores.")
  + P("<b>e) Dos lecturas.</b> (1) Si la cátedra la toma como distractor de \"inversión de variables\" (el ajuste es de Y sobre X, \"a mayor colesterol, mayor glucemia\", y e lo enuncia al revés), entonces la única respuesta es c. (2) Si se lee como asociación (la correlación es simétrica: con r > 0, a mayor glucemia también se observa en promedio mayor colesterol), e también es verdadera; pero no hay opción \"c y e\", así que la opción más segura sigue siendo c. Escribí primero c.")
  + AGR("La lectura de la recta como valor esperado y la simetría de la correlación son conceptos de estadística, no de las diapositivas."),
  practica=True)

d(4, "EV-P7",
  "De una muestra de 8 observaciones conjuntas de valores de dos variables X e Y se obtiene la siguiente información: Σxᵢ = 24, Σxᵢyᵢ = 64, Σyᵢ = 40, s_y² = 12, s_x² = 6. Opciones: Ninguna de las anteriores / La recta de regresión vale Y = 8.5 − 1.167x / La predicción para x = 4 es 4.83 / Las observaciones se muestran como una nube de puntos creciente / La recta de regresión vale Y = −1.167 + 8.5x.",
  P("Convención: varianzas y covarianza dividiendo por n (la que reproduce las opciones).")
  + CMD(f"x̄ = 24/8 = {n(e7['xm'],0)}        ȳ = 40/8 = {n(e7['ym'],0)}\ns_xy = Σxy/n − x̄·ȳ = 64/8 − 3·5 = {n(e7['sxy'],0)}\na₁ = s_xy / s_x² = −7/6 = {n(e7['b'],4)}\na₀ = ȳ − a₁·x̄ = 5 + (7/6)·3 = {n(e7['a'],1)}\n\nY = 8.5 − 1.167x        Y(4) = 8.5 − 4.667 = {n(e7['p4'],3)}\nr = s_xy/√(s_x² s_y²) = −7/√72 = {n(e7['r'],3)}")
  + P("<b>Respuesta: \"La recta de regresión vale Y = 8.5 − 1.167x\".</b> La predicción en x = 4 es 3.83, no 4.83 (corrida en 1). La covarianza negativa indica nube <b>decreciente</b>. Y = −1.167 + 8.5x intercambia ordenada y pendiente.")
  + P(f"Con n − 1 en la covarianza (s_xy = −56/7 = {n(e7['n1']['sxy'],0)}) saldría Y = {n(e7['n1']['a'],0)} − {n(-e7['n1']['b'],3)}x y predicción {n(e7['n1']['p4'],3)}: no coincide con ninguna opción, así que la convención de la cátedra es n. s_y² sólo hace falta para r.")
  + AGR("La recta desde medias y varianzas es la misma recta de MC de la p8 escrita en forma estadística; la fórmula no está en las diapositivas."),
  practica=True)

d(11, None,
  "TP2 ej. 1: Investigue otras métricas comúnmente utilizadas para evaluar la bondad del ajuste. Explique cómo se interpretan y dé ejemplos.",
  P("Las diapositivas dan el residuo, E, E_RMS y r (p18–19). Otras métricas (todo agregado):")
  + TBL(["Métrica", "Fórmula", "Interpretación"], [
      ["R² (coef. de determinación)", "1 − S_R/S_T = r²", "Fracción de la dispersión de y explicada por el modelo (0 a 1)."],
      ["R² ajustado", "1 − (1 − R²)(n − 1)/(n − p − 1)", "Penaliza agregar parámetros: sirve para comparar polinomios de distinto grado."],
      ["MSE / RMSE", "Σ(yᵢ − Fᵢ)²/n y su raíz", "Error típico en unidades de y (la cátedra usa n + 1)."],
      ["MAE", "Σ|yᵢ − Fᵢ|/n", "Error medio absoluto; menos sensible a outliers que RMSE."],
      ["MAPE", "(100/n)Σ|(yᵢ − Fᵢ)/yᵢ|", "Error porcentual; no sirve si hay yᵢ cerca de 0."],
      ["AIC / BIC", "n·ln(S_R/n) + 2p (AIC)", "Comparan modelos balanceando ajuste y número de parámetros; menor es mejor."],
      ["Gráfico de residuos", "Rᵢ vs xᵢ", "Patrón ⇒ modelo inadecuado; abanico ⇒ dispersión variable."]])
  + P(f"Ejemplo con el ajuste de la p9 (F = −0.734 + 0.739x, n = 6, p = 1): R² = {n(ej['R2'],4)}, R² ajustado = {n(ej['R2adj'],4)}, MAE = {n(ej['MAE'],4)}, E_RMS(cátedra) = {n(ej['Erms_cat'],4)}, RMSE(n) = {n(ej['Erms_n'],4)}. Lectura: la recta explica el {n(100*ej['R2'],1)} % de la dispersión y se equivoca en promedio {n(ej['MAE'],2)} unidades de y (MAE).")
  + AGR("Todo el contenido de este ejercicio (salvo E, E_RMS y r) es agregado."),
  practica=True)

g = R["tp2e2"]
d(4, None,
  "TP2 ej. 2: Para estimar g se lanza una pelota desde lo alto de una torre de 30 m y se mide la velocidad en distintos pisos: X (m) = 0, 5, 10, 15, 20, 25; V (m/s) = 0, 9.85, 14.32, 17.63, 19.34, 22.41. Con el eje positivo hacia abajo vale v² = 2gx. Usando regresión por mínimos cuadrados, estime g.",
  P("<b>Linealización:</b> llamando Y = v², el modelo es Y = (2g)·x: una recta en (x, Y). Valores de Y = v²:")
  + TBL(["x", "0", "5", "10", "15", "20", "25"], [["v²"] + [n(t_, 2) for t_ in g['v2']]])
  + P("<b>Opción 1 (modelo exacto, recta por el origen):</b> Y = k·x con una sola incógnita k = 2g. E = Σ(Yᵢ − kxᵢ)²; dE/dk = −2Σ(Yᵢ − kxᵢ)xᵢ = 0 → k = ΣxᵢYᵢ / Σxᵢ².")
  + CMD(f"ΣxᵢYᵢ = {n(g['SxY'],4)}     Σxᵢ² = {n(g['Sxx'],0)}\nk = {n(g['k0'],4)}   →   g = k/2 = {n(g['g0'],3)} m/s²     (R² = {n(g['R2_0'],4)})")
  + P(f"<b>Opción 2 (recta completa de la p8):</b> Y = a₀ + a₁x con N = 6, Sx = {n(g['lin']['Sx'],0)}, Sy = {n(g['lin']['Sy'],4)}, Sxx = {n(g['lin']['Sxx'],0)}, Sxy = {n(g['lin']['Sxy'],4)} → a₁ = {n(g['a1'],4)}, a₀ = {n(g['a0'],4)}; g = a₁/2 = {n(g['g1'],3)} m/s² (R² = {n(g['R2_1'],4)}). La ordenada debería ser 0 por la física; el valor chico que sale es error experimental.")
  + P(f"<b>Control con potencia</b> (p13): ln v = ln c + p ln x sin el punto x = 0 da p = {n(g['pot']['p'],4)} (≈ 1/2, confirma v ∝ √x) y g = c²/2 = {n(g['pot']['g'],3)} m/s².")
  + P(f"<b>Resultado:</b> g ≈ {n(g['g0'],2)} m/s² (opción 1).")
  + AGR("La recta forzada por el origen no está en las diapositivas (es MC con una sola incógnita). El valor de referencia 9.81 m/s² es agregado: la estimación queda apenas por encima."),
  practica=True)

p3 = R["tp2e3"]
d(6, None,
  "TP2 ej. 3: Población mundial (billones): 1850 → 1.3; 1900 → 1.6; 1950 → 3; 1980 → 4.4; 2000 → 6. Asumiendo p = b·e^{mx} (x: año): a. Linealice la ecuación y, usando mínimos cuadrados, calcule m y b. b. Estime la población en 1970.",
  P("<b>a)</b> ln p = ln b + m·x: recta en (x, ln p) con pendiente m y ordenada ln b (p12).")
  + TBL(["x (año)", "1850", "1900", "1950", "1980", "2000"], [["ln p"] + [n(t_, 4) for t_ in p3['lnp']]])
  + CMD(f"N = 5   Sx = {n(p3['lin']['Sx'],0)}   Sy = {n(p3['lin']['Sy'],4)}   Sxx = {n(p3['lin']['Sxx'],0)}   Sxy = {n(p3['lin']['Sxy'],4)}\n\nm    = (N·Sxy − Sx·Sy)/(N·Sxx − Sx²)   = {n(p3['m'],6)} 1/año\nln b = (Sxx·Sy − Sxy·Sx)/(N·Sxx − Sx²) = {n(p3['lnb'],4)}\nb    = e^(ln b) = {p3['b']:.4e} billones")
  + P(f"R² en ln p = {n(p3['R2_log'],4)}.")
  + P(f"<b>b)</b> p(1970) = b·e^{{m·1970}} = <b>{n(p3['p1970'],3)} billones</b>.")
  + P(f"Observación numérica: con x = año, Sxx y Sx² son números enormes y casi iguales (mal condicionamiento) y b sale ridículamente chico. Centrando t = año − 1850 se obtiene el mismo m = {n(p3['t']['m'],6)}, b = {n(p3['t']['b'],4)} (población estimada en 1850) y la misma estimación {n(p3['t']['p1970'],3)} billones para 1970. Con este m la población se duplica cada ln 2/m = {n(p3['duplica'],1)} años.")
  + AGR("El cambio de variable t = año − 1850 y el tiempo de duplicación no están en el enunciado ni en las diapositivas."),
  practica=True)

tb = R["tp2e4"]["tab"]
d(7, None,
  "TP2 ej. 4: Lea los datos de datos.csv (101 puntos, x de 0 a 10). a. Grafique los datos. b. Realice una aproximación polinomial simple; elija el orden que considere más adecuado. c. ¿Qué coeficientes obtuvo? Grafique el ajuste sobre los datos originales.",
  P("<b>a)</b> El gráfico muestra una nube creciente, curva (concavidad positiva) y con ruido: no es una recta.")
  + P("<b>b)</b> Se ajustan polinomios de grado 1 a 6 (np.polyfit) y se comparan:")
  + TBL(["Grado", "R²", "E_RMS (cátedra)", "S_R"], [[str(t_['g']), n(t_['R2'], 4), n(t_['Erms_cat'], 3), n(t_['Sr'], 1)] for t_ in tb])
  + P("De grado 1 a 2 el error cae a menos de la mitad y R² salta; de 2 en adelante la mejora es mínima (cuarto decimal) y cada grado suma inestabilidad (p15). <b>Grado 2.</b>")
  + P("<b>c)</b> Coeficientes (np.polyfit devuelve de mayor a menor potencia):")
  + CMD(f"y ≈ {n(tb[1]['c'][0],4)} x² + {n(tb[1]['c'][1],4)} x {n(tb[1]['c'][2],4).replace('−','− ')}\n\nnp.polyfit(x, y, 2) → [{n(tb[1]['c'][0],6)}, {n(tb[1]['c'][1],6)}, {n(tb[1]['c'][2],6)}]")
  + P(f"En el gráfico la parábola sigue la nube y la recta (grado 1, y ≈ {n(tb[0]['c'][0],3)}x {n(tb[0]['c'][1],3).replace('−','− ')}) queda por encima en el medio y por debajo en los extremos.")
  + TRAMPA("Dar la tabla de R² por grado y preguntar qué grado elegir; el distractor es \"el de mayor R²\".", real=False),
  practica=True, fig=fig_tp4())

d(4, None,
  "Con los datos del ejemplo de la cátedra, x = 1, 2, 3, 4, 5, 6 e y = 0, 0.6, 1.77, 1.92, 3.31, 3.52: a) calcule la recta de MC desde las sumatorias; b) calcule los residuos, E, E_RMS (fórmula de la cátedra) y r; c) prediga y para x = 7.",
  P(f"<b>a)</b> N = 6, Sx = {n(ej['Sx'],0)}, Sy = {n(ej['Sy'],2)}, Sxx = {n(ej['Sxx'],0)}, Sxy = {n(ej['Sxy'],2)}.")
  + CMD(f"N·Sxx − Sx² = 6·91 − 21² = {n(6*91-21**2,0)}\na₀ = (91·11.12 − 51.86·21)/105 = {n(ej['a0'],4)}\na₁ = (6·51.86 − 21·11.12)/105 = {n(ej['a1'],4)}\nF(x) = {n(ej['a0'],3)} + {n(ej['a1'],3)}x      (la diapositiva: −0.734 + 0.739x)")
  + P("<b>b)</b>")
  + TBL(["x", "1", "2", "3", "4", "5", "6"], [["Rᵢ"] + [n(t_, 4) for t_ in ej['res']]])
  + CMD(f"E = ΣRᵢ² = {n(ej['E'],4)}\nE_RMS = (E/(n+1))^(1/2) = (E/7)^(1/2) = {n(ej['Erms_cat'],4)}\nS_T = Σ(yᵢ − ȳ)² = {n(ej['St'],4)}      S_R = {n(ej['Sr'],4)}\nr = √((S_T − S_R)/S_T) = {n(ej['r'],4)}")
  + P(f"<b>c)</b> F(7) = {n(ej['a0'],4)} + {n(ej['a1'],4)}·7 = {n(ej['pred7'],3)} (extrapolación, un paso fuera del rango).")
  + TRAMPA("Dar las sumatorias y pedir a₀, a₁ o la predicción con opciones que intercambian a₀ y a₁, o E_RMS dividiendo por n en lugar de n + 1.", real=False),
  practica=True)

v = R["evp7v"]
d(4, None,
  f"Variante de EV-P7: de n = {v['n']} observaciones se obtiene Σxᵢ = {v['Sx']}, Σyᵢ = {v['Sy']}, Σxᵢyᵢ = {v['Sxy']}, s_x² = {v['sx2']}, s_y² = {v['sy2']} (dividiendo por n). Halle la recta de regresión, la predicción para x = 7 y el coeficiente de correlación. ¿La nube es creciente o decreciente?",
  CMD(f"x̄ = {v['Sx']}/{v['n']} = {n(v['xm'],0)}      ȳ = {v['Sy']}/{v['n']} = {n(v['ym'],0)}\ns_xy = {v['Sxy']}/{v['n']} − {n(v['xm'],0)}·{n(v['ym'],0)} = {n(v['sxy'],0)}\na₁ = s_xy/s_x² = {n(v['b'],0)}      a₀ = ȳ − a₁x̄ = {n(v['a'],0)}\nY = {n(v['a'],0)} + {n(v['b'],0)}x      Y(7) = {n(v['p7'],0)}\nr = s_xy/√(s_x²·s_y²) = {n(v['sxy'],0)}/√{n(v['sx2']*v['sy2'],0)} = {n(v['r'],3)}")
  + P("Covarianza positiva: nube <b>creciente</b>. r = 0.667: asociación moderada.")
  + TRAMPA(f"Así lo tomaría la cátedra: opción con los coeficientes intercambiados (Y = {n(v['b'],0)} + {n(v['a'],0)}x, que en x = 7 da {n(v['p7_swap'],0)}) y la predicción corrida en 1.", real=False),
  practica=True)

cq_ = R["cuad"]
d(7, None,
  "Ajuste una parábola F = a₀ + a₁x + a₂x² por mínimos cuadrados a x = 0, 1, 2, 3, 4; y = 2.1, 7.7, 13.6, 27.2, 40.9 armando las ecuaciones normales. Compare con la recta y prediga y(5).",
  P("Sumatorias: N = 5, Σx = 10, Σx² = 30, Σx³ = 100, Σx⁴ = 354, Σy = 91.5, Σxy = 280.1, Σx²y = 961.3. Ecuaciones normales (p15):")
  + CMD(" 5 a₀ +  10 a₁ +  30 a₂ =  91.5\n10 a₀ +  30 a₁ + 100 a₂ = 280.1\n30 a₀ + 100 a₁ + 354 a₂ = 961.3")
  + CMD(f"a₀ = {n(cq_['a'][0],4)}   a₁ = {n(cq_['a'][1],4)}   a₂ = {n(cq_['a'][2],4)}\nF(x) = {n(cq_['a'][0],3)} + {n(cq_['a'][1],3)}x + {n(cq_['a'][2],3)}x²\nF(5) = {n(cq_['p5'],2)}")
  + P(f"R² cuadrático = {n(cq_['R2'],4)} frente a {n(cq_['R2lin'],4)} de la recta: la curvatura es real (la recta deja residuos con patrón). np.polyfit da los mismos coeficientes.")
  + TRAMPA("Pedir el sistema de ecuaciones normales con las sumatorias correctas (el error típico es poner Σx² donde va Σx³) o cuántas ecuaciones hay.", real=False),
  practica=True)

pw = R["pot"]
d(6, None,
  "Los datos x = 1, 2, 3, 4, 5; y = 0.5, 1.7, 3.4, 5.7, 8.4 parecen seguir y = a₂x^{b₂}. Linealice, calcule a₂ y b₂ por mínimos cuadrados y estime y(6).",
  P("ln y = ln a₂ + b₂ ln x: recta en (ln x, ln y) (p13).")
  + TBL(["ln x"] + [n(t_, 4) for t_ in pw['lnx']], [["ln y"] + [n(t_, 4) for t_ in pw['lny']]])
  + CMD(f"b₂ = pendiente = {n(pw['b2'],4)}\nln a₂ = ordenada = {n(pw['lna2'],4)}   →   a₂ = {n(pw['a2'],4)}\ny ≈ {n(pw['a2'],3)}·x^{n(pw['b2'],3)}        y(6) = {n(pw['p6'],2)}")
  + P(f"R² en (ln x, ln y) = {n(pw['R2_log'],5)}: el patrón potencia describe los datos casi perfecto.")
  + TRAMPA("Dar la recta ln y vs ln x y preguntar el modelo; distractores con la exponencial (ln y vs x) o dejando a₂ = ordenada sin exponenciar.", real=False),
  practica=True)

mm = R["mult"]
d(8, None,
  "Ajuste el plano y = a₀ + a₁x₁ + a₂x₂ a los datos (x₁, x₂, y): (0, 1, 4.1), (1, 0, 4.9), (1, 2, 8.8), (2, 1, 9.2), (2, 3, 13.1), (3, 2, 12.8). Prediga y para x₁ = 2, x₂ = 2.",
  P("Sistema de ecuaciones normales (p16, con las erratas corregidas):")
  + CMD(f"⎡ N      Σx₁    Σx₂   ⎤   ⎡ {int(mm['M'][0,0])}  {int(mm['M'][0,1])}  {int(mm['M'][0,2])} ⎤\n⎢ Σx₁    Σx₁²   Σx₁x₂ ⎥ = ⎢ {int(mm['M'][1,0])} {int(mm['M'][1,1])} {int(mm['M'][1,2])} ⎥\n⎣ Σx₂    Σx₁x₂  Σx₂²  ⎦   ⎣ {int(mm['M'][2,0])} {int(mm['M'][2,1])} {int(mm['M'][2,2])} ⎦\n\nlado derecho: Σy = {n(mm['b'][0],1)}, Σx₁y = {n(mm['b'][1],1)}, Σx₂y = {n(mm['b'][2],1)}")
  + CMD(f"a₀ = {n(mm['a'][0],4)}   a₁ = {n(mm['a'][1],4)}   a₂ = {n(mm['a'][2],4)}\ny(2, 2) = {n(mm['p'],3)}        R² = {n(mm['R2'],4)}")
  + TRAMPA("Como EV-P1, dar el plano ya ajustado y pedir una predicción; o pedir la matriz de ecuaciones normales (ojo con las erratas de la p16).", real=False),
  practica=True)

d(7, None,
  "TP2 ej. 4 (gráfico del ajuste): con el polinomio elegido para datos.csv, compare gráficamente la recta y la parábola sobre los datos. ¿Qué se ve en los residuos de cada uno?",
  P(f"Ajustes (np.polyfit): grado 1: y ≈ {n(tb[0]['c'][0],3)}x {n(tb[0]['c'][1],3).replace('−','− ')}; grado 2: y ≈ {n(tb[1]['c'][0],3)}x² + {n(tb[1]['c'][1],3)}x {n(tb[1]['c'][2],3).replace('−','− ')}.")
  + P(f"La recta queda por encima de los datos en el centro y por debajo en los extremos: sus residuos forman una U (patrón), E_RMS(cátedra) = {n(tb[0]['Erms_cat'],2)}. La parábola pasa por el medio de la nube en todo el rango y sus residuos quedan repartidos al azar, E_RMS = {n(tb[1]['Erms_cat'],2)}: lo que queda es el ruido de los datos, que ningún grado mayor elimina (grado 6: {n(tb[5]['Erms_cat'],2)}).")
  + TRAMPA("Mostrar el gráfico de los dos ajustes y preguntar cuál elegir y por qué; la justificación que se espera es la forma de los residuos, no sólo el número.", real=False),
  practica=True, fig=fig_tp4_ajuste())

unit = {
    "id": 2, "titulo": "Ajuste de curvas", "corto": "Ajuste de curvas",
    "fuente": "ajuste.pdf (23 diapositivas) + TP2", "peso": 2,
    "intro": intro, "secciones": S, "conexiones": conexiones, "fidelidad": fidelidad,
    "auditoria": auditoria, "glosario": glosario, "quiz": Q, "desarrollo": dev,
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "m2.js")
with open(out, "w", encoding="utf-8") as f:
    f.write("ASO.addUnit(" + json.dumps(unit, ensure_ascii=False, indent=1) + ");\n")
print("ok", len(Q), "preguntas,", len(S), "secciones,", len(dev), "desarrollo")
