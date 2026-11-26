#!/usr/bin/env python3
# Genera /home/claude/hub/data/m4.js — U4: Ecuaciones diferenciales ordinarias (EDO)
# Todos los números salen de verif_u4.py (se importa: R = resultados).
import json, os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verif_u4 as V
from figs import nueva, svg

R = V.R
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "m4.js")

# ---------------------------------------------------------------- helpers HTML
def D(term, txt, src):
    return f'<div class="def"><b>{term}</b><p>{txt}</p><span class="src">{src}</span></div>'
def AGR(txt):
    return f'<div class="agregado"><b>Agregado</b><p>{txt}</p></div>'
def TRAMPA(txt, real=True):
    rot = "Así lo toma la cátedra" if real else "Así lo tomaría la cátedra"
    return f'<div class="trampa"><b>{rot}</b><p>{txt}</p></div>'
def EJ(txt):
    return f'<div class="ej"><b>Ejemplo</b><p>{txt}</p></div>'
def TBL(head, rows):
    hh = "".join(f"<th>{x}</th>" for x in head)
    b = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<div class="tbl"><table><thead><tr>{hh}</tr></thead><tbody>{b}</tbody></table></div>'
def UL(*items): return "<ul>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>"
def OL(*items): return "<ol>" + "".join(f"<li>{i}</li>" for i in items) + "</ol>"
P = lambda t: f"<p>{t}</p>"
H4 = lambda t: f"<h4>{t}</h4>"
CMD = lambda t: f'<pre class="cmd">{t}</pre>'
FIG = lambda s: f'<figure class="fig">{s}</figure>'

def c(x, nd=4):
    """número con nd decimales y coma decimal."""
    return (f"{x:.{nd}f}").replace(".", ",").replace("-", "−")
def g3(x):
    """3 cifras significativas, sin notación científica, coma decimal."""
    if x == 0: return "0"
    d = 2 - int(math.floor(math.log10(abs(x))))
    s = f"{x:.{d}f}" if d > 0 else str(int(round(x, d)))
    return s.replace(".", ",").replace("-", "−")
def gs(x, k=3):
    """k cifras significativas con notación ×10^ para chicos."""
    if x == 0: return "0"
    e = int(math.floor(math.log10(abs(x))))
    if -3 <= e <= 4:
        d = max(k - 1 - e, 0)
        return f"{x:.{d}f}".replace(".", ",").replace("-", "−")
    m = x / 10 ** e
    sup = str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{m:.{k-1}f}".replace(".", ",").replace("-", "−") + "×10" + sup

def tabla_err(filas, dec_real=6, aprox=g3, err=g3):
    rows = []
    for t, yr, ya, el, eg in filas:
        rows.append([c(t, 1), c(yr, dec_real), aprox(ya), "—" if el is None else err(el), "—" if t == filas[0][0] else err(eg)])
    return TBL(["t", "y real", "y aprox", "E local = |yᵢ₊₁ − yᵢ|", "E global = |y(tᵢ) − yᵢ|"], rows)

M = R["multi"]; RK = R["rk"]; EMP = R["emp"]

import sympy as _sp
_DER = {1: "y′", 2: "y″", 3: "y‴", 4: "y⁗", 5: "y⁽⁵⁾", 6: "y⁽⁶⁾"}
_SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
def TY(nombre):
    """Primer término del error local de un multipaso, desde verif_u4 (p. ej. '−(2/9)·h³·y‴')."""
    m = M[nombre]; e = _sp.sympify(m["c0"]); sym = list(e.free_symbols)[0]
    k = int(str(sym)[1:]); co = _sp.Rational(e / sym)
    sg = "−" if co < 0 else ""; co = abs(co)
    cs = "" if co == 1 else (str(co) + "·" if co.q == 1 else f"({co.p}/{co.q})·")
    pot = "" if m["k0"] == 1 else str(m["k0"]).translate(_SUP)
    return f"{sg}{cs}h{pot}·{_DER[k]}"

# ---------------------------------------------------------------- figuras
def fig_stiff_cat():
    fig, ax = nueva()
    t = R["stiff_cat_t"]
    ax.plot(t, R["stiff_cat_ex"], "k-", lw=2.2, label="solución real")
    for nom, st in (("Euler", "o-"), ("RK2", "s-"), ("RK4", "^-")):
        ax.plot(t, R["stiff_cat"][nom], st, ms=4, label=f"{nom}, h = 0,1")
    ax.set_xlabel("t"); ax.set_ylabel("y")
    ax.set_title("y′ = −21y + e^(−t), y(0) = 0, h = 0,1 (verif_u4.py)")
    ax.legend(fontsize=8, ncol=2)
    return svg(fig)

def fig_tp_b():
    fig, ax = nueva(1, 2, ancho=8.4, alto=3.4)
    ts = [-1 + 0.2 * i for i in range(6)]
    tf = [-1 + 0.01 * i for i in range(101)]
    ax[0].plot(tf, [t - math.exp(-5 * t) for t in tf], "k-", lw=2.2, label="exacta")
    for nom, st in (("Euler", "o-"), ("RK2", "s-"), ("RK4", "^-"), ("Heun", "d-"), ("AB2", "x-")):
        ax[0].plot(ts, R["tp_b"][0.2][nom][2], st, ms=4, label=nom)
    ax[0].set_ylim(-160, 40); ax[0].set_xlabel("t"); ax[0].set_title("h = 0,2 (AB2 se va de escala)")
    ax[0].legend(fontsize=7, ncol=2)
    for hh in (0.2, 0.01):
        n = int(round(1 / hh))
        tt = [-1 + hh * i for i in range(n + 1)]
        e = [abs(y - (t - math.exp(-5 * t))) for y, t in zip(R["tp_b"][hh]["RK4"][2], tt)]
        ax[1].semilogy(tt, [max(v, 1e-12) for v in e], "o-" if hh == 0.2 else "-", ms=4, label=f"RK4, h = {str(hh).replace('.', ',')}")
        e = [abs(y - (t - math.exp(-5 * t))) for y, t in zip(R["tp_b"][hh]["Euler"][2], tt)]
        ax[1].semilogy(tt, [max(v, 1e-12) for v in e], "s--" if hh == 0.2 else "--", ms=4, label=f"Euler, h = {str(hh).replace('.', ',')}")
    ax[1].set_xlabel("t"); ax[1].set_title("|E global| (escala log)"); ax[1].legend(fontsize=7)
    return svg(fig)

def fig_q_stiff():
    fig, ax = nueva()
    t = R["stiff_cat_t"]
    ax.plot(t, R["stiff_cat_ex"], "k-", lw=2.2, label="solución real")
    ax.plot(t, R["stiff_cat"]["Euler"], "o-", ms=4, label="método A")
    ax.plot(t, R["stiff_cat"]["RK4"], "^-", ms=4, label="método B")
    ax.set_xlabel("t"); ax.set_ylabel("y"); ax.legend(fontsize=8)
    ax.set_title("y′ = −21y + e^(−t), y(0) = 0, h = 0,1")
    return svg(fig)

def fig_q_t4():
    fig, ax = nueva()
    tf = [i / 100 for i in range(101)]
    ax.plot(tf, [math.exp(t) for t in tf], "k-", lw=2.2, label="solución del PVI: eᵗ")
    ax.plot(tf, [math.exp(t / 2) for t in tf], "k:", lw=1.6, label="e^(t/2)")
    for n, st in ((10, "o"), (100, "-")):
        tt, yy = R["fig_t4"][n]
        ax.plot(tt, yy, st, ms=4, label=f"método, h = {str(1/n).replace('.', ',')}")
    ax.set_xlabel("t"); ax.set_ylabel("y"); ax.legend(fontsize=8)
    ax.set_title("yₙ₊₁ = yₙ₋₁ + h·y′ₙ aplicado a y′ = y, y(0) = 1")
    return svg(fig)

def fig_q_loglog():
    fig, ax = nueva()
    hs = R["fig_h"]["h"]
    for nom, st, lab in (("Euler", "o-", "curva 1"), ("RK2", "s-", "curva 2"), ("RK4", "^-", "curva 3")):
        ax.loglog(hs, R["fig_h"][nom], st, ms=5, label=lab)
    ax.minorticks_off(); ax.set_xticks(hs); ax.set_xticklabels([str(v).replace(".", ",") for v in hs])
    ax.set_xlabel("h"); ax.set_ylabel("|E global| en t = 2")
    ax.set_title("y′ = 1 + y/t, y(1) = 2: error final según h"); ax.legend(fontsize=8)
    return svg(fig)

def fig_q_eg():
    fig, ax = nueva()
    t = R["fig_eg"]["t"]
    ax.semilogy(t[1:], R["fig_eg"]["euler"][1:], "o-", ms=4, label="curva P")
    ax.semilogy(t[1:], R["fig_eg"]["rk4"][1:], "^-", ms=4, label="curva Q")
    ax.set_xlabel("t"); ax.set_ylabel("|E global| (log)")
    ax.set_title("y′ = 1 + y/t, y(1) = 2, h = 0,1, dos métodos de un paso"); ax.legend(fontsize=8)
    return svg(fig)

F_STIFF_CAT = fig_stiff_cat()

# ---------------------------------------------------------------- SECCIONES
S = []
def sec(i, titulo, estrellas, pags, *html):
    S.append({"id": f"m4-s{i}", "titulo": titulo, "estrellas": estrellas, "pags": pags, "html": "".join(html)})

sec(1, "Ecuaciones diferenciales: definiciones y Ley de Malthus", 1, "EDO1 p2–6",
    P("Muchos sistemas se describen por cómo <b>cambia</b> una variable respecto de otra (una tasa), no por la variable misma. Cuando la hipótesis es una tasa de cambio, el modelo es una ecuación diferencial y la incógnita es una función."),
    P("Ejemplos de la cátedra: dinámica poblacional <code>dP(t)/dt = rP(t)</code> (la población cambia respecto del tiempo, variable independiente; la incógnita es P(t)) y la ecuación de la onda <code>∂²u/∂t² = c²·∂²u/∂x²</code> (derivadas parciales: no es EDO)."),
    D("Ecuación Diferencial", "Ecuación que implica una o más derivadas de una función desconocida (incógnita).", "EDO1 p4"),
    D("Ecuación Diferencial Ordinaria (EDO)", "Sólo poseen derivadas de una sola variable independiente.", "EDO1 p4"),
    D("EDO de Primer Orden", "Es una EDO que sólo tiene derivadas de primer orden y tiene la forma: y′ = f(x, y).", "EDO1 p4"),
    D("EDO de Orden Superior", "Es una EDO en la que se encuentran derivadas de orden superior (y′, y″, y‴, etc).", "EDO1 p5"),
    CMD("y⁽ᵐ⁾(t) = f(x, y, y′, y″, …, y⁽ᵐ⁻¹⁾),   a ≤ t ≤ b"),
    P("Ejemplos de la diapositiva: <code>y″ − 2y′ + y = t·eᵗ − t</code> (orden 2, f(t, y, y′)) y <code>x²y‴ − 2xy′ + y − sin(x) = x·eˣ − ln t</code> (orden 3). Las preguntas que deja la cátedra: ¿cuál es el orden?, ¿cuál es la incógnita?, ¿cómo es f?"),
    H4("Ley de Malthus (p6)"),
    P("<code>dP/dt = rP</code>. Suponiendo r = 1: <code>P(t) = e^(t·c₁)</code>; haciendo variar c₁ hay <b>infinitas curvas</b>. Con la condición <code>P(t₀) = P₀</code> se despeja c₁ y queda <b>una curva particular</b>. Esa es la idea del problema de valor inicial: la EDO sola da una familia; la condición inicial elige una."),
    AGR("Con la forma de la diapositiva, P(t) = e^(t·c₁) no es la solución general de P′ = P (que es P = C·eᵗ); la idea que la cátedra quiere mostrar es la misma: una constante libre → familia de curvas; la condición inicial la fija."),
)

sec(2, "Problema de valor inicial (PVI): existencia, unicidad y qué es una solución", 1, "EDO1 p7–10",
    P("Antes de aproximar hay que saber si hay algo que aproximar: el PVI fija la curva y el teorema de existencia y unicidad garantiza que esa curva exista y sea una sola."),
    CMD("PVI:   y′ = f(x, y)\n       y(x₀) = y₀"),
    P("La diapositiva también muestra el PVI de un <b>sistema</b> (<code>u′ = f₁(t,u,v)</code>, <code>v′ = f₂(t,u,v)</code>, <code>u(t₀) = u₀</code>, <code>v(t₀) = v₀</code>) y el de una EDO de <b>orden 2</b> (<code>y″ − 2y′ + y = t·eᵗ − t</code>, <code>y(x₀) = y₀</code>, <code>y′(x₀) = y₁</code>): una condición inicial por cada ecuación / por cada derivada hasta el orden − 1."),
    D("Teorema de existencia y unicidad de la solución", "Si la función f definida en una región D satisface la condición de Lipschitz, para algún L&gt;0 (cte de Lipschitz), entonces se puede definir una región en la que el PVI tiene una solución única. ‖f(t, y) − f(t, ŷ)‖ ≤ L‖y − ŷ‖  ∀(t, y), (t, ŷ) ∈ D", "EDO1 p8"),
    D("¿Cuándo y es solución del PVI?", "La función y(x) es solución del PVI sobre un intervalo I si, ∀x ∈ I: 1) (x₀, y₀) ∈ al dominio; 2) y(x₀) = y₀; 3) y′(x) existe y además y′ = f(x, y).", "EDO1 p9"),
    H4("Solución analítica vs numérica (p10)"),
    TBL(["", "Analítica", "Numérica"], [
        ["Resultado", "Una función y (una curva)", "Una secuencia de puntos (xₙ, yₙ) que aproximan el valor de la función"],
        ["Cómo", "Se resuelve la EDO", "Partición (en gral. regular) de [a, b] con a = x₀; a partir de la condición inicial se estiman los valores en los siguientes puntos de la grilla"],
    ]),
)

sec(3, "Solución numérica: la grilla, h y n", 2, "EDO1 p11–14; EDO2 p2–4",
    P("Un método numérico no devuelve una función: devuelve valores en puntos. Hay que decidir dónde (la grilla), hasta dónde y con qué paso."),
    TBL(["x", "y(x)", "yᵢ (valor aproximado, método EDO)"], [
        ["a = x₀", "y(x₀) = y₀", ""], ["x₁ = x₀ + h", "—", "y₁ = MÉTODO"], ["x₂ = x₁ + h", "—", "y₂ = MÉTODO"]]),
    P("Se obtiene la secuencia <code>(x₀, y(x₀)), (x₁, y₁), (x₂, y₂), (x₃, y₃), …</code>. Las tres preguntas de la cátedra:"),
    UL("<b>¿Hasta cuándo sigo?</b> La partición se hace en [a, b]: se sigue hasta xₙ = b. Sabiendo h y b se tienen todas las abscisas, <code>xᵢ = x₀ + i·h</code>, y la cantidad de subintervalos <code>n = (b − a)/h</code>.",
       "<b>¿Qué tamaño debe tener h?</b> n ≫ → mejor aproximación (h ≪); n ≫ → más cálculos (costo computacional; ≫ error de redondeo).",
       "<b>¿Qué tan buena es la aproximación?</b> Depende de la ESTABILIDAD (cómo los errores de los datos influyen en los resultados): del PVI y del método (§4)."),
    EJ("Ejemplo de la cátedra (EDO2 p17): [0, 1] con h = 0,2 → n = (1 − 0)/0,2 = 5; x₀ = 0, x₁ = 0,2, x₂ = 0,4, …, x₅ = 1."),
)

sec(4, "Estabilidad del PVI y estabilidad del método", 3, "EDO1 p15–17",
    P("Los datos siempre tienen errores (redondeo, medición). La pregunta es si esos errores se quedan chicos o crecen. La cátedra separa dos estabilidades: la del <b>problema</b> y la del <b>método</b>."),
    H4("a) Estabilidad del PVI — condición intrínseca del problema"),
    P("Gráficos de la p16: <b>EDO estable</b> y′ = −y (las curvas con distinta condición inicial se juntan), <b>EDO inestable</b> y′ = y (se separan), <b>EDO marginalmente estable</b> y′ = 1/2 (rectas paralelas: la distancia se mantiene)."),
    D("PVI bien planteado (bien condicionado)", "∃ una solución y(x) única para el problema; ∃ k, ξ &gt; 0 (ctes) tal que existe una única solución z(t) al PVI z′ = f(x, y) + a(x), z(x₀) = y₀ + b₀, tal que |z(x) − y(x)| &lt; kξ si |b₀| &lt; ξ ∧ |a(x)| &lt; ξ. Condición intrínseca del problema.", "EDO1 p16"),
    H4("b) Estabilidad del método — condición intrínseca del método"),
    D("Método numérico estable", "Un método numérico es estable si pequeñas perturbaciones en las condiciones iniciales producen cambios igualmente pequeños en las aproximaciones posteriores.", "EDO1 p17; EDO3 p4"),
    P("El diagrama de la p17: la <b>calidad de las aproximaciones</b> depende de tres cosas juntas: estabilidad del método, <b>tamaño del paso h</b> y estabilidad del PVI. Un método estable no arregla un PVI inestable, y un PVI estable puede dar basura con un método usado fuera de su rango de h (stiff, §18)."),
    TBL(["", "Estabilidad del PVI", "Estabilidad del método"], [
        ["De qué depende", "Del problema (f y la condición inicial)", "Del método (y, en la práctica, de h)"],
        ["Qué mide", "Si perturbar la condición inicial o f cambia poco la solución exacta", "Si perturbar la condición inicial cambia poco las aproximaciones"],
        ["Ejemplo de la cátedra", "y′ = −y estable; y′ = y inestable; y′ = 1/2 marginal", "Trapecio: “el más estable (A-estable)”; Euler explícito: “estable (aunque no el más estable)”"],
    ]),
    TRAMPA("IT25-T2 P5: <i>“Defina la estabilidad de un método para resolver EDOs. Mencione al menos dos métodos estables y explique.”</i> Se pide la del <b>método</b> (no la del PVI): definición textual de la p17 + dos métodos con su razón (desarrollo IT25-T2 P5)."),
)

sec(5, "Errores: truncamiento y redondeo; error local y error global", 3, "EDO2 p5–7",
    P("Como la solución es numérica, siempre hay error. Hay que saber de dónde viene (método o computadora), cómo se acumula paso a paso y cómo se mide en la tabla."),
    D("Error de truncamiento", "Asociado fundamentalmente al método, puede crecer si el método no es estable. Se pueden establecer estrategias como bajar el tamaño del paso (h), o cambiar por un método más estable.", "EDO2 p5"),
    D("Error de redondeo", "Asociado al uso del computador. Cuanto más cálculos se realicen, los errores de redondeo serán más significativos.", "EDO2 p5"),
    P("“La inestabilidad del método está asociada a la propagación del error. Al aproximar el valor yᵢ a partir de un único punto (x₀, y₀) se cometen errores de truncamiento. Pueden observarse dos tipos de errores” (p6):"),
    D("E. Local", "Se genera en cada paso, cada vez que se realiza un salto en la grilla. Es la diferencia que existe entre el valor calculado xₙ₊₁ y el valor de la solución en el punto anterior xₙ.", "EDO2 p6"),
    D("E. Global", "Teóricamente corresponde a la sumatoria de todos los errores locales. Tiene el sentido de analizar la diferencia entre la solución aproximada (yᵢ) y la verdadera (y(xᵢ)).", "EDO2 p6"),
    D("Resultado general", "Si el error de truncamiento local de un método es O(hⁿ⁺¹) entonces el error de truncamiento global es O(hⁿ).", "EDO2 p6"),
    H4("Cómo se llenan las columnas en las tablas de la cátedra (p7)"),
    TBL(["x", "y(x)", "yᵢ", "Elocal", "Eglobal"], [
        ["a = x₀", "y(x₀) = y₀", "", "", ""], ["x₁ = x₀ + h", "—", "y₁ = MÉTODO", "|y₁ − y₀|", "|y₁ − y(x₁)|"], ["x₂ = x₁ + h", "—", "y₂ = MÉTODO", "|y₂ − y₁|", "|y₂ − y(x₂)|"]]),
    P("Ésa es la convención que se usa en todos los ejemplos (Euler, AB2, Heun) y en el TP4: <b>E local = |yᵢ₊₁ − yᵢ|</b> (salto entre aproximaciones consecutivas) y <b>E global = |yᵢ − y(xᵢ)|</b> (contra la solución exacta). El EG “no siempre se puede calcular” (EDO2 p21): hace falta la solución exacta."),
    AGR("En la bibliografía el <i>error de truncamiento local</i> es otra cosa: el error que comete <b>un</b> paso arrancando desde el valor exacto, τ = y(xₙ₊₁) − [y(xₙ) + h·F(…)]; es el que se usa para el orden (O(hᵖ⁺¹)) y para la consistencia (§13). La columna “E local” de las tablas de la cátedra (|yᵢ₊₁ − yᵢ|) mide cuánto avanzó la aproximación, no ese τ; por eso en sus tablas el EL es grande (0,3; 0,35…) aunque el método sea bueno. En el TP4 usá la convención de la tabla; en teoría (orden, consistencia), la de Taylor."),
    TRAMPA("IT25-T3 P5 y IT25-T4 P5 piden la relación de consistencia y convergencia con “los errores”: consistencia ↔ <b>error local</b>; convergencia ↔ <b>error global</b> (EDO3 p3). Mezclarlas es el error típico."),
)

sec(6, "Clasificación de métodos: explícito/implícito, un paso/multipaso", 3, "EDO2 p8",
    P("Las preguntas 4 de los cuatro temas del integral 2025 dan una fórmula y piden clasificarla. La clasificación sale de mirar <b>qué valores aparecen del lado derecho</b>."),
    TBL(["Métodos Explícitos", "Métodos Implícitos"], [["yₙ₊₁ = MÉTODO(xₙ, yₙ)", "yₙ₊₁ = MÉTODO(xₙ, yₙ, xₙ₊₁, yₙ₊₁)"]]),
    TBL(["Métodos de un paso", "Métodos multipaso"], [["Solo utilizan información del punto anterior xₙ₋₁ para obtener la aproximación de y(xₙ)", "Usan más de un punto anterior para estimar el próximo valor de la función (necesitan xₙ₋₁, xₙ₋₂, …, etc)"]]),
    H4("Receta para los 4 incisos (IT25 P4)"),
    OL("<b>¿Es Runge-Kutta?</b> Sí si se puede escribir como <code>yₙ₊₁ = yₙ + h·F(h, xₙ, yₙ, f)</code> (expresión general RK, EDO2 p11), con F un promedio de pendientes evaluadas a partir de (xₙ, yₙ) únicamente. Los RK son “de un solo paso” (EDO2 p12): si la fórmula usa yₙ₋₁, yₙ₋₂… <b>no</b> es RK.",
       "<b>¿De qué orden?</b> Desarrollar por Taylor: el orden p es la potencia hasta la que el método “coincide con la ST” (EDO2 p12); error local O(hᵖ⁺¹) → global O(hᵖ) (resultado general). Si el desarrollo no coincide ni en el término h·y′, el método <b>no es consistente</b> (§13).",
       "<b>¿Un paso o multipaso?</b> Contar los puntos anteriores: sólo (xₙ, yₙ) → un paso; yₙ₋₁, yₙ₋₂… → multipaso de k pasos (y entonces “no se autoinicia”, EDO3 p10).",
       "<b>¿Implícito o explícito?</b> Si a la derecha aparece yₙ₊₁ (o y′ₙ₊₁ = f(xₙ₊₁, yₙ₊₁), o una k que depende de sí misma) → implícito; si sólo aparecen valores ya conocidos → explícito."),
    TRAMPA("Mecanismo de distractor de F18-P1: atribuirle a un método el rasgo de otra familia (“es de un solo paso”, “es implícito”, “pertenece a un P-C” para AB3, que es multipaso explícito). En cada inciso escribí el rasgo <b>y la razón</b> (qué valores usa la fórmula): la nota parcial se pierde por no justificar."),
)

sec(7, "Serie de Taylor como método", 1, "EDO2 p9–10",
    P("Es el punto de partida de todos los métodos: si conociera todas las derivadas de y en xᵢ, podría avanzar con su desarrollo."),
    CMD("xᵢ = x₀ + i·h\nyᵢ₊₁ = y(xᵢ) + h·y′(xᵢ) + h²·y″(xᵢ)/2 + h³·y‴(xᵢ)/3! + … + hᵏ⁺¹·y⁽ᵏ⁺¹⁾(ξ)/(k+1)!,   ξ ∈ (xᵢ, xᵢ₊₁)\n                                                    └── error al truncar la serie"),
    P("Observaciones de la cátedra:"),
    UL("Método de <b>un solo paso</b>: para estimar yᵢ₊₁ sólo necesito (xᵢ, yᵢ).",
       "Necesito conocer las <b>derivadas</b> de la función y.",
       "Es una serie infinita que en la práctica será truncada (e. truncamiento del orden del término truncado)."),
    P("<b>Ventaja:</b> “Es el más estable de todos los métodos”. <b>Desventaja:</b> “Necesita poder calcular derivadas (costoso), por lo tanto no se utiliza en la práctica”."),
)

sec(8, "Métodos de Runge-Kutta: expresión general", 3, "EDO2 p11–12",
    P("Los RK resuelven la desventaja de Taylor: en lugar de calcular derivadas de y, combinan evaluaciones de f (que es dato) para reproducir la serie hasta cierto término."),
    D("Expresión general de los métodos RK", "yₙ₊₁ = yₙ + h·F(h, xₙ, yₙ, f), donde h es el paso y F la función de incremento.", "EDO2 p11"),
    P("Observaciones (p12):"),
    UL("Método de <b>un solo paso</b>.", "Modificación de la serie de Taylor.",
       "Aproximan el valor de la derivada mediante una función de incremento <code>F(h, xₙ, yₙ, f) ≈ f(xₙ, yₙ)</code>.",
       "<b>Coinciden con la ST hasta el término de hᵖ, donde p es el orden del método.</b>",
       "Error de truncamiento global es O(hᵖ)."),
    CMD("yᵢ₊₁ = y(xᵢ) + h·[ y′(xᵢ) + h·y″(xᵢ)/2 + h²·y‴(xᵢ)/3! + … ] + hᵏ⁺¹·y⁽ᵏ⁺¹⁾(ξ)/(k+1)!\n             └────────── F aproxima este corchete ──────────┘"),
    P("Conjunto de métodos desarrollado alrededor del año 1900 por los matemáticos alemanes C. Runge y M. W. Kutta."),
    TRAMPA("IT25-T1 P4 pregunta “¿Es un método Runge-Kutta? ¿por qué?”. La justificación que espera la cátedra es la forma <code>yₙ₊₁ = yₙ + h·F(h, xₙ, yₙ, f)</code> con F = combinación de evaluaciones de f, y que sea de un solo paso."),
)

eu = R["cat_euler"]
sec(9, "Método de Euler explícito (RK1)", 3, "EDO2 p13–21",
    P("Es el método más simple: avanzar por la recta tangente. Sirve de base (y de método de arranque) para todo lo demás."),
    CMD("Δy/Δx ≈ f(x, y)   →   (yₙ₊₁ − yₙ)/(xₙ₊₁ − xₙ) ≈ f(xₙ, yₙ),  con xₙ₊₁ − xₙ = h\n\nyₙ₊₁ = yₙ + h·f(xₙ, yₙ)     n = 0, 1, …\nVal. Nuevo = Val. Anterior + Paso · Pendiente"),
    P("<b>Función de incremento:</b> F = f(xₙ, yₙ). <b>Orden:</b> coincide con 2 términos de la ST (y(xₙ) + h·y′(xₙ)); el término del error es h²·y″/2 → <b>error trunc local O(h²), orden del método O(h)</b> (p15)."),
    H4("Observaciones, ventajas y desventajas (p16)"),
    UL("De un solo paso.", "De orden p = 1 (E. truncamiento = p + 1 = 2).", "Coincide con dos términos de la ST.",
       "Aproxima la y′ con f(x, y) → ventaja sobre ST.", "La función de incremento F = f(x, y).",
       "Es explícito porque para calcular yᵢ₊₁ sólo necesita el valor de yᵢ."),
    P("<b>Ventajas:</b> método estable (aunque no es el más estable de todos); fácil de implementar. <b>Desventajas:</b> orden de convergencia lineal O(h), p = 1 → E. truncamiento muy alto; el error crece cuando crece el paso, para una buena aproximación se necesita un paso muy chico."),
    H4("Ejemplo de la cátedra (p17–20): y′ = y − t² + 1, 0 ≤ t ≤ 1, y(0) = 0,5, h = 0,2, exacta y(t) = (t + 1)² − 0,5eᵗ"),
    P(f"Primer paso: f(0; 0,5) = 0,5 − 0 + 1 = 1,5 → y₁ = 0,5 + 0,2·1,5 = {c(eu[1][2],1)}. EL = |0,8 − 0,5| = {c(eu[1][3],1)}; EG = |{c(eu[1][1],8)} − 0,8| = {c(eu[1][4],8)}."),
    tabla_err(eu, 8, aprox=lambda v: c(v, 6), err=lambda v: c(v, 8)),
    P("(Tabla recalculada en <code>verif_u4.py</code>: coincide con la p20. La p17 tiene la columna y(real) equivocada y el título dice “y(2)” aunque el intervalo es [0, 1].)"),
    P("Preguntas de la p21: ¿es un buen resultado? (el EG crece paso a paso hasta " + c(eu[-1][4], 4) + " en t = 1); ¿qué me dice el EG? (qué tan lejos estoy de la solución, sólo si la conozco); ¿qué me dice el EL? (cuánto cambia la aproximación en cada salto); ¿cómo impacta h? (h menor → menor error, más cálculos); ¿cómo impacta el orden? (orden mayor → error menor con el mismo h: ver RK2, RK4, AB2)."),
    AGR(f"Orden comprobado en <code>verif_u4.py</code>: por Taylor el error local es (h²/2)·y″ (orden 1) y, al dividir h por 2, el error global en t = 1 se divide por ≈ 2 (órdenes empíricos {', '.join(c(o,2) for o in EMP['Euler explícito']['ords'])})."),
)

ei = R["cat_eimpl"]
sec(10, "Método de Euler implícito (RK1)", 3, "EDO2 p22–25",
    P("Es Euler pero tomando la pendiente en el punto de llegada. Cuesta más (hay que despejar yₙ₊₁), y a cambio es más estable."),
    CMD("yₙ₊₁ = yₙ + h·f(xₙ₊₁, yₙ₊₁)     n = 0, 1, …"),
    P("<b>Función de incremento:</b> comparando con la expresión general RK, <code>F_I = f(xₙ₊₁, yₙ₊₁)</code>. <b>Orden:</b> “el orden viene dado por la cantidad de aproximaciones de derivadas de la serie de Taylor. En este caso el orden es 1”. <b>¿Cómo calculamos yₙ₊₁?</b> Despejando (aparece de los dos lados)."),
    H4("Ejemplo de la cátedra: y′ = y − t² + 1, y(0) = 0,5, h = 0,2"),
    CMD("yₙ₊₁ = yₙ + h(yₙ₊₁ − t²ₙ₊₁ + 1)\nyₙ₊₁ − h·yₙ₊₁ = yₙ − h·t²ₙ₊₁ + h\nyₙ₊₁ = (yₙ − h·t²ₙ₊₁ + h) / (1 − h)"),
    P(f"y₁ = (0,5 − 0,2·0,2² + 0,2)/(1 − 0,2) = {c(ei[1][2],3)} (la diapositiva escribe 0,864999); EG = {c(ei[1][4],4)}; EL = {c(ei[1][3],3)}. La diapositiva dice “COMPLETAR LA TABLA!”: está completa en el ejercicio de práctica correspondiente."),
    P("Preguntas de la cátedra: <b>¿siempre se podrá “explicitar” la ecuación?</b> No: acá se pudo porque f es lineal en y; con f no lineal (p. ej. y² o e^y) no hay despeje algebraico. <b>¿De qué otra manera se podría usar?</b> Como corrector, iterando a partir de una predicción explícita (eso es el esquema predictor-corrector, §17). <b>¿Cómo es la aproximación comparada con Euler explícito?</b> Del mismo orden; en este ejemplo el EG final es " + c(ei[-1][4], 4) + " contra " + c(eu[-1][4], 4) + " del explícito (uno se pasa por arriba y el otro queda por abajo)."),
    AGR("La ventaja del implícito aparece en problemas stiff: aplicado a y′ = λy con λ &lt; 0, el factor por paso es 1/(1 − hλ), que tiene módulo &lt; 1 para cualquier h &gt; 0; el explícito tiene factor 1 + hλ y sólo es estable si |1 + hλ| ≤ 1. Es lo que la cátedra resume como “los métodos implícitos generalmente son más estables (región de estabilidad más amplia)” (EDO3 p4)."),
)

sec(11, "RK2: familia, Trapecio / Euler modificado y punto medio", 3, "EDO2 p26–28",
    P("Con dos evaluaciones de f por paso se puede coincidir con la serie de Taylor hasta h²: orden 2. Distintas elecciones de las constantes dan distintos métodos de la familia."),
    D("Métodos Runge-Kutta 2 (RK2)", "Es una familia de métodos Runge-Kutta de orden 2 cuya función de incremento es: F(h, xₙ, yₙ, f) = γ₁f(x, y) + γ₂f(x + αh, y + βh·f(x, y)). Asignando diferentes valores a las constantes se obtienen los diferentes métodos.", "EDO2 p26"),
    H4("Método del Trapecio o Euler Modificado (p27)"),
    CMD("yₙ₊₁ = yₙ + ½·h·( f(xₙ, yₙ) + f(xₙ₊₁, yₙ + h·f(xₙ, yₙ)) )\n\nOtra manera:  k₁ = f(xₙ, yₙ)\n              k₂ = f(xₙ₊₁, yₙ + h·f(xₙ, yₙ))\n              yₙ₊₁ = yₙ + ½·h·(k₁ + k₂)"),
    UL("“Coinciden 3 términos con la ST” → “Error trunc local O(h³). Orden método es 2 y el error de trunc global es O(h²)”.",
       "La función de incremento es F(h, x, y, f) = ½(f(x, y) + f(x + h, y + h·f(x, y))) → γ₁ = 1 − γ₂, γ₂ = ½ = γ₁, α = 1, β = 2 (así en la diapositiva).",
       "“Es el método más estable (A-estable)”; la diapositiva lo rotula “(IMPLÍCITO)”."),
    AGR(f"Dos aclaraciones verificadas en <code>verif_u4.py</code>. (1) <b>β = 2 es una errata</b>: con la fórmula escrita, y + h·f(x, y) corresponde a <b>β = 1</b>; con β = 2 el método baja a orden {RK['RK2 con β=2 (errata de la diapositiva)']['butcher']}. (2) La fórmula de la p27 es <b>explícita</b> (k₂ usa la predicción de Euler, no yₙ₊₁): es el “Euler modificado”, que coincide con un Heun de una sola corrección. El <b>Trapecio implícito</b> es <code>yₙ₊₁ = yₙ + (h/2)(f(xₙ, yₙ) + f(xₙ₊₁, yₙ₊₁))</code>, que la cátedra usa como corrector en Heun (EDO4 p4); ése es el A-estable. Los dos son de orden 2 (Taylor: error local h³ en ambos). En un examen: si te piden “Trapecio, implícito, A-estable”, escribí la fórmula con yₙ₊₁ adentro."),
    H4("Método del punto medio (p28)"),
    CMD("yₙ₊₁ = yₙ + h·f( xₙ + h/2,  yₙ + (h/2)·f(xₙ, yₙ) )     n = 0, 1, …"),
    P("Las preguntas de la cátedra, respondidas (verificado por Taylor y por órdenes empíricos en <code>verif_u4.py</code>):"),
    UL(f"¿De qué orden es? <b>2</b> (error local h³; órdenes empíricos {', '.join(c(o,2) for o in EMP['Punto medio']['ords'])}).",
       "¿De qué orden es el error de truncamiento global? <b>O(h²)</b>.",
       "¿Implícito o explícito? <b>Explícito</b>: la pendiente en el medio se evalúa con una predicción de Euler a medio paso, todo con valores conocidos.",
       "¿Constantes? <b>γ₁ = 0, γ₂ = 1, α = ½, β = ½</b>."),
    TBL(["Método", "γ₁", "γ₂", "α", "β", "Tipo"], [
        ["Euler modificado (fórmula de la p27)", "½", "½", "1", "1 (diapositiva: 2)", "explícito, orden 2"],
        ["Punto medio", "0", "1", "½", "½", "explícito, orden 2"]]),
)

sec(12, "RK4 y la relación entre evaluaciones y orden", 3, "EDO2 p29–30",
    P("Es el RK más usado: cuatro evaluaciones de f por paso para orden 4."),
    CMD("yₙ₊₁ = yₙ + (h/6)(k₁ + 2k₂ + 2k₃ + k₄)\nk₁ = f(xₙ, yₙ)\nk₂ = f(xₙ + h/2, yₙ + h·k₁/2)\nk₃ = f(xₙ + h/2, yₙ + h·k₂/2)\nk₄ = f(xₙ + h, yₙ + h·k₃)"),
    P(f"Preguntas de la cátedra: <b>¿de qué orden es? ¿qué significa?</b> Orden 4: coincide con la serie de Taylor hasta h⁴ (error local O(h⁵)); al dividir h por 2 el error global se divide por ≈ 2⁴ = 16 (empírico en <code>verif_u4.py</code>: órdenes {', '.join(c(o,2) for o in EMP['RK4']['ords'])}). <b>¿Error global?</b> O(h⁴). <b>¿Implícito o explícito?</b> Explícito: cada k usa sólo las anteriores."),
    H4("Métodos Runge-Kutta: costo vs orden (p30)"),
    UL("Si aumenta el orden del método, h no necesariamente debe ser muy chico.",
       "Si el orden aumenta, aumenta el número de evaluaciones. Esto incrementa el costo (mayor esfuerzo de cálculo)."),
    TBL(["#Eval", "2", "3", "4", "5", "6", "7", "8 o +"], [["Mejor ELoc posible", "O(h²)", "O(h³)", "O(h⁴)", "O(h⁴)", "O(h⁵)", "O(h⁶)", "O(hⁿ⁻²)"]]),
    P("4 evaluaciones = “+ difundido”; 5 evaluaciones = “&gt; trabajo computacional que beneficio respecto al error” (con 5 no se gana orden respecto de 4)."),
    AGR("En esa tabla lo que la cátedra llama “ELoc” es el orden del error <b>global</b> (con 4 evaluaciones, RK4 tiene error global O(h⁴) y local O(h⁵)). Para escribir: RK4 → error local O(h⁵), global O(h⁴), orden 4."),
)

sec(13, "Evaluación de métodos: consistencia, convergencia y estabilidad", 3, "EDO3 p2–4",
    P("Saber el orden no alcanza: hay que saber si el método aproxima la EDO correcta (consistencia), si la solución numérica llega a la verdadera al achicar h (convergencia) y si los errores no se disparan (estabilidad). Las preguntas 5 del integral 2025 salen de esta diapositiva."),
    P("El repaso de la p2 ordena lo visto: métodos de un solo paso (Serie de Taylor; RK1 Euler explícito e implícito; RK2 Trapecio o Euler modificado y punto medio; RK4) y estabilidad (del PVI y del método)."),
    D("Consistencia", "Un método de diferencias con error de truncamiento local Tn en el n-ésimo paso es consistente con la EDO que aproxima si cumple que lim_{h→0} max_{1≤i≤n} |Tn| = 0. Obs: depende de h; relaciona el ELocal con el método.", "EDO3 p3"),
    D("Convergencia", "Un método de diferencias es convergente con respecto a la EDO que aproxima si verifica que lim_{h→0} max_{1≤i≤n} |yᵢ − y(xᵢ)| = 0. Obs: no depende de h; para cualquier h, la sol numérica va a tender a la real; relaciona el EGlobal con el método.", "EDO3 p3"),
    D("Estabilidad", "Un método numérico es estable si pequeñas perturbaciones en las condiciones iniciales producen cambios igualmente pequeños en las aproximaciones posteriores.", "EDO3 p4"),
    P("Observaciones de la cátedra (p4):"),
    UL("<b>Si un método es consistente y estable entonces es convergente.</b>",
       "En gral, si el problema es bien condicionado entonces los métodos RK son estables aunque pueden ser lentos.",
       "Los métodos implícitos generalmente son más estables (región de estabilidad más amplia)."),
    TBL(["", "Consistencia", "Convergencia", "Estabilidad"], [
        ["Qué pide", "El error de truncamiento local tiende a 0 con h → 0", "El error global máximo tiende a 0 con h → 0", "Perturbaciones chicas → cambios chicos"],
        ["Error con el que se relaciona", "Error <b>local</b>", "Error <b>global</b>", "Propagación del error (EDO2 p6)"],
        ["Depende de h (según la cátedra)", "Sí", "No", "—"],
        ["Cómo se verifica", "Taylor: el método coincide con la ST al menos hasta h·y′ (orden p ≥ 1)", "Consistente + estable ⇒ convergente; o comparando con la exacta al achicar h", "Perturbar y₀ / mirar si el error crece"],
    ]),
    AGR(f"(1) La cátedra enuncia sólo “consistente + estable ⇒ convergente”. El recíproco completo (teorema de equivalencia de Dahlquist, para métodos multipaso lineales: convergente ⇔ consistente y cero-estable) no está en el material. Ojo: esa “estabilidad” es para h → 0; un método convergente igual puede ser inestable con un h fijo grande (Euler en el ejemplo stiff de la cátedra, §18). (2) Cómo verificar consistencia por Taylor: el método debe reproducir y(xₙ₊₁) = y(xₙ) + h·y′(xₙ) + O(h²); el primer término que no coincide da el orden. Ejemplo (verif_u4): <code>yₙ₊₁ = yₙ₋₁ + h·y′ₙ</code> deja un residuo {M['IT25-T4 (h, como dice el enunciado)']['c0']}·h, es decir h·y′: <b>no es consistente</b> (aproxima y′ = f/2)."),
    TRAMPA("IT25-T3 P5: <i>“Defina consistencia de un método para resolver EDOs. Explicar su relación con los errores.”</i> IT25-T4 P5: <i>“Defina convergencia… Explique cómo se relaciona con los errores.”</i> F20-P2: <i>“Definir consistencia…”</i>. WA25 V/F: <i>“método convergente → es estable”</i>, <i>“método no consistente para PVI, se recomienda cambiar método”</i>. Respuestas completas en Desarrollo."),
)

sec(14, "Elección del tamaño de paso h", 2, "EDO3 p5–6",
    P("h controla a la vez precisión, costo y estabilidad, y no se pueden maximizar las tres."),
    TBL(["Si…", "Ventaja", "Desventaja"], [
        ["h ≪", "&gt; precisión", "&gt; costo computacional; &lt; estabilidad"],
        ["h ≫", "&gt; estabilidad; &lt; costo computacional", "&lt; precisión"]]),
    P("“<b>h debe ser suficientemente grande para minimizar el costo computacional conservando estabilidad y precisión</b>” (p5)."),
    P("Ejemplo de la p6 (señal sismológica: un tramo tranquilo y después la onda P y la onda S con mucha variabilidad): ¿orden alto o bajo?, ¿implícito o explícito?, ¿qué h? Respuesta de la diapositiva: <b>≪ variabilidad → h ≫; ≫ variabilidad → h ≪</b>. Un h fijo obliga a usar el h chico en todo el intervalo; de ahí la motivación del paso variable (§15)."),
    AGR("El “h ≪ → &lt; estabilidad” de la cátedra se entiende por el redondeo: más pasos → más operaciones → más error de redondeo acumulado (EDO1 p14, EDO2 p5). En problemas stiff pasa lo contrario: con h grande el método explícito se vuelve inestable (§18). Por eso la frase final de la diapositiva pide un h “suficientemente grande” pero conservando estabilidad."),
)

sec(15, "Métodos de paso variable: Runge-Kutta-Fehlberg (RKF45)", 3, "EDO3 p7–9",
    P("Con h fijo se paga el paso chico en todo el intervalo. El paso variable achica h sólo donde hace falta y lo agranda donde la solución es tranquila."),
    P("<b>OBJETIVO</b> (p7): ajustar el tamaño del paso (h) automáticamente (<b>paso adaptativo</b>), de modo que el error de truncamiento esté acotado (no supere un valor aceptable). Pero no conozco el EG porque no conozco la solución → “<b>Aproximo el EG con el EL mediante el uso de dos métodos de RK de diferentes órdenes en cada punto</b> para determinar cuál es el valor aceptable de h (que controle el EG). Así encuentro la cantidad y posición de los nodos para usar en la aproximación que garantizan el E. de truncamiento acotado”."),
    H4("RKF45 (p8; la diapositiva escribe “Felberhg”)"),
    CMD("RK4:  yₙ₊₁ = yₙ + 23/216·k₁ + 1408/2565·k₃ + 219/4104·k₄ − 1/5·k₅        (así en la diapositiva)\nRK5:  ȳₙ₊₁ = yₙ + 16/135·k₁ + 6656/12825·k₃ + 28561/56430·k₄ − 9/50·k₅ + 2/55·k₆\n\nk₁ = h·f(t, y)\nk₂ = h·f(t + h/4, y + k₁/4)\nk₃ = h·f(t + 3h/8, y + 3k₁/32 + 9k₂/32)\nk₄ = h·f(t + 12h/13, y + 1932k₁/2197 − 7200k₂/2197 + 7296k₃/2197)\nk₅ = h·f(t + h, y + 439k₁/216 − 8k₂ + 3680k₃/513 − 845k₄/4104)\nk₆ = h·f(t + h/2, y − 8k₁/27 + 2k₂ − 3544k₃/2565 + 1859k₄/4104 − 11k₅/40)"),
    UL("“Para cada cálculo de yᵢ₊₁ se usan dos aproximaciones diferentes. Si la precisión esperada no se alcanza con estos métodos, se reduce el tamaño del paso y se vuelve a realizar el cálculo. Si por el contrario la precisión en ambas aproximaciones es mejor que la esperada, el tamaño del paso se incrementa.”",
       "“Se usan los métodos de RK4 y RK5 (ambos métodos tienen Elocal O(h⁴)). Se calcula el EL en ambos casos y a partir de ahí tomar la decisión de agrandar o achicar h.”"),
    AGR(f"Verificado en <code>verif_u4.py</code> con las condiciones de orden: (1) los coeficientes de la fórmula “RK4” de la diapositiva tienen una <b>errata</b>: son <b>25/216</b> y <b>2197/4104</b> (con 23/216 y 219/4104 los pesos ni siquiera suman 1: orden {RK[list(RK)[-1]]['butcher']}). (2) Con los coeficientes correctos, la fórmula de 4º orden da orden {RK[list(RK)[-3]]['butcher']} y la de 5º, orden {RK[list(RK)[-2]]['butcher']}: el nombre 45 viene de ahí. Las seis k son las mismas para las dos fórmulas: por eso son 6 evaluaciones por paso y no 4 + 6. RKF45 es <b>de un paso, explícito</b> y <b>no</b> es predictor-corrector (no itera un corrector implícito: compara dos fórmulas explícitas)."),
    H4("Ejemplo (p9): y′ = 1 + y², y(0) = 0, exacta tan(t)"),
    TBL(["", "RKF45 (h variable)", "RK4 (h = 0,1 constante)"], [
        ["Nodos usados hasta t = 1,4", "11 (t = 0; 0,2; …; 1,2; 1,3; 1,35; 1,4: el paso se achica cerca de la asíntota de tan t)", "15"],
        ["Error en t = 1,4", "−0,0006208", "0,0059089"]]),
    P("Conclusiones de la diapositiva: “El error de truncamiento en ambos métodos es de orden 4”; “RKF45 ≫ cálculos y evaluaciones de la función”."),
    TRAMPA("WA25: <i>“RKF: objetivo, en qué consiste, cuándo se usa, por qué. Orden, ¿multipaso? ¿implícito? ¿P-C? ¿Con qué otro método lo comparo?”</i> F18-P2 lo pide en un cuadro junto a Trapecio, RK4 y Heun. Respuesta completa en Desarrollo."),
)

ab = R["cat_ab2"]
sec(16, "Métodos multipaso: Adams-Bashforth y Adams-Moulton", 3, "EDO3 p10–18",
    P("Los RK tiran la información de los pasos anteriores y evalúan f varias veces por paso. Los multipaso reutilizan los valores ya calculados: más orden con menos evaluaciones nuevas, a cambio de no poder arrancar solos."),
    H4("Características (p10)"),
    UL("Utilizan más de un paso anterior para el cálculo de la aproximación siguiente.",
       "Ej.: un método multipaso de 2 pasos explícito tiene la forma <code>yᵢ₊₁ = MÉTODO(h, xᵢ₋₁, yᵢ₋₁, xᵢ, yᵢ)</code>.",
       "Se pueden ver los métodos de un solo paso como un caso particular de los multipasos.",
       "Pueden ser explícitos o implícitos.", "<b>No se autoinician!</b>"),
    H4("Métodos de Adams (p11)"),
    CMD("yₙ₊₁ = yₙ + ∫ₓₙ^ₓₙ₊₁ P(x) dx        (se integra un polinomio interpolante que pasa por los puntos anteriores)\n\nyₙ₊₁ = yₙ + h·Σ bⱼ·f(xₙ₋ⱼ, yₙ₋ⱼ)       forma general de los métodos Adams\n\nExplícitos  b₋₁ = 0  → Adams-Bashforth        Implícitos  b₋₁ ≠ 0  → Adams-Moulton"),
    P("(b₋₁ es el coeficiente de f(xₙ₊₁, yₙ₊₁): si está, el método es implícito.)"),
    H4("Adams-Bashforth (p12, p17)"),
    CMD("AB2 (2 pasos):  yₙ₊₁ = yₙ + (h/2)[3f(xₙ, yₙ) − f(xₙ₋₁, yₙ₋₁)]\nAB3 (3 pasos):  yₙ₊₁ = yₙ + (h/12)[23f(xₙ, yₙ) − 16f(xₙ₋₁, yₙ₋₁) + 5f(xₙ₋₂, yₙ₋₂)]\nAB4 (4 pasos):  yₙ₊₁ = yₙ + (h/24)[55f(xₙ, yₙ) − 59f(xₙ₋₁, yₙ₋₁) + 37f(xₙ₋₂, yₙ₋₂) − 9f(xₙ₋₃, yₙ₋₃)]"),
    P("AB2 según la cátedra: necesita dos puntos anteriores (método de 2 pasos); dado el PVI sólo tengo (x₀, y₀) y necesito (x₁, y₁) para iniciar AB2 → <b>(x₁, y₁): uso un método de un orden menor (x ejemplo Euler)</b>; error de truncamiento de orden O(h²); “&gt; número de evaluaciones de la función f”."),
    AGR(f"Órdenes verificados por Taylor (<code>verif_u4.py</code>): AB2 error local {TY('AB2')} → orden {M['AB2']['p']}; AB3 → orden {M['AB3 (F18-P1)']['p']}; AB4 → orden {M['AB4']['p']}. En general ABk es de orden k. Sobre las evaluaciones: guardando f(xₙ₋₁, yₙ₋₁) del paso anterior, AB2 hace <b>una</b> evaluación nueva de f por paso (RK2 hace dos); en el ejemplo de la cátedra AB2 evalúa f dos veces en cada fila porque recalcula la anterior."),
    H4("Ejemplo de la cátedra (p13–15): y′ = y − x² + 1, y(0) = 0,5, h = 0,2"),
    P(f"y₁ = y₀ + h·f(x₀, y₀) = {c(ab[1][2],1)} (<b>método de inicio</b>: Euler). y₂ = y₁ + (h/2)[3f(x₁, y₁) − f(x₀, y₀)] = 0,8 + 0,1·(3·1,76 − 1,5) = {c(ab[2][2],3)}; y₃, y₄, y₅ con AB2."),
    tabla_err(ab, 8, aprox=lambda v: c(v, 6), err=lambda v: c(v, 6)),
    P(f"(La diapositiva deja vacías las columnas EL y EG; se completaron con la convención de la cátedra. Los yₐₚᵣₒₓ coinciden con la p15.) Comparando con Euler (p16): mismo arranque, pero EG final {c(ab[-1][4],4)} contra {c(eu[-1][4],4)} de Euler: AB2 es de orden 2. En común: los dos son explícitos y usan f de puntos ya calculados."),
    H4("Adams-Moulton (p18)"),
    CMD("AM2 (“de 2 pasos”):  yₙ₊₁ = yₙ + (h/12)[5f(xₙ₊₁, yₙ₊₁) + 8f(xₙ, yₙ) − f(xₙ₋₁, yₙ₋₁)]\nAM3 (“de 3 pasos”):  yₙ₊₁ = yₙ + (h/24)[9f(xₙ₊₁, yₙ₊₁) + 19f(xₙ, yₙ) − 5f(xₙ₋₁, yₙ₋₁) + f(xₙ₋₂, yₙ₋₂)]"),
    UL("“Los métodos implícitos (AM) tienen el problema de que primero se debe convertir algebraicamente a una manera explícita, aunque en general dan mejores resultados que los métodos explícitos (AB).”",
       "“En la práctica se utilizan para mejorar las aproximaciones obtenidas con otro método y dan el nombre a un nuevo tipo de métodos denominados predictor-corrector.”"),
    AGR(f"Taylor (<code>verif_u4.py</code>): el que la cátedra llama AM2 es de orden {M['AM2 (cátedra: \'de 2 pasos\')']['p']} y el AM3, de orden {M['AM3 (cátedra: \'de 3 pasos\')']['p']}. En EDO4 p4 la misma fórmula (h/24)(9, 19, −5, 1) aparece como “AM4”: la numeración cambia entre clases (por cantidad de puntos o por orden). En el examen escribí la fórmula y aclarás su orden."),
    TRAMPA("F18-P1: <i>“Y_{n+1} = Y_n + (1/12)H[23F(n) − 16F(n−1) + 5F(n−2)]. Responder V o F, justificar: a) Es un método de un solo paso. b) Es un método implícito. c) Pertenece a un esquema Predictor-Corrector.”</i> Es AB3: los tres incisos son falsos tal como están escritos (ver Desarrollo). WA25: <i>“Los métodos multipaso explícitos se basan en…”</i> → en integrar un polinomio interpolante que pasa por los puntos anteriores (Adams-Bashforth)."),
)

hn = R["cat_heun"]; tr = R["cat_heun_traza"]
sec(17, "Métodos predictor-corrector: Heun, Milne, AB4+AM4", 3, "EDO4 p2–8",
    P("Los implícitos son más estables pero hay que despejar yₙ₊₁. El predictor-corrector evita el despeje: predice con un explícito y corrige iterando el implícito."),
    P("<b>OBJETIVO</b> (p2): “Mejorar una aproximación inicial mediante correcciones sucesivas de acuerdo al error.”"),
    UL("Se transforma en un método <b>iterativo</b>.",
       "Trabaja en dos etapas: 1) <b>Predicción</b> de la aproximación (con un método <b>explícito</b>, en gral de orden n) y 2) <b>Corrección</b> del valor predicho usando un método <b>implícito</b> de orden n+1.",
       "<b>El orden del método es el orden del corrector.</b>", "&gt; costo."),
    CMD("y₁ = MetExplícito                         ← predictor\nMIENTRAS |error| &gt; E && it &lt; itMax         ← corrector\n    Yant ← y₁\n    y₁ ← MetImplícito(usa y₁ en f(x₁, y₁))\n    error = |Yant − y₁|\ny₂ = MetExplícito …"),
    H4("Esquemas de la p4"),
    CMD("HEUN    P: Euler explícito   yₙ₊₁ = yₙ + h·f(xₙ, yₙ)\n        C: Trapecio          yₙ₊₁ = yₙ + (h/2)(f(xₙ, yₙ) + f(xₙ₊₁, yₙ₊₁))\n\nMILNE   P:  yₙ₊₁ = yₙ₋₃ + (4/3)h(2fₙ − fₙ₋₁ + 2fₙ₋₂)\n        C:  yₙ₊₁ = yₙ₋₁ + (h/3)(fₙ₊₁ + 4fₙ + fₙ₋₁)\n\nAB4+AM4 P:  yₙ₊₁ = yₙ + (h/24)(55fₙ − 59fₙ₋₁ + 37fₙ₋₂ − 9fₙ₋₃)\n        C:  yₙ₊₁ = yₙ + (h/24)(9fₙ₊₁ + 19fₙ − 5fₙ₋₁ + fₙ₋₂)"),
    P("Preguntas de la cátedra sobre Heun, respondidas:"),
    UL("<b>¿Es un método RK? ¿por qué?</b> Con una sola corrección es exactamente el RK2 Euler modificado (yₙ₊₁ = yₙ + ½h(k₁ + k₂)): tiene la forma yₙ₊₁ = yₙ + h·F y es de un paso. Iterado como P-C deja de ser un RK “puro” (la cantidad de evaluaciones depende de la tolerancia), pero sigue siendo de un paso.",
       "<b>¿De qué orden es?</b> 2: el orden del corrector (Trapecio).",
       "<b>¿Cuántas veces se evalúa la función por punto de la grilla?</b> 1 para el predictor (f(xₙ, yₙ), que el corrector reutiliza) + 1 por cada corrección (f(xₙ₊₁, y₁)): en el ejemplo, 2 correcciones → 3 evaluaciones.",
       f"<b>¿De qué orden es Milne?</b> Orden 4 (predictor y corrector de orden {M['Milne predictor']['p']} y {M['Milne corrector']['p']} por Taylor, <code>verif_u4.py</code>). Es multipaso (4 puntos anteriores en el predictor).",
       "<b>¿Cuándo usaría P-C?</b> Cuando quiero la estabilidad/precisión de un implícito sin despejarlo (f no lineal en y) y puedo controlar el error con una tolerancia.",
       "<b>¿Por qué el predictor es como mínimo de orden n si el corrector es n+1?</b> Porque el corrector parte de la predicción: si es muy mala hacen falta muchas iteraciones (más costo) o no converge; con orden n alcanza con pocas correcciones para llegar al orden n+1."),
    H4("Ejemplo de Heun (p5–7): y′ = y − x² + 1, y(0) = 0,5, h = 0,2, EL &lt; 10⁻², ItMax = 5"),
    P(f"Primer punto: predictor y₁ = 0,5 + 0,2·1,5 = {c(tr[0][1],1)}. Corrección 1: y₁ = 0,5 + 0,1·(f(0; 0,5) + f(0,2; 0,8)) = {c(tr[1][1],3)}, error = |{c(tr[1][1],3)} − 0,8| = {c(tr[1][2],3)} ≥ 10⁻² → sigue. Corrección 2: y₁ = {c(tr[2][1],4)}, error = {c(tr[2][2],4)} &lt; 10⁻² → corta."),
    tabla_err(hn, 5, aprox=lambda v: c(v, 4), err=lambda v: c(v, 4)),
    P("La p7 muestra lo mismo redondeado (EL 0,3; 0,4; 0,4; 0,5; 0,5 y EG 0,0007; 0,002; 0,003; 0,004; 0,006). Erratas menores de las p5–6: la tolerancia aparece como 10⁻² y como “error &gt;= 10⁻¹”; el primer error figura 0,028 y la cuenta da 0,026; y(real) en 0,2 figura 0,82923 (es 0,82930)."),
    P(f"Comparando AB2 y Heun (p8): Heun tiene un EG final {c(hn[-1][4],4)} contra {c(ab[-1][4],4)} de AB2 (los dos son de orden 2, pero Heun corrige con el Trapecio implícito y evalúa f más veces). En común: orden 2, arrancan con Euler. Preguntas de la p7: ¿asegura convergencia? Sólo si la iteración del corrector converge (en stiff con h grande puede no hacerlo, §18). ¿Costo? Una evaluación por corrección. ¿Consistente? Sí (Euler y Trapecio lo son)."),
    TRAMPA("F20-P1 y F18-P2 piden Heun en el cuadro (orden, ecuación, tipo, ventajas/desventajas). WA25 pregunta si RKF es P-C: no. F18-P1 pregunta si AB3 “pertenece a un esquema P-C”: puede <b>usarse</b> como predictor de uno, pero por sí solo es un AB explícito."),
)

sc = R["stiff_cat"]; amp = R["amp"]
sec(18, "Ecuaciones stiff (rígidas)", 3, "EDO4 p9–12",
    P("Hay ecuaciones donde un método que funciona bien en general da basura salvo que h sea ridículamente chico. Es el caso de las stiff: lo que obliga a achicar h es un término que dura muy poco."),
    D("Sistema rígido (stiff)", "Un sistema rígido (stiff) es aquel que tiene componentes que cambian rápidamente, junto con componentes de cambio lento. En muchos casos los componentes de variación rápida son efímeros o transitorios, que desaparecen, después de lo cual la solución es dominada por componentes de variación lenta. Aunque los fenómenos transitorios existen sólo en una pequeña parte del intervalo de integración, pueden determinar el tiempo en toda la solución. En general son ecuaciones cuyas soluciones exactas tienen un término de la forma exp(−Cx), donde C es una cte positiva grande y es denominado término transitorio (decae rápidamente a cero al aumentar t) lo que hace inestable la ecuación.", "EDO4 p9"),
    P("El gráfico de la p9 (corriente de un rotor): variación muy rápida e inestable al principio (término transitorio) y después variación lenta."),
    H4("Ejemplo de la cátedra (p10–11): y′ = −21y + e^(−t), y(0) = 0"),
    P("Exacta (sympy en <code>verif_u4.py</code>): y = (e^(−t) − e^(−21t))/20: el término e^(−21t) es el transitorio (C = 21)."),
    FIG(F_STIFF_CAT),
    UL(f"Con <b>h = 0,1</b>: Euler oscila con amplitud creciente (en t = 1 da {c(sc['Euler'][-1],4)}, la exacta es {c(R['stiff_cat_ex'][-1],4)}); RK2 se va para abajo ({c(sc['RK2'][-1],4)}); RK4 sigue a la solución ({c(sc['RK4'][-1],4)}).",
       "Con <b>h variable</b> (RK23 y RK45) los dos siguen a la solución.",
       "Heun con h = 0,1 (11 puntos) se va a −4×10⁹ en la diapositiva; con h = 0,001 (1001 puntos) da la curva correcta.",
       "Pregunta de la p11: ¿qué métodos elegirían? → ver las observaciones."),
    AGR(f"Por qué pasa (verif_u4): para y′ = λy con λ = −21 y h = 0,1, cada paso multiplica el error por un factor: Euler 1 + hλ = {c(amp['euler'],2)} (|·| &gt; 1 → oscila y crece), RK2 1 + hλ + (hλ)²/2 = {c(amp['rk2'],3)} (crece), RK4 = {c(amp['rk4'],3)} (decrece: estable), Euler implícito 1/(1 − hλ) = {c(amp['eimpl'],3)}, Trapecio (1 + hλ/2)/(1 − hλ/2) = {c(amp['trap'],3)} (estables para cualquier h). En Heun iterado, la corrección multiplica la diferencia entre iteraciones por |hλ/2| = {c(amp['heun_iter'],2)} &gt; 1: la iteración del corrector diverge. Con h = 0,001 Euler da un error de {gs(R['stiff_cat_euler_h0001'])} en t = 1."),
    H4("Observaciones de la cátedra (p12) = las estrategias"),
    UL("Para ecuaciones stiff, la elección de un tamaño de paso grande hace que en general el método no converja a la solución exacta. Una estrategia → <b>achicar el tamaño del paso</b> (métodos de paso fijo).",
       "Los <b>métodos de paso variable</b> suelen andar mejor ya que ajustan el tamaño del paso.",
       "En todos los casos de ecuaciones stiff lo recomendable es tener un <b>método con una región de estabilidad grande</b> que permite mejores aproximaciones. Los <b>métodos implícitos</b> en general tienen mayor región de estabilidad.",
       "<b>El único método A-estable es el método del trapecio.</b> Aunque este no proporciona aproximaciones exactas con tamaños de pasos grandes, los errores no crecen exponencialmente.",
       "En general las ecuaciones stiff son fáciles de reconocer."),
    AGR("“El único A-estable es el trapecio” es la frase de la cátedra; en la bibliografía, Euler implícito también es A-estable, y lo que es único es otra cosa: entre los métodos multipaso lineales ningún A-estable supera el orden 2, y el Trapecio es el de orden 2 con menor error (segunda barrera de Dahlquist). En el examen citá la frase de la cátedra y, si querés, agregá que Euler implícito también tiene región de estabilidad que contiene todo el semiplano izquierdo."),
    TRAMPA("IT25-T1 P5: <i>“Explicar qué posibles estrategias recomendaría para una mejor convergencia si se intenta resolver numéricamente una EDO stiff. Explicar en cada caso las ventajas y desventajas. Defina ecuación diferencial stiff.”</i> Las tres estrategias son las tres primeras viñetas de la p12; las ventajas/desventajas salen de §14, §15 y §16–17. F20-P1 agrega la columna “para ecuación stiff” en el cuadro."),
)

sis = R["sis_cat"]
sec(19, "Sistemas de EDOs", 2, "EDO4 p13–16",
    P("Muchos modelos tienen varias incógnitas acopladas (presa-depredador, circuitos). Los métodos para una EDO se aplican igual, componente a componente."),
    CMD("y₁′ = f₁(x, y₁, y₂, …, yₙ)          y₁(x₀) = y₁⁰\ny₂′ = f₂(x, y₁, y₂, …, yₙ)          y₂(x₀) = y₂⁰\n…                                   …\nyₙ′ = fₙ(x, y₁, y₂, …, yₙ)          yₙ(x₀) = yₙ⁰\n\nVectorialmente:   Y′(x) = F(x, Y),   Y(x₀) = Y₀"),
    P("Observaciones de la cátedra:"),
    UL("Se tiene un valor inicial por cada ecuación.", "Se resuelven de manera <b>sincrónica (o simultánea)</b> todas las ecuaciones.",
       "Al finalizar se tiene un vector solución.", "Cálculo del error se hace mediante <b>norma infinito</b> (el mayor elemento en valor absoluto).",
       "Se generalizan los métodos para una EDO."),
    H4("Ejemplo (p14–16)"),
    CMD("u₁′ = 3u₁ + 2u₂ − (2t² + 1)e^(2t)\nu₂′ = 4u₁ + u₂ + (t² + 2t − 4)e^(2t)\nu₁(0) = 1,  u₂(0) = 1        Calcular U = (u₁, u₂) para t = 1 con h = 0,2"),
    P("N = 2 (orden del sistema); f₁ y f₂ son los lados derechos. En la tabla cada fila calcula <code>u₁ⁱ⁺¹ = Met(t, u₁ⁱ, u₂ⁱ, f₁)</code> y <code>u₂ⁱ⁺¹ = Met(t, u₁ⁱ, u₂ⁱ, f₂)</code> con los valores <b>del paso anterior de las dos</b>, y EL = ‖(u₁ⁱ − u₁ⁱ⁺¹, u₂ⁱ − u₂ⁱ⁺¹)‖∞."),
    P(f"Con Euler, primer paso (verif_u4): f₁(0, 1, 1) = 3 + 2 − 1 = {c(sis['f0'][0],0)}, f₂(0, 1, 1) = 4 + 1 − 4 = {c(sis['f0'][1],0)} → u₁¹ = 1 + 0,2·4 = {c(sis['euler1'][0],1)}, u₂¹ = 1 + 0,2·1 = {c(sis['euler1'][1],1)}; EL = ‖(0,8; 0,2)‖∞ = {c(sis['el'],1)}."),
    AGR("La p16 escribe Euler como u₁¹ = u₁⁰ <b>−</b> h·(…): es una errata, Euler suma: u¹ = u⁰ <b>+</b> h·f(t, u⁰)."),
)

sec(20, "EDO de orden superior → sistema de primer orden", 2, "EDO4 p17–18",
    P("Los métodos están escritos para y′ = f(x, y). Una EDO de orden m se transforma en un sistema de m ecuaciones de primer orden y se resuelve con §19."),
    CMD("y⁽ᵐ⁾(x) = f(x, y, y′, y″, …, y⁽ᵐ⁻¹⁾)\ny(x₀) = α₀,  y′(x₀) = α₁,  …,  y⁽ᵐ⁻¹⁾(x₀) = α₍ₘ₋₁₎          Ecuación diferencial de orden m"),
    D("Idea", "La idea es llevar a la ecuación diferencial de orden superior a un sistema de ecuaciones diferenciales ordinarias de primer orden mediante una serie de reemplazos de variables (no siempre es posible).", "EDO4 p17"),
    H4("Ejemplo de la cátedra (p18): y‴ = 3y″ − 2√y′ + e^(xy), y(0) = 3, y′(0) = 4, y″(0) = 5"),
    P("EDO de 3er orden, 3 valores iniciales. Cambio de variables <code>y₁ = y, y₂ = y′, y₃ = y″</code>:"),
    CMD("y₁′ = y₂\ny₂′ = y₃\ny₃′ = 3y₃ − 2√y₂ + e^(x·y₁)\ny₁(0) = 3,  y₂(0) = 4,  y₃(0) = 5"),
    AGR("La diapositiva escribe y₃′ = 3<b>y₂</b> − 2√y₂ + e^(x·y₁); como y″ = y₃, el primer término es 3<b>y₃</b>."),
    H4("Proceso (WA25: “resolver una EDO de orden superior usando cualquier método”)"),
    OL("Despejar la derivada de mayor orden: y⁽ᵐ⁾ = f(x, y, …, y⁽ᵐ⁻¹⁾).",
       "Definir m variables: y₁ = y, y₂ = y′, …, yₘ = y⁽ᵐ⁻¹⁾.",
       "Escribir el sistema: y₁′ = y₂, …, yₘ₋₁′ = yₘ, yₘ′ = f(x, y₁, …, yₘ), con las m condiciones iniciales.",
       "Aplicar el método elegido a las m ecuaciones <b>simultáneamente</b> en cada paso (cada k o cada pendiente es un vector).",
       "El resultado buscado es la primera componente (y₁ = y); el error se mide con norma infinito."),
)

sec(21, "Cuadro comparativo de métodos", 3, "Síntesis EDO2–EDO4",
    P("F20-P1 y F18-P2 piden exactamente esto: orden, ecuación, tipo, ventajas, desventajas y comportamiento con stiff. Todo lo que figura sale de las diapositivas citadas; los órdenes están verificados por Taylor en <code>verif_u4.py</code>."),
    TBL(["Método", "Ecuación", "Orden (EG)", "Tipo", "Ventajas", "Desventajas", "Stiff"], [
        ["Serie de Taylor", "yᵢ₊₁ = y + hy′ + h²y″/2 + …", "el del término truncado", "1 paso, explícito", "“El más estable de todos”", "Necesita derivadas de y (costoso), no se usa en la práctica", "—"],
        ["Euler explícito (RK1)", "yₙ₊₁ = yₙ + h·f(xₙ, yₙ)", "O(h)", "1 paso, explícito, RK", "Fácil; estable (no el más)", "Error alto; necesita h muy chico", "Malo con h grande (oscila y crece, EDO4 p10)"],
        ["Euler implícito (RK1)", "yₙ₊₁ = yₙ + h·f(xₙ₊₁, yₙ₊₁)", "O(h)", "1 paso, implícito, RK", "Región de estabilidad amplia", "Hay que despejar (no siempre se puede)", "Bueno (implícito)"],
        ["Trapecio", "yₙ₊₁ = yₙ + (h/2)(f(xₙ, yₙ) + f(xₙ₊₁, yₙ₊₁))", "O(h²)", "1 paso, implícito, RK2", "“El más estable (A-estable)”", "Hay que despejar o iterar; con h grande no es exacto", "El recomendado: “los errores no crecen exponencialmente”"],
        ["Euler modificado (RK2)", "yₙ₊₁ = yₙ + (h/2)(k₁ + k₂), k₂ = f(xₙ₊₁, yₙ + hk₁)", "O(h²)", "1 paso, explícito, RK2", "Orden 2 con 2 evaluaciones", "2 evaluaciones por paso", "Inestable con h grande (deriva, EDO4 p10)"],
        ["Punto medio (RK2)", "yₙ₊₁ = yₙ + h·f(xₙ + h/2, yₙ + (h/2)f(xₙ, yₙ))", "O(h²)", "1 paso, explícito, RK2", "Orden 2 con 2 evaluaciones", "2 evaluaciones por paso", "Como RK2"],
        ["Heun (P-C)", "P: Euler explícito; C: Trapecio, iterado hasta EL &lt; tol", "O(h²) (el del corrector)", "1 paso, P-C (P explícito, C implícito)", "Mejora la predicción; no hay que despejar", "&gt; costo: 1 evaluación por corrección", "Con h grande la iteración puede divergir (EDO4 p11)"],
        ["RK4", "yₙ₊₁ = yₙ + (h/6)(k₁ + 2k₂ + 2k₃ + k₄)", "O(h⁴)", "1 paso, explícito, RK", "“+ difundido”; h no tiene que ser muy chico", "4 evaluaciones por paso", "Anduvo con h = 0,1 en el ejemplo; explícito: región acotada"],
        ["RKF45", "RK4 y RK5 con las mismas 6 k", "4 (fórmulas de orden 4 y 5)", "1 paso, explícito, paso variable", "Ajusta h; error acotado sin conocer la exacta", "≫ cálculos y evaluaciones (6 por paso)", "“Suele andar mejor” (paso variable)"],
        ["AB2", "yₙ₊₁ = yₙ + (h/2)[3fₙ − fₙ₋₁]", "O(h²)", "2 pasos, explícito, Adams", "Reutiliza f anteriores", "No se autoinicia (arranque con Euler)", "Explícito: región chica"],
        ["AB3", "yₙ₊₁ = yₙ + (h/12)[23fₙ − 16fₙ₋₁ + 5fₙ₋₂]", "O(h³)", "3 pasos, explícito, Adams", "Orden 3", "Necesita 3 valores de arranque", "Explícito: región chica"],
        ["AM (AM2/AM3 de la cátedra)", "yₙ₊₁ = yₙ + (h/12)[5fₙ₊₁ + 8fₙ − fₙ₋₁] / (h/24)[9fₙ₊₁ + 19fₙ − 5fₙ₋₁ + fₙ₋₂]", "O(h³) / O(h⁴)", "multipaso, implícito, Adams", "“En general dan mejores resultados que AB”", "Hay que explicitar; se usan como corrector", "Implícito: más estable que AB"],
        ["Milne (P-C)", "P: yₙ₋₃ + (4/3)h(2fₙ − fₙ₋₁ + 2fₙ₋₂); C: yₙ₋₁ + (h/3)(fₙ₊₁ + 4fₙ + fₙ₋₁)", "O(h⁴)", "multipaso, P-C", "Orden 4", "Necesita 4 valores de arranque; &gt; costo", "—"],
    ]),
    AGR("Las columnas “Stiff” de RK2/RK4 salen del ejemplo de la cátedra (EDO4 p10) y la de AB, de la regla “los implícitos tienen mayor región de estabilidad” (EDO3 p4, EDO4 p12): la diapositiva no analiza AB con stiff. En el TP4 1b (stiff) AB2 con h = 0,2 explota (ver Práctica)."),
    TRAMPA("Rasgo de otra familia (mecanismo 9): “Heun es multipaso”, “RKF es predictor-corrector”, “AB2 es de un paso”, “RK4 es implícito porque k₄ usa xₙ + h”. Cada celda “tipo” se justifica mirando qué valores usa la fórmula (§6).", real=False),
)

sec(22, "TP4: cálculo a mano con aritmética de 3 dígitos", 2, "TP4",
    P("No hay examen práctico de EDO en el material; lo único es el TP4 (predicción). Hay que saber hacer un paso de cada método a mano, completar el cuadro y sacar conclusiones."),
    P("Enunciado: resolver a lápiz con <b>h = 0,2</b> los PVI a) y′ = 1 + y/t, 1 ≤ t ≤ 2, y(1) = 2 (exacta t ln t + 2t); b) y′ = 5e^(5t)(y − t)² + 1, −1 ≤ t ≤ 0, y(−1) = −1 − e⁵ ≈ −149,413 (exacta t − e^(−5t)); c) y′ = 3t/y, 1 ≤ t ≤ 2, y(1) = 2 (exacta √(3t² + 1)); d) bacterias que crecen a una tasa de 3000e^(2t/5), con 7.500 en t = 0, hasta t = 1. Métodos: <b>Runge-Kutta 2, Runge-Kutta 4, Adams-Bashforth 2, Heun</b>. “Realizar las operaciones utilizando aritmética de tres dígitos. Calcular el error local y global en cada caso que sea posible.” Cuadro: <code>t | y real | y aprox | E Local | E Global</code>. Después: implementar en Python con h = 10⁻², 10⁻³, 10⁻⁴ y concluir sobre estabilidad, tamaño de paso, consistencia y convergencia."),
    H4("Convenciones usadas en este hub (todas en verif_u4.py)"),
    UL("<b>Aritmética de 3 dígitos:</b> se redondea (no se trunca) a 3 cifras significativas el resultado de <b>cada</b> operación (+, −, ×, ÷ y cada e^x, √, ln), en el orden en que está escrita la fórmula. h/2 = 0,1 es exacto; en RK4 se hace h·(k₁ + 2k₂ + 2k₃ + k₄) y después ÷6 (en vez de multiplicar por h/6 = 0,0333, que mete un redondeo extra). La consigna no explicita si es redondeo o truncamiento: se eligió redondeo.",
       "<b>RK2:</b> la variante que desarrolla la cátedra (Trapecio/Euler modificado, EDO2 p27): k₁ = f(tₙ, yₙ), k₂ = f(tₙ + h, yₙ + h·k₁), yₙ₊₁ = yₙ + (h/2)(k₁ + k₂).",
       "<b>Heun:</b> el de la cátedra (EDO4 p4–7): predictor Euler, corrector Trapecio <b>iterado</b> mientras |Yant − y| ≥ 10⁻² y it &lt; 5.",
       "<b>AB2:</b> arranque con Euler, como en el ejemplo de la cátedra (EDO3 p13: “método de inicio”).",
       "<b>Errores:</b> E local = |yᵢ₊₁ − yᵢ| y E global = |y(tᵢ) − yᵢ| (EDO2 p7). La y real y los errores se calculan con precisión completa (si se redondeara la y real a 3 cifras, el EG quedaría tapado por el redondeo) y se muestran con 3 cifras significativas."),
    P(f"Lo que muestra la cuenta del 1a: con 3 dígitos, RK2, RK4 y Heun dan <b>los mismos</b> valores (y(2) ≈ {g3(R['tp_a']['RK2']['filas'][-1][2])}), porque la diferencia entre ellos (con precisión completa: RK2 {c(R['tp_a_full']['RK2'][-1],4)}, Heun {c(R['tp_a_full']['Heun'][-1],4)}, RK4 {c(R['tp_a_full']['RK4'][-1],4)}, exacta {c(V.ya_ex(2.0),4)}) es menor que el redondeo de la tercera cifra. El error que domina es el de <b>redondeo</b>, no el de truncamiento: por eso el TP después pide implementarlo en Python."),
    TRAMPA("Predicción (no hay examen práctico de EDO en el material): un paso de RK2/RK4/Heun/AB2 a mano con aritmética de 3 dígitos y el cuadro con EL y EG con la convención de la cátedra. Ejercicios resueltos en Desarrollo → Práctica.", real=False),
)

# ---------------------------------------------------------------- AUDITORÍA
AUD = [
    ("Definiciones: ecuación diferencial, EDO, primer orden, orden superior (EDO1 p2–5)", "m4-s1"),
    ("Ley de Malthus (EDO1 p6)", "m4-s1"),
    ("Problema de valor inicial (PVI) (EDO1 p7, p9–10)", "m4-s2"),
    ("Teorema de existencia y unicidad / Lipschitz (EDO1 p8)", "m4-s2"),
    ("PVI: grilla, ¿hasta cuándo?, tamaño de h, calidad (EDO1 p11–14; EDO2 p2–4)", "m4-s3"),
    ("PVI: estabilidad del PVI y del método (EDO1 p15–17)", "m4-s4"),
    ("PVI: errores de truncamiento y redondeo, E local, E global, resultado general (EDO2 p5–7)", "m4-s5"),
    ("Clasificación de métodos (EDO2 p8)", "m4-s6"),
    ("Serie de Taylor (EDO2 p9–10)", "m4-s7"),
    ("Métodos de Runge-Kutta (RK) (EDO2 p11–12)", "m4-s8"),
    ("Método de Euler Explícito (RK1) y su ejemplo (EDO2 p13–21)", "m4-s9"),
    ("Método de Euler Implícito (RK1) (EDO2 p22–25)", "m4-s10"),
    ("Métodos Runge-Kutta 2 (RK2): Trapecio / Euler modificado, punto medio (EDO2 p26–28)", "m4-s11"),
    ("Métodos Runge-Kutta 4 (RK4) (EDO2 p29)", "m4-s12"),
    ("Métodos Runge-Kutta: #eval vs orden (EDO2 p30)", "m4-s12"),
    ("Repaso (EDO3 p2)", "m4-s13"),
    ("Evaluación de métodos: consistencia, convergencia, estabilidad (EDO3 p3–4)", "m4-s13"),
    ("Elección de h y ejemplo sismológico (EDO3 p5–6)", "m4-s14"),
    ("Métodos de paso variable / RKF45 y ejemplo (EDO3 p7–9)", "m4-s15"),
    ("Métodos Multipasos / Métodos de Adams (EDO3 p10–11)", "m4-s16"),
    ("Método de Adams-Bashforth 2 (AB2), ejemplo y comparación (EDO3 p12–16)", "m4-s16"),
    ("Métodos de Adams-Bashforth AB3, AB4 (EDO3 p17)", "m4-s16"),
    ("Método de Adams-Moulton (EDO3 p18)", "m4-s16"),
    ("Métodos Predictor-Corrector (EDO4 p2–4)", "m4-s17"),
    ("Método de Heun: ejemplo y comparación con AB2 (EDO4 p5–8)", "m4-s17"),
    ("Ecuaciones Stiff (EDO4 p9–12)", "m4-s18"),
    ("Sistemas de EDOs (EDO4 p13–16)", "m4-s19"),
    ("EDO de orden superior (EDO4 p17–18)", "m4-s20"),
    ("TP4: ej. 1 a–d con RK2, RK4, AB2, Heun, 3 dígitos (TP4)", "m4-s22"),
    ("Cuadro comparativo (F20-P1, F18-P2)", "m4-s21"),
]
auditoria = [{"tema": t, "seccion": s, "estado": "cubierto"} for t, s in AUD]

# ---------------------------------------------------------------- GLOSARIO
glosario = [
    {"termino": "PVI", "def": "Problema de valor inicial: y′ = f(x, y) con y(x₀) = y₀. La condición inicial elige una curva de la familia de soluciones.", "seccion": "m4-s2"},
    {"termino": "Condición de Lipschitz", "def": "‖f(t, y) − f(t, ŷ)‖ ≤ L‖y − ŷ‖ en D; si f la cumple, el PVI tiene solución única en una región (EDO1 p8).", "seccion": "m4-s2"},
    {"termino": "Estabilidad del PVI", "def": "Condición intrínseca del problema: el PVI está bien planteado si perturbaciones chicas en f o en y₀ cambian poco la solución exacta.", "seccion": "m4-s4"},
    {"termino": "Estabilidad del método", "def": "Un método es estable si pequeñas perturbaciones en las condiciones iniciales producen cambios igualmente pequeños en las aproximaciones posteriores (condición intrínseca del método).", "seccion": "m4-s4"},
    {"termino": "Error de truncamiento", "def": "Asociado al método (cortar la serie de Taylor); se combate bajando h o cambiando a un método más estable / de mayor orden.", "seccion": "m4-s5"},
    {"termino": "Error de redondeo", "def": "Asociado al uso del computador: crece con la cantidad de cálculos (más pasos, h más chico).", "seccion": "m4-s5"},
    {"termino": "Error local (cátedra)", "def": "Se genera en cada salto de la grilla; en las tablas de la cátedra, |yᵢ₊₁ − yᵢ|. Se relaciona con la consistencia.", "seccion": "m4-s5"},
    {"termino": "Error global", "def": "Suma teórica de los errores locales; |yᵢ − y(xᵢ)|. Sólo se calcula si se conoce la exacta. Se relaciona con la convergencia.", "seccion": "m4-s5"},
    {"termino": "Orden de un método", "def": "p tal que el método coincide con la serie de Taylor hasta hᵖ: error local O(hᵖ⁺¹), global O(hᵖ).", "seccion": "m4-s8"},
    {"termino": "Función de incremento", "def": "F(h, xₙ, yₙ, f) en yₙ₊₁ = yₙ + h·F: la pendiente promedio que usa un RK (Euler: F = f(xₙ, yₙ)).", "seccion": "m4-s8"},
    {"termino": "Método explícito", "def": "yₙ₊₁ = MÉTODO(xₙ, yₙ): el lado derecho sólo usa valores ya conocidos.", "seccion": "m4-s6"},
    {"termino": "Método implícito", "def": "yₙ₊₁ = MÉTODO(xₙ, yₙ, xₙ₊₁, yₙ₊₁): yₙ₊₁ aparece a la derecha (en f); hay que despejar o iterar.", "seccion": "m4-s6"},
    {"termino": "Método de un paso", "def": "Sólo usa el punto anterior (xₙ, yₙ). Todos los RK, Taylor, Heun.", "seccion": "m4-s6"},
    {"termino": "Método multipaso", "def": "Usa más de un punto anterior (yₙ₋₁, yₙ₋₂…). No se autoinicia: necesita un método de arranque.", "seccion": "m4-s16"},
    {"termino": "Consistencia", "def": "lim_{h→0} max|Tₙ| = 0 con Tₙ el error de truncamiento local: el método aproxima la EDO correcta. Depende de h; relaciona el E local con el método.", "seccion": "m4-s13"},
    {"termino": "Convergencia", "def": "lim_{h→0} max|yᵢ − y(xᵢ)| = 0: la solución numérica tiende a la real. Relaciona el E global con el método. Consistente + estable ⇒ convergente.", "seccion": "m4-s13"},
    {"termino": "A-estable", "def": "Método cuyo error no crece para y′ = λy con Re λ < 0 sin importar h. Según la cátedra, el único A-estable es el Trapecio (el más estable).", "seccion": "m4-s18"},
    {"termino": "Ecuación stiff", "def": "Tiene componentes que cambian rápido (transitorio exp(−Cx), C grande) junto con componentes lentas; obliga a h muy chico con métodos explícitos.", "seccion": "m4-s18"},
    {"termino": "Término transitorio", "def": "exp(−Cx) con C positiva grande: decae rápido a cero pero hace inestable la ecuación para los métodos con h grande.", "seccion": "m4-s18"},
    {"termino": "Paso variable (adaptativo)", "def": "Ajustar h automáticamente para que el error de truncamiento quede acotado; RKF45 compara RK4 y RK5 en cada paso.", "seccion": "m4-s15"},
    {"termino": "Adams-Bashforth / Adams-Moulton", "def": "Multipaso de Adams (integran un polinomio interpolante): AB explícito (b₋₁ = 0), AM implícito (b₋₁ ≠ 0).", "seccion": "m4-s16"},
    {"termino": "Predictor-corrector", "def": "Predice con un explícito (orden n) y corrige iterando un implícito (orden n+1); el orden es el del corrector. Ej.: Heun (Euler + Trapecio).", "seccion": "m4-s17"},
    {"termino": "Norma infinito", "def": "Mayor componente en valor absoluto; con ella se mide el error en sistemas de EDO.", "seccion": "m4-s19"},
]

# ---------------------------------------------------------------- QUIZ
Q = []
def q(sec_, tema, pide, sobre, opciones, expl, fig=None):
    it = {"id": f"m4-q{len(Q)+1:02d}", "seccion": f"m4-s{sec_}", "tema": tema, "pide": pide, "sobre": sobre,
          "opciones": [{"t": t, "v": v} for t, v in (sorted(opciones, key=lambda o: o[1]) if pide == "una" else opciones)], "explicacion": expl}
    if fig: it["fig"] = fig
    Q.append(it)

qe = R["q_eimpl"]
q(6, "Clasificar Euler implícito", "una",
  "El método yₙ₊₁ = yₙ + h·f(xₙ₊₁, yₙ₊₁) es…",
  [("De un paso e implícito: la pendiente se evalúa en el punto que se quiere calcular", True),
   ("De un paso y explícito: para avanzar sólo usa el punto anterior (xₙ, yₙ) de la grilla", False),
   ("Multipaso e implícito: usa dos abscisas de la grilla, xₙ y xₙ₊₁, en cada paso", False),
   ("Multipaso y explícito de orden 2: combina dos valores de y en cada paso", False)],
  "Es Euler implícito (EDO2 p22). Es de un paso porque el único punto anterior que usa es (xₙ, yₙ); es implícito porque yₙ₊₁ aparece dentro de f, igual que en la definición yₙ₊₁ = MÉTODO(xₙ, yₙ, xₙ₊₁, yₙ₊₁). Usar xₙ₊₁ no lo hace multipaso (multipaso = usar puntos anteriores como yₙ₋₁). Decir que es explícito confunde ‘usa un solo punto anterior’ con ‘no usa la incógnita’. Su orden es 1, no 2.")

q(16, "AB3 (F18-P1)", "verdaderas",
  "el método Yₙ₊₁ = Yₙ + (1/12)H[23F(n) − 16F(n−1) + 5F(n−2)]",
  [("Es un método de un solo paso", False),
   ("Es un método implícito", False),
   ("Necesita valores de arranque para los primeros puntos, porque no se autoinicia", True),
   ("Su error local es O(h⁴) y su error global O(h³)", True)],
  f"Es AB3 (EDO3 p17). Usa F en n, n−1 y n−2: es multipaso de 3 pasos, no de un paso. No aparece F(n+1): b₋₁ = 0, es explícito (Adams-Bashforth). Como todo multipaso, ‘no se autoinicia’ (EDO3 p10): y₁ e y₂ salen de un método de un paso. Por Taylor (verif_u4.py) el error local es {TY('AB3 (F18-P1)')}: orden {M['AB3 (F18-P1)']['p']}, global O(h³).")

q(13, "Consistencia (definición de la cátedra)", "una",
  "Según la cátedra, la consistencia de un método para resolver EDO…",
  [("Pide que el error de truncamiento local tienda a 0 con h y relaciona el E local con el método", True),
   ("Pide que |yᵢ − y(xᵢ)| tienda a 0 con h y relaciona el E global con el método", False),
   ("No depende de h y garantiza que para cualquier h la solución numérica tiende a la real", False),
   ("Garantiza por sí sola que el error global tienda a 0, aunque el método sea inestable", False)],
  "EDO3 p3: consistente si lim_{h→0} max|Tₙ| = 0 con Tₙ el error de truncamiento local; ‘depende de h’ y ‘relaciona el ELocal con el método’. La que habla de |yᵢ − y(xᵢ)| y del E global es la convergencia, y ‘no depende de h… para cualquier h tiende a la real’ es la observación de la cátedra sobre convergencia. La consistencia sola no alcanza: hace falta además estabilidad para tener convergencia (EDO3 p4).")

q(13, "V/F de WA25 sobre consistencia, convergencia y estabilidad", "verdaderas",
  "consistencia, convergencia y estabilidad de un método para un PVI",
  [("Si un método no es consistente para un PVI, conviene cambiar de método en lugar de achicar h", True),
   ("Si un método es convergente, entonces es estable para cualquier tamaño de paso h", False),
   ("Si un método es estable, puedo aumentar h libremente y bajar los cálculos sin perder precisión", False),
   ("Si un método es consistente y estable, entonces es convergente", True)],
  "Un método inconsistente no aproxima la EDO correcta: achicar h no lo arregla (ej.: yₙ₊₁ = yₙ₋₁ + h·y′ₙ converge a la solución de y′ = f/2), así que se cambia de método. ‘Consistente y estable ⇒ convergente’ es la observación textual de EDO3 p4. ‘Convergente ⇒ estable para cualquier h’ es falsa: Euler explícito es convergente y con h = 0,1 en y′ = −21y + e^(−t) es inestable (EDO4 p10). ‘Estable ⇒ aumentar h sin perder precisión’ es falsa: h grande → menor precisión (EDO3 p5) porque el error global es O(hᵖ).")

q(18, "Definición de stiff", "una",
  "¿Cuál es la definición de ecuación stiff que da la cátedra?",
  [("Tiene componentes que cambian rápido junto con otras de cambio lento; el transitorio exp(−Cx) con C grande la hace inestable", True),
   ("Es una EDO cuya solución crece exponencialmente, como y′ = y, de modo que cualquier error inicial se amplifica", False),
   ("Es una EDO de orden superior que no se puede llevar a un sistema de primer orden con cambio de variables", False),
   ("Es una EDO que sólo se puede resolver con métodos implícitos porque su f no cumple la condición de Lipschitz", False)],
  "EDO4 p9: sistema rígido = componentes de cambio rápido junto con componentes de cambio lento; el término exp(−Cx) con C grande es el transitorio. y′ = y es el ejemplo de EDO inestable de EDO1 p16 (otro concepto: el PVI mismo amplifica). La conversión a sistema (no siempre posible) es el tema de orden superior. Que ‘sólo se puede resolver con implícitos’ es un absoluto falso: la cátedra también achica h o usa paso variable.")

q(18, "Stiff con gráfico", "una",
  "La figura muestra y′ = −21y + e^(−t), y(0) = 0, resuelta con h = 0,1 por dos métodos (A y B). ¿Qué explica la curva del método A?",
  [("Es Euler explícito con un h fuera de su región de estabilidad: cada paso amplifica el error y la curva oscila creciendo", True),
   ("Es un error de redondeo: con tantos puntos en [0, 1] se acumulan los redondeos y la curva se separa de la real", False),
   ("El método A no es consistente con la EDO; por eso, aunque se achique h, nunca se acercará a la solución real", False),
   ("La EDO es inestable (como y′ = y): cualquier método amplifica las perturbaciones, B también lo hará más adelante", False)],
  f"Es el ejemplo de EDO4 p10 recalculado en verif_u4.py: con hλ = {c(amp['z'],1)} el factor de Euler es 1 + hλ = {c(amp['euler'],1)}, de módulo mayor que 1, y el error cambia de signo y crece; B (RK4, factor {c(amp['rk4'],3)}) sigue a la solución. No es redondeo: son sólo 10 pasos. Euler es consistente: con h = 0,001 el error en t = 1 es {gs(R['stiff_cat_euler_h0001'])}. La EDO no es inestable: la exacta (e^(−t) − e^(−21t))/20 decae; lo que la hace difícil es el transitorio.", fig=fig_q_stiff())

q(18, "Estrategias para stiff", "verdaderas",
  "las estrategias para resolver una EDO stiff según la cátedra",
  [("Achicar el tamaño del paso es una estrategia (métodos de paso fijo)", True),
   ("Los métodos de paso variable suelen andar mejor porque ajustan el tamaño del paso", True),
   ("Un método explícito de orden alto como RK4 resuelve cualquier ecuación stiff con cualquier h", False),
   ("Conviene un método con región de estabilidad grande; los implícitos en general la tienen", True)],
  "Las tres verdaderas son las observaciones de EDO4 p12: achicar h (paso fijo), paso variable (ajusta h) y métodos con región de estabilidad grande (implícitos; el Trapecio es A-estable). La de RK4 es un absoluto falso: RK4 es explícito y su región de estabilidad es acotada; en el ejemplo anduvo con h = 0,1 pero con un h más grande (o un C más grande) también falla. En el TP4 1b con h = 0,2 RK4 tiene errores de decenas en los primeros pasos.")

q(9, "Orden de Euler explícito", "una",
  "¿Por qué Euler explícito es de orden 1?",
  [("Coincide con 2 términos de la serie de Taylor: error local O(h²) y global O(h)", True),
   ("Porque usa un solo punto anterior para calcular el siguiente, es decir, es de un paso", False),
   ("Coincide con 1 término de la serie de Taylor: error local O(h) y global O(1)", False),
   ("Porque evalúa f una sola vez por paso y cada evaluación aporta un orden", False)],
  "EDO2 p15–16: yₙ₊₁ = y(xₙ) + h·y′(xₙ) coincide con dos términos de la ST; el primer término que falta es h²y″/2 → error trunc local O(h²), orden del método 1, error global O(h). Ser de un paso es otra clasificación (Taylor y RK4 también son de un paso y tienen otro orden). Si coincidiera con un solo término no sería consistente. La relación evaluaciones–orden de EDO2 p30 no es ‘una evaluación por orden’ (con 5 evaluaciones se sigue en orden 4).")

q(11, "Constantes del punto medio", "una",
  "Con F = γ₁f(x, y) + γ₂f(x + αh, y + βh·f(x, y)), ¿qué constantes dan el método del punto medio yₙ₊₁ = yₙ + h·f(xₙ + h/2, yₙ + (h/2)f(xₙ, yₙ))?",
  [("γ₁ = 0, γ₂ = 1, α = ½, β = ½", True),
   ("γ₁ = ½, γ₂ = ½, α = 1, β = 1", False),
   ("γ₁ = 1, γ₂ = 0, α = ½, β = ½", False),
   ("γ₁ = 0, γ₂ = 1, α = 1, β = ½", False)],
  "En el punto medio toda la pendiente se toma en el segundo término (γ₁ = 0, γ₂ = 1), evaluado en x + h/2 (α = ½) con y + (h/2)f (β = ½). γ₁ = γ₂ = ½, α = β = 1 es el Euler modificado de EDO2 p27 (la diapositiva escribe β = 2: errata). γ₁ = 1, γ₂ = 0 deja F = f(x, y): es Euler. α = 1 con β = ½ evalúa en x + h con medio paso en y: no es ninguno de los métodos de la cátedra y baja a orden 1.")

q(12, "RK4: evaluaciones y orden", "una",
  "Sobre RK4 (yₙ₊₁ = yₙ + (h/6)(k₁ + 2k₂ + 2k₃ + k₄)), ¿qué es correcto?",
  [("Evalúa f 4 veces por paso; error local O(h⁵), global O(h⁴)", True),
   ("Evalúa f 4 veces por paso; error local O(h⁴), global O(h⁵)", False),
   ("Evalúa f 6 veces por paso porque combina una fórmula de orden 4 y otra de orden 5", False),
   ("Es implícito porque k₄ se evalúa en xₙ + h, el punto que se quiere calcular", False)],
  "RK4 tiene cuatro k (cuatro evaluaciones) y orden 4: coincide con la ST hasta h⁴, error local O(h⁵) y global O(h⁴) (resultado general de EDO2 p6; órdenes empíricos ≈ 3,99 en verif_u4.py). Decir local O(h⁴) y global O(h⁵) invierte local y global. Seis evaluaciones con fórmulas de orden 4 y 5 es RKF45. k₄ usa xₙ + h pero con yₙ + h·k₃, que ya se conoce: es explícito.")

q(16, "Adams y multipaso", "verdaderas",
  "los métodos multipaso de Adams",
  [("Se basan en integrar un polinomio interpolante que pasa por los puntos anteriores", True),
   ("En la forma general yₙ₊₁ = yₙ + hΣbⱼf(xₙ₋ⱼ, yₙ₋ⱼ), b₋₁ ≠ 0 corresponde a Adams-Bashforth", False),
   ("No se autoinician: necesitan otro método para los primeros puntos", True),
   ("Adams-Moulton es explícito porque sus coeficientes se conocen de antemano", False)],
  "EDO3 p11: yₙ₊₁ = yₙ + ∫P(x)dx con P interpolante por los puntos anteriores (es la respuesta a la V/F incompleta de WA25 ‘los métodos multipaso explícitos se basan en…’). b₋₁ = 0 es Adams-Bashforth (explícito) y b₋₁ ≠ 0, Adams-Moulton (implícito): la afirmación que asigna b₋₁ ≠ 0 a Adams-Bashforth está invertida. ‘No se autoinician’ es textual de EDO3 p10. Que los coeficientes sean conocidos no tiene nada que ver: AM es implícito porque aparece f(xₙ₊₁, yₙ₊₁).")

q(16, "Arranque de AB2", "una",
  "En el ejemplo de la cátedra con AB2, ¿cómo se obtiene y₁?",
  [("Con un método de un paso de un orden menor, por ejemplo Euler: y₁ = y₀ + h·f(x₀, y₀)", True),
   ("Con el propio AB2, tomando f(x₋₁, y₋₁) = f(x₀, y₀) porque todavía no hay punto anterior", False),
   ("Con Adams-Moulton como corrector, iterando hasta que el error local sea menor que 10⁻²", False),
   ("No hace falta: AB2 es explícito y por eso se autoinicia desde la condición inicial", False)],
  "EDO3 p12–13: ‘dado el PVI, sólo tengo (x₀, y₀) y necesito (x₁, y₁) para iniciar AB2 → uso un método de un orden menor (x ejemplo Euler)’; en la tabla y₁ = 0,8 figura como ‘método de inicio’. Inventar f(x₋₁, y₋₁) no está en el material. AM como corrector es el esquema P-C (EDO4), no el arranque. Ser explícito no lo autoinicia: ningún multipaso se autoinicia (EDO3 p10).")

q(17, "Predictor-corrector", "verdaderas",
  "los métodos predictor-corrector",
  [("El orden del método es el orden del corrector", True),
   ("El predictor es un método explícito y el corrector un método implícito", True),
   ("Tienen menor costo que aplicar sólo el predictor, porque convergen en menos pasos", False),
   ("Heun es un P-C con Euler explícito como predictor y Trapecio como corrector", True)],
  "EDO4 p2: predicción con un método explícito (en gral de orden n), corrección con un implícito (orden n+1), ‘el orden del método es el orden del corrector’ y ‘> costo’: cada corrección es una evaluación más de f, así que no cuestan menos que el predictor solo. EDO4 p4: Heun = P Euler explícito + C Trapecio.")

q(15, "RKF45", "una",
  "¿En qué consiste Runge-Kutta-Fehlberg (RKF45)?",
  [("En cada paso calcula dos aproximaciones RK de distinto orden, estima el error y agranda o achica h", True),
   ("Es un multipaso: usa los cuatro puntos anteriores para predecir y el quinto para corregir", False),
   ("Es un predictor-corrector: RK4 predice y RK5 corrige iterando hasta que EL < ε", False),
   ("Es un RK implícito de orden 5 con paso fijo, más estable que RK4 para ecuaciones stiff", False)],
  "EDO3 p7–8: el objetivo es el paso adaptativo; ‘aproximo el EG con el EL mediante el uso de dos métodos de RK de diferentes órdenes en cada punto’; si no se alcanza la precisión se reduce h, si sobra se incrementa. Es de un paso y explícito (las k usan sólo valores conocidos). No es P-C: no itera un corrector implícito; las dos fórmulas son explícitas y usan las mismas seis k. Paso fijo es justamente lo que no es.")

q(15, "RKF45 V/F", "verdaderas",
  "el método RKF45 de la cátedra",
  [("Es un método de paso variable (adaptativo)", True),
   ("Requiere más evaluaciones de f por paso que RK4", True),
   ("Su objetivo es que el error de truncamiento quede acotado", True),
   ("Necesita la solución exacta para calcular el error global en cada paso", False)],
  "EDO3 p7–9: ajusta h automáticamente (paso adaptativo) para que el error de truncamiento no supere un valor aceptable; ‘RKF45 ≫ cálculos y evaluaciones de la función’ (6 k por paso contra 4). No usa la exacta: justamente ‘no conozco el EG porque no conozco la solución’, por eso estima el error comparando dos fórmulas de distinto orden.")

q(5, "Error global", "una",
  "Según la cátedra, el error global en xᵢ es…",
  [("|yᵢ − y(xᵢ)|: teóricamente, la sumatoria de todos los errores locales hasta ese punto", True),
   ("|yᵢ₊₁ − yᵢ|: la diferencia entre dos aproximaciones consecutivas de la grilla", False),
   ("El error que se genera sólo en el último paso, al llegar a xₙ = b", False),
   ("El error de redondeo acumulado, que no depende del método sino del computador", False)],
  "EDO2 p6–7: E. Global ‘teóricamente corresponde a la sumatoria de todos los errores locales’ y es la diferencia entre la aproximada yᵢ y la verdadera y(xᵢ). |yᵢ₊₁ − yᵢ| es la columna E local de las tablas de la cátedra. El error de un solo paso es el local. El redondeo es otro tipo de error (del computador); el global de truncamiento depende del método (O(hᵖ)).")

q(5, "Resultado general local → global", "una",
  "Si el error de truncamiento local de un método es O(h⁴), entonces…",
  [("El error global es O(h³) y el método es de orden 3", True),
   ("El error global es O(h⁴) y el método es de orden 4", False),
   ("El error global es O(h⁵) y el método es de orden 5", False),
   ("El error global es O(h⁴) y el método es de orden 3", False)],
  "Resultado general (EDO2 p6): local O(hⁿ⁺¹) ⇒ global O(hⁿ); con n + 1 = 4, global O(h³) y orden 3 (así AB3 o RK3). Confundir local con global da ‘orden 4’; sumar en lugar de restar da O(h⁵). ‘Global O(h⁴) y orden 3’ mezcla las dos cosas: si el método es de orden 3 su global es O(h³), no O(h⁴).")

q(6, "Clasificar IT25-T2", "verdaderas",
  "el método yₙ₊₁ = (4/3)yₙ − (1/3)yₙ₋₁ + (2/3)h·y′ₙ₊₁",
  [("No es Runge-Kutta, porque usa yₙ₋₁ y los RK son de un solo paso", True),
   ("Es implícito, porque y′ₙ₊₁ = f(xₙ₊₁, yₙ₊₁) depende de la incógnita", True),
   ("Es de Adams-Moulton, porque es un multipaso implícito", False),
   ("Es de orden 3, porque combina tres términos (yₙ, yₙ₋₁ y y′ₙ₊₁)", False)],
  f"Usa yₙ₋₁ → multipaso, no RK (los RK son de un paso, EDO2 p12). y′ₙ₊₁ es f(xₙ₊₁, yₙ₊₁) → implícito. No es de Adams: la forma de Adams es yₙ₊₁ = yₙ + hΣbⱼf (EDO3 p11) y acá yₙ y yₙ₋₁ tienen coeficientes 4/3 y −1/3 (es la BDF2, nombre que no da la cátedra). El orden no se cuenta por términos: por Taylor (verif_u4.py) el error local es {TY('IT25-T2')} → orden {M['IT25-T2']['p']}.")

q(6, "IT25-T4 tal como está escrito", "una",
  "La figura aplica yₙ₊₁ = yₙ₋₁ + h·y′ₙ (con h, como en el enunciado de IT25-T4) a y′ = y, y(0) = 1. ¿Qué muestra?",
  [("Que el método no es consistente: al achicar h converge a e^(t/2), la solución de y′ = y/2, y no a eᵗ", True),
   ("Que es el método del punto medio de 2 pasos, de orden 2, y que con h = 0,01 ya coincide con eᵗ", False),
   ("Que es consistente pero inestable: el error crece con t y se corrige achicando el paso h", False),
   ("Que es Adams-Bashforth de 2 pasos con un error de redondeo que se acumula en cada paso", False)],
  f"Por Taylor (verif_u4.py): y(xₙ₊₁) − y(xₙ₋₁) − h·y′(xₙ) = {TY('IT25-T4 (h, como dice el enunciado)')} + O(h³): falta coincidir en el término h·y′, así que no es consistente. En la figura, con h = 0,1 y h = 0,01 la aproximación se pega a e^(t/2) (en t = 1 da {c(R['fig_t4'][10][1][-1],3)} y {c(R['fig_t4'][100][1][-1],3)}, frente a e^(1/2) = {c(R['t4h_sqrt_e'],3)}), no a eᵗ: achicar h no lo arregla. Con 2h sería el punto medio de 2 pasos (orden 2). No es AB2 (que usa yₙ y 3fₙ − fₙ₋₁) y el desvío no es redondeo.", fig=fig_q_t4())

q(21, "Orden desde gráfico log-log", "una",
  "La figura muestra el error global en t = 2 en función de h (ejes logarítmicos) para Euler, RK2 y RK4 en y′ = 1 + y/t. ¿Qué asignación es correcta?",
  [("Curva 1: Euler (pendiente ≈ 1); curva 2: RK2 (≈ 2); curva 3: RK4 (≈ 4)", True),
   ("Curva 1: RK4, porque es la que tiene el error más grande con h = 0,2", False),
   ("Curva 3: Euler, porque su error es el más chico al tener menos evaluaciones", False),
   ("Las tres tienen pendiente 1, porque el error global es O(h) en todo método de un paso", False)],
  f"Si el error global es O(hᵖ), en log-log es una recta de pendiente p. Las pendientes medidas en verif_u4.py son {c(R['pend']['Euler'],2)}, {c(R['pend']['RK2'],2)} y {c(R['pend']['RK4'],2)}: Euler (orden 1), RK2 (orden 2) y RK4 (orden 4). RK4 tiene el error más chico, no el más grande; menos evaluaciones no significa menos error (EDO2 p30: más orden ⇒ más evaluaciones). El orden no es 1 para todo método de un paso: depende de hasta qué término coincide con Taylor.", fig=fig_q_loglog())

q(21, "Error global vs t", "una",
  "La figura muestra |E global| (escala log) en función de t para dos métodos de un paso con h = 0,1 en y′ = 1 + y/t, y(1) = 2. Si son Euler explícito y RK4, ¿qué afirmación es correcta?",
  [("Q es RK4: su error global es O(h⁴), del orden de 10⁻⁶, mientras P (Euler) es O(h), del orden de 10⁻²", True),
   ("P es RK4: al evaluar f cuatro veces acumula más redondeo y su error global termina siendo mayor", False),
   ("Q es Euler: su error local es menor porque sólo coincide con dos términos de la serie de Taylor", False),
   ("No se puede saber: con el mismo h los dos métodos tienen el mismo error, sólo cambia el costo", False)],
  f"Los datos de la figura (verif_u4.py): en t = 2, Euler {gs(R['fig_eg']['euler'][-1])} y RK4 {gs(R['fig_eg']['rk4'][-1])}. RK4 coincide con Taylor hasta h⁴ y Euler sólo hasta h: con h = 0,1 eso da varios órdenes de magnitud de diferencia. El redondeo es despreciable frente a esos valores con precisión doble. ‘Coincide con dos términos’ describe a Euler, pero eso le da más error, no menos. El orden sí cambia el error con el mismo h (EDO2 p21: ‘cómo impacta el orden del método en la calidad’).", fig=fig_q_eg())

q(21, "Cuadro F20-P1", "verdaderas",
  "Trapecio, Euler explícito, AB2 y Heun",
  [("Según la cátedra, el Trapecio es el único método A-estable", True),
   ("AB2 es de un paso porque sólo agrega una evaluación nueva de f por paso", False),
   ("Heun es de orden 2 porque ese es el orden de su corrector, el Trapecio", True),
   ("Euler explícito es el recomendado para ecuaciones stiff porque es el más fácil de implementar", False)],
  "EDO4 p12: ‘El único método A-estable es el método del trapecio’. Heun es P-C y ‘el orden del método es el orden del corrector’ (EDO4 p2): Trapecio, orden 2. AB2 usa fₙ y fₙ₋₁: es multipaso de 2 pasos aunque reutilice la evaluación anterior. Euler explícito es fácil pero con stiff y h grande oscila y diverge (EDO4 p10): para stiff se recomiendan métodos con región de estabilidad grande.")

q(19, "Sistemas y orden superior", "verdaderas",
  "los sistemas de EDO y las EDO de orden superior",
  [("Una EDO de orden m se lleva a un sistema de m ecuaciones de primer orden con m condiciones iniciales", True),
   ("En un sistema, el error se calcula con la norma infinito", True),
   ("Las ecuaciones se resuelven una por una: primero u₁ en todo [a, b] y después u₂", False),
   ("El cambio de variables para llevar una EDO de orden superior a un sistema siempre es posible", False)],
  "EDO4 p13 y p17: un valor inicial por ecuación, error con norma infinito, y la EDO de orden m se lleva a un sistema de primer orden (y₁ = y, y₂ = y′, …). Las ecuaciones se resuelven ‘de manera sincrónica (o simultánea)’: cada paso usa u₁ y u₂ del paso anterior, no se puede terminar una antes de empezar la otra. La cátedra aclara ‘(no siempre es posible)’: el ‘siempre’ es el absoluto falso.")

q(20, "Cambio de variables", "una",
  "Para y‴ = 3y″ − 2√y′ + e^(xy), con y₁ = y, y₂ = y′, y₃ = y″, la ecuación de y₃′ es…",
  [("y₃′ = 3y₃ − 2√y₂ + e^(x·y₁)", True),
   ("y₃′ = 3y₃ − 2√y₃ + e^(x·y₁)", False),
   ("y₃′ = 3y₂ − 2√y₁ + e^(x·y₃)", False),
   ("y₃′ = 3y″ − 2√y′ + e^(xy), con y₃ = y‴", False)],
  "y₃′ = y‴ y se reemplaza cada derivada por su variable: y″ = y₃, y′ = y₂, y = y₁ (EDO4 p18; la diapositiva escribe 3y₂ en el primer término, errata: es 3y₃). Poner √y₃ confunde y′ con y″ (√y′ = √y₂); la opción con 3y₂ − 2√y₁ + e^(x·y₃) corre todas las variables un lugar; la que deja 3y″ − 2√y′ no hace el cambio de variables (el sistema tiene que quedar sólo en y₁, y₂, y₃ y sus primeras derivadas).")

q(10, "Euler implícito a mano", "una",
  f"En el ejemplo de la cátedra con Euler implícito (y′ = y − t² + 1, y(0) = 0,5, h = 0,2), ¿cuánto vale y₁?",
  [(f"y₁ = (0,5 − 0,2·0,2² + 0,2)/(1 − 0,2) ≈ {c(qe['ok'],3)}", True),
   (f"y₁ = 0,5 + 0,2·(0,5 − 0 + 1) = {c(qe['expl'],1)}", False),
   (f"y₁ = 0,5 − 0,2·0,2² + 0,2 = {c(qe['sin_div'],3)}", False),
   (f"y₁ = (0,5 − 0,2·0,2² + 0,2)/(1 + 0,2) ≈ {c(qe['div_mas'],3)}", False)],
  f"Despeje de EDO2 p24: yₙ₊₁ − h·yₙ₊₁ = yₙ − h·t²ₙ₊₁ + h ⇒ yₙ₊₁ = (yₙ − h·t²ₙ₊₁ + h)/(1 − h) → {c(qe['ok'],3)} (la diapositiva: 0,864999). {c(qe['expl'],1)} es Euler explícito (pendiente en el punto de partida). {c(qe['sin_div'],3)} se olvida de pasar h·yₙ₊₁ al otro lado. Dividir por 1 + h es un error de signo al despejar. Calculado en verif_u4.py.")

q(9, "EL y EG en la tabla de la cátedra", "una",
  f"En el ejemplo de Euler de la cátedra, en t = 0,4: y aprox = {c(eu[2][2],3)}, y anterior = {c(eu[1][2],1)}, y real = {c(eu[2][1],8)}. Con la convención de la cátedra:",
  [(f"EL = {c(eu[2][3],3)} y EG = {c(eu[2][4],8)}", True),
   (f"EL = {c(eu[2][4],8)} y EG = {c(eu[2][3],3)}", False),
   (f"EL = {c(eu[2][3],3)} y EG = {c(eu[1][4]+eu[2][4],8)}", False),
   (f"EL = {c(abs(eu[2][1]-eu[1][1]),8)} y EG = {c(eu[2][4],8)}", False)],
  f"EDO2 p7: EL = |yᵢ₊₁ − yᵢ| = |1,152 − 0,8| = {c(eu[2][3],3)}; EG = |y(xᵢ) − yᵢ| = |1,21408765 − 1,152| = {c(eu[2][4],8)} (p19). La opción con EL = {c(eu[2][4],8)} y EG = {c(eu[2][3],3)} los intercambia. Un EG de {c(eu[1][4]+eu[2][4],8)} sale de sumar los EG de dos filas (el EG ya es acumulado, no se suma). Un EL de {c(abs(eu[2][1]-eu[1][1]),8)} es la diferencia entre dos valores reales: eso no es ningún error del método.")

q(4, "Estabilidad del PVI y del método", "verdaderas",
  "la estabilidad según la cátedra",
  [("La estabilidad del PVI es una condición intrínseca del problema", True),
   ("La estabilidad del método es una condición intrínseca del método", True),
   ("En los gráficos de la cátedra, y′ = y es el ejemplo de EDO estable", False),
   ("Si el PVI es inestable, usar un método estable lo vuelve un problema estable", False)],
  "EDO1 p16–17: la estabilidad del PVI es ‘condición intrínseca del problema’ y la del método ‘condición intrínseca del método’. y′ = y es el ejemplo de EDO inestable (las curvas se separan); la estable es y′ = −y. Un método estable no cambia la naturaleza del PVI: la calidad depende de ambas estabilidades y de h (diagrama de EDO1 p17).")

q(14, "Elección de h", "verdaderas",
  "el tamaño de paso h",
  [("Un h muy chico implica más costo computacional y más error de redondeo", True),
   ("Un h grande implica menor precisión", True),
   ("Si la solución tiene mucha variabilidad, conviene un h grande para ahorrar cálculos", False),
   ("Si aumenta el orden del método, h necesariamente debe ser muy chico", False)],
  "EDO1 p14 y EDO3 p5: n ≫ (h ≪) → más cálculos y ≫ error de redondeo; h ≫ → menor precisión. EDO3 p6: ‘≫ variabilidad → h ≪’ (la afirmación de h grande con mucha variabilidad lo invierte). EDO2 p30: ‘si aumenta el orden del método, h no necesariamente debe ser muy chico’; la afirmación con ‘necesariamente’ lo niega con un absoluto.")

q(10, "Métodos implícitos", "verdaderas",
  "los métodos implícitos",
  [("Euler implícito y Euler explícito tienen el mismo orden", True),
   ("Un método implícito siempre se puede explicitar algebraicamente", False),
   ("Un método implícito puede usarse como corrector de un esquema predictor-corrector", True),
   ("Según la cátedra, los AM en general dan mejores resultados que los AB", True)],
  "EDO2 p25: ‘ambos son del mismo orden’ (1). La cátedra pregunta ‘¿siempre se podrá explicitar la ecuación?’: no, sólo cuando f permite despejar yₙ₊₁ (en el ejemplo es lineal en y). EDO3 p18: los AM ‘en general dan mejores resultados que los métodos explícitos (AB)’ y ‘en la práctica se utilizan para mejorar las aproximaciones obtenidas con otro método’ → predictor-corrector.")

q(16, "Clasificar AM3 de la cátedra", "una",
  "El método yₙ₊₁ = yₙ + (h/24)[9f(xₙ₊₁, yₙ₊₁) + 19f(xₙ, yₙ) − 5f(xₙ₋₁, yₙ₋₁) + f(xₙ₋₂, yₙ₋₂)] es…",
  [("Multipaso implícito de Adams-Moulton, de orden 4 por Taylor", True),
   ("Multipaso explícito de Adams-Bashforth, de orden 4 por Taylor", False),
   ("De un paso e implícito, porque es un Runge-Kutta de 4 etapas", False),
   ("Multipaso implícito de Adams-Moulton, de orden 3 porque es de 3 pasos", False)],
  f"Usa fₙ₋₁ y fₙ₋₂ → multipaso; aparece f(xₙ₊₁, yₙ₊₁) (b₋₁ = 9/24 ≠ 0) → implícito, Adams-Moulton (EDO3 p11 y p18). Por Taylor (verif_u4.py) el error local es {TY('AM3 (cátedra: \'de 3 pasos\')')} → orden {M['AM3 (cátedra: \'de 3 pasos\')']['p']}. La cátedra lo llama ‘AM3, de 3 pasos’ y en EDO4 ‘AM4’: el número de pasos no es el orden. No es RK: usa puntos anteriores.")

# ---------------------------------------------------------------- DESARROLLO
DES = []
def d(sec_, real, preg, modelo, practica=False, fig=None):
    it = {"id": f"m4-d{len(DES)+1:02d}", "seccion": f"m4-s{sec_}", "real": real, "pregunta": preg, "modelo": modelo}
    if practica: it["practica"] = True
    if fig: it["fig"] = fig
    DES.append(it)

kA = RK["IT25-T1 lectura A: RK3 de Kutta (k2 en x+h/2 con k1; k3 en x+h con -k1+2k2)"]
kB = RK["IT25-T1 lectura B: k3 = f(x+h, y+h·k2)"]
kD = RK["IT25-T1 lectura D: implícito (Lobatto IIIA, etapas acopladas)"]
d(6, "IT25-T1 P4",
  "Dado el siguiente método para resolver EDO, yₙ₊₁ = yₙ + h·[(1/6)(k₁ + 4k₂ + k₃)]: a. ¿Es un método Runge Kutta? ¿por qué? b. ¿De qué orden es? ¿por qué? c. ¿Es de un solo paso o es multipaso? ¿por qué? d. ¿Es un método implícito o explícito? ¿por qué?",
  P("<b>Antes de responder:</b> el enunciado no define k₁, k₂, k₃. El orden y el carácter implícito/explícito dependen de cómo se calculan esas etapas; hay que decir qué se supone. La lectura habitual (escribila primero) es el RK de 3 etapas con pesos 1/6, 4/6, 1/6:")
  + CMD("k₁ = f(xₙ, yₙ)\nk₂ = f(xₙ + h/2, yₙ + (h/2)·k₁)\nk₃ = f(xₙ + h, yₙ − h·k₁ + 2h·k₂)")
  + AGR("Estas etapas (RK3 “de Kutta”) no están en las diapositivas: la cátedra sólo da RK1, RK2, RK4 y RKF45. Se eligen porque son las que dan orden 3 con esos pesos.")
  + P("<b>a. Sí, es un método Runge-Kutta</b>, porque tiene la expresión general de los RK, <code>yₙ₊₁ = yₙ + h·F(h, xₙ, yₙ, f)</code> (EDO2 p11), con función de incremento F = (1/6)(k₁ + 4k₂ + k₃): un promedio ponderado de pendientes (los pesos suman 1/6 + 4/6 + 1/6 = 1) evaluadas a partir de (xₙ, yₙ). Es una modificación de la serie de Taylor que aproxima la derivada con evaluaciones de f, sin calcular derivadas de y.")
  + P(f"<b>b. Orden 3</b> (con las etapas de arriba), porque coincide con la serie de Taylor hasta el término h³: error de truncamiento local O(h⁴) y, por el resultado general, error global O(h³). Verificado en <code>verif_u4.py</code> por Taylor (primer término que no coincide: h^{kA['kmin']}) y con las condiciones de orden (orden {kA['butcher']}); al dividir h por 2 el error global se divide por ≈ 8 (órdenes empíricos {', '.join(c(o,2) for o in EMP['RK3 Kutta (IT25-T1 lect. A)']['ords'])}). Coincide con la tabla de la cátedra (EDO2 p30): 3 evaluaciones → mejor orden posible 3.")
  + P(f"Si las etapas fueran otras el orden cambia, por eso hay que aclararlo: con k₃ = f(xₙ + h, yₙ + h·k₂) baja a orden {kB['butcher']}; con etapas implícitas acopladas (Lobatto IIIA, agregado) sube a orden {kD['butcher']}. Lo que sí vale para cualquier lectura: como los pesos suman 1, el método es <b>consistente</b> (orden ≥ 1).")
  + P("<b>c. De un solo paso</b>, porque para calcular yₙ₊₁ sólo usa el punto anterior (xₙ, yₙ); las k se calculan dentro del paso, no son puntos anteriores de la grilla (como todo RK, EDO2 p12).")
  + P("<b>d. Explícito</b> (con las etapas de arriba), porque cada kᵢ usa sólo valores ya conocidos (xₙ, yₙ y las k anteriores): yₙ₊₁ = MÉTODO(xₙ, yₙ). Sería implícito si alguna k dependiera de sí misma o de yₙ₊₁.")
  + AGR("Si f depende sólo de x, este método es la regla de Simpson: yₙ₊₁ = yₙ + (h/6)(f(xₙ) + 4f(xₙ + h/2) + f(xₙ₊₁)). Es lo que pasa en el problema de bacterias del TP4 (1d)."))

t2 = M["IT25-T2"]; t3 = M["IT25-T3"]
d(6, "IT25-T2 P4",
  "Dado el método yₙ₊₁ = (4/3)yₙ − (1/3)yₙ₋₁ + (2/3)h·y′ₙ₊₁: a. ¿Es un método Runge Kutta? ¿por qué? b. ¿De qué orden es? ¿por qué? c. ¿Es de un solo paso o es multipaso? ¿por qué? d. ¿Es un método implícito o explícito? ¿por qué?",
  P("<b>a. No es un método Runge-Kutta</b>, porque usa un valor de un punto anterior al último (yₙ₋₁): los RK son de un solo paso y tienen la forma yₙ₊₁ = yₙ + h·F(h, xₙ, yₙ, f) (EDO2 p11–12). Tampoco es de Adams: la forma de Adams es yₙ₊₁ = yₙ + h·Σbⱼf(xₙ₋ⱼ, yₙ₋ⱼ) (EDO3 p11), y acá combina dos valores de y (4/3 y −1/3) con una sola f.")
  + AGR("Su nombre en la bibliografía es BDF2 (fórmula de diferencias hacia atrás de 2 pasos). La cátedra no lo nombra.")
  + P(f"<b>b. Orden 2</b>, porque al desarrollar por Taylor alrededor de xₙ coincide con la serie hasta h²: el error de truncamiento local es {t2['c0'].replace('*Y3','·y‴').replace('-','−')}·h³ + … = O(h³), y por el resultado general el error global es O(h²). Cuenta (verif_u4.py): 4/3 − 1/3 = 1 (coincide el término y), y el término en h·y′ también coincide, así que es consistente; el primer término que no se cancela es el de h³.")
  + CMD("y(xₙ₊₁) − (4/3)y(xₙ) + (1/3)y(xₙ₋₁) − (2/3)h·y′(xₙ₊₁) = " + TY("IT25-T2") + "(xₙ) + O(h⁴)")
  + P(f"Verificación numérica en y′ = y − t² + 1: al dividir h por 2 el error global se divide por ≈ 4 (órdenes {', '.join(c(o,2) for o in EMP['IT25-T2 (BDF2)']['ords'])}).")
  + P("<b>c. Multipaso (de 2 pasos)</b>, porque usa dos puntos anteriores, yₙ y yₙ₋₁. Por eso no se autoinicia: y₁ hay que obtenerlo con un método de un paso (EDO3 p10).")
  + P("<b>d. Implícito</b>, porque y′ₙ₊₁ = f(xₙ₊₁, yₙ₊₁) depende de la incógnita yₙ₊₁: tiene la forma yₙ₊₁ = MÉTODO(xₙ, yₙ, xₙ₊₁, yₙ₊₁) (EDO2 p8). Hay que despejar yₙ₊₁ o usarlo como corrector de un predictor-corrector.")
  + AGR("Como implícito tiene región de estabilidad amplia (la BDF2 es A-estable): es una opción para ecuaciones stiff."))

d(6, "IT25-T3 P4",
  "Dado el método yₙ₊₁ = (18/11)yₙ − (9/11)yₙ₋₁ + (2/11)yₙ₋₂ + (6/11)h·y′ₙ₊₁: a. ¿Es un método Runge Kutta? ¿por qué? b. ¿De qué orden es? ¿por qué? c. ¿Es de un solo paso o es multipaso? ¿por qué? d. ¿Es un método implícito o explícito? ¿por qué?",
  P("<b>a. No es Runge-Kutta</b>, porque usa yₙ₋₁ e yₙ₋₂: los RK son de un solo paso (EDO2 p12). Tampoco tiene la forma de Adams (yₙ₊₁ = yₙ + hΣbⱼf): combina tres valores anteriores de y con una sola f.")
  + AGR("Es la BDF3 (diferencias hacia atrás de 3 pasos). La cátedra no la nombra.")
  + P(f"<b>b. Orden 3</b>, porque por Taylor coincide con la serie hasta h³: error local {t3['c0'].replace('*Y4','·y⁗').replace('-','−')}·h⁴ = O(h⁴) → global O(h³) (resultado general EDO2 p6). Es consistente: 18/11 − 9/11 + 2/11 = 1 y el término en h·y′ también coincide (verif_u4.py). Numéricamente, al dividir h por 2 el error se divide por ≈ 8 (órdenes {', '.join(c(o,2) for o in EMP['IT25-T3 (BDF3)']['ords'])}).")
  + CMD("y(xₙ₊₁) − (18/11)y(xₙ) + (9/11)y(xₙ₋₁) − (2/11)y(xₙ₋₂) − (6/11)h·y′(xₙ₊₁) = −(3/22)h⁴·y⁗(xₙ) + O(h⁵)")
  + P("<b>c. Multipaso de 3 pasos</b>, porque necesita yₙ, yₙ₋₁ e yₙ₋₂. No se autoinicia: y₁ e y₂ salen de un método de un paso.")
  + P("<b>d. Implícito</b>, porque y′ₙ₊₁ = f(xₙ₊₁, yₙ₊₁) contiene la incógnita (forma yₙ₊₁ = MÉTODO(xₙ, yₙ, xₙ₊₁, yₙ₊₁)). Hay que despejar o iterar (predictor-corrector)."))

t4 = M["IT25-T4 (h, como dice el enunciado)"]; t4b = M["IT25-T4 (lectura 2h)"]
d(6, "IT25-T4 P4",
  "Dado el método yₙ₊₁ = yₙ₋₁ + h·y′ₙ: a. ¿Es un método Runge Kutta? ¿por qué? b. ¿De qué orden es? ¿por qué? c. ¿Es de un solo paso o es multipaso? ¿por qué? d. ¿Es un método implícito o explícito? ¿por qué?",
  P("<b>Ojo con el enunciado:</b> dice h·y′ₙ. Escrito así, el método <b>no es consistente</b>. El método conocido con esta forma (punto medio de 2 pasos, o “leapfrog”) lleva <b>2h</b>. Puede ser una errata del enunciado o una trampa. Los incisos a, c y d valen igual para las dos lecturas; en b escribí primero el análisis literal y después la lectura 2h.")
  + P("<b>a. No es Runge-Kutta</b>, porque usa yₙ₋₁ (un punto anterior a xₙ): los RK son de un solo paso con la forma yₙ₊₁ = yₙ + h·F(h, xₙ, yₙ, f) (EDO2 p11–12). Tampoco es de Adams, porque no parte de yₙ sino de yₙ₋₁ (la forma de Adams es yₙ₊₁ = yₙ + hΣbⱼf).")
  + P("<b>b. Lectura literal (con h): orden 0, no es consistente.</b> Por Taylor alrededor de xₙ:")
  + CMD("y(xₙ₊₁) − y(xₙ₋₁) = 2h·y′(xₙ) + (h³/3)·y‴(xₙ) + …\ny(xₙ₊₁) − y(xₙ₋₁) − h·y′(xₙ) = h·y′(xₙ) + (h³/3)·y‴(xₙ) + …     (verif_u4.py)")
  + P(f"El error local es O(h) = O(h⁰⁺¹): el método ni siquiera coincide con la ST en el término h·y′, así que su orden es 0 y el error de truncamiento (por unidad de paso) no tiende a 0 → no cumple lim_{{h→0}} max|Tₙ| = 0: <b>no es consistente</b>, y por lo tanto no es convergente. Numéricamente (verif_u4.py, y′ = y, y(0) = 1): en t = 1 da {c(R['t4h'][0][1],5)} con h = 0,1 y {c(R['t4h'][-1][1],5)} con h = 10⁻⁴: converge a e^(1/2) = {c(R['t4h_sqrt_e'],5)}, la solución de y′ = y/2, y no a e = {c(R['t4h_e'],5)}. Achicar h no lo arregla: hay que cambiar de método.")
  + P(f"<b>Lectura con 2h</b> (yₙ₊₁ = yₙ₋₁ + 2h·y′ₙ): <b>orden 2</b>, porque el residuo de Taylor empieza en {TY('IT25-T4 (lectura 2h)')} → error local O(h³), global O(h²) (órdenes empíricos {', '.join(c(o,2) for o in EMP['IT25-T4 lectura 2h (punto medio 2 pasos)']['ords'])}). Es el punto medio de 2 pasos: la pendiente en el medio del intervalo [xₙ₋₁, xₙ₊₁].")
  + AGR("El nombre “punto medio de 2 pasos / leapfrog” no está en el material (el “punto medio” de la cátedra es el RK2 de un paso, EDO2 p28). Además es débilmente estable: en problemas con solución decreciente (y′ = λy, λ &lt; 0) aparece una oscilación espuria que crece, por eso no se recomienda para stiff.")
  + P("<b>c. Multipaso de 2 pasos</b>, porque usa yₙ₋₁ y yₙ (en y′ₙ = f(xₙ, yₙ)). No se autoinicia: y₁ se obtiene con un método de un paso.")
  + P("<b>d. Explícito</b>, porque el lado derecho sólo tiene valores conocidos (yₙ₋₁ y f(xₙ, yₙ)); no aparece yₙ₊₁: yₙ₊₁ = MÉTODO(xₙ₋₁, yₙ₋₁, xₙ, yₙ).")
  + TRAMPA("En el parcial corregido de IT25-T4 la pregunta 4 vale ~2 puntos. Si sólo escribís “orden 2, es el punto medio” y el enunciado decía h, la respuesta es incorrecta para el método escrito; si sólo escribís “no consistente” sin justificar con Taylor, perdés la razón. Escribí las dos lecturas con su desarrollo."))

d(18, "IT25-T1 P5",
  "Explicar qué posibles estrategias recomendaría para una mejor convergencia si se intenta resolver numéricamente una EDO stiff. Explicar en cada caso las ventajas y desventajas. Defina ecuación diferencial stiff.",
  P("<b>Definición</b> (EDO4 p9): una ecuación (sistema) stiff o rígida es aquella que tiene <b>componentes que cambian rápidamente, junto con componentes de cambio lento</b>. Los componentes rápidos suelen ser transitorios que desaparecen, después de lo cual dominan los lentos; aunque los transitorios existen sólo en una pequeña parte del intervalo, pueden determinar el tiempo (el paso) en toda la solución. En general la solución exacta tiene un término exp(−Cx) con C positiva grande (término transitorio), que decae rápidamente a cero pero hace inestable la ecuación para los métodos numéricos. Ejemplo de la cátedra: y′ = −21y + e^(−t), y(0) = 0, con exacta (e^(−t) − e^(−21t))/20.")
  + P("<b>Estrategias</b> (EDO4 p12), con ventajas y desventajas:")
  + TBL(["Estrategia", "Ventajas", "Desventajas"], [
      ["1. <b>Achicar el tamaño del paso h</b> (métodos de paso fijo)", "Sirve con cualquier método, incluso explícitos simples; con h suficientemente chico el método vuelve a ser estable y converge (Heun con h = 0,001 da la curva correcta, EDO4 p11)", "h chico en <b>todo</b> el intervalo, aunque el transitorio dure poco: muchos pasos → mayor costo computacional y más error de redondeo (EDO1 p14, EDO3 p5)"],
      ["2. <b>Métodos de paso variable</b> (RKF45, RK23)", "“Suelen andar mejor ya que ajustan el tamaño del paso”: h chico sólo en el transitorio y grande donde la solución es lenta; controlan el error sin conocer la exacta", "≫ cálculos y evaluaciones de f por paso (6 en RKF45); si el método base es explícito, en el tramo lento h sigue limitado por la estabilidad"],
      ["3. <b>Métodos con región de estabilidad grande: implícitos</b> (Trapecio, Euler implícito, AM)", "“Los métodos implícitos en general tienen mayor región de estabilidad”; el Trapecio es A-estable: “aunque no proporciona aproximaciones exactas con pasos grandes, los errores no crecen exponencialmente”", "Hay que despejar yₙ₊₁ (no siempre se puede) o iterar un corrector (predictor-corrector: &gt; costo); con h grande pierden precisión"],
  ])
  + P(f"Por qué falla un explícito (verif_u4.py, ejemplo de la cátedra con h = 0,1): Euler multiplica el error por 1 + hλ = {c(amp['euler'],2)} en cada paso (|·| &gt; 1: oscila y crece) y RK2 por {c(amp['rk2'],3)}; RK4 ({c(amp['rk4'],3)}), Euler implícito ({c(amp['eimpl'],3)}) y Trapecio ({c(amp['trap'],3)}) quedan por debajo de 1.")
  + AGR(f"Cuidado con el predictor-corrector con h grande: en Heun cada corrección multiplica la diferencia por |hλ/2| = {c(amp['heun_iter'],2)} &gt; 1 en el ejemplo, la iteración del corrector diverge (es lo que muestra la diapositiva de Heun con h = 0,1). Para stiff, el implícito hay que resolverlo exactamente (despeje o Newton), no por iteración de punto fijo."))

d(4, "IT25-T2 P5",
  "Defina la estabilidad de un método para resolver EDOs. Mencione al menos dos métodos estables y explique.",
  D("Estabilidad del método", "Un método numérico es estable si pequeñas perturbaciones en las condiciones iniciales producen cambios igualmente pequeños en las aproximaciones posteriores.", "EDO1 p17; EDO3 p4")
  + P("Es una <b>condición intrínseca del método</b> (distinta de la estabilidad del PVI, que es del problema). La inestabilidad del método está asociada a la <b>propagación del error</b> (EDO2 p6): un método inestable amplifica en cada paso los errores de truncamiento y redondeo. Importa porque, junto con la consistencia, garantiza la convergencia: “si un método es consistente y estable entonces es convergente” (EDO3 p4).")
  + P("<b>Métodos estables:</b>")
  + UL("<b>Trapecio</b> (yₙ₊₁ = yₙ + (h/2)(f(xₙ, yₙ) + f(xₙ₊₁, yₙ₊₁))): para la cátedra es “el método más estable (A-estable)” y “el único método A-estable”: aunque con h grande no da aproximaciones exactas, los errores no crecen exponencialmente. Es implícito, de orden 2.",
       "<b>Euler implícito</b> (yₙ₊₁ = yₙ + h·f(xₙ₊₁, yₙ₊₁)) y en general los implícitos como Adams-Moulton: “los métodos implícitos generalmente son más estables (región de estabilidad más amplia)”.",
       "<b>Euler explícito</b>: la cátedra lo da como “método estable (aunque no es el más estable de todos)”: para problemas bien condicionados y h chico el error no se dispara; con h grande en un stiff sí (EDO4 p10).",
       "<b>Serie de Taylor</b>: “es el más estable de todos los métodos”, pero no se usa porque necesita las derivadas de y.",
       "<b>Runge-Kutta</b> en general: “si el problema es bien condicionado entonces los métodos RK son estables aunque pueden ser lentos”.")
  + P("Para la respuesta alcanzan dos con su explicación; los más seguros son <b>Trapecio</b> (A-estable) y <b>Euler implícito</b> (implícito → región de estabilidad amplia).")
  + AGR(f"Cómo se ve la estabilidad (verif_u4.py): para y′ = λy con λ = −21 y h = 0,1, el factor por paso es Euler explícito {c(amp['euler'],2)} (inestable), Euler implícito 1/(1 − hλ) = {c(amp['eimpl'],3)} y Trapecio (1 + hλ/2)/(1 − hλ/2) = {c(amp['trap'],3)}: los dos implícitos tienen |factor| &lt; 1 para cualquier h &gt; 0 cuando λ &lt; 0."))

d(13, "IT25-T3 P5",
  "Defina consistencia de un método para resolver EDOs. Explicar su relación con los errores.",
  D("Consistencia", "Un método de diferencias con error de truncamiento local Tn en el n-ésimo paso es consistente con la EDO que aproxima si cumple que lim_{h→0} max_{1≤i≤n} |Tn| = 0.", "EDO3 p3")
  + P("Es decir: al achicar el paso, el error que el método comete en cada paso tiende a cero; el método aproxima la EDO correcta. Observaciones de la cátedra: <b>depende de h</b> y <b>relaciona el ELocal con el método</b>.")
  + P("<b>Relación con los errores:</b>")
  + UL("Con el <b>error local (de truncamiento)</b>: la consistencia es una condición sobre él. Por Taylor, si el método coincide con la serie hasta hᵖ, el error local es O(hᵖ⁺¹); el método es consistente si p ≥ 1 (coincide al menos en y(xₙ) + h·y′(xₙ)).",
       "Con el <b>error global</b>: por el resultado general (EDO2 p6), local O(hᵖ⁺¹) ⇒ global O(hᵖ). La consistencia sola no asegura que el error global vaya a 0: hace falta además <b>estabilidad</b> (que los errores locales no se amplifiquen al propagarse). “Si un método es consistente y estable entonces es convergente” (EDO3 p4), y la convergencia es la que se refiere al error global.",
       "Con el <b>redondeo</b>: no interviene en la definición; un h muy chico baja el truncamiento pero sube el redondeo.")
  + EJ(f"Euler explícito: error local (h²/2)y″ → consistente (orden 1). En cambio yₙ₊₁ = yₙ₋₁ + h·y′ₙ deja un residuo h·y′ que no se va: no es consistente y converge a la solución de otra ecuación (y′ = f/2) aunque se achique h (verif_u4.py). Por eso, si un método no es consistente, se cambia de método.")
  + AGR("En la bibliografía (Burden-Faires) el Tₙ de esta definición es el error de un paso dividido por h; con esa convención consistente significa que coincide con la serie al menos en el término h·y′, que es lo que se escribió arriba."))

d(13, "IT25-T4 P5",
  "Defina convergencia de un método para EDOs. Explique cómo se relaciona con los errores.",
  D("Convergencia", "Un método de diferencias es convergente con respecto a la EDO que aproxima si verifica que lim_{h→0} max_{1≤i≤n} |yᵢ − y(xᵢ)| = 0.", "EDO3 p3")
  + P("Es decir: al achicar el paso, la máxima diferencia entre la solución numérica y la verdadera en los puntos de la grilla tiende a cero. Observaciones de la cátedra: no depende de h; para cualquier h la solución numérica va a tender a la real; <b>relaciona el EGlobal con el método</b>.")
  + P("<b>Relación con los errores:</b>")
  + UL("Con el <b>error global</b>: |yᵢ − y(xᵢ)| es exactamente el error global (EDO2 p6–7); converger es que el error global máximo tienda a 0. Para un método de orden p, el error global es O(hᵖ): cuanto mayor el orden, más rápido converge al achicar h.",
       "Con el <b>error local</b>: el global es “teóricamente la sumatoria de todos los errores locales”. Si el error local es O(hᵖ⁺¹) (método consistente, p ≥ 1), el global es O(hᵖ). Pero para que la suma no se dispare los errores no deben amplificarse: <b>consistente + estable ⇒ convergente</b> (EDO3 p4).",
       "Con el <b>redondeo</b>: en la práctica h no puede ir a 0: más pasos → más error de redondeo (EDO2 p5).")
  + EJ(f"Euler en y′ = y − t² + 1 (verif_u4.py): el error global en t = 1 es {gs(EMP['Euler explícito']['errs'][0])} con h = 0,1 y {gs(EMP['Euler explícito']['errs'][1])} con h = 0,05: se divide por ≈ 2 al dividir h por 2 (orden 1). RK4: {gs(EMP['RK4']['errs'][0])} → {gs(EMP['RK4']['errs'][1])} (se divide por ≈ 16, orden 4).")
  + AGR("Ojo con el V/F de WA25 “método convergente → es estable”: la cátedra sólo enuncia consistente + estable ⇒ convergente. En la teoría (Dahlquist), un multipaso lineal convergente es cero-estable, pero eso es para h → 0; con un h fijo grande un método convergente puede ser inestable (Euler en el stiff de EDO4 p10)."))

d(21, "F20-P1",
  "Cuadro comparativo de Trapecio, Euler explícito, AB2 y Heun: orden, tipo, ecuación, ventajas/desventajas, para ecuación stiff.",
  TBL(["", "Trapecio", "Euler explícito", "AB2", "Heun"], [
      ["Ecuación", "yₙ₊₁ = yₙ + (h/2)(f(xₙ, yₙ) + f(xₙ₊₁, yₙ₊₁))", "yₙ₊₁ = yₙ + h·f(xₙ, yₙ)", "yₙ₊₁ = yₙ + (h/2)[3f(xₙ, yₙ) − f(xₙ₋₁, yₙ₋₁)]", "P: yₙ₊₁ = yₙ + h·f(xₙ, yₙ); C: yₙ₊₁ = yₙ + (h/2)(f(xₙ, yₙ) + f(xₙ₊₁, yₙ₊₁)) iterado"],
      ["Orden", "2 (local O(h³), global O(h²))", "1 (local O(h²), global O(h))", "2 (EG O(h²))", "2 (el del corrector)"],
      ["Tipo", "Un paso, implícito, RK2", "Un paso, explícito, RK1", "Multipaso (2 pasos), explícito, Adams-Bashforth", "Un paso, predictor-corrector (P explícito, C implícito)"],
      ["Ventajas", "“El más estable (A-estable)”", "Fácil de implementar; estable (no el más)", "Orden 2 reutilizando f del paso anterior", "Mejora la predicción de Euler hasta una tolerancia; no hay que despejar"],
      ["Desventajas", "Implícito: despejar yₙ₊₁ (no siempre se puede) o iterar; con h grande no es exacto", "Orden lineal: error alto; necesita h muy chico", "No se autoinicia (arranque con Euler); explícito", "&gt; costo: una evaluación de f por corrección"],
      ["Para stiff", "El recomendado: único A-estable según la cátedra, los errores no crecen exponencialmente", "Malo con h grande: oscila y crece (EDO4 p10)", "Explícito → región de estabilidad chica (en TP4 1b con h = 0,2 explota)", "Con h grande la iteración del corrector puede divergir (EDO4 p11, h = 0,1); con h chico funciona"],
  ])
  + P(f"Comparación numérica de la cátedra (y′ = y − t² + 1, h = 0,2, EG en t = 1, verif_u4.py): Euler {c(eu[-1][4],4)}; AB2 {c(ab[-1][4],4)}; Heun {c(hn[-1][4],4)}.")
  + AGR("La fila “AB2 para stiff” y el resultado del TP4 1b no están en las diapositivas (la cátedra no aplica AB2 a stiff); se apoyan en la regla de EDO3 p4 / EDO4 p12 y en la corrida de verif_u4.py."))

d(13, "F20-P2", "Definir consistencia de un método numérico para solución de EDOs.",
  D("Consistencia", "Un método de diferencias con error de truncamiento local Tn en el n-ésimo paso es consistente con la EDO que aproxima si cumple que lim_{h→0} max_{1≤i≤n} |Tn| = 0.", "EDO3 p3")
  + P("Depende de h y relaciona el error local con el método: cuando h → 0, el error de truncamiento que se comete en cada paso tiende a 0. En la práctica se verifica por Taylor: el método tiene que coincidir con la serie al menos en y(xₙ) + h·y′(xₙ) (orden p ≥ 1). Junto con la estabilidad da convergencia (EDO3 p4).")
  + EJ("Todos los métodos de la cátedra (Euler, RK2, RK4, AB, AM, Heun, RKF) son consistentes. yₙ₊₁ = yₙ₋₁ + h·y′ₙ no lo es (residuo h·y′)."))

ab3 = M["AB3 (F18-P1)"]
d(16, "F18-P1",
  "Yₙ₊₁ = Yₙ + (1/12)H[23F(n) − 16F(n−1) + 5F(n−2)]. Responder si es verdadero o falso, justificar: a) Es un método de un solo paso. b) Es un método implícito. c) Pertenece a un esquema Predictor-Corrector.",
  P("Es <b>Adams-Bashforth de 3 pasos (AB3)</b>, EDO3 p17.")
  + P("<b>a) Falso.</b> Usa F(n), F(n−1) y F(n−2), es decir, información de xₙ, xₙ₋₁ y xₙ₋₂: es <b>multipaso de 3 pasos</b>. Un método de un paso sólo usa el punto anterior (EDO2 p8). Consecuencia: no se autoinicia; Y₁ e Y₂ se calculan con otro método.")
  + P("<b>b) Falso.</b> Es <b>explícito</b>: no aparece F(n+1) = f(xₙ₊₁, Yₙ₊₁); en la forma general de Adams b₋₁ = 0 → Adams-Bashforth (EDO3 p11). Yₙ₊₁ = MÉTODO(xₙ, Yₙ, xₙ₋₁, Yₙ₋₁, xₙ₋₂, Yₙ₋₂).")
  + P("<b>c) Falso tal como está</b> (justificar): el método en sí no es predictor-corrector, porque no tiene una etapa de corrección con un método implícito ni itera; es un AB explícito. Lo que sí es cierto es que <b>puede usarse como predictor</b> de un esquema P-C (por ejemplo, con un Adams-Moulton como corrector, como AB4+AM4 en EDO4 p4). Si la cátedra quiso preguntar “¿puede formar parte de un esquema P-C?”, la respuesta es Verdadero como predictor; escribí primero Falso con esta aclaración.")
  + P(f"Extra que suma: orden 3 (Taylor, verif_u4.py: error local {ab3['c0'].replace('*Y4','·y⁗')}·h⁴ → global O(h³))."))

d(21, "F18-P2",
  "Cuadro comparativo dando orden, ecuación, tipo, ventajas y desventajas: Trapecio, Runge Kutta 4, Runge Kutta Fehlberg y Heun.",
  TBL(["", "Trapecio", "RK4", "RKF45", "Heun"], [
      ["Ecuación", "yₙ₊₁ = yₙ + (h/2)(f(xₙ, yₙ) + f(xₙ₊₁, yₙ₊₁))", "yₙ₊₁ = yₙ + (h/6)(k₁ + 2k₂ + 2k₃ + k₄); k₁ = f(xₙ, yₙ), k₂ = f(xₙ + h/2, yₙ + hk₁/2), k₃ = f(xₙ + h/2, yₙ + hk₂/2), k₄ = f(xₙ + h, yₙ + hk₃)", "Dos fórmulas con las mismas k₁…k₆: yₙ₊₁ = yₙ + 25/216k₁ + 1408/2565k₃ + 2197/4104k₄ − 1/5k₅ y ȳₙ₊₁ = yₙ + 16/135k₁ + 6656/12825k₃ + 28561/56430k₄ − 9/50k₅ + 2/55k₆", "P: Euler explícito; C: Trapecio, iterado hasta |Yant − y| &lt; tol"],
      ["Orden", "2", "4 (local O(h⁵), global O(h⁴))", "4 (la cátedra: “el error de truncamiento es de orden 4”; las fórmulas son de orden 4 y 5)", "2 (el del corrector)"],
      ["Tipo", "Un paso, implícito", "Un paso, explícito, RK", "Un paso, explícito, paso variable (adaptativo)", "Un paso, predictor-corrector"],
      ["Ventajas", "A-estable, “el más estable”", "“+ difundido”; h no tiene que ser muy chico", "Ajusta h solo; error acotado sin conocer la exacta; mejor error con menos nodos (ej. EDO3 p9)", "Mejora iterativamente la predicción; no hay que despejar"],
      ["Desventajas", "Hay que despejar o iterar", "4 evaluaciones de f por paso", "≫ cálculos y evaluaciones (6 por paso)", "&gt; costo por las correcciones"],
  ])
  + AGR("Los coeficientes 25/216 y 2197/4104 corrigen una errata de EDO3 p8 (23/216 y 219/4104), verificada en verif_u4.py."))

d(15, "WA25",
  "Runge-Kutta-Fehlberg: objetivo, en qué consiste, cuándo se usa, por qué. Orden, ¿multipaso? ¿implícito? ¿P-C (predictor-corrector)? ¿Con qué otro método lo comparo?",
  UL("<b>Objetivo:</b> ajustar el tamaño del paso h automáticamente (paso adaptativo) para que el error de truncamiento quede acotado, sin superar un valor aceptable (EDO3 p7).",
     "<b>En qué consiste:</b> en cada paso calcula dos aproximaciones de distinto orden, RK4 y RK5, con las mismas seis k (k₁ = h·f(t, y), …, k₆). La diferencia estima el error local; si la precisión esperada no se alcanza, se reduce h y se vuelve a calcular; si es mejor que la esperada, h se incrementa (EDO3 p8).",
     "<b>Cuándo se usa:</b> cuando la solución tiene zonas de distinta variabilidad (tramos tranquilos y tramos con cambios rápidos), en ecuaciones stiff (“los métodos de paso variable suelen andar mejor”, EDO4 p12) y cuando no se conoce la solución exacta para controlar el error.",
     "<b>Por qué:</b> porque no conozco el EG (no conozco la solución): “aproximo el EG con el EL mediante el uso de dos métodos de RK de diferentes órdenes en cada punto”; así se determina la cantidad y posición de los nodos que garantizan el error de truncamiento acotado.",
     "<b>Orden:</b> 4. La cátedra dice “ambos métodos tienen Elocal O(h⁴)” y “el error de truncamiento en ambos métodos es de orden 4” (EDO3 p8–9).",
     "<b>¿Multipaso?</b> No: es de un paso, cada yᵢ₊₁ sale sólo de (tᵢ, yᵢ).",
     "<b>¿Implícito?</b> No: es explícito, cada k usa sólo las anteriores.",
     "<b>¿Predictor-corrector?</b> No: no predice con un explícito para corregir iterando un implícito; compara dos fórmulas explícitas para decidir el paso.",
     "<b>¿Con qué se compara?</b> Con RK4 de paso constante (EDO3 p9): para y′ = 1 + y², y(0) = 0, RKF45 llega a t = 1,4 con error −0,0006208 usando 11 nodos (achica h cerca de la asíntota de tan t), y RK4 con h = 0,1 tiene error 0,0059089; a cambio, RKF45 hace “≫ cálculos y evaluaciones de la función”.")
  + AGR(f"Por las condiciones de orden (verif_u4.py) la fórmula de 4º orden es de orden {RK[list(RK)[-3]]['butcher']} y la de 5º de orden {RK[list(RK)[-2]]['butcher']} (de ahí el 45); el valor que se avanza tiene orden 4."))

d(13, "WA25",
  "Verdadero o falso, justificar: a) Si un método no es consistente para un PVI, se recomienda cambiar de método. b) Si un método es convergente, entonces es estable. c) Si un método es estable, puedo aumentar h y bajar la cantidad de cálculos. d) Los métodos multipaso explícitos se basan en… (completar).",
  P("<b>a) Verdadero.</b> Si no es consistente, el error de truncamiento local no tiende a 0 con h: el método aproxima otra ecuación y no converge (consistencia es condición necesaria). Achicar h no lo arregla (a diferencia de un problema de estabilidad, donde la cátedra propone “bajar el tamaño del paso o cambiar por un método más estable”, EDO2 p5). Ej.: yₙ₊₁ = yₙ₋₁ + h·y′ₙ converge a la solución de y′ = f/2 (verif_u4.py).")
  + P("<b>b) Falso como afirmación general</b> (con la justificación): la cátedra enuncia sólo la implicación “consistente y estable ⇒ convergente” (EDO3 p4), no su recíproca. Además la estabilidad en la práctica depende de h: Euler explícito es convergente (consistente y estable para h chico) y aun así con h = 0,1 en la ecuación stiff y′ = −21y + e^(−t) es inestable (oscila y crece, EDO4 p10).")
  + AGR("Matiz de teoría: para métodos multipaso lineales, convergente ⇒ cero-estable (teorema de Dahlquist); esa “estabilidad” es para h → 0. Si te piden V/F sin más, escribí Falso con el contraejemplo del stiff y nombrá que la cátedra sólo da la implicación consistente + estable ⇒ convergente.")
  + P("<b>c) Falso.</b> La estabilidad sola no garantiza precisión: con h grande el error global O(hᵖ) crece (“h ≫ → &lt; precisión”, EDO3 p5), y en problemas stiff un h grande hace que el método “en general no converja” (EDO4 p12). La cátedra pide un h “suficientemente grande para minimizar el costo computacional <b>conservando estabilidad y precisión</b>”.")
  + P("<b>d) Completar:</b> “…en integrar un polinomio interpolante que pasa por los puntos anteriores” (métodos de Adams-Bashforth: yₙ₊₁ = yₙ + ∫ₓₙ^ₓₙ₊₁ P(x)dx = yₙ + hΣbⱼf(xₙ₋ⱼ, yₙ₋ⱼ) con b₋₁ = 0, EDO3 p11). Por usar puntos anteriores no se autoinician."))

sis2 = R["sis"]
d(20, "WA25",
  "Explicar el proceso para resolver una EDO de orden superior usando cualquier método.",
  P("La idea (EDO4 p17): “llevar a la ecuación diferencial de orden superior a un sistema de ecuaciones diferenciales ordinarias de primer orden mediante una serie de reemplazos de variables (no siempre es posible)”, y resolver el sistema con el método elegido.")
  + OL("Escribir la EDO despejando la derivada de mayor orden: y⁽ᵐ⁾ = f(x, y, y′, …, y⁽ᵐ⁻¹⁾), con las m condiciones iniciales y(x₀) = α₀, …, y⁽ᵐ⁻¹⁾(x₀) = α₍ₘ₋₁₎.",
       "Cambio de variables: y₁ = y, y₂ = y′, …, yₘ = y⁽ᵐ⁻¹⁾.",
       "Sistema de primer orden: y₁′ = y₂, y₂′ = y₃, …, yₘ′ = f(x, y₁, …, yₘ); condiciones y₁(x₀) = α₀, …, yₘ(x₀) = α₍ₘ₋₁₎. Vectorialmente Y′ = F(x, Y), Y(x₀) = Y₀.",
       "Aplicar el método (Euler, RK4, AB2…) a todas las ecuaciones de manera <b>sincrónica</b>: en cada paso se calculan todas las componentes con los valores del paso anterior (en RK, cada k es un vector).",
       "La solución buscada es la primera componente y₁; el error se mide con norma infinito.")
  + EJ("Ejemplo de la cátedra (EDO4 p18): y‴ = 3y″ − 2√y′ + e^(xy), y(0) = 3, y′(0) = 4, y″(0) = 5 → y₁′ = y₂, y₂′ = y₃, y₃′ = 3y₃ − 2√y₂ + e^(x·y₁), y₁(0) = 3, y₂(0) = 4, y₃(0) = 5 (la diapositiva escribe 3y₂: errata).")
  + EJ(f"TP4 ej. 2d: y″ − 2y′ + 2y = e^(2t)·sen t, y(0) = −0,4, y′(0) = −0,6 → y₁′ = y₂, y₂′ = e^(2t)·sen t − 2y₁ + 2y₂. Un paso de Euler con h = 0,2: y₁ = −0,4 + 0,2·(−0,6) = {c(sis2['euler1'][0],2)}, y₂ = −0,6 + 0,2·({c(sis2['f0'][1],1)}) = {c(sis2['euler1'][1],2)} (verif_u4.py; resuelto completo en Práctica)."))

# ---- variantes rotadas
d(6, None, "Clasificá el método yₙ₊₁ = yₙ + h·f(xₙ₊₁, yₙ₊₁): a. ¿Es RK? ¿por qué? b. ¿Orden? ¿por qué? c. ¿Un paso o multipaso? ¿por qué? d. ¿Implícito o explícito? ¿por qué?",
  P("Es <b>Euler implícito</b> (EDO2 p22).")
  + P("<b>a. Sí, es RK (RK1)</b>, porque tiene la forma yₙ₊₁ = yₙ + h·F(h, xₙ, yₙ, f) con función de incremento F_I = f(xₙ₊₁, yₙ₊₁) (EDO2 p23) y es de un solo paso.")
  + P(f"<b>b. Orden 1</b>, porque coincide con dos términos de la serie de Taylor: error local O(h²) (Taylor en verif_u4.py: {TY('Euler implícito')}), global O(h). La cátedra: “el orden viene dado por la cantidad de aproximaciones de derivadas de la serie de Taylor. En este caso el orden es 1”. Órdenes empíricos {', '.join(c(o,2) for o in EMP['Euler implícito']['ords'])}.")
  + P("<b>c. Un paso</b>, porque sólo usa (xₙ, yₙ) como punto anterior.")
  + P("<b>d. Implícito</b>, porque yₙ₊₁ aparece dentro de f: yₙ₊₁ = MÉTODO(xₙ, yₙ, xₙ₊₁, yₙ₊₁). Hay que despejarlo (no siempre se puede) o iterarlo como corrector. A cambio es más estable que Euler explícito (región de estabilidad más amplia)."))

d(11, None, "Clasificá el método del Trapecio yₙ₊₁ = yₙ + (h/2)(f(xₙ, yₙ) + f(xₙ₊₁, yₙ₊₁)) (los 4 incisos con su porqué) y decí en qué se diferencia de la fórmula “Euler modificado” de la cátedra.",
  P("<b>a. Sí, es RK (RK2)</b>, porque se escribe yₙ₊₁ = yₙ + h·F con F = ½(f(xₙ, yₙ) + f(xₙ₊₁, yₙ₊₁)), promedio de dos pendientes, y es de un paso (EDO2 p26–27).")
  + P(f"<b>b. Orden 2</b>, porque coincide con 3 términos de la ST: error local O(h³), global O(h²) (EDO2 p27; Taylor {TY('Trapecio (implícito)')} en verif_u4.py; órdenes empíricos {', '.join(c(o,2) for o in EMP['Trapecio implícito']['ords'])}).")
  + P("<b>c. Un paso</b>, porque sólo usa el punto anterior (xₙ, yₙ).")
  + P("<b>d. Implícito</b>, porque yₙ₊₁ aparece en f(xₙ₊₁, yₙ₊₁). Es “el más estable (A-estable)”.")
  + P("<b>Diferencia con Euler modificado:</b> en la fórmula de EDO2 p27, yₙ₊₁ = yₙ + ½h(f(xₙ, yₙ) + f(xₙ₊₁, yₙ + h·f(xₙ, yₙ))), la yₙ₊₁ de adentro se reemplaza por la predicción de Euler: queda <b>explícito</b>, también de orden 2, y es el Heun con una sola corrección. La diapositiva lo rotula “(IMPLÍCITO)” porque nace del Trapecio; la fórmula escrita es explícita.")
  + AGR("La diapositiva da β = 2 para esa fórmula: es una errata, con y + h·f(x, y) es β = 1 (con β = 2 el método baja a orden 1, verif_u4.py)."))

d(17, None, "Respecto al método de Heun (P: Euler explícito; C: Trapecio): ¿es un método RK? ¿por qué? ¿De qué orden es? ¿Cuántas veces se evalúa la función para aproximar un punto de la grilla? ¿Es multipaso? ¿Es implícito?",
  P("Son las preguntas de EDO4 p4.")
  + P("<b>¿RK?</b> Con una sola corrección es exactamente el RK2 Euler modificado: yₙ₊₁ = yₙ + (h/2)(k₁ + k₂), k₁ = f(xₙ, yₙ), k₂ = f(xₙ₊₁, yₙ + h·k₁), con la forma yₙ₊₁ = yₙ + h·F y de un paso. En la versión de la cátedra el corrector se itera hasta que |Yant − y| &lt; tolerancia: es un <b>predictor-corrector</b> de un paso; deja de ser un RK con un número fijo de etapas.")
  + P("<b>Orden 2</b>, porque en un P-C “el orden del método es el orden del corrector” (EDO4 p2) y el Trapecio es de orden 2 (error local O(h³)).")
  + P(f"<b>Evaluaciones:</b> f(xₙ, yₙ) una vez (la usan el predictor y todas las correcciones) + una f(xₙ₊₁, y) por corrección. En el ejemplo de la cátedra el primer punto necesitó 2 correcciones (0,8 → {c(tr[1][1],3)} → {c(tr[2][1],4)}): 3 evaluaciones.")
  + P("<b>¿Multipaso?</b> No: sólo usa (xₙ, yₙ). <b>¿Implícito?</b> El corrector (Trapecio) es implícito, pero no se despeja: se itera a partir del valor predicho; el esquema completo se calcula con valores conocidos."))

d(11, None, "Clasificá el método del punto medio yₙ₊₁ = yₙ + h·f(xₙ + h/2, yₙ + (h/2)f(xₙ, yₙ)): a. ¿Es RK? b. ¿Orden y orden del error global? c. ¿Un paso o multipaso? d. ¿Implícito o explícito? e. ¿Valores de las constantes de la función de incremento? (cada uno con su porqué)",
  P("Son las preguntas de EDO2 p28.")
  + P("<b>a. Sí, RK2</b>, porque es de un paso con F = f(xₙ + h/2, yₙ + (h/2)f(xₙ, yₙ)), un caso de la familia F = γ₁f(x, y) + γ₂f(x + αh, y + βh·f(x, y)).")
  + P(f"<b>b. Orden 2, error global O(h²)</b>, porque coincide con la ST hasta h² (error local O(h³)): verif_u4.py por Taylor y por condiciones de orden (orden {RK['Punto medio (RK2, γ1=0, γ2=1, α=β=1/2)']['butcher']}); órdenes empíricos {', '.join(c(o,2) for o in EMP['Punto medio']['ords'])}.")
  + P("<b>c. Un paso</b>: sólo usa (xₙ, yₙ). <b>d. Explícito</b>: la pendiente del medio se evalúa con una predicción de Euler a medio paso, con valores conocidos.")
  + P("<b>e. γ₁ = 0, γ₂ = 1, α = ½, β = ½</b>: toda la pendiente es la del medio, evaluada en x + h/2 con y + (h/2)f(x, y)."))

d(16, None, "Clasificá AB2: yₙ₊₁ = yₙ + (h/2)[3f(xₙ, yₙ) − f(xₙ₋₁, yₙ₋₁)] (RK sí/no, orden, uno/multipaso, implícito/explícito, con su porqué) y explicá cómo se arranca.",
  P("<b>a. No es RK</b>, porque usa f(xₙ₋₁, yₙ₋₁): es multipaso; los RK son de un paso. Es de <b>Adams-Bashforth</b> (forma yₙ₊₁ = yₙ + hΣbⱼf con b₋₁ = 0).")
  + P(f"<b>b. Orden 2</b> (“error de truncamiento es de orden O(h²)”, EDO3 p12), porque integra un polinomio interpolante lineal por (xₙ₋₁, fₙ₋₁), (xₙ, fₙ): por Taylor el error local es {TY('AB2')} (verif_u4.py) → global O(h²).")
  + P("<b>c. Multipaso de 2 pasos</b>: necesita dos puntos anteriores. <b>d. Explícito</b>: no aparece f(xₙ₊₁, yₙ₊₁).")
  + P(f"<b>Arranque:</b> dado el PVI sólo se tiene (x₀, y₀); y₁ se calcula con un método de un paso “de un orden menor (x ejemplo Euler)”: y₁ = y₀ + h·f(x₀, y₀). En el ejemplo de la cátedra: y₁ = {c(ab[1][2],1)}, y₂ = {c(ab[2][2],3)} (EDO3 p13–15)."))

d(16, None, "Clasificá el método yₙ₊₁ = yₙ + (h/12)[5f(xₙ₊₁, yₙ₊₁) + 8f(xₙ, yₙ) − f(xₙ₋₁, yₙ₋₁)]: a. ¿Es RK? b. ¿Orden? c. ¿Un paso o multipaso? d. ¿Implícito o explícito? e. ¿Para qué se usa en la práctica? (con su porqué)",
  P("Es el que la cátedra llama <b>Adams-Moulton 2 (AM2, “de 2 pasos”)</b>, EDO3 p18.")
  + P("<b>a. No es RK</b>, porque usa f(xₙ₋₁, yₙ₋₁) (un punto anterior a xₙ); es un método de Adams (yₙ₊₁ = yₙ + hΣbⱼf).")
  + P(f"<b>b. Orden {M['AM2 (cátedra: \'de 2 pasos\')']['p']}</b>, porque por Taylor el error local es {TY('AM2 (cátedra: \'de 2 pasos\')')} (verif_u4.py) → global O(h³). La cátedra no da el orden de AM; el número “2” es el de pasos, no el orden.")
  + P("<b>c. Multipaso</b> (usa xₙ y xₙ₋₁). <b>d. Implícito</b>: aparece f(xₙ₊₁, yₙ₊₁) (b₋₁ = 5/12 ≠ 0, Adams-Moulton).")
  + P("<b>e.</b> “Los métodos implícitos (AM) tienen el problema de que primero se debe convertir algebraicamente a una manera explícita, aunque en general dan mejores resultados que los métodos explícitos (AB). En la práctica se utilizan para mejorar las aproximaciones obtenidas con otro método” → como <b>corrector</b> en un esquema predictor-corrector (con un AB como predictor)."))

d(17, None, "El método de Milne usa P: yₙ₊₁ = yₙ₋₃ + (4/3)h(2fₙ − fₙ₋₁ + 2fₙ₋₂) y C: yₙ₊₁ = yₙ₋₁ + (h/3)(fₙ₊₁ + 4fₙ + fₙ₋₁). ¿De qué tipo es? ¿De qué orden? ¿Por qué el predictor tiene que ser al menos de orden n si el corrector es n+1? ¿Cuándo usarías un esquema P-C?",
  P("<b>Tipo:</b> predictor-corrector <b>multipaso</b>: el predictor es explícito y usa yₙ₋₃ y f en xₙ, xₙ₋₁, xₙ₋₂ (4 puntos anteriores); el corrector es implícito (aparece fₙ₊₁) y usa yₙ₋₁. No es RK y no se autoinicia: y₁, y₂, y₃ salen de un método de un paso (p. ej. RK4).")
  + P(f"<b>Orden 4</b>, porque en un P-C el orden es el del corrector (EDO4 p2); por Taylor (verif_u4.py) el corrector tiene error local {M['Milne corrector']['c0'].replace('*Y5','·y⁽⁵⁾').replace('-','−')}·h⁵ (orden {M['Milne corrector']['p']}) y el predictor {M['Milne predictor']['c0'].replace('*Y5','·y⁽⁵⁾')}·h⁵ (orden {M['Milne predictor']['p']}).")
  + P("<b>¿Por qué el predictor al menos de orden n?</b> Porque el corrector arranca desde la predicción: si el predictor es mucho peor, hacen falta muchas iteraciones del corrector (más evaluaciones, &gt; costo) o la iteración puede no converger; con un predictor de orden n, pocas correcciones alcanzan para llegar al orden n+1 del corrector.")
  + P("<b>¿Cuándo usar P-C?</b> Cuando se quiere la precisión/estabilidad de un método implícito sin despejar yₙ₊₁ (f no lineal en y) y se puede controlar el error con una tolerancia (el algoritmo de EDO4 p3: mientras |error| &gt; E y it &lt; itMax).")
  + AGR("La cátedra no da el orden de Milne (lo pregunta en EDO4 p4); el 4 sale de la verificación por Taylor."))

# ---- PRÁCTICA
A = R["tp_a"]
def det_rk2(dd, t, y):
    return f"k₁ = f({t}; {y}) = {g3(float(dd['k1']))}; k₂ = f(t + h; y + h·k₁) = {g3(float(dd['k2']))}; y = {y} + 0,1·({g3(float(dd['k1']))} + {g3(float(dd['k2']))}) = {g3(float(dd['y']))}"

a2 = A["RK2"]["det"]
d(22, None, "TP4 ej. 1a con RK2: y′ = 1 + y/t, 1 ≤ t ≤ 2, y(1) = 2 (exacta y = t ln t + 2t), h = 0,2, aritmética de 3 dígitos. Usá el RK2 de la cátedra (Euler modificado: k₁ = f(tₙ, yₙ), k₂ = f(tₙ + h, yₙ + h·k₁), yₙ₊₁ = yₙ + (h/2)(k₁ + k₂)) y completá el cuadro t | y real | y aprox | E local | E global.",
  P("Convención: se redondea a 3 cifras significativas el resultado de cada operación; E local = |yᵢ₊₁ − yᵢ|, E global = |y(tᵢ) − yᵢ| (EDO2 p7), con la y real en precisión completa (verif_u4.py).")
  + P(f"<b>Paso 1</b> (t = 1 → 1,2): k₁ = 1 + 2/1 = {g3(float(a2[0]['k1']))}; y + h·k₁ = 2 + 0,2·3 = 2,6; k₂ = 1 + 2,6/1,2 = 1 + 2,17 = {g3(float(a2[0]['k2']))}; y₁ = 2 + 0,1·(3 + 3,17) = 2 + 0,1·6,17 = 2 + 0,617 = {g3(float(a2[0]['y']))}.")
  + P(f"<b>Paso 2</b> (1,2 → 1,4): k₁ = 1 + 2,62/1,2 = {g3(float(a2[1]['k1']))}; k₂ = f(1,4; 2,62 + 0,2·3,18) = {g3(float(a2[1]['k2']))}; y₂ = {g3(float(a2[1]['y']))}.")
  + tabla_err(A["RK2"]["filas"])
  + P(f"Con precisión completa RK2 da y(2) = {c(R['tp_a_full']['RK2'][-1],4)} (exacta {c(V.ya_ex(2.0),4)}); con 3 dígitos {g3(A['RK2']['filas'][-1][2])}. El EG crece de a poco (se acumula), el EL es el salto entre filas (≈ 0,6–0,7, no mide error del método).")
  + TRAMPA("Predicción: un paso de RK2 a mano con la aritmética de 3 dígitos y la fila del cuadro con EL y EG. Errores típicos: usar h·(k₁ + k₂) en lugar de (h/2)(k₁ + k₂), evaluar k₂ en tₙ en vez de tₙ + h, o calcular el EL contra la exacta.", real=False),
  practica=True)

a4 = A["RK4"]["det"]
d(22, None, "TP4 ej. 1a con RK4: y′ = 1 + y/t, y(1) = 2, h = 0,2, aritmética de 3 dígitos. Calculá k₁…k₄ del primer paso y completá el cuadro hasta t = 2.",
  P("RK4 de la cátedra (EDO2 p29). Convención de 3 dígitos: cada operación se redondea; h/2 = 0,1; se calcula h·(k₁ + 2k₂ + 2k₃ + k₄) y después ÷ 6.")
  + P(f"<b>Paso 1</b> (t = 1, y = 2): k₁ = 1 + 2/1 = {g3(float(a4[0]['k1']))}; k₂ = f(1,1; 2 + 0,1·3 = 2,3) = 1 + 2,3/1,1 = {g3(float(a4[0]['k2']))}; k₃ = f(1,1; 2 + 0,1·{g3(float(a4[0]['k2']))}) = {g3(float(a4[0]['k3']))}; k₄ = f(1,2; 2 + 0,2·{g3(float(a4[0]['k3']))}) = {g3(float(a4[0]['k4']))}; suma k₁ + 2k₂ + 2k₃ + k₄ = {g3(float(a4[0]['s']))}; y₁ = 2 + 0,2·{g3(float(a4[0]['s']))}/6 = {g3(float(a4[0]['y']))}.")
  + tabla_err(A["RK4"]["filas"])
  + P(f"Con precisión completa RK4 da y(2) = {c(R['tp_a_full']['RK4'][-1],6)} contra la exacta {c(V.ya_ex(2.0),6)} (error {gs(abs(R['tp_a_full']['RK4'][-1] - V.ya_ex(2.0)))}). Con 3 dígitos los valores son <b>los mismos que con RK2</b>: la ventaja de orden 4 queda tapada por el redondeo de la tercera cifra (el EG que se ve, ≈ 10⁻³, es redondeo).")
  + TRAMPA("Predicción: calcular las cuatro k de un paso de RK4 con 3 dígitos. Errores típicos: evaluar k₂ y k₃ en tₙ + h en lugar de tₙ + h/2, usar k₁ en vez de k₂ para k₃, olvidar los factores 2 de k₂ y k₃.", real=False),
  practica=True)

aab = A["AB2"]["det"]
d(22, None, "TP4 ej. 1a con Adams-Bashforth 2: y′ = 1 + y/t, y(1) = 2, h = 0,2, aritmética de 3 dígitos. ¿Cómo se obtiene y₁? Completá el cuadro.",
  P("AB2 necesita dos puntos; como en el ejemplo de la cátedra (EDO3 p13), y₁ sale de <b>Euler</b> (“método de inicio”).")
  + P(f"<b>Arranque:</b> f(1; 2) = {g3(float(aab[0]['f']))}; y₁ = 2 + 0,2·3 = {g3(A['AB2']['filas'][1][2])}.")
  + P(f"<b>Paso 2 (AB2):</b> f(1,2; 2,6) = 1 + 2,6/1,2 = {g3(float(aab[1]['f']))}; 3f₁ = {g3(float(aab[1]['tres']))}; 3f₁ − f₀ = {g3(float(aab[1]['dif']))}; y₂ = 2,6 + 0,1·{g3(float(aab[1]['dif']))} = {g3(A['AB2']['filas'][2][2])}.")
  + tabla_err(A["AB2"]["filas"])
  + P(f"El EG de AB2 ({g3(A['AB2']['filas'][-1][4])} en t = 2) es mayor que el de RK2/RK4/Heun ({g3(A['RK2']['filas'][-1][4])}) porque arrastra el error del arranque con Euler (EG = {g3(A['AB2']['filas'][1][4])} ya en t = 1,2). Con precisión completa AB2 da {c(R['tp_a_full']['AB2'][-1],4)}.")
  + TRAMPA("Predicción: “¿AB2 se autoinicia?” (no: arranque con Euler) y el primer paso AB2 con 3f(xₙ, yₙ) − f(xₙ₋₁, yₙ₋₁). Error típico: invertir los coeficientes (3fₙ₋₁ − fₙ) o multiplicar por h en vez de h/2.", real=False),
  practica=True)

ah = A["Heun"]
d(22, None, "TP4 ej. 1a con Heun: y′ = 1 + y/t, y(1) = 2, h = 0,2, aritmética de 3 dígitos, con el Heun de la cátedra (P: Euler; C: Trapecio iterado mientras |Yant − y| ≥ 10⁻², itMax = 5). Completá el cuadro y compará los cuatro métodos.",
  P("Algoritmo de EDO4 p3–6: predictor yₚ = yₙ + h·f(tₙ, yₙ); repetir y ← yₙ + (h/2)(f(tₙ, yₙ) + f(tₙ₊₁, y)) hasta que el cambio sea menor que 10⁻².")
  + P(f"<b>Paso 1:</b> P: y = 2 + 0,2·3 = {g3(float(ah['det'][0][0][1]))}. C1: y = 2 + 0,1·(3 + f(1,2; 2,6)) = 2 + 0,1·(3 + 3,17) = {g3(float(ah['det'][0][1][1]))} (cambio {g3(ah['det'][0][1][2])}, no menor que 10⁻²). C2: y = 2 + 0,1·(3 + f(1,2; 2,62)) = {g3(float(ah['det'][0][2][1]))} (cambio 0) → corta. Correcciones por paso: {ah['its']}.")
  + tabla_err(ah["filas"])
  + TBL(["Método", "y(2) con 3 dígitos", "EG", "y(2) precisión completa", "EG"], [
      [m, g3(A[m]["filas"][-1][2]), g3(A[m]["filas"][-1][4]), c(R["tp_a_full"][m][-1], 6), gs(abs(R["tp_a_full"][m][-1] - V.ya_ex(2.0)))] for m in ("RK2", "RK4", "AB2", "Heun")])
  + P("<b>Conclusiones</b> (las que pide el TP): (1) con precisión completa el orden se nota: RK4 ≪ Heun &lt; RK2 &lt; AB2 (AB2 arrastra el arranque con Euler); (2) con 3 dígitos RK2, RK4 y Heun coinciden: domina el error de redondeo, no el de truncamiento; (3) los cuatro son consistentes (coinciden con la ST al menos en h·y′) y el PVI es estable (las soluciones no se separan), así que todos convergen al achicar h; (4) el EL de la tabla (|yᵢ₊₁ − yᵢ|) no sirve para comparar métodos: es casi igual en todos.")
  + TRAMPA("Predicción: comparar los métodos en el cuadro y justificar con orden, arranque y redondeo. Error típico: concluir que RK4 “no mejora” a RK2 mirando sólo la cuenta con 3 dígitos.", real=False),
  practica=True)

C_ = R["tp_c"]; c2 = C_["RK2"]["det"][0]; c4 = C_["RK4"]["det"][0]
d(22, None, "TP4 ej. 1c, un paso a mano: y′ = 3t/y, y(1) = 2 (exacta √(3t² + 1)), h = 0,2, aritmética de 3 dígitos. Calculá y(1,2) con RK2 (Euler modificado) y con RK4, con sus E local y E global.",
  P(f"<b>RK2:</b> k₁ = 3·1/2 = {g3(float(c2['k1']))}; y + h·k₁ = 2 + 0,3 = 2,3; k₂ = 3·1,2/2,3 = 3,6/2,3 = {g3(float(c2['k2']))}; y₁ = 2 + 0,1·(1,5 + {g3(float(c2['k2']))}) = {g3(float(c2['y']))}.")
  + P(f"<b>RK4:</b> k₁ = {g3(float(c4['k1']))}; k₂ = f(1,1; 2 + 0,1·1,5 = 2,15) = 3,3/2,15 = {g3(float(c4['k2']))}; k₃ = f(1,1; 2 + 0,1·{g3(float(c4['k2']))}) = {g3(float(c4['k3']))}; k₄ = f(1,2; 2 + 0,2·{g3(float(c4['k3']))}) = {g3(float(c4['k4']))}; suma = {g3(float(c4['s']))}; y₁ = 2 + 0,2·{g3(float(c4['s']))}/6 = {g3(float(c4['y']))}.")
  + TBL(["Método", "y(1,2) aprox", "y real", "E local", "E global"], [
      [m, g3(C_[m]["filas"][1][2]), c(C_[m]["filas"][1][1], 6), g3(C_[m]["filas"][1][3]), g3(C_[m]["filas"][1][4])] for m in ("RK2", "RK4")])
  + P(f"Los dos dan {g3(C_['RK2']['filas'][1][2])}: con 3 dígitos no se distinguen. Cuadro completo hasta t = 2 (verif_u4.py): RK2/RK4/Heun terminan en {g3(C_['RK2']['filas'][-1][2])} (exacta {c(C_['RK2']['filas'][-1][1],4)}); AB2 en {g3(C_['AB2']['filas'][-1][2])}.")
  + TRAMPA("Predicción: un paso de RK2 y RK4 con 3 dígitos para una f no lineal en y. Error típico: en RK4 evaluar k₂ con y + h·k₁ (paso entero) en vez de y + (h/2)·k₁.", real=False),
  practica=True)

hc = C_["Heun"]["det"]; abc = C_["AB2"]["det"]
d(22, None, "TP4 ej. 1c, Heun y AB2: y′ = 3t/y, y(1) = 2, h = 0,2, aritmética de 3 dígitos. a) Un paso de Heun (P: Euler; C: Trapecio iterado con tolerancia 10⁻²). b) Los dos primeros pasos de AB2 (arranque con Euler).",
  P(f"<b>a) Heun.</b> P: y = 2 + 0,2·(3·1/2) = 2 + 0,3 = {g3(float(hc[0][0][1]))}. C1: y = 2 + 0,1·(1,5 + 3·1,2/2,3) = 2 + 0,1·(1,5 + 1,57) = {g3(float(hc[0][1][1]))} (cambio {g3(hc[0][1][2])}: no es menor que 10⁻², sigue). C2: y = 2 + 0,1·(1,5 + 3,6/2,31) = {g3(float(hc[0][2][1]))} (cambio 0) → y(1,2) ≈ {g3(C_['Heun']['filas'][1][2])}; EL = {g3(C_['Heun']['filas'][1][3])}, EG = {g3(C_['Heun']['filas'][1][4])}.")
  + P(f"<b>b) AB2.</b> Arranque (Euler): y₁ = 2 + 0,2·1,5 = {g3(C_['AB2']['filas'][1][2])} (EG = {g3(C_['AB2']['filas'][1][4])}). Paso AB2: f(1,2; 2,3) = {g3(float(abc[1]['f']))}; 3f₁ − f₀ = {g3(float(abc[1]['tres']))} − 1,5 = {g3(float(abc[1]['dif']))}; y₂ = 2,3 + 0,1·{g3(float(abc[1]['dif']))} = {g3(C_['AB2']['filas'][2][2])} (exacta {c(C_['AB2']['filas'][2][1],4)}, EG = {g3(C_['AB2']['filas'][2][4])}).")
  + TRAMPA("Predicción: mostrar el criterio de corte del corrector (|Yant − y| contra la tolerancia) y el arranque de AB2. Error típico: cortar Heun después del predictor (eso es Euler) o arrancar AB2 con el valor exacto.", real=False),
  practica=True)

Dd = R["tp_d"]; d2_ = Dd["RK2"]["det"][0]; d4_ = Dd["RK4"]["det"][0]
d(22, None, "TP4 ej. 1d (bacterias): el número de bacterias crece a una tasa de 3000e^(2t/5) por unidad de tiempo; en t = 0 hay 7.500. Planteá el PVI, obtené la exacta y hacé un paso (h = 0,2, 3 dígitos) con RK2 y con RK4. ¿Qué tiene de particular este PVI?",
  P(f"<b>PVI:</b> P′(t) = 3000e^(2t/5), P(0) = 7500, 0 ≤ t ≤ 1. <b>Exacta</b> (sympy en verif_u4.py): P(t) = {R['tp_d_sol'].replace('*exp(2*t/5)','·e^(2t/5)')}; P(1) = {c(R['tp_d_exacta_1'],2)}.")
  + P("<b>Particularidad:</b> f depende sólo de t (no de P). Entonces RK2 Euler modificado es la regla del trapecio aplicada a ∫f dt y RK4 es la regla de Simpson (k₂ = k₃); el problema es una integral, no hay propagación de error a través de f.")
  + P(f"<b>RK2:</b> k₁ = 3000·e⁰ = {g3(float(d2_['k1']))}; k₂ = 3000·e^(0,4·0,2) = 3000·e^0,08 = 3000·1,08 = {g3(float(d2_['k2']))}; P₁ = 7500 + 0,1·(3000 + 3240) = {g3(float(d2_['y']))}.")
  + P(f"<b>RK4:</b> k₁ = {g3(float(d4_['k1']))}; k₂ = k₃ = 3000·e^0,04 = 3000·1,04 = {g3(float(d4_['k2']))}; k₄ = {g3(float(d4_['k4']))}; suma = {g3(float(d4_['s']))}; P₁ = 7500 + 0,2·{g3(float(d4_['s']))}/6 = {g3(float(d4_['y']))}.")
  + P(f"Exacta en 0,2: {c(Dd['RK2']['filas'][1][1],2)} → EG = {g3(Dd['RK2']['filas'][1][4])} en los dos (EL = {g3(Dd['RK2']['filas'][1][3])}). Con 3 cifras significativas, a partir de 10.000 la resolución es de 100 bacterias: en t = 1 todos dan {g3(Dd['RK4']['filas'][-1][2])} contra {c(R['tp_d_exacta_1'],1)} (EG {g3(Dd['RK4']['filas'][-1][4])}), error dominado por el redondeo.")
  + TRAMPA("Predicción: plantear el PVI a partir del enunciado en palabras (tasa = P′, condición inicial) y reconocer que con f(t) los RK se reducen a reglas de integración.", real=False),
  practica=True)

B = R["tp_b"]
d(18, None, "TP4 ej. 1b (stiff): y′ = 5e^(5t)(y − t)² + 1, −1 ≤ t ≤ 0, y(−1) = −1 − e⁵ ≈ −149,413, exacta y = t − e^(−5t). ¿Qué pasa con h = 0,2 y con h chico? ¿Qué se concluye sobre la estabilidad de la ecuación y de los métodos?",
  P(f"Corridas con precisión completa (verif_u4.py; con 3 dígitos es peor: el primer paso de Euler da {g3(R['tp_b_d_euler1'])} contra la exacta {c(V.yb_ex(-0.8),2)}). Valor final y(0) (exacta −1):")
  + TBL(["Método", "h = 0,2", "h = 0,1", "h = 0,01", "h = 0,001"], [
      [m] + [gs(B[hh][m][0], 4) for hh in (0.2, 0.1, 0.01, 0.001)] for m in ("Euler", "RK2", "RK4", "AB2", "Heun")])
  + P(f"Con h = 0,2, en la grilla (exacta: {', '.join(c(v,1) for v in R['tp_b_ex'])}): Euler da {', '.join(c(v,1) for v in B[0.2]['Euler'][2][:3])}…: el primer paso salta de {c(B[0.2]['Euler'][2][0],1)} a {c(B[0.2]['Euler'][2][1],1)} (f(−1, y₀) = 5e⁵ + 1 = {c(V.fb(-1.0, -1 - math.exp(5)),1)}, la pendiente inicial es enorme) y el error en t = −0,8 es {c(abs(B[0.2]['Euler'][2][1] - V.yb_ex(-0.8)),1)}. RK4: {', '.join(c(v,1) for v in B[0.2]['RK4'][2][:4])}… (cerca). AB2 explota ({gs(B[0.2]['AB2'][0],3)} en t = 0).")
  + P("<b>Conclusiones:</b>")
  + UL("La ecuación es <b>stiff</b>: la exacta tiene el término transitorio e^(−5t), que en [−1, 0] va de e⁵ = " + c(math.exp(5),1) + " a 1: cambia muy rápido al principio y lento al final.",
       f"Con <b>h = 0,2</b> los métodos explícitos no siguen la solución: los errores en los primeros pasos son de decenas; AB2 es inestable y diverge. Sobre la solución exacta, ∂f/∂y = {c(float(R['tp_b_dfdy']),0)}: con h = 0,2 queda h·∂f/∂y = {c(0.2*float(R['tp_b_dfdy']),0)}, justo en el borde de la región de estabilidad de Euler y RK2 y fuera de la de AB2 (agregado).",
       "Con <b>h chico</b> (0,01; 0,001) todos convergen a −1: los métodos son consistentes y, con h suficientemente chico, estables → convergen (EDO3 p4).",
       "Estrategias (EDO4 p12): achicar h, paso variable, o un método implícito (región de estabilidad grande).")
  + AGR("La región de estabilidad (para y′ = λy: Euler |1 + hλ| ≤ 1, RK2 |1 + hλ + (hλ)²/2| ≤ 1, AB2 sólo −1 &lt; hλ &lt; 0) no está en las diapositivas; se usa para explicar por qué AB2 explota con hλ = −2 y RK4 (región hasta ≈ −2,8) no.")
  + TRAMPA("Predicción: el TP pide “¿qué puede concluir respecto a la estabilidad de las ecuaciones y a los métodos aplicados?” y “¿converge? ¿por qué?”. La respuesta tiene que separar estabilidad del PVI (stiff, transitorio) y del método (depende de h), y decir que achicar h la recupera.", real=False),
  practica=True, fig=fig_tp_b())

d(10, None, "Euler implícito, completar la tabla de la cátedra (EDO2 p24, “COMPLETAR LA TABLA!”): y′ = y − t² + 1, 0 ≤ t ≤ 1, y(0) = 0,5, h = 0,2, exacta (t + 1)² − 0,5eᵗ. Usá el despeje yₙ₊₁ = (yₙ − h·t²ₙ₊₁ + h)/(1 − h) y compará con Euler explícito.",
  P("Fórmula (despeje de la cátedra): <code>yₙ₊₁ = (yₙ − h·t²ₙ₊₁ + h)/(1 − h)</code>. Primer paso: y₁ = (0,5 − 0,2·0,04 + 0,2)/0,8 = 0,692/0,8 = " + c(ei[1][2], 3) + ".")
  + tabla_err(ei, 8, aprox=lambda v: c(v, 6), err=lambda v: c(v, 6))
  + P(f"Comparación con Euler explícito (EG en t = 1): implícito {c(ei[-1][4],4)}, explícito {c(eu[-1][4],4)}. Son del mismo orden (1) y errores parecidos: el implícito sobreestima (usa la pendiente al final del intervalo, que es mayor) y el explícito subestima. Su ventaja no es la precisión sino la estabilidad (región más amplia, §10).")
  + TRAMPA("Predicción: despejar yₙ₊₁ de un implícito con f lineal en y y llenar la tabla. Error típico: olvidarse de dividir por (1 − h).", real=False),
  practica=True)

d(20, None, "TP4 ej. 2d (orden superior): y″ − 2y′ + 2y = e^(2t)·sen t, 0 ≤ t ≤ 1, y(0) = −0,4, y′(0) = −0,6, exacta y = (1/5)e^(2t)(sen t − 2cos t). Llevala a un sistema de primer orden y hacé un paso con h = 0,2 con Euler y con RK2 (Euler modificado). Calculá el EL con norma infinito y el EG de y.",
  P("<b>Sistema</b> (y₁ = y, y₂ = y′): y₁′ = y₂; y₂′ = e^(2t)·sen t − 2y₁ + 2y₂; y₁(0) = −0,4, y₂(0) = −0,6.")
  + P(f"F(0, u₀) = (−0,6; e⁰·0 − 2(−0,4) + 2(−0,6)) = ({c(sis2['f0'][0],1)}; {c(sis2['f0'][1],1)}).")
  + P(f"<b>Euler:</b> u₁ = u₀ + 0,2·F = ({c(sis2['euler1'][0],2)}; {c(sis2['euler1'][1],2)}). EL = ‖u₁ − u₀‖∞ = {c(sis2['el_euler'],3)}.")
  + P(f"<b>RK2:</b> k₁ = F(0, u₀) = ({c(sis2['k1'][0],1)}; {c(sis2['k1'][1],1)}); u₀ + h·k₁ = ({c(sis2['euler1'][0],2)}; {c(sis2['euler1'][1],2)}); k₂ = F(0,2; esa predicción) = ({c(sis2['k2'][0],2)}; {c(sis2['k2'][1],4)}); u₁ = u₀ + 0,1·(k₁ + k₂) = ({c(sis2['rk2_1'][0],3)}; {c(sis2['rk2_1'][1],4)}). EL = {c(sis2['el_rk2'],3)}.")
  + P(f"<b>Exacta</b> en t = 0,2: y = {c(sis2['ex1'],4)}, y′ = {c(sis2['exp1'],4)} (verificada con sympy). EG de y: Euler {c(abs(sis2['euler1'][0]-sis2['ex1']),4)}; RK2 {c(abs(sis2['rk2_1'][0]-sis2['ex1']),4)}.")
  + P("(Cuentas en precisión completa: el ejercicio 2 del TP es para scipy.integrate y no pide 3 dígitos.)")
  + TRAMPA("Predicción: el proceso de WA25 aplicado a un caso concreto. Error típico: actualizar y₁ y usar ese valor nuevo para y₂ en el mismo paso (las ecuaciones se resuelven sincrónicamente con los valores del paso anterior).", real=False),
  practica=True)

d(6, None, "Verificá por Taylor el orden de AB2, yₙ₊₁ = yₙ + (h/2)[3f(xₙ, yₙ) − f(xₙ₋₁, yₙ₋₁)], y decí si es consistente.",
  P("Se reemplaza la solución exacta en la fórmula y se desarrolla alrededor de xₙ (y⁽ᵏ⁾ = y⁽ᵏ⁾(xₙ), f = y′):")
  + CMD("y(xₙ₊₁)  = y + h·y′ + (h²/2)·y″ + (h³/6)·y‴ + …\ny′(xₙ₋₁) = y′ − h·y″ + (h²/2)·y‴ − …\n(h/2)[3y′(xₙ) − y′(xₙ₋₁)] = (h/2)[2y′ + h·y″ − (h²/2)·y‴ + …] = h·y′ + (h²/2)·y″ − (h³/4)·y‴ + …\n\nLTE = y(xₙ₊₁) − y(xₙ) − (h/2)[3y′ₙ − y′ₙ₋₁] = (1/6 + 1/4)·h³·y‴ + … = (5/12)·h³·y‴ + …")
  + P(f"Coinciden los términos en h y h² → error local O(h³) → <b>orden 2</b>, error global O(h²) (como dice la cátedra en EDO3 p12). Es <b>consistente</b> porque coincide con la ST al menos en y + h·y′. Coincide con verif_u4.py (coeficiente {M['AB2']['c0'].replace('*Y3','·y‴')}).")
  + TRAMPA("Predicción: justificar el inciso b) “¿de qué orden es? ¿por qué?” con el desarrollo de Taylor, no sólo con el nombre del método.", real=False))

# ---------------------------------------------------------------- UNIDAD
unidad = {
    "id": 4,
    "titulo": "Resolución de Ecuaciones Diferenciales Ordinarias",
    "corto": "EDO",
    "fuente": "edo1.pdf (17 diap.), edo2.pdf (30), edo3.pdf (18), edo4.pdf (18) y tp4.pdf (3 pág.)",
    "peso": 3,
    "intro": P("Un PVI y′ = f(x, y), y(x₀) = y₀ casi nunca se resuelve analíticamente: los métodos numéricos dan una secuencia de puntos (xₙ, yₙ) sobre una grilla de paso h. La unidad recorre los métodos de un paso (Taylor, Euler explícito e implícito, RK2, RK4, RKF45), los multipaso (Adams-Bashforth y Adams-Moulton), los predictor-corrector (Heun, Milne) y los criterios para evaluarlos (orden, consistencia, convergencia, estabilidad), más stiff, sistemas y orden superior.")
             + P("<b>Es la unidad que más pesa:</b> en los cuatro temas del integral 2025 la pregunta 4 da un método y pide los mismos 4 incisos (¿RK? ¿orden? ¿un paso o multipaso? ¿implícito o explícito?, cada uno con su porqué) y la 5 pide stiff, estabilidad, consistencia o convergencia y su relación con los errores. Los finales piden cuadros comparativos (Trapecio, Euler, AB2, Heun, RK4, RKF). La práctica de EDO no tiene examen en el material: sólo el TP4 (h = 0,2, 3 dígitos, RK2, RK4, AB2, Heun)."),
    "secciones": S,
    "conexiones": UL("Con U2 (ajuste de curvas): los métodos de Adams integran un <b>polinomio interpolante</b> por los puntos anteriores; el TP4 compara aproximación y valor real con errores, como el ECM del ajuste.",
                     "Con U3 (Fourier): la señal sismológica de EDO3 p6 (variabilidad alta → h chico) es la misma idea que Nyquist: muestrear más fino donde la señal cambia rápido.",
                     "Con integración numérica: si f depende sólo de x, RK2 Euler modificado es la regla del trapecio y RK4 es Simpson (TP4 1d)."),
    "fidelidad": {
        "fuentes": "Las 83 diapositivas de edo1–edo4 renderizadas en hojas 2×2 y leídas (zoom en EDO2 p27, EDO3 p8 y EDO4 p7); el TP4 se renderizó y se leyó (el texto extraído tenía las fórmulas rotas). Inventario en calibracion/01-inventario-u4.md. Todos los números (tablas, órdenes por Taylor, órdenes empíricos, cuentas de 3 dígitos, stiff, sistemas, figuras) salen de tools/gen/verif_u4.py.",
        "agregado": "Nombres BDF2/BDF3 (IT25-T2/T3), punto medio de 2 pasos / leapfrog (IT25-T4 lectura 2h) y etapas RK3 de Kutta y Lobatto IIIA (IT25-T1); órdenes de AM2, AM3, AB3, AB4 y Milne por Taylor (la cátedra no los da); región de estabilidad y factores de amplificación para y′ = λy; teorema de Dahlquist (convergente ⇒ cero-estable) y la aclaración de que Euler implícito también es A-estable; diferencia entre el E local de las tablas (|yᵢ₊₁ − yᵢ|) y el error de truncamiento local de la bibliografía; correcciones de erratas (β = 2 → 1, coeficientes de RKF45 25/216 y 2197/4104, signo de Euler en el sistema, 3y₂ → 3y₃); evaluaciones de f por paso en AB2.",
        "huecos": "Ninguna diapositiva ilegible. Sin objetivos oficiales (no hay diapositiva de objetivos). EDO3 p9 (ejemplo RKF45) es una tabla tomada de la bibliografía: se transcribieron sólo los valores finales. No hay examen práctico de EDO en el material: toda la práctica es predicción sobre el TP4.",
    },
    "auditoria": auditoria,
    "glosario": glosario,
    "quiz": Q,
    "desarrollo": DES,
}

# Ampliación: ejercicios de práctica del Burden (cap. 5), ver burden_u4.py
from burden_u4 import items as _burden_items
unidad["desarrollo"] = unidad["desarrollo"] + _burden_items()
src = "ASO.addUnit(" + json.dumps(unidad, ensure_ascii=False, indent=1) + ");\n"
open(OUT, "w", encoding="utf-8").write(src)

# métricas rápidas de balance
import statistics
larga = 0; props = []; lc = []; li = []
for x in Q:
    ops = x["opciones"]; corr = [o for o in ops if o["v"]]
    m = max(ops, key=lambda o: len(o["t"])); larga += m in corr
    props.append(len(corr) / len(ops))
    lc += [len(o["t"]) for o in corr]; li += [len(o["t"]) for o in ops if o not in corr]
print(f"escrito {OUT}: {len(S)} secciones, {len(Q)} quiz, {len(DES)} desarrollo")
print(f"balance: más larga correcta {100*larga/len(Q):.0f}% (esperado {100*statistics.mean(props):.0f}%), ratio {statistics.mean(lc)/statistics.mean(li):.2f}, verdaderas {sum(x['pide']=='verdaderas' for x in Q)}")
for x in Q:
    ops = x["opciones"]; m = max(ops, key=lambda o: len(o["t"]))
    print(f"  {x['id']} {x['pide'][:4]} larga_ok={m['v']} lens={[len(o['t']) for o in ops]} v={[int(o['v']) for o in ops]}")
