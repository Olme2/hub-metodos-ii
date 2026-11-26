"""Ejercicios de práctica de U4 tomados del Burden-Faires, Análisis numérico, 10.ª ed., cap. 5.
Amplía la unidad sin tocar lo aprobado: gen_u4.py llama a items() y los agrega al final de `desarrollo`.
Todo número sale de estas funciones (precisión doble). Las soluciones exactas están verificadas con sympy
en burden_verif.py. Convenciones de la cátedra (EDO2 p7, EDO3 p12, EDO4 p4–6):
  E local = |y_{i+1} − y_i|, E global = |y(t_i) − y_i|; RK2 = Euler modificado; AB2 arranca con Euler;
  Heun = predictor Euler + corrector Trapecio iterado mientras |dif| >= 1e-2, itMax 5.
Uso suelto: python3 burden_u4.py   (imprime todas las tablas para revisar)
"""
import math
from figs import nueva, svg

def fmt(x, d=6):
    return f"{x:.{d}f}".replace(".", ",")

# ---------- métodos (convención de la cátedra) ----------
def euler(f, t0, y0, h, n):
    T, Y, det = [t0], [y0], []
    for i in range(n):
        k = f(T[-1], Y[-1]); Y.append(Y[-1] + h*k); T.append(T[-1] + h); det.append({"f": k})
    return T, Y, det

def rk2(f, t0, y0, h, n):
    T, Y, det = [t0], [y0], []
    for i in range(n):
        t, y = T[-1], Y[-1]
        k1 = f(t, y); k2 = f(t + h, y + h*k1)
        Y.append(y + h/2*(k1 + k2)); T.append(t + h); det.append({"k1": k1, "k2": k2})
    return T, Y, det

def punto_medio(f, t0, y0, h, n):
    T, Y, det = [t0], [y0], []
    for i in range(n):
        t, y = T[-1], Y[-1]
        k1 = f(t, y); k2 = f(t + h/2, y + h/2*k1)
        Y.append(y + h*k2); T.append(t + h); det.append({"k1": k1, "k2": k2})
    return T, Y, det

def rk4(f, t0, y0, h, n):
    T, Y, det = [t0], [y0], []
    for i in range(n):
        t, y = T[-1], Y[-1]
        k1 = f(t, y); k2 = f(t + h/2, y + h/2*k1); k3 = f(t + h/2, y + h/2*k2); k4 = f(t + h, y + h*k3)
        Y.append(y + h/6*(k1 + 2*k2 + 2*k3 + k4)); T.append(t + h); det.append({"k1": k1, "k2": k2, "k3": k3, "k4": k4})
    return T, Y, det

def ab2(f, t0, y0, h, n, arranque="euler"):
    T, Y = [t0], [y0]
    y1 = y0 + h*f(t0, y0) if arranque == "euler" else rk4(f, t0, y0, h, 1)[1][1]
    T.append(t0 + h); Y.append(y1); det = [{"arranque": arranque}]
    for i in range(1, n):
        fn, fm = f(T[i], Y[i]), f(T[i-1], Y[i-1])
        Y.append(Y[i] + h/2*(3*fn - fm)); T.append(T[i] + h); det.append({"fn": fn, "fn-1": fm})
    return T, Y, det

def heun(f, t0, y0, h, n, tol=1e-2, itmax=5):
    T, Y, det = [t0], [y0], []
    for i in range(n):
        t, y = T[-1], Y[-1]
        fn = f(t, y); yp = y + h*fn; tr = [("P", yp)]
        ya = yp
        for it in range(itmax):
            yc = y + h/2*(fn + f(t + h, ya)); d = abs(yc - ya); tr.append(("C", yc, d)); ya = yc
            if d < tol: break
        Y.append(ya); T.append(t + h); det.append({"traza": tr})
    return T, Y, det

# ---------- presentación ----------
def tabla(T, Y, real):
    filas = []
    for i, (t, y) in enumerate(zip(T, Y)):
        yr = real(t)
        el = "—" if i == 0 else fmt(abs(Y[i] - Y[i-1]))
        eg = "—" if i == 0 else fmt(abs(yr - y))
        filas.append(f"<tr><td>{fmt(t, 2)}</td><td>{fmt(yr)}</td><td>{fmt(y)}</td><td>{el}</td><td>{eg}</td></tr>")
    return ('<div class="tbl"><table><thead><tr><th>t</th><th>y real</th><th>y aprox</th>'
            '<th>E local = |yᵢ₊₁ − yᵢ|</th><th>E global = |y(tᵢ) − yᵢ|</th></tr></thead><tbody>'
            + "".join(filas) + "</tbody></table></div>")

FUENTE = ("Burden y Faires, <i>Análisis numérico</i>, 10.ª ed., cap. 5 (el libro que pasó la cátedra). "
          "El Burden publica respuestas sólo de los ejercicios impares y sin las columnas de error; "
          "la tabla de acá está calculada en precisión doble con la convención de la cátedra.")

def PRED(txt):
    return f'<div class="trampa"><b>Así lo tomaría la cátedra</b><p>{txt}</p></div>'

def AGR(txt):
    return f'<div class="agregado"><b>Agregado</b><p>{txt}</p></div>'

def item(n, seccion, pregunta, modelo, fig=None):
    d = {"id": f"m4-b{n:02d}", "seccion": seccion, "real": None, "practica": True,
         "pregunta": pregunta, "modelo": modelo}
    if fig: d["fig"] = fig
    return d

# ---------- PVI (soluciones exactas verificadas en burden_verif.py) ----------
E, L, C, S, Q = math.exp, math.log, math.cos, math.sin, math.sqrt
pvi = {
 "a": (lambda t, y: t*E(3*t) - 2*y, lambda t: t*E(3*t)/5 - E(3*t)/25 + E(-2*t)/25,
       "y′ = t·e³ᵗ − 2y, 0 ≤ t ≤ 1, y(0) = 0", "y(t) = (1/5)·t·e³ᵗ − (1/25)·e³ᵗ + (1/25)·e⁻²ᵗ"),
 "b": (lambda t, y: 1 + (t - y)**2, lambda t: t + 1/(1 - t),
       "y′ = 1 + (t − y)², 2 ≤ t ≤ 3, y(2) = 1", "y(t) = t + 1/(1 − t)"),
 "c": (lambda t, y: 1 + y/t, lambda t: t*L(t) + 2*t,
       "y′ = 1 + y/t, 1 ≤ t ≤ 2, y(1) = 2", "y(t) = t·ln t + 2t"),
 "d": (lambda t, y: C(2*t) + S(3*t), lambda t: S(2*t)/2 - C(3*t)/3 + 4/3,
       "y′ = cos 2t + sen 3t, 0 ≤ t ≤ 1, y(0) = 1", "y(t) = (1/2)·sen 2t − (1/3)·cos 3t + 4/3"),
 "3d": (lambda t, y: -t*y + 4*t/y, lambda t: Q(4 - 3*E(-t*t)),
        "y′ = −t·y + 4t/y, 0 ≤ t ≤ 1, y(0) = 1", "y(t) = √(4 − 3·e^(−t²))"),
}

def items():
    out = []
    # B1 Euler 5.2-1a
    f, r, enun, sol = pvi["a"]; T, Y, det = euler(f, 0, 0, 0.5, 2)
    out.append(item(1, "m4-s9",
        f"Burden §5.2 ej. 1a: resolvé con Euler explícito {enun}, h = 0,5. Exacta: {sol}. Completá el cuadro t | y real | y aprox | E local | E global.",
        f"<p><b>Método (EDO2 p13):</b> yₙ₊₁ = yₙ + h·f(tₙ, yₙ).</p>"
        f"<p><b>Paso 1</b> (t = 0 → 0,5): f(0; 0) = 0·e⁰ − 2·0 = {fmt(det[0]['f'])} → y₁ = 0 + 0,5·0 = {fmt(Y[1])}.</p>"
        f"<p><b>Paso 2</b> (0,5 → 1): f(0,5; {fmt(Y[1])}) = 0,5·e^1,5 − 2·{fmt(Y[1])} = {fmt(det[1]['f'])} → y₂ = {fmt(Y[1])} + 0,5·{fmt(det[1]['f'])} = {fmt(Y[2])}.</p>"
        + tabla(T, Y, r) +
        f"<p>El error global en t = 1 es grande ({fmt(abs(r(1) - Y[2]), 4)} sobre {fmt(r(1), 4)}): orden 1 con h = 0,5. Achicando h el EG baja aproximadamente en proporción a h.</p>"
        + PRED("Dos pasos de Euler con h grande y la pregunta \"¿es un buen resultado? ¿qué harías?\": la respuesta esperada es achicar h o cambiar a un método de mayor orden, mirando el EG (EDO2 p21).")
        + f"<p class=\"src\">{FUENTE}</p>"))
    # B2 Euler 5.2-1d
    f, r, enun, sol = pvi["d"]; T, Y, det = euler(f, 0, 1, 0.25, 4)
    out.append(item(2, "m4-s9",
        f"Burden §5.2 ej. 1d: Euler explícito para {enun}, h = 0,25. Exacta: {sol}. Cuadro completo.",
        f"<p>Acá f depende sólo de t: Euler es la regla del rectángulo izquierdo aplicada a ∫ f(t) dt.</p>"
        f"<p><b>Paso 1:</b> f(0) = cos 0 + sen 0 = {fmt(det[0]['f'])} → y₁ = 1 + 0,25·1 = {fmt(Y[1])}.</p>"
        f"<p><b>Paso 2:</b> f(0,25) = cos 0,5 + sen 0,75 = {fmt(det[1]['f'])} → y₂ = {fmt(Y[2])}.</p>"
        + tabla(T, Y, r) + PRED("Mismo esquema con una f que depende sólo de t: ojo con evaluar en radianes.")
        + f"<p class=\"src\">{FUENTE}</p>"))
    # B3 RK2 5.4-1c
    f, r, enun, sol = pvi["c"]; T, Y, det = rk2(f, 1, 2, 0.25, 4)
    out.append(item(3, "m4-s11",
        f"Burden §5.4 ej. 1c: RK2 (Euler modificado de la cátedra) para {enun}, h = 0,25. Exacta: {sol}. Es el mismo PVI del TP4 1a con otro h.",
        f"<p><b>Método (EDO2 p27):</b> k₁ = f(tₙ, yₙ), k₂ = f(tₙ + h, yₙ + h·k₁), yₙ₊₁ = yₙ + (h/2)(k₁ + k₂). El Burden lo llama \"método de Euler modificado\": es el mismo.</p>"
        f"<p><b>Paso 1:</b> k₁ = 1 + 2/1 = {fmt(det[0]['k1'])}; yₙ + h·k₁ = 2 + 0,25·3 = 2,75; k₂ = 1 + 2,75/1,25 = {fmt(det[0]['k2'])}; y₁ = 2 + 0,125·({fmt(det[0]['k1'])} + {fmt(det[0]['k2'])}) = {fmt(Y[1])}.</p>"
        + tabla(T, Y, r) +
        f"<p>Compará con TP4 1a (h = 0,2): el EG en t = 2 de acá ({fmt(abs(r(2) - Y[-1]), 5)}) es mayor porque h es mayor; con orden 2 el EG escala aprox. con h².</p>"
        + PRED("Un paso de RK2 escribiendo k₁ y k₂ explícitamente: se corrige que uses tₙ + h (no tₙ + h/2, que es punto medio).")
        + f"<p class=\"src\">{FUENTE}</p>"))
    # B4 RK2 5.4-1b
    f, r, enun, sol = pvi["b"]; T, Y, det = rk2(f, 2, 1, 0.5, 2)
    out.append(item(4, "m4-s11",
        f"Burden §5.4 ej. 1b: RK2 (Euler modificado) para {enun}, h = 0,5. Exacta: {sol}.",
        f"<p><b>Paso 1</b> (t = 2 → 2,5): k₁ = 1 + (2 − 1)² = {fmt(det[0]['k1'])}; y + h·k₁ = 1 + 0,5·2 = 2; k₂ = f(2,5; 2) = 1 + 0,5² = {fmt(det[0]['k2'])}; y₁ = 1 + 0,25·(2 + 1,25) = {fmt(Y[1])}.</p>"
        f"<p><b>Paso 2</b> (2,5 → 3): k₁ = {fmt(det[1]['k1'])}; k₂ = {fmt(det[1]['k2'])}; y₂ = {fmt(Y[2])}.</p>"
        + tabla(T, Y, r) + PRED("Dos pasos con números redondos: es el tipo de cuenta que entra en un parcial a mano.")
        + f"<p class=\"src\">{FUENTE}</p>"))
    # B5 punto medio 5.4-5 (sobre 1d)
    f, r, enun, sol = pvi["d"]; T, Y, det = punto_medio(f, 0, 1, 0.25, 4)
    T2, Y2, _ = euler(f, 0, 1, 0.25, 4); T3, Y3, _ = rk2(f, 0, 1, 0.25, 4)
    out.append(item(5, "m4-s11",
        f"Burden §5.4 ej. 5 (repetir el ej. 1 con punto medio), caso d: {enun}, h = 0,25, con el método del punto medio. Exacta: {sol}. Compará el EG final con Euler y con Euler modificado.",
        f"<p><b>Método (EDO2 p28):</b> yₙ₊₁ = yₙ + h·f(tₙ + h/2, yₙ + (h/2)·f(tₙ, yₙ)).</p>"
        f"<p><b>Paso 1:</b> k₁ = f(0) = {fmt(det[0]['k1'])}; k₂ = f(0,125) = cos 0,25 + sen 0,375 = {fmt(det[0]['k2'])}; y₁ = 1 + 0,25·{fmt(det[0]['k2'])} = {fmt(Y[1])}.</p>"
        + tabla(T, Y, r) +
        f"<p><b>EG en t = 1:</b> Euler {fmt(abs(r(1) - Y2[-1]), 5)} · Euler modificado {fmt(abs(r(1) - Y3[-1]), 5)} · punto medio {fmt(abs(r(1) - Y[-1]), 5)}. Los dos RK2 son de orden 2: con el mismo h su EG final es entre 3 y 6 veces menor que el de Euler.</p>"
        + PRED("\"¿Punto medio y Euler modificado son el mismo método?\": no; ambos son RK2 (orden 2), con distintas constantes γ₁, γ₂, α, β (EDO2 p26–28).")
        + f"<p class=\"src\">{FUENTE}</p>"))
    # B6 RK4 sobre 5.4-1a
    f, r, enun, sol = pvi["a"]; T, Y, det = rk4(f, 0, 0, 0.5, 2); _, Ye, _ = euler(f, 0, 0, 0.5, 2)
    k = det[0]
    out.append(item(6, "m4-s12",
        f"PVI del Burden §5.4 ej. 1a resuelto con RK4: {enun}, h = 0,5. Exacta: {sol}. Calculá k₁…k₄ del primer paso y el cuadro.",
        f"<p><b>Método (EDO2 p29):</b> k₁ = f(tₙ, yₙ), k₂ = f(tₙ + h/2, yₙ + (h/2)k₁), k₃ = f(tₙ + h/2, yₙ + (h/2)k₂), k₄ = f(tₙ + h, yₙ + h·k₃), yₙ₊₁ = yₙ + (h/6)(k₁ + 2k₂ + 2k₃ + k₄).</p>"
        f"<p><b>Paso 1:</b> k₁ = {fmt(k['k1'])}; k₂ = {fmt(k['k2'])}; k₃ = {fmt(k['k3'])}; k₄ = {fmt(k['k4'])}; y₁ = 0 + (0,5/6)·({fmt(k['k1'])} + 2·{fmt(k['k2'])} + 2·{fmt(k['k3'])} + {fmt(k['k4'])}) = {fmt(Y[1])}.</p>"
        + tabla(T, Y, r) +
        f"<p>Con el mismo h = 0,5, Euler da EG = {fmt(abs(r(1) - Ye[-1]), 4)} en t = 1 y RK4 {fmt(abs(r(1) - Y[-1]), 4)}: con un h grande RK4 todavía sirve (EDO2 p30: \"si aumenta el orden, h no necesariamente debe ser muy chico\"), a costa de 4 evaluaciones de f por paso.</p>"
        + PRED("Un paso de RK4 con los cuatro k: el error típico es usar yₙ + h·k₂ en k₃ (va h/2).")
        + f"<p class=\"src\">{FUENTE}</p>"))
    # B7 RK4 5.4-3d
    f, r, enun, sol = pvi["3d"]; T, Y, det = rk4(f, 0, 1, 0.1, 2); k = det[0]
    out.append(item(7, "m4-s12",
        f"PVI del Burden §5.4 ej. 3d con RK4: {enun}, h = 0,1, dos pasos. Exacta: {sol}.",
        f"<p><b>Paso 1:</b> k₁ = f(0; 1) = 0 → {fmt(k['k1'])}; k₂ = f(0,05; 1) = {fmt(k['k2'])}; k₃ = {fmt(k['k3'])}; k₄ = {fmt(k['k4'])}; y₁ = {fmt(Y[1])}.</p>"
        + tabla(T, Y, r) +
        "<p>El EG queda por debajo de 10⁻⁵: RK4 con h = 0,1 es prácticamente exacto para una solución suave.</p>"
        + AGR("La extracción del PDF del Burden muestra la exacta como 4 − 3e^(−t²) sin la raíz; con y(0) = 1 sólo cumple √(4 − 3e^(−t²)) (verificado con sympy en burden_verif.py).")
        + PRED("Mismo pedido con una f no lineal en y: hay que llevar bien los paréntesis al evaluar 4t/y.")
        + f"<p class=\"src\">{FUENTE}</p>"))
    # B8 AB2 5.6-1a
    f, r, enun, sol = pvi["a"]; T, Y, det = ab2(f, 0, 0, 0.2, 5); Tb, Yb, _ = ab2(f, 0, 0, 0.2, 5, "rk4")
    out.append(item(8, "m4-s16",
        f"Burden §5.6 ej. 1a con Adams-Bashforth de 2 pasos: {enun}, h = 0,2. Exacta: {sol}. Arrancá como la cátedra (y₁ con Euler).",
        f"<p><b>Método (EDO3 p12):</b> yₙ₊₁ = yₙ + (h/2)(3f(tₙ, yₙ) − f(tₙ₋₁, yₙ₋₁)). Es multipaso: para y₁ no hay yₙ₋₁, así que se usa un método de inicio; la cátedra usa Euler (EDO3 p13).</p>"
        f"<p><b>Inicio:</b> y₁ = 0 + 0,2·f(0; 0) = {fmt(Y[1])}.</p>"
        f"<p><b>Paso 2:</b> f(0,2; y₁) = {fmt(det[1]['fn'])}, f(0; 0) = {fmt(det[1]['fn-1'])} → y₂ = {fmt(Y[1])} + 0,1·(3·{fmt(det[1]['fn'])} − {fmt(det[1]['fn-1'])}) = {fmt(Y[2])}.</p>"
        + tabla(T, Y, r) +
        f"<p>El Burden arranca con valores de RK4; con ese arranque y(1) = {fmt(Yb[-1])} (EG {fmt(abs(r(1) - Yb[-1]), 5)}) contra {fmt(Y[-1])} (EG {fmt(abs(r(1) - Y[-1]), 5)}) arrancando con Euler: el error del arranque se arrastra, aunque acá pesa poco: domina el error de AB2 con h = 0,2 sobre una solución que crece como e³ᵗ. En general conviene un arranque de orden igual o mayor que el del método.</p>"
        + PRED("\"¿Por qué AB2 necesita otro método para arrancar?\" y un par de pasos a mano: decir que es multipaso y qué método de inicio usaste.")
        + f"<p class=\"src\">{FUENTE}</p>"))
    # B9 Heun 5.6-1d
    f, r, enun, sol = pvi["d"]; T, Y, det = heun(f, 0, 1, 0.2, 5)
    tr = det[0]["traza"]
    traza = " → ".join([f"P: {fmt(tr[0][1])}"] + [f"C: {fmt(x[1])} (dif {fmt(x[2], 6)})" for x in tr[1:]])
    out.append(item(9, "m4-s17",
        f"PVI del Burden §5.6 ej. 1d con Heun predictor-corrector (como la cátedra: corrector iterado mientras la diferencia sea ≥ 10⁻², itMax 5): {enun}, h = 0,2. Exacta: {sol}.",
        f"<p><b>Método (EDO4 p4):</b> predictor Euler yₙ₊₁⁽⁰⁾ = yₙ + h·f(tₙ, yₙ); corrector Trapecio yₙ₊₁ = yₙ + (h/2)(f(tₙ, yₙ) + f(tₙ₊₁, yₙ₊₁)), repetido con el último valor.</p>"
        f"<p><b>Paso 1:</b> {traza}. Como f no depende de y, el corrector da lo mismo en la segunda vuelta: converge en una corrección.</p>"
        + tabla(T, Y, r) + PRED("Traza de un paso P-C (valor predicho, cada corrección y la diferencia contra la tolerancia) y cuántas iteraciones hicieron falta.")
        + f"<p class=\"src\">{FUENTE}</p>"))
    # B10 stiff 5.11-1b
    fs = lambda t, y: -20*(y - t*t) + 2*t; rs = lambda t: t*t + E(-20*t)/3
    T, Y, _ = euler(fs, 0, 1/3, 0.1, 10); T4, Y4, _ = rk4(fs, 0, 1/3, 0.1, 10)
    Tf, Yf, _ = euler(fs, 0, 1/3, 0.01, 100)
    fig, ax = nueva()
    import numpy as np
    tt = np.linspace(0, 1, 300); ax.plot(tt, [rs(x) for x in tt], color="k", label="exacta")
    ax.plot(T, Y, "o-", label="Euler h = 0,1"); ax.plot(T4, Y4, "s--", label="RK4 h = 0,1")
    ax.plot(Tf[::10], Yf[::10], "^:", label="Euler h = 0,01"); ax.set_xlabel("t"); ax.set_ylabel("y"); ax.legend(fontsize=8)
    filas = "".join(f"<tr><td>{fmt(T[i], 1)}</td><td>{fmt(rs(T[i]))}</td><td>{fmt(Y[i])}</td><td>{fmt(Y4[i])}</td></tr>" for i in range(0, 11, 2))
    out.append(item(10, "m4-s18",
        "Burden §5.11 ej. 1b (stiff): y′ = −20(y − t²) + 2t, 0 ≤ t ≤ 1, y(0) = 1/3, h = 0,1; exacta y = t² + (1/3)e^(−20t). Resolvé con Euler y con RK4, compará con Euler h = 0,01 y explicá qué pasa.",
        "<p><b>Por qué es stiff (EDO4 p9):</b> la solución tiene el término transitorio (1/3)e^(−20t), que muere enseguida, y el resto (t²) cambia lento. Ese término exige h chico aunque ya no aporte nada.</p>"
        '<div class="tbl"><table><thead><tr><th>t</th><th>y real</th><th>Euler h = 0,1</th><th>RK4 h = 0,1</th></tr></thead><tbody>' + filas + "</tbody></table></div>"
        f"<p>La tabla muestra un paso de cada dos. Euler con h = 0,1 queda por encima de la exacta en los pasos pares y por debajo en los impares, con un error de alrededor de 1/3 en valor absoluto que no achica: oscila alrededor de la solución sin amortiguarse; con h = 0,01 sigue bien a la exacta (EG en t = 1: {fmt(abs(rs(1) - Yf[-1]), 5)}). RK4 con h = 0,1 se mantiene acotado (EG en t = 1: {fmt(abs(rs(1) - Y4[-1]), 5)}).</p>"
        + AGR("Explicación con y′ = λy y λ = −20: Euler multiplica el error por (1 + hλ) = 1 − 2 = −1 en cada paso, así que el error del transitorio no decae, cambia de signo. Hace falta |1 + hλ| &lt; 1, o sea h &lt; 0,1. RK4 tiene una región de estabilidad más grande (hλ = −2 cae adentro). Un método implícito (Euler implícito, Trapecio) no tiene esa restricción.")
        + PRED("Con un gráfico así: \"¿qué está pasando? ¿qué estrategia recomendaría?\" — achicar h, paso variable, o un método implícito/de mayor región de estabilidad (EDO4 p12), con sus costos."),
        fig=svg(fig)))
    # B11 sistema 5.9-1a con Euler
    u = lambda t: (E(5*t)/3 - E(-t)/3 + E(2*t), E(5*t)/3 + 2*E(-t)/3 + t*t*E(2*t))
    F = lambda t, a, b: (3*a + 2*b - (2*t*t + 1)*E(2*t), 4*a + b + (t*t + 2*t - 4)*E(2*t))
    U = [(0.0, 1.0, 1.0)]; h = 0.2
    for i in range(2):
        t, a, b = U[-1]; fa, fb = F(t, a, b); U.append((t + h, a + h*fa, b + h*fb))
    f0 = F(0, 1, 1)
    filas = "".join(f"<tr><td>{fmt(t, 1)}</td><td>{fmt(u(t)[0])}</td><td>{fmt(a)}</td><td>{fmt(u(t)[1])}</td><td>{fmt(b)}</td><td>{'—' if i == 0 else fmt(max(abs(a - U[i-1][1]), abs(b - U[i-1][2])))}</td></tr>" for i, (t, a, b) in enumerate(U))
    out.append(item(11, "m4-s19",
        "Burden §5.9 ej. 1a (es el ejemplo de sistemas de la cátedra, EDO4 p14): u₁′ = 3u₁ + 2u₂ − (2t² + 1)e²ᵗ, u₂′ = 4u₁ + u₂ + (t² + 2t − 4)e²ᵗ, u₁(0) = u₂(0) = 1, h = 0,2. Dos pasos de Euler; EL con norma infinito.",
        f"<p>Euler se aplica a cada componente con los valores del paso anterior: uⱼ⁽ⁿ⁺¹⁾ = uⱼ⁽ⁿ⁾ + h·fⱼ(tₙ, u₁⁽ⁿ⁾, u₂⁽ⁿ⁾).</p>"
        f"<p><b>Paso 1:</b> f₁(0; 1, 1) = 3 + 2 − 1 = {fmt(f0[0])}; f₂(0; 1, 1) = 4 + 1 − 4 = {fmt(f0[1])} → u₁ = 1 + 0,2·4 = {fmt(U[1][1])}; u₂ = 1 + 0,2·1 = {fmt(U[1][2])}.</p>"
        '<div class="tbl"><table><thead><tr><th>t</th><th>u₁ real</th><th>u₁ aprox</th><th>u₂ real</th><th>u₂ aprox</th><th>EL = ‖uⁿ⁺¹ − uⁿ‖∞</th></tr></thead><tbody>' + filas + "</tbody></table></div>"
        + AGR("Exactas (verificadas con sympy): u₁ = (1/3)e⁵ᵗ − (1/3)e⁻ᵗ + e²ᵗ, u₂ = (1/3)e⁵ᵗ + (2/3)e⁻ᵗ + t²e²ᵗ. En EDO4 p16 la cátedra escribe u₁¹ = u₁⁰ − h·(…): el signo correcto es +.")
        + PRED("Un paso de un método sobre un sistema de 2 ecuaciones, con la norma infinito para el EL.")))
    # B12 orden superior 5.9-3a
    ry = lambda t: t**3*E(t)/6 - t*E(t) + 2*E(t) - t - 2
    V = [(0.0, 0.0, 0.0)]; h = 0.1
    for i in range(2):
        t, a, b = V[-1]; V.append((t + h, a + h*b, b + h*(2*b - a + t*E(t) - t)))
    filas = "".join(f"<tr><td>{fmt(t, 1)}</td><td>{fmt(ry(t))}</td><td>{fmt(a)}</td><td>{fmt(b)}</td></tr>" for (t, a, b) in V)
    out.append(item(12, "m4-s20",
        "Burden §5.9 ej. 3a: y″ − 2y′ + y = t·eᵗ − t, 0 ≤ t ≤ 1, y(0) = y′(0) = 0, h = 0,1. Escribilo como sistema de primer orden y hacé dos pasos de Euler.",
        "<p><b>Reducción (EDO4 p17):</b> u₁ = y, u₂ = y′ → u₁′ = u₂; u₂′ = y″ = 2u₂ − u₁ + t·eᵗ − t; u₁(0) = 0, u₂(0) = 0.</p>"
        f"<p><b>Paso 1:</b> u₁ = 0 + 0,1·0 = 0; u₂ = 0 + 0,1·(0 − 0 + 0 − 0) = 0. <b>Paso 2:</b> u₁ = {fmt(V[2][1])}; u₂ = 0 + 0,1·(0,1·e^0,1 − 0,1) = {fmt(V[2][2])}.</p>"
        '<div class="tbl"><table><thead><tr><th>t</th><th>y real</th><th>u₁ = y aprox</th><th>u₂ = y′ aprox</th></tr></thead><tbody>' + filas + "</tbody></table></div>"
        "<p>Con h = 0,1 los primeros pasos de Euler todavía no se mueven (la solución arranca muy plana): para ver algo hay que avanzar más o usar RK4. Lo que se evalúa acá es la reducción a sistema.</p>"
        + AGR("Exacta (verificada con sympy): y = (1/6)t³eᵗ − t·eᵗ + 2eᵗ − t − 2.")
        + PRED("WA25: \"proceso para resolver una EDO de orden superior usando cualquier método\": cambio de variables, sistema de primer orden con sus condiciones iniciales, y aplicar el método componente a componente.")))
    return out

if __name__ == "__main__":
    import re
    for d in items():
        print("=" * 70); print(d["id"], d["pregunta"][:110])
        txt = re.sub(r"<tr>", "\n", d["modelo"]); txt = re.sub(r"<[^>]+>", " ", txt)
        print(re.sub(r"[ \t]+", " ", txt)[:1400])
