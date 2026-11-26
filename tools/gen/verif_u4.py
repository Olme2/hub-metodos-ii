#!/usr/bin/env python3
"""Verificación numérica y simbólica de la U4 (EDO).
Todo número que aparece en data/m4.js sale de acá.
Correr:  python3 verif_u4.py            (imprime todo)
El generador gen_u4.py importa este módulo y usa R (diccionario de resultados).
"""
import math
from fractions import Fraction as Fr
import sympy as sp

R = {}          # resultados que usa el generador
SALIDA = []     # líneas para imprimir

def out(*a):
    SALIDA.append(" ".join(str(x) for x in a))

# ---------------------------------------------------------------------------
# 1) ORDEN DE MÉTODOS MULTIPASO (y de un paso escritos como multipaso) POR TAYLOR
#    LTE = y(x_{n+1}) - sum a_j y(x_{n-j}) - h sum b_j y'(x_{n-j})   (b_{-1} = término implícito)
#    Se desarrolla alrededor de x_n con y^{(k)}(x_n) = Yk.
# ---------------------------------------------------------------------------
h = sp.symbols('h', positive=True)
Y = sp.symbols('Y0:9')   # Y[k] = y^{(k)}(x_n)

def taylor_y(s, orden=8):
    """y(x_n + s h) como serie en h."""
    return sum(Y[k] * (s * h) ** k / sp.factorial(k) for k in range(orden))

def taylor_yp(s, orden=8):
    """y'(x_n + s h)."""
    return sum(Y[k + 1] * (s * h) ** k / sp.factorial(k) for k in range(orden - 1))

def lte_multipaso(a, b, bimpl=0, escala=1):
    """a: dict {j: coef de y_{n-j}}, b: dict {j: coef de h*y'_{n-j}}, bimpl: coef de h*y'_{n+1}.
    escala multiplica h en el término de derivadas (para la lectura '2h')."""
    expr = taylor_y(1) - sum(c * taylor_y(-j) for j, c in a.items()) \
        - escala * h * (bimpl * taylor_yp(1) + sum(c * taylor_yp(-j) for j, c in b.items()))
    expr = sp.expand(expr)
    pol = sp.Poly(expr, h)
    coefs = {m[0]: sp.simplify(c) for m, c in zip(pol.monoms(), pol.coeffs())}
    # primer término no nulo
    k0 = min(coefs)
    return k0, coefs[k0], coefs

MULTI = {
    "IT25-T2": dict(a={0: Fr(4, 3), 1: Fr(-1, 3)}, b={}, bimpl=Fr(2, 3)),
    "IT25-T3": dict(a={0: Fr(18, 11), 1: Fr(-9, 11), 2: Fr(2, 11)}, b={}, bimpl=Fr(6, 11)),
    "IT25-T4 (h, como dice el enunciado)": dict(a={1: 1}, b={0: 1}),
    "IT25-T4 (lectura 2h)": dict(a={1: 1}, b={0: 2}),
    "Euler explícito": dict(a={0: 1}, b={0: 1}),
    "Euler implícito": dict(a={0: 1}, b={}, bimpl=1),
    "Trapecio (implícito)": dict(a={0: 1}, b={0: Fr(1, 2)}, bimpl=Fr(1, 2)),
    "AB2": dict(a={0: 1}, b={0: Fr(3, 2), 1: Fr(-1, 2)}),
    "AB3 (F18-P1)": dict(a={0: 1}, b={0: Fr(23, 12), 1: Fr(-16, 12), 2: Fr(5, 12)}),
    "AB4": dict(a={0: 1}, b={0: Fr(55, 24), 1: Fr(-59, 24), 2: Fr(37, 24), 3: Fr(-9, 24)}),
    "AM2 (cátedra: 'de 2 pasos')": dict(a={0: 1}, b={0: Fr(8, 12), 1: Fr(-1, 12)}, bimpl=Fr(5, 12)),
    "AM3 (cátedra: 'de 3 pasos')": dict(a={0: 1}, b={0: Fr(19, 24), 1: Fr(-5, 24), 2: Fr(1, 24)}, bimpl=Fr(9, 24)),
    "Milne predictor": dict(a={3: 1}, b={0: Fr(8, 3), 1: Fr(-4, 3), 2: Fr(8, 3)}),
    "Milne corrector": dict(a={1: 1}, b={0: Fr(4, 3), 1: Fr(1, 3)}, bimpl=Fr(1, 3)),
}

out("=" * 78)
out("1) ORDEN POR TAYLOR — métodos lineales multipaso (LTE = C·h^(p+1)·y^(p+1) + ...)")
R["multi"] = {}
for nom, d in MULTI.items():
    d2 = {k: ({j: sp.Rational(c.numerator, c.denominator) if isinstance(c, Fr) else c for j, c in v.items()}
              if isinstance(v, dict) else (sp.Rational(v.numerator, v.denominator) if isinstance(v, Fr) else v))
          for k, v in d.items()}
    k0, c0, coefs = lte_multipaso(**d2)
    consistente = k0 >= 2
    p = k0 - 1
    out(f"  {nom:40s} primer término de LTE: h^{k0} · ({c0})   → orden p = {p}   consistente: {consistente}")
    R["multi"][nom] = dict(k0=k0, c0=str(c0), p=p, consistente=consistente)

# ---------------------------------------------------------------------------
# 2) ORDEN DE MÉTODOS RK — (a) condiciones de Butcher exactas hasta orden 5
#                          (b) Taylor directo con una f(x,y) concreta no lineal
# ---------------------------------------------------------------------------
Q = sp.Rational
TABLAS = {
    "Euler explícito (RK1)": ([[0]], [1]),
    "Euler implícito (RK1)": ([[1]], [1]),
    "Trapecio implícito (y_{n+1}=y_n+h/2(f_n+f_{n+1}))": ([[0, 0], [Q(1, 2), Q(1, 2)]], [Q(1, 2), Q(1, 2)]),
    "Euler modificado / Heun 1 corrección (RK2, γ1=γ2=1/2, α=β=1)": ([[0, 0], [1, 0]], [Q(1, 2), Q(1, 2)]),
    "RK2 con β=2 (errata de la diapositiva)": ([[0, 0], [2, 0]], [Q(1, 2), Q(1, 2)]),
    "Punto medio (RK2, γ1=0, γ2=1, α=β=1/2)": ([[0, 0], [Q(1, 2), 0]], [0, 1]),
    "RK4 clásico": ([[0, 0, 0, 0], [Q(1, 2), 0, 0, 0], [0, Q(1, 2), 0, 0], [0, 0, 1, 0]],
                    [Q(1, 6), Q(1, 3), Q(1, 3), Q(1, 6)]),
    "IT25-T1 lectura A: RK3 de Kutta (k2 en x+h/2 con k1; k3 en x+h con -k1+2k2)":
        ([[0, 0, 0], [Q(1, 2), 0, 0], [-1, 2, 0]], [Q(1, 6), Q(2, 3), Q(1, 6)]),
    "IT25-T1 lectura B: k3 = f(x+h, y+h·k2)": ([[0, 0, 0], [Q(1, 2), 0, 0], [0, 1, 0]], [Q(1, 6), Q(2, 3), Q(1, 6)]),
    "IT25-T1 lectura C: k3 = f(x+h, y+h·k1)": ([[0, 0, 0], [Q(1, 2), 0, 0], [1, 0, 0]], [Q(1, 6), Q(2, 3), Q(1, 6)]),
    "IT25-T1 lectura D: implícito (Lobatto IIIA, etapas acopladas)":
        ([[0, 0, 0], [Q(5, 24), Q(1, 3), Q(-1, 24)], [Q(1, 6), Q(2, 3), Q(1, 6)]], [Q(1, 6), Q(2, 3), Q(1, 6)]),
}
A5 = [[0] * 6 for _ in range(6)]
A5[1][0] = Q(1, 4)
A5[2][:2] = [Q(3, 32), Q(9, 32)]
A5[3][:3] = [Q(1932, 2197), Q(-7200, 2197), Q(7296, 2197)]
A5[4][:4] = [Q(439, 216), -8, Q(3680, 513), Q(-845, 4104)]
A5[5][:5] = [Q(-8, 27), 2, Q(-3544, 2565), Q(1859, 4104), Q(-11, 40)]
TABLAS["RKF45: fórmula 'RK4' (pesos 25/216, 1408/2565, 2197/4104, -1/5)"] = (A5, [Q(25, 216), 0, Q(1408, 2565), Q(2197, 4104), Q(-1, 5), 0])
TABLAS["RKF45: fórmula 'RK5' (pesos 16/135, 6656/12825, 28561/56430, -9/50, 2/55)"] = (A5, [Q(16, 135), 0, Q(6656, 12825), Q(28561, 56430), Q(-9, 50), Q(2, 55)])
TABLAS["RKF45 'RK4' tal como lo escribe la diapositiva (23/216, 1408/2565, 219/4104, -1/5)"] = (A5, [Q(23, 216), 0, Q(1408, 2565), Q(219, 4104), Q(-1, 5), 0])

def orden_butcher(A, b):
    A = sp.Matrix(A); b = sp.Matrix(b); s = len(b)
    one = sp.ones(s, 1)
    c = A * one
    C = sp.diag(*c)
    dot = lambda u, v: (u.T * v)[0]
    conds = {
        1: [(dot(b, one), Q(1))],
        2: [(dot(b, c), Q(1, 2))],
        3: [(dot(b, C * c), Q(1, 3)), (dot(b, A * c), Q(1, 6))],
        4: [(dot(b, C * C * c), Q(1, 4)), (dot(b, C * A * c), Q(1, 8)), (dot(b, A * C * c), Q(1, 12)), (dot(b, A * A * c), Q(1, 24))],
        5: [(dot(b, C ** 3 * c), Q(1, 5)), (dot(b, C * C * A * c), Q(1, 10)), (dot(b, C * A * C * c), Q(1, 15)),
            (dot(b, C * A * A * c), Q(1, 30)), (dot(b, sp.Matrix([(A * c)[i] ** 2 for i in range(s)])), Q(1, 20)),
            (dot(b, A * C * C * c), Q(1, 20)), (dot(b, A * C * A * c), Q(1, 40)), (dot(b, A * A * C * c), Q(1, 60)),
            (dot(b, A * A * A * c), Q(1, 120))],
    }
    p = 0
    for k in range(1, 6):
        if all(sp.simplify(v - w) == 0 for v, w in conds[k]):
            p = k
        else:
            break
    return p

# Taylor directo con f concreta (no lineal en y y dependiente de x)
xs, ys = sp.symbols('x y')
F = xs * ys ** 2 + ys + xs ** 2
NT = 6

def derivs_exactas(Fexpr, n):
    """y', y'', ... en (x,y) por derivada total."""
    d = [Fexpr]
    for _ in range(n - 1):
        g = d[-1]
        d.append(sp.expand(sp.diff(g, xs) + sp.diff(g, ys) * Fexpr))
    return d

DER = derivs_exactas(F, NT)
x0v, y0v = Q(1, 3), Q(1, 2)

def serie_exacta():
    return sp.expand(y0v + sum(DER[k].subs({xs: x0v, ys: y0v}) * h ** (k + 1) / sp.factorial(k + 1) for k in range(NT)))

def trunc(e, n=NT + 1):
    e = sp.expand(e)
    return sum(e.coeff(h, k) * h ** k for k in range(n + 1))

def fser(xv, yv):
    """f evaluada con argumentos que son series en h, truncada."""
    return trunc(F.subs({xs: xv, ys: yv}))

def paso_rk(A, b):
    s = len(b)
    c = [sum(A[i]) for i in range(s)]
    k = [sp.Integer(0)] * s
    explicito = all(A[i][j] == 0 for i in range(s) for j in range(i, s))
    iters = 1 if explicito else NT + 2
    for _ in range(iters):
        knew = []
        for i in range(s):
            ks = k if not explicito else knew + k[len(knew):]
            yi = y0v + h * sum(A[i][j] * ks[j] for j in range(s))
            knew.append(fser(x0v + c[i] * h, trunc(yi)))
        k = knew
    return trunc(y0v + h * sum(b[i] * k[i] for i in range(s)))

ex = serie_exacta()
out("=" * 78)
out("2) ORDEN DE MÉTODOS RK — condiciones de Butcher (exactas) y Taylor directo con f(x,y)=x·y²+y+x²")
R["rk"] = {}
for nom, (A, b) in TABLAS.items():
    pB = orden_butcher(A, b)
    dif = sp.expand(paso_rk(A, b) - ex)
    kmin = min(k for k in range(NT + 2) if dif.coeff(h, k) != 0)
    pT = kmin - 1
    out(f"  {nom}\n      Butcher: orden {pB}   |  Taylor: error local h^{kmin} → orden {pT}")
    R["rk"][nom] = dict(butcher=pB, taylor=pT, kmin=kmin)

# ---------------------------------------------------------------------------
# 3) ORDEN EMPÍRICO (convergencia) en y' = y - t² + 1, y(0)=0.5, [0,1]
# ---------------------------------------------------------------------------
fEj = lambda t, y: y - t * t + 1
yEj = lambda t: (t + 1) ** 2 - 0.5 * math.exp(t)

def int_unpaso(paso, f, t0, y0, hh, n):
    t, y = t0, y0
    for _ in range(n):
        y = paso(f, t, y, hh); t += hh
    return y

def euler(f, t, y, hh): return y + hh * f(t, y)
def rk2(f, t, y, hh):
    k1 = f(t, y); k2 = f(t + hh, y + hh * k1); return y + hh / 2 * (k1 + k2)
def pmedio(f, t, y, hh):
    return y + hh * f(t + hh / 2, y + hh / 2 * f(t, y))
def rk4(f, t, y, hh):
    k1 = f(t, y); k2 = f(t + hh / 2, y + hh * k1 / 2); k3 = f(t + hh / 2, y + hh * k2 / 2); k4 = f(t + hh, y + hh * k3)
    return y + hh / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
def kutta3(f, t, y, hh):
    k1 = f(t, y); k2 = f(t + hh / 2, y + hh * k1 / 2); k3 = f(t + hh, y - hh * k1 + 2 * hh * k2)
    return y + hh / 6 * (k1 + 4 * k2 + k3)

def multipaso_lineal(f, coefy, coefbf, bimpl, t0, y0, hh, n, escala=1):
    """Para y' = y - t² + 1 (lineal en y): el término implícito se despeja exacto.
    Arranque con la solución exacta (para medir sólo el orden del método)."""
    k = max(list(coefy) + list(coefbf) + [0]) + 1
    ts = [t0 + i * hh for i in range(n + 1)]
    ys = [yEj(ts[i]) for i in range(k)]
    for m in range(k - 1, n):
        rhs = sum(c * ys[m - j] for j, c in coefy.items()) + escala * hh * sum(c * f(ts[m - j], ys[m - j]) for j, c in coefbf.items())
        tn1 = ts[m + 1]
        # y_{n+1} = rhs + escala*h*bimpl*(y_{n+1} - t² + 1)
        g = escala * hh * bimpl
        ys.append((rhs + g * (-tn1 * tn1 + 1)) / (1 - g))
    return ys[-1]

out("=" * 78)
out("3) ORDEN EMPÍRICO: error global en t=1 para y'=y-t²+1, y(0)=0.5 (log2 del cociente al dividir h por 2)")
R["emp"] = {}
def emp(nom, fun):
    errs = []
    for n in (10, 20, 40, 80):
        errs.append(abs(fun(1.0 / n, n) - yEj(1.0)))
    ords = [math.log2(errs[i] / errs[i + 1]) for i in range(3)]
    out(f"  {nom:38s} errores {['%.2e' % e for e in errs]}  órdenes {['%.2f' % o for o in ords]}")
    R["emp"][nom] = dict(errs=errs, ords=ords)

for nom, p in [("Euler explícito", euler), ("RK2 Euler modificado", rk2), ("Punto medio", pmedio), ("RK3 Kutta (IT25-T1 lect. A)", kutta3), ("RK4", rk4)]:
    emp(nom, lambda hh, n, p=p: int_unpaso(p, fEj, 0.0, 0.5, hh, n))
emp("Euler implícito", lambda hh, n: multipaso_lineal(fEj, {0: 1}, {}, 1, 0.0, 0.5, hh, n))
emp("Trapecio implícito", lambda hh, n: multipaso_lineal(fEj, {0: 1}, {0: 0.5}, 0.5, 0.0, 0.5, hh, n))
emp("AB2 (arranque exacto)", lambda hh, n: multipaso_lineal(fEj, {0: 1}, {0: 1.5, 1: -0.5}, 0, 0.0, 0.5, hh, n))
emp("AB3 (arranque exacto)", lambda hh, n: multipaso_lineal(fEj, {0: 1}, {0: 23 / 12, 1: -16 / 12, 2: 5 / 12}, 0, 0.0, 0.5, hh, n))
emp("IT25-T2 (BDF2)", lambda hh, n: multipaso_lineal(fEj, {0: 4 / 3, 1: -1 / 3}, {}, 2 / 3, 0.0, 0.5, hh, n))
emp("IT25-T3 (BDF3)", lambda hh, n: multipaso_lineal(fEj, {0: 18 / 11, 1: -9 / 11, 2: 2 / 11}, {}, 6 / 11, 0.0, 0.5, hh, n))
emp("IT25-T4 lectura 2h (punto medio 2 pasos)", lambda hh, n: multipaso_lineal(fEj, {1: 1}, {0: 1}, 0, 0.0, 0.5, hh, n, escala=2))
# IT25-T4 con h: no converge
out("  IT25-T4 con h (y_{n+1}=y_{n-1}+h·f_n), y'=y, y(0)=1, valor en t=1 (exacto e=2.71828; e^{1/2}=1.64872):")
R["t4h"] = []
for n in (10, 100, 1000, 10000):
    hh = 1.0 / n
    ys_ = [1.0, math.exp(hh)]
    for m in range(1, n):
        ys_.append(ys_[m - 1] + hh * ys_[m])
    out(f"     h=1/{n:<6d} y_n(1) = {ys_[-1]:.5f}")
    R["t4h"].append((n, ys_[-1]))
R["t4h_e"] = math.e
R["t4h_sqrt_e"] = math.exp(0.5)

# ---------------------------------------------------------------------------
# 4) EJEMPLOS DE LA CÁTEDRA (y' = y - t² + 1, y(0)=0.5, h=0.2) — precisión completa
#    E local (cátedra) = |y_{i+1} - y_i| ; E global = |y(t_i) - y_i|
# ---------------------------------------------------------------------------
out("=" * 78)
out("4) EJEMPLOS DE LA CÁTEDRA, y'=y-t²+1, y(0)=0.5, h=0.2 (precisión completa)")
H = 0.2
T = [round(i * H, 10) for i in range(6)]

def tabla(ys_, yreal=yEj, ts=T):
    filas = []
    for i, (t, y) in enumerate(zip(ts, ys_)):
        el = None if i == 0 else abs(ys_[i] - ys_[i - 1])
        eg = abs(yreal(t) - y)
        filas.append((t, yreal(t), y, el, eg))
    return filas

def imprime(nom, filas):
    out(f"  {nom}")
    out("     t     y real        y aprox       E local       E global")
    for t, yr, ya, el, eg in filas:
        out(f"    {t:4.1f}  {yr:11.8f}  {ya:11.8f}  {'' if el is None else '%11.8f' % el:>11s}  {eg:11.8f}")

ys_e = [0.5]
for i in range(5): ys_e.append(euler(fEj, T[i], ys_e[-1], H))
R["cat_euler"] = tabla(ys_e); imprime("Euler explícito", R["cat_euler"])
ys_ei = [0.5]
for i in range(5): ys_ei.append((ys_ei[-1] - H * T[i + 1] ** 2 + H) / (1 - H))
R["cat_eimpl"] = tabla(ys_ei); imprime("Euler implícito (despeje y_{n+1}=(y_n-h t²_{n+1}+h)/(1-h))", R["cat_eimpl"])
ys_ab = [0.5, euler(fEj, 0, 0.5, H)]
for i in range(1, 5): ys_ab.append(ys_ab[i] + H / 2 * (3 * fEj(T[i], ys_ab[i]) - fEj(T[i - 1], ys_ab[i - 1])))
R["cat_ab2"] = tabla(ys_ab); imprime("AB2 con arranque Euler (cátedra)", R["cat_ab2"])

def heun_pc(f, t, y, hh, tol=1e-2, itmax=5, rd=lambda v: v, traza=None):
    """Heun de la cátedra: P = Euler explícito, C = Trapecio, se itera el corrector
    MIENTRAS error >= tol && it < itMax, con error = |Yant - y1|."""
    fy = rd(f(t, y))
    y1 = rd(y + rd(hh * fy))
    if traza is not None: traza.append(("P", y1))
    err = abs(y1 - y); it = 0
    while err >= tol and it < itmax:
        yant = y1
        y1 = rd(y + rd(rd(hh / 2) * rd(fy + rd(f(t + hh, y1)))))
        err = abs(yant - y1); it += 1
        if traza is not None: traza.append(("C", y1, err))
    return y1, it

ys_h = [0.5]; tr0 = []
for i in range(5):
    yv, _ = heun_pc(fEj, T[i], ys_h[-1], H, traza=tr0 if i == 0 else None)
    ys_h.append(yv)
R["cat_heun"] = tabla(ys_h); R["cat_heun_traza"] = tr0
imprime("Heun P-C (Euler + Trapecio iterado, tol 1e-2, itMax 5)", R["cat_heun"])
out("     traza del primer paso:", tr0)
ys_r2 = [0.5]
for i in range(5): ys_r2.append(rk2(fEj, T[i], ys_r2[-1], H))
R["cat_rk2"] = tabla(ys_r2); imprime("RK2 Euler modificado (sin iterar)", R["cat_rk2"])
ys_r4 = [0.5]
for i in range(5): ys_r4.append(rk4(fEj, T[i], ys_r4[-1], H))
R["cat_rk4"] = tabla(ys_r4); imprime("RK4", R["cat_rk4"])

# ---------------------------------------------------------------------------
# 5) ARITMÉTICA DE 3 DÍGITOS (TP4): redondeo a 3 cifras significativas después de CADA operación
# ---------------------------------------------------------------------------
def r3(v):
    if v == 0 or not math.isfinite(v): return v
    return float(f"{v:.3g}")

class D:
    """Número con aritmética de 3 cifras significativas: cada +, -, ×, ÷ redondea el resultado."""
    __slots__ = ("v",)
    def __init__(self, v): self.v = r3(float(v.v if isinstance(v, D) else v))
    def _o(self, o): return o.v if isinstance(o, D) else float(o)
    def __add__(self, o): return D(self.v + self._o(o))
    __radd__ = __add__
    def __sub__(self, o): return D(self.v - self._o(o))
    def __rsub__(self, o): return D(self._o(o) - self.v)
    def __mul__(self, o): return D(self.v * self._o(o))
    __rmul__ = __mul__
    def __truediv__(self, o): return D(self.v / self._o(o))
    def __rtruediv__(self, o): return D(self._o(o) / self.v)
    def __neg__(self): return D(-self.v)
    def __abs__(self): return D(abs(self.v))
    def __float__(self): return self.v
    def __repr__(self): return f"{self.v:g}"

def dexp(x): return D(math.exp(float(x)))
def dsqrt(x): return D(math.sqrt(float(x)))
def dlog(x): return D(math.log(float(x)))

def rk2_d(f, t, y, hh, pasos=None):
    k1 = f(t, y)
    k2 = f(t + hh, y + hh * k1)
    yn = y + (hh / 2) * (k1 + k2)
    if pasos is not None: pasos.append(dict(k1=k1, k2=k2, y=yn))
    return yn

def rk4_d(f, t, y, hh, pasos=None):
    m = hh / 2
    k1 = f(t, y)
    k2 = f(t + m, y + m * k1)
    k3 = f(t + m, y + m * k2)
    k4 = f(t + hh, y + hh * k3)
    s = k1 + 2 * k2 + 2 * k3 + k4
    yn = y + (hh * s) / 6
    if pasos is not None: pasos.append(dict(k1=k1, k2=k2, k3=k3, k4=k4, s=s, y=yn))
    return yn

def euler_d(f, t, y, hh): return y + hh * f(t, y)

def ab2_d(f, ts, ys_, hh, i):
    return ys_[i] + (hh / 2) * (3 * f(ts[i], ys_[i]) - f(ts[i - 1], ys_[i - 1]))

def heun_d(f, t, y, hh, traza=None, tol=1e-2, itmax=5):
    fy = f(t, y)
    y1 = y + hh * fy
    if traza is not None: traza.append(("P", y1))
    err = abs(float(y1) - float(y)); it = 0
    while err >= tol and it < itmax:
        yant = y1
        y1 = y + (hh / 2) * (fy + f(t + hh, y1))
        err = abs(float(yant) - float(y1)); it += 1
        if traza is not None: traza.append(("C", y1, err))
    return y1, it

def corre_tp(nom, f, yreal, t0, y0, hh, n, metodos=("RK2", "RK4", "AB2", "Heun"), fexacta=None):
    ts = [D(t0 + i * hh) for i in range(n + 1)]
    res = {}
    for m in metodos:
        ys_ = [D(y0)]; det = []; its = []
        for i in range(n):
            if m == "RK2": ys_.append(rk2_d(f, ts[i], ys_[i], D(hh), det))
            elif m == "RK4": ys_.append(rk4_d(f, ts[i], ys_[i], D(hh), det))
            elif m == "AB2":
                if i == 0:
                    fi = f(ts[0], ys_[0]); det.append(dict(f=fi)); ys_.append(ys_[0] + D(hh) * fi)
                else:
                    fi = f(ts[i], ys_[i]); fm = f(ts[i - 1], ys_[i - 1]); tres = 3 * fi; dif = tres - fm
                    yn = ys_[i] + (D(hh) / 2) * dif
                    det.append(dict(f=fi, fm=fm, tres=tres, dif=dif)); ys_.append(yn)
            elif m == "Heun":
                tr = []
                yv, it = heun_d(f, ts[i], ys_[i], D(hh), traza=tr)
                ys_.append(yv); its.append(it); det.append(tr)
        filas = []
        for i in range(n + 1):
            t = float(ts[i]); yr = yreal(t); ya = float(ys_[i])
            el = None if i == 0 else abs(ya - float(ys_[i - 1]))
            filas.append((t, yr, ya, el, abs(yr - ya)))
        res[m] = dict(filas=filas, det=det, its=its)
        out(f"  [{nom}] {m} (3 dígitos)")
        out("     t     y real       y aprox     E local     E global")
        for t, yr, ya, el, eg in filas:
            out(f"    {t:4.1f}  {yr:10.6f}  {ya:10.4g}  {'' if el is None else '%10.3g' % el:>10s}  {eg:10.3g}")
        if its: out("     iteraciones del corrector por paso:", its)
        out("     detalle 1er paso:", det[0], ("| 2do paso: " + str(det[1])) if len(det) > 1 else "")
    return res

out("=" * 78)
out("5) TP4 ej. 1 — aritmética de 3 cifras significativas (redondeo en cada operación)")
# a) y' = 1 + y/t, 1<=t<=2, y(1)=2, exacta t ln t + 2t
fa = lambda t, y: 1 + y / t
ya_ex = lambda t: t * math.log(t) + 2 * t
R["tp_a"] = corre_tp("1a", fa, ya_ex, 1.0, 2.0, 0.2, 5)
# misma cuenta con precisión completa, para comparar el efecto del redondeo
fa_f = lambda t, y: 1 + y / t
full = {}
for nom, p in [("RK2", rk2), ("RK4", rk4)]:
    yy = [2.0]
    for i in range(5): yy.append(p(fa_f, 1 + 0.2 * i, yy[-1], 0.2))
    full[nom] = yy
yy = [2.0, euler(fa_f, 1.0, 2.0, 0.2)]
for i in range(1, 5): yy.append(yy[i] + 0.1 * (3 * fa_f(1 + 0.2 * i, yy[i]) - fa_f(1 + 0.2 * (i - 1), yy[i - 1])))
full["AB2"] = yy
yy = [2.0]
for i in range(5): yy.append(heun_pc(fa_f, 1 + 0.2 * i, yy[-1], 0.2)[0])
full["Heun"] = yy
R["tp_a_full"] = full
out("  1a precisión completa y(2):", {k: round(v[-1], 6) for k, v in full.items()}, " exacto:", round(ya_ex(2.0), 6))

# c) y' = 3t/y, y(1)=2, exacta sqrt(3t²+1)
fc = lambda t, y: (3 * t) / y
yc_ex = lambda t: math.sqrt(3 * t * t + 1)
R["tp_c"] = corre_tp("1c", fc, yc_ex, 1.0, 2.0, 0.2, 5)
# d) bacterias: P' = 3000 e^{2t/5}, P(0)=7500 → exacta 7500 e^{0.4 t}
fd = lambda t, y: 3000 * dexp((2 * t) / 5)
yd_ex = lambda t: 7500 * math.exp(0.4 * t)
R["tp_d"] = corre_tp("1d", fd, yd_ex, 0.0, 7500.0, 0.2, 5)
R["tp_d_exacta_1"] = yd_ex(1.0)
# chequeo simbólico de la exacta de d)
tt = sp.symbols('t'); Pf = sp.Function('P')
sol = sp.dsolve(sp.Eq(Pf(tt).diff(tt), 3000 * sp.exp(sp.Rational(2, 5) * tt)), Pf(tt), ics={Pf(0): 7500})
R["tp_d_sol"] = str(sol.rhs)
out("  1d solución exacta (sympy):", sol.rhs, " P(1) =", round(yd_ex(1.0), 4))
# chequeo de las exactas de a) y c) con sympy
ya_s = tt * sp.log(tt) + 2 * tt
out("  chequeo exacta 1a:", sp.simplify(ya_s.diff(tt) - (1 + ya_s / tt)), ", y(1) =", ya_s.subs(tt, 1))
yc_s = sp.sqrt(3 * tt ** 2 + 1)
out("  chequeo exacta 1c:", sp.simplify(yc_s.diff(tt) - 3 * tt / yc_s), ", y(1) =", yc_s.subs(tt, 1))
yb_s = tt - sp.exp(-5 * tt)
out("  chequeo exacta 1b:", sp.simplify(yb_s.diff(tt) - (5 * sp.exp(5 * tt) * (yb_s - tt) ** 2 + 1)), ", y(-1) =", yb_s.subs(tt, -1), "=", float(yb_s.subs(tt, -1)))

# ---------------------------------------------------------------------------
# 6) TP4 1b (stiff): y' = 5e^{5t}(y-t)²+1, -1<=t<=0, y(-1) = -1-e^5, exacta t - e^{-5t}
# ---------------------------------------------------------------------------
out("=" * 78)
out("6) TP4 1b (stiff) — precisión completa (con 3 dígitos el resultado es todavía peor)")
fb = lambda t, y: 5 * math.exp(5 * t) * (y - t) ** 2 + 1
yb_ex = lambda t: t - math.exp(-5 * t)
def integra(p, f, t0, y0, hh, n, multi=False):
    ts = [t0 + i * hh for i in range(n + 1)]; ys_ = [y0]
    for i in range(n):
        if multi and i > 0:
            ys_.append(ys_[i] + hh / 2 * (3 * f(ts[i], ys_[i]) - f(ts[i - 1], ys_[i - 1])))
        else:
            try:
                ys_.append(p(f, ts[i], ys_[i], hh))
            except OverflowError:
                ys_.append(float('inf'))
        if not math.isfinite(ys_[-1]) or abs(ys_[-1]) > 1e300:
            ys_[-1] = float('inf')
    return ts, ys_
def heun_f(f, t, y, hh): return heun_pc(f, t, y, hh)[0]
R["tp_b"] = {}
for hh, n in ((0.2, 5), (0.1, 10), (0.01, 100), (0.001, 1000)):
    fila = {}
    for nom, p, mlt in (("Euler", euler, False), ("RK2", rk2, False), ("RK4", rk4, False), ("AB2", euler, True), ("Heun", heun_f, False)):
        with_err = True
        try:
            ts, ys_ = integra(p, fb, -1.0, -1 - math.exp(5), hh, n, multi=mlt)
        except OverflowError:
            ts, ys_ = [0], [float('inf')]
        fila[nom] = (ys_[-1], abs(ys_[-1] - yb_ex(0.0)) if math.isfinite(ys_[-1]) else float('inf'), ys_)
    R["tp_b"][hh] = fila
    out(f"  h={hh}: y(0) aprox = " + ", ".join(f"{k}: {v[0]:.6g}" for k, v in fila.items()) + f"   (exacto {yb_ex(0.0)})")
# primeros pasos con h=0.2 en RK4 y Euler
for nom in ("Euler", "RK2", "RK4", "AB2", "Heun"):
    out(f"    h=0.2 {nom}: ", ["%.4g" % v for v in R["tp_b"][0.2][nom][2]])
out("    exactos en la grilla h=0.2:", ["%.4g" % yb_ex(-1 + 0.2 * i) for i in range(6)])
R["tp_b_ex"] = [yb_ex(-1 + 0.2 * i) for i in range(6)]
# derivada parcial ∂f/∂y sobre la solución exacta: 10 e^{5t}(y-t) = -10 (constante)
yy_, ttt = sp.symbols('y t')
fbs = 5 * sp.exp(5 * ttt) * (yy_ - ttt) ** 2 + 1
dfdy = sp.simplify(sp.diff(fbs, yy_).subs(yy_, ttt - sp.exp(-5 * ttt)))
R["tp_b_dfdy"] = str(dfdy)
out("  ∂f/∂y evaluada sobre la solución exacta:", dfdy, " → h·λ con h=0.2:", float(dfdy) * 0.2)
# 3 dígitos, primer paso de Euler con h=0.2
fb_d = lambda t, y: 5 * dexp(5 * t) * ((y - t) * (y - t)) + 1
y0d = D(-1 - math.exp(5))
R["tp_b_d_euler1"] = float(euler_d(fb_d, D(-1), y0d, D(0.2)))
out("  1b con 3 dígitos: y(-1) =", y0d, "  Euler 1 paso → y(-0.8) =", R["tp_b_d_euler1"], " exacto", yb_ex(-0.8))

# ---------------------------------------------------------------------------
# 7) EJEMPLO STIFF DE LA CÁTEDRA: y' = -21y + e^{-t}, y(0)=0
# ---------------------------------------------------------------------------
out("=" * 78)
out("7) Stiff de la cátedra y'=-21y+e^{-t}, y(0)=0, h=0.1 en [0,1]")
fs = lambda t, y: -21 * y + math.exp(-t)
ys_ex = lambda t: (math.exp(-t) - math.exp(-21 * t)) / 20
ts_s = sp.symbols('t'); ysf = sp.Function('y')
sol_s = sp.dsolve(sp.Eq(ysf(ts_s).diff(ts_s), -21 * ysf(ts_s) + sp.exp(-ts_s)), ysf(ts_s), ics={ysf(0): 0})
out("  exacta (sympy):", sol_s.rhs)
R["stiff_cat"] = {}
for nom, p in (("Euler", euler), ("RK2", rk2), ("RK4", rk4), ("Heun P-C", heun_f)):
    ts, yy = integra(p, fs, 0.0, 0.0, 0.1, 10)
    R["stiff_cat"][nom] = yy
    out(f"  {nom:9s}", ["%.4g" % v for v in yy])
R["stiff_cat_t"] = [0.1 * i for i in range(11)]
R["stiff_cat_ex"] = [ys_ex(0.1 * i) for i in range(11)]
out("  exacta   ", ["%.4g" % v for v in R["stiff_cat_ex"]])
# factores de amplificación para y' = λy con hλ = -2.1
z = -2.1
R["amp"] = dict(z=z, euler=1 + z, rk2=1 + z + z * z / 2, rk4=1 + z + z ** 2 / 2 + z ** 3 / 6 + z ** 4 / 24,
                trap=(1 + z / 2) / (1 - z / 2), eimpl=1 / (1 - z), heun_iter=abs(-21 * 0.1 / 2))
out("  factor de amplificación por paso con hλ=-2.1:", {k: round(v, 4) for k, v in R["amp"].items()})
out("  (heun_iter = |h·λ/2|: si >1 la iteración del corrector diverge)")
ts_, yy_h = integra(euler, fs, 0.0, 0.0, 0.001, 1000)
R["stiff_cat_euler_h0001"] = abs(yy_h[-1] - ys_ex(1.0))
out("  Euler h=0.001: error en t=1 =", "%.3g" % R["stiff_cat_euler_h0001"])

# ---------------------------------------------------------------------------
# 8) SISTEMAS / ORDEN SUPERIOR — TP4 ej. 2d: y'' - 2y' + 2y = e^{2t} sen t, y(0)=-0.4, y'(0)=-0.6
#    y1 = y, y2 = y'  →  y1' = y2 ; y2' = e^{2t} sen t - 2 y1 + 2 y2
# ---------------------------------------------------------------------------
out("=" * 78)
out("8) TP4 ej. 2d como sistema, h=0.2 (precisión completa)")
def F2(t, u): return [u[1], math.exp(2 * t) * math.sin(t) - 2 * u[0] + 2 * u[1]]
y2d_ex = lambda t: 0.2 * math.exp(2 * t) * (math.sin(t) - 2 * math.cos(t))
y2d_exp = lambda t: 0.2 * math.exp(2 * t) * (4 * math.sin(t) - 3 * math.cos(t))
tq = sp.symbols('t'); ysol = sp.Rational(1, 5) * sp.exp(2 * tq) * (sp.sin(tq) - 2 * sp.cos(tq))
out("  chequeo exacta 2d:", sp.simplify(ysol.diff(tq, 2) - 2 * ysol.diff(tq) + 2 * ysol - sp.exp(2 * tq) * sp.sin(tq)),
    " y(0) =", ysol.subs(tq, 0), " y'(0) =", ysol.diff(tq).subs(tq, 0), " y'(t) =", sp.simplify(ysol.diff(tq)))
u0 = [-0.4, -0.6]; hh = 0.2
fu0 = F2(0.0, u0)
ue1 = [u0[i] + hh * fu0[i] for i in range(2)]
k1 = F2(0.0, u0); k2 = F2(hh, [u0[i] + hh * k1[i] for i in range(2)])
ur1 = [u0[i] + hh / 2 * (k1[i] + k2[i]) for i in range(2)]
R["sis"] = dict(f0=fu0, euler1=ue1, k1=k1, k2=k2, rk2_1=ur1, ex1=y2d_ex(0.2), exp1=y2d_exp(0.2),
                el_euler=max(abs(ue1[i] - u0[i]) for i in range(2)), el_rk2=max(abs(ur1[i] - u0[i]) for i in range(2)))
out("  f(0,u0) =", fu0, "\n  Euler 1 paso:", ue1, "\n  RK2: k1 =", k1, " k2 =", k2, " → u1 =", ur1,
    "\n  exacta en 0.2: y =", y2d_ex(0.2), " y' =", y2d_exp(0.2),
    "\n  E local (norma infinito |u1-u0|): Euler", R["sis"]["el_euler"], " RK2", R["sis"]["el_rk2"])
# sistema de la cátedra (edo4): u1' = 3u1 + 2u2 - (2t²+1)e^{2t}; u2' = 4u1 + u2 + (t²+2t-4)e^{2t}
def Fc(t, u): return [3 * u[0] + 2 * u[1] - (2 * t * t + 1) * math.exp(2 * t), 4 * u[0] + u[1] + (t * t + 2 * t - 4) * math.exp(2 * t)]
fc0 = Fc(0.0, [1.0, 1.0]); uc1 = [1 + 0.2 * fc0[0], 1 + 0.2 * fc0[1]]
R["sis_cat"] = dict(f0=fc0, euler1=uc1, el=max(abs(uc1[0] - 1), abs(uc1[1] - 1)))
out("  sistema de la cátedra, Euler 1 paso (con +h·f): f(0,u0) =", fc0, " u(0.2) =", uc1, " EL(norma inf) =", R["sis_cat"]["el"])

# ---------------------------------------------------------------------------
# 9) Datos para figuras del quiz: error global vs t (1a, Euler vs RK4) y orden con h
# ---------------------------------------------------------------------------
out("=" * 78)
out("9) Datos para figuras")
ts9 = [1 + 0.1 * i for i in range(11)]
def corre(p, f, t0, y0, hh, n):
    yy = [y0]
    for i in range(n): yy.append(p(f, t0 + i * hh, yy[-1], hh))
    return yy
eE = [abs(a - ya_ex(t)) for a, t in zip(corre(euler, fa_f, 1.0, 2.0, 0.1, 10), ts9)]
eR = [abs(a - ya_ex(t)) for a, t in zip(corre(rk4, fa_f, 1.0, 2.0, 0.1, 10), ts9)]
R["fig_eg"] = dict(t=ts9, euler=eE, rk4=eR)
out("  1a h=0.1: EG final Euler = %.3g, RK4 = %.3g" % (eE[-1], eR[-1]))
# error final vs h (log-log) para Euler, RK2, RK4 en 1a
hs = [0.2, 0.1, 0.05, 0.025, 0.0125]
R["fig_h"] = dict(h=hs)
for nom, p in (("Euler", euler), ("RK2", rk2), ("RK4", rk4)):
    R["fig_h"][nom] = [abs(corre(p, fa_f, 1.0, 2.0, hh_, int(round(1 / hh_)))[-1] - ya_ex(2.0)) for hh_ in hs]
    out(f"  1a error en t=2 vs h, {nom}:", ["%.3g" % e for e in R["fig_h"][nom]])


# ---------------------------------------------------------------------------
# 10) Números para preguntas del quiz
# ---------------------------------------------------------------------------
out("=" * 78)
out("10) Quiz")
# Euler implícito, primer paso del ejemplo de la cátedra y distractores con errores típicos
num = 0.5 - 0.2 * 0.2 ** 2 + 0.2
R["q_eimpl"] = dict(ok=num / (1 - 0.2), sin_div=num, div_mas=num / (1 + 0.2), expl=0.5 + 0.2 * fEj(0, 0.5))
out("  Euler implícito y1:", R["q_eimpl"])
# IT25-T4 con h vs exacta para y'=y, y(0)=1 en [0,1]
def t4h_traj(n):
    hh = 1.0 / n; yy = [1.0, math.exp(hh)]
    for m in range(1, n): yy.append(yy[m - 1] + hh * yy[m])
    return [i * hh for i in range(n + 1)], yy
R["fig_t4"] = {n: t4h_traj(n) for n in (10, 100)}
out("  IT25-T4 (h) en t=1: n=10 →", R["fig_t4"][10][1][-1], " n=100 →", R["fig_t4"][100][1][-1])
# pendientes log-log de fig_h
for nom in ("Euler", "RK2", "RK4"):
    e = R["fig_h"][nom]
    out(f"  pendiente log-log {nom}:", round(math.log(e[0] / e[-1]) / math.log(0.2 / 0.0125), 2))
R["pend"] = {nom: math.log(R["fig_h"][nom][0] / R["fig_h"][nom][-1]) / math.log(0.2 / 0.0125) for nom in ("Euler", "RK2", "RK4")}

if __name__ == "__main__":
    print("\n".join(SALIDA))
