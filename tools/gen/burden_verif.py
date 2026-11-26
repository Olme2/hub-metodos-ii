"""Verifica con sympy que las soluciones exactas de los ejercicios del Burden (10.ª ed., cap. 5)
extraídas del PDF realmente resuelven cada PVI (la extracción del PDF puede tener errores)."""
import sympy as sp
t = sp.symbols('t', positive=True)
Y = sp.Function('y')
P = {
 "5.2-1a": (lambda y: t*sp.exp(3*t) - 2*y, sp.Rational(1,5)*t*sp.exp(3*t) - sp.Rational(1,25)*sp.exp(3*t) + sp.Rational(1,25)*sp.exp(-2*t), 0, 0),
 "5.2-1b": (lambda y: 1 + (t - y)**2, t + 1/(1 - t), 2, 1),
 "5.2-1c": (lambda y: 1 + y/t, t*sp.log(t) + 2*t, 1, 2),
 "5.2-1d": (lambda y: sp.cos(2*t) + sp.sin(3*t), sp.sin(2*t)/2 - sp.cos(3*t)/3 + sp.Rational(4,3), 0, 1),
 "5.4-3a": (lambda y: 2 - 2*t*y/(t**2 + 1), (2*t + 1)/(t**2 + 1), 0, 1),
 "5.4-3b": (lambda y: y**2/(1 + t), -1/sp.log(t + 1), 1, -1/sp.log(2)),
 "5.4-3c": (lambda y: (y**2 + y)/t, 2*t/(1 - 2*t), 1, -2),
 "5.4-3d": (lambda y: -t*y + 4*t/y, sp.sqrt(4 - 3*sp.exp(-t**2)), 0, 1),
 "5.11-1a": (lambda y: -9*y, sp.exp(1 - 9*t), 0, sp.E),
 "5.11-1b": (lambda y: -20*(y - t**2) + 2*t, t**2 + sp.exp(-20*t)/3, 0, sp.Rational(1,3)),
 "5.11-1c": (lambda y: -20*y + 20*sp.sin(t) + sp.cos(t), sp.sin(t) + sp.exp(-20*t), 0, 1),
 "5.11-1d": (lambda y: 50/y - 50*y, sp.sqrt(1 + sp.exp(-100*t)), 0, sp.sqrt(2)),
}
ok = True
for k, (f, sol, t0, y0) in P.items():
    res = sp.simplify(sp.diff(sol, t) - f(sol))
    ci = sp.simplify(sol.subs(t, t0) - y0) if t0 != 0 else sp.simplify(sp.limit(sol, t, 0) - y0)
    print(f"{k:8s} residuo EDO = {res}   condición inicial = {ci}")
    ok &= (res == 0 and ci == 0)
print("NOTA 5.4-3d: el PDF dice y = 4 − 3e^(−t²); con y(0)=1 sólo cumple sqrt(4 − 3e^(−t²)) — se usa la raíz (extracción del PDF sin la raíz).")
# sistemas y orden superior
u1, u2 = sp.Function('u1'), sp.Function('u2')
s1 = sp.dsolve([sp.Eq(u1(t).diff(t), 3*u1(t) + 2*u2(t) - (2*t**2 + 1)*sp.exp(2*t)),
                sp.Eq(u2(t).diff(t), 4*u1(t) + u2(t) + (t**2 + 2*t - 4)*sp.exp(2*t))], ics={u1(0): 1, u2(0): 1})
print("5.9-1a:", [sp.simplify(e.rhs) for e in s1])
ya = t**3*sp.exp(t)/6 - t*sp.exp(t) + 2*sp.exp(t) - t - 2
print("5.9-3a residuo:", sp.simplify(ya.diff(t, 2) - 2*ya.diff(t) + ya - (t*sp.exp(t) - t)), " y(0) =", ya.subs(t, 0), " y'(0) =", ya.diff(t).subs(t, 0))
print("TODO OK" if ok else "HAY SOLUCIONES QUE NO CUMPLEN")
