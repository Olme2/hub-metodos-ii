#!/usr/bin/env python3
"""Verificación numérica U2 — Ajuste de curvas. Todo número que usa m2.js sale de acá.
Convenciones:
 - E_RMS de la cátedra (U2 p18): sqrt(E/(n+1)). Se calcula también la usual sqrt(E/n).
 - r de la cátedra (U2 p19): r = sqrt((St-Sr)/St)  (0 <= r <= 1). Se calcula también Pearson con signo.
 - EV-P7: varianza/covarianza dividiendo por n (la que reproduce las opciones del Moodle);
   se calcula también con n-1 para mostrar que no coincide.
Uso: python3 verif_u2.py   (imprime todo). gen_u2.py hace `from verif_u2 import R, datos`.
"""
import numpy as np

R = {}

def lin(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); N = len(x)
    Sx, Sy, Sxy, Sxx = x.sum(), y.sum(), (x * y).sum(), (x * x).sum()
    den = N * Sxx - Sx ** 2
    a0 = (Sxx * Sy - Sxy * Sx) / den
    a1 = (N * Sxy - Sx * Sy) / den
    return dict(N=N, Sx=Sx, Sy=Sy, Sxy=Sxy, Sxx=Sxx, a0=a0, a1=a1)

def bondad(y, yap):
    y = np.asarray(y, float); yap = np.asarray(yap, float); n = len(y)
    St = ((y - y.mean()) ** 2).sum(); Sr = ((y - yap) ** 2).sum()
    return dict(St=St, Sr=Sr, R2=1 - Sr / St, r=np.sqrt(max(0, (St - Sr) / St)),
                E=Sr, Erms_cat=np.sqrt(Sr / (n + 1)), Erms_n=np.sqrt(Sr / n))

# ---------- 1. Ejemplo de la diapositiva (p9-10)
x = np.arange(1, 7); y = np.array([0, 0.6, 1.77, 1.92, 3.31, 3.52])
e = lin(x, y); e.update(bondad(y, e["a0"] + e["a1"] * x))
e["res"] = y - (e["a0"] + e["a1"] * x)
e["pred7"] = e["a0"] + e["a1"] * 7
e["pearson"] = np.corrcoef(x, y)[0, 1]
e["MAE"] = np.abs(e["res"]).mean()
e["R2adj"] = 1 - (1 - e["R2"]) * (len(x) - 1) / (len(x) - 2)
R["ej"] = e

# ---------- 2. COVID (p20-21): x=1..7
xc = np.arange(1, 8); yc = np.array([159, 156, 187, 150, 219, 198, 218])
c = lin(xc, yc); c.update(bondad(yc, c["a0"] + c["a1"] * xc)); R["covid"] = c

# ---------- 3. ET25-T1 ej 2
xe = np.array([3, 6, 7, 8, 9.]); ye = np.array([47, 493, 1215, 3132, 8292.])
t = lin(xe, ye); t.update(bondad(ye, t["a0"] + t["a1"] * xe))
t["coef_examen"] = 1146.59578596; t["int_examen"] = -4931.29430665
t["d_coef"] = t["a1"] - t["coef_examen"]; t["d_int"] = t["a0"] - t["int_examen"]
t["res"] = ye - (t["a0"] + t["a1"] * xe)
t["yap"] = t["a0"] + t["a1"] * xe
# R2 con los coeficientes del examen
t["R2_examen"] = bondad(ye, t["int_examen"] + t["coef_examen"] * xe)["R2"]
# Exponencial linealizado ln y = ln a + b x
le = lin(xe, np.log(ye)); a_exp = np.exp(le["a0"]); b_exp = le["a1"]
yexp = a_exp * np.exp(b_exp * xe)
t["exp"] = dict(lna=le["a0"], a=a_exp, b=b_exp, R2_log=bondad(np.log(ye), le["a0"] + le["a1"] * xe)["R2"],
                R2_orig=bondad(ye, yexp)["R2"], yap=yexp, factor=np.exp(b_exp),
                Erms_cat=bondad(ye, yexp)["Erms_cat"])
t["Erms_cat_lin"] = t["Erms_cat"]
# cuadrático, para comparar
cq = np.polyfit(xe, ye, 2); t["cuad"] = dict(c=cq, R2=bondad(ye, np.polyval(cq, xe))["R2"], yap=np.polyval(cq, xe))
# predicción en x=10 con cada modelo
t["p10_lin"] = t["a0"] + t["a1"] * 10; t["p10_exp"] = a_exp * np.exp(b_exp * 10); t["p10_cuad"] = np.polyval(cq, 10)
R["et1"] = t

# ---------- 4. EV-P1 acero
u, v = 0.05, 1125
R["evp1"] = dict(dur=161.3365 + 32.9688 * u - 0.0855 * v, t_u=32.9688 * u, t_v=0.0855 * v)
cu = [0.02, 0.10, 0.18]; T = [1000, 1100, 1200, 1300]
tab = {0.02: [78.9, 65.1, 55.2, 56.4], 0.10: [80.9, 69.7, 57.4, 55.4], 0.18: [85.3, 71.8, 60.7, 58.9]}
U, V, Y = [], [], []
for c_ in cu:
    for j, Tj in enumerate(T):
        U.append(c_); V.append(Tj); Y.append(tab[c_][j])
A = np.column_stack([np.ones(12), U, V]); M = A.T @ A; bvec = A.T @ np.array(Y)
R["evp1"]["coef_tabla"] = np.linalg.solve(M, bvec)
Y2 = list(Y); Y2[U.index(0.10) + 2] = 57.9
R["evp1"]["coef_tabla_579"] = np.linalg.lstsq(A, np.array(Y2), rcond=None)[0]
R["evp1"]["M"] = M; R["evp1"]["b"] = bvec
R["evp1"]["y_tabla"] = R["evp1"]["coef_tabla"] @ [1, u, v]

# ---------- 5. EV-P4
R["evp4"] = dict(y100=20 + 100 / 4, pend=1 / 4)

# ---------- 6. EV-P7
n = 8; Sx = 24; Sxy = 64; Sy = 40; sy2 = 12; sx2 = 6
xm, ym = Sx / n, Sy / n
sxy = Sxy / n - xm * ym
b = sxy / sx2; a = ym - b * xm
R["evp7"] = dict(xm=xm, ym=ym, sxy=sxy, b=b, a=a, p4=a + b * 4, r=sxy / np.sqrt(sx2 * sy2),
                 r2=(sxy / np.sqrt(sx2 * sy2)) ** 2, p4_swap=b + a * 4)
# con n-1 en la covarianza (y sx2 tomado como cuasivarianza)
sxy1 = (Sxy - n * xm * ym) / (n - 1); b1 = sxy1 / sx2; a1 = ym - b1 * xm
R["evp7"]["n1"] = dict(sxy=sxy1, b=b1, a=a1, p4=a1 + b1 * 4)
# Variante rotada (práctica)
n_, Sx_, Sy_, Sxy_, sx2_, sy2_ = 10, 50, 80, 440, 4, 9
xm_, ym_ = Sx_ / n_, Sy_ / n_; sxy_ = Sxy_ / n_ - xm_ * ym_; b_ = sxy_ / sx2_; a_ = ym_ - b_ * xm_
R["evp7v"] = dict(n=n_, Sx=Sx_, Sy=Sy_, Sxy=Sxy_, sx2=sx2_, sy2=sy2_, xm=xm_, ym=ym_, sxy=sxy_, b=b_, a=a_,
                  p7=a_ + b_ * 7, r=sxy_ / np.sqrt(sx2_ * sy2_), p7_swap=b_ + a_ * 7)
# Variante para el quiz
n_, Sx_, Sy_, Sxy_, sx2_, sy2_ = 5, 15, 40, 108, 2, 10
xm_, ym_ = Sx_ / n_, Sy_ / n_; sxy_ = Sxy_ / n_ - xm_ * ym_; b_ = sxy_ / sx2_; a_ = ym_ - b_ * xm_
R["evp7q"] = dict(n=n_, Sx=Sx_, Sy=Sy_, Sxy=Sxy_, sx2=sx2_, sy2=sy2_, xm=xm_, ym=ym_, sxy=sxy_, b=b_, a=a_,
                  p5=a_ + b_ * 5, r=sxy_ / np.sqrt(sx2_ * sy2_), p5_swap=b_ + a_ * 5)

# ---------- 7. TP2 ej 2: v^2 = 2 g x
xg = np.array([0, 5, 10, 15, 20, 25.]); vg = np.array([0, 9.85, 14.32, 17.63, 19.34, 22.41])
Yg = vg ** 2
k0 = (xg * Yg).sum() / (xg * xg).sum()               # recta por el origen: Y = k x
lg = lin(xg, Yg)                                      # recta con ordenada
R["tp2e2"] = dict(v2=Yg, Sxx=(xg * xg).sum(), SxY=(xg * Yg).sum(), k0=k0, g0=k0 / 2,
                  a0=lg["a0"], a1=lg["a1"], g1=lg["a1"] / 2,
                  R2_0=bondad(Yg, k0 * xg)["R2"], R2_1=bondad(Yg, lg["a0"] + lg["a1"] * xg)["R2"],
                  lin=lg)
# variante potencia ln v = ln c + p ln x (sin x=0)
lp = lin(np.log(xg[1:]), np.log(vg[1:]))
R["tp2e2"]["pot"] = dict(p=lp["a1"], c=np.exp(lp["a0"]), g=np.exp(lp["a0"]) ** 2 / 2)

# ---------- 8. TP2 ej 3: p = b e^{m x}
xa = np.array([1850, 1900, 1950, 1980, 2000.]); pa = np.array([1.3, 1.6, 3, 4.4, 6])
la = lin(xa, np.log(pa)); m = la["a1"]; lnb = la["a0"]; bb = np.exp(lnb)
R["tp2e3"] = dict(lnp=np.log(pa), m=m, lnb=lnb, b=bb, p1970=bb * np.exp(m * 1970), lin=la,
                  R2_log=bondad(np.log(pa), lnb + m * xa)["R2"], yap=bb * np.exp(m * xa))
# misma cosa con t = año - 1850
lt = lin(xa - 1850, np.log(pa))
R["tp2e3"]["t"] = dict(m=lt["a1"], b=np.exp(lt["a0"]), p1970=np.exp(lt["a0"]) * np.exp(lt["a1"] * 120))
R["tp2e3"]["duplica"] = np.log(2) / m

# ---------- 9. TP2 ej 4: datos.csv
d = np.loadtxt("/home/claude/fuentes/tp/datos.csv", delimiter=",", skiprows=1)
xd, yd = d.T
tab4 = []
for g in range(1, 7):
    cg = np.polyfit(xd, yd, g); bg = bondad(yd, np.polyval(cg, xd))
    tab4.append(dict(g=g, c=cg, R2=bg["R2"], Erms_cat=bg["Erms_cat"], Sr=bg["Sr"]))
R["tp2e4"] = dict(n=len(xd), xmin=xd.min(), xmax=xd.max(), tab=tab4, x=xd, y=yd)

# ---------- 10. Práctica: cuadrática por ecuaciones normales
xq = np.array([0, 1, 2, 3, 4.]); yq = np.array([2.1, 7.7, 13.6, 27.2, 40.9])
Mq = np.array([[len(xq), xq.sum(), (xq**2).sum()], [xq.sum(), (xq**2).sum(), (xq**3).sum()],
               [(xq**2).sum(), (xq**3).sum(), (xq**4).sum()]])
bq = np.array([yq.sum(), (xq * yq).sum(), (xq**2 * yq).sum()])
aq = np.linalg.solve(Mq, bq)
R["cuad"] = dict(x=xq, y=yq, M=Mq, b=bq, a=aq, R2=bondad(yq, np.polyval(aq[::-1], xq))["R2"],
                 polyfit=np.polyfit(xq, yq, 2), p5=np.polyval(aq[::-1], 5),
                 R2lin=bondad(yq, lin(xq, yq)["a0"] + lin(xq, yq)["a1"] * xq)["R2"])

# ---------- 11. Práctica: potencia y = a2 x^b2
xp = np.array([1, 2, 3, 4, 5.]); yp = np.array([0.5, 1.7, 3.4, 5.7, 8.4])
lpw = lin(np.log(xp), np.log(yp))
R["pot"] = dict(x=xp, y=yp, lnx=np.log(xp), lny=np.log(yp), b2=lpw["a1"], lna2=lpw["a0"], a2=np.exp(lpw["a0"]),
                p6=np.exp(lpw["a0"]) * 6 ** lpw["a1"], R2_log=bondad(np.log(yp), lpw["a0"] + lpw["a1"] * np.log(xp))["R2"])

# ---------- 12. Práctica: regresión múltiple chica
x1 = np.array([0, 1, 1, 2, 2, 3.]); x2 = np.array([1, 0, 2, 1, 3, 2.]); ym3 = np.array([4.1, 4.9, 8.8, 9.2, 13.1, 12.8])
Am = np.column_stack([np.ones(6), x1, x2]); Mm = Am.T @ Am; bm = Am.T @ ym3; am = np.linalg.solve(Mm, bm)
R["mult"] = dict(x1=x1, x2=x2, y=ym3, M=Mm, b=bm, a=am, p=am @ [1, 2, 2],
                 R2=bondad(ym3, Am @ am)["R2"])

# ---------- 13. Quiz: valores para preguntas numéricas
# r desde St y Sr
R["q_r"] = dict(St=200.0, Sr=50.0, r=np.sqrt((200 - 50) / 200), R2=(200 - 50) / 200)
R["q_r2"] = dict(St=80.0, Sr=20.0, r=np.sqrt((80 - 20) / 80), wrong=20 / 80, wrong2=np.sqrt(20 / 80))
# ERMS con 4 residuos
res = np.array([1, -2, 2, -1.]); E_ = (res ** 2).sum()
R["q_rms"] = dict(res=res, E=E_, n=4, cat=np.sqrt(E_ / 5), usual=np.sqrt(E_ / 4), sin_raiz=E_ / 5)
# exponencial linealizada: ln y = 0.7 + 0.3 x
R["q_exp"] = dict(lna=0.7, b=0.3, a=np.exp(0.7), y5=np.exp(0.7) * np.exp(0.3 * 5), y5_mal=0.7 + 0.3 * 5)
# potencia: ln y = 1.2 + 1.5 ln x
R["q_pot"] = dict(lna=1.2, b=1.5, a=np.exp(1.2))
# predicción con el ejemplo de la diapositiva en x = 8
R["q_ej8"] = dict(p=-0.734 + 0.739 * 8, swap=0.739 - 0.734 * 8)
# predicción múltiple quiz: y = 161.3365 + 32.9688u - 0.0855v con u=0.10, v=1200
R["q_acero"] = dict(dur=161.3365 + 32.9688 * 0.10 - 0.0855 * 1200, sin_u=161.3365 - 0.0855 * 1200,
                    signo=161.3365 + 32.9688 * 0.10 + 0.0855 * 1200)

# ---------- 14. Datasets sintéticos para figuras (fijos)
def datos():
    D = {}
    rng = np.random.default_rng(7)
    x = np.sort(rng.uniform(0.05, 10, 160))
    y = -300 + 700 * np.log(x + 0.25) + rng.normal(0, 1, x.size) * (25 + 45 * x)
    D["it3"] = (x, y)
    rng = np.random.default_rng(3)
    x = np.linspace(0, 10, 40); D["expo"] = (x, 2 * np.exp(0.55 * x) * (1 + rng.normal(0, 0.06, x.size)))
    rng = np.random.default_rng(4)
    x = np.linspace(0, 10, 30); D["parab"] = (x, 100 - 0.9 * x ** 2 + rng.normal(0, 2.5, x.size))
    rng = np.random.default_rng(5)
    x = np.linspace(1, 30, 30); D["ruido"] = (x, rng.normal(0, 0.45, x.size))
    rng = np.random.default_rng(6)
    x = np.linspace(0, 10, 25); D["curva"] = (x, 3 + 0.4 * x ** 2 + rng.normal(0, 1.5, x.size))
    rng = np.random.default_rng(8)
    x = np.linspace(1, 20, 25); D["pot"] = (x, 1.5 * x ** 1.8 * (1 + rng.normal(0, 0.05, x.size)))
    rng = np.random.default_rng(9)
    x = np.linspace(0.5, 10, 25); D["recip"] = (x, 1 / (0.4 * x + 0.5) * (1 + rng.normal(0, 0.03, x.size)))
    rng = np.random.default_rng(10)
    x = np.linspace(0, 10, 30); D["lin"] = (x, 2 + 1.5 * x + rng.normal(0, 1.2, x.size))
    rng = np.random.default_rng(11)
    x = np.linspace(0, 6.28 * 2, 40); D["seno"] = (x, np.sin(x) + rng.normal(0, 0.1, x.size))
    return D

D = datos()
F = {}
for k, (xx, yy) in D.items():
    l = lin(xx, yy); bo = bondad(yy, l["a0"] + l["a1"] * xx)
    F[k] = dict(a0=l["a0"], a1=l["a1"], r=bo["r"], R2=bo["R2"], pearson=np.corrcoef(xx, yy)[0, 1])
# residuos de la recta sobre 'curva' y R2 cuadrático
cc = np.polyfit(*D["curva"], 2); F["curva"]["R2cuad"] = bondad(D["curva"][1], np.polyval(cc, D["curva"][0]))["R2"]
cl = np.polyfit(np.log(D["it3"][0]), D["it3"][1], 1)
F["it3"]["R2log"] = bondad(D["it3"][1], np.polyval(cl, np.log(D["it3"][0])))["R2"]
cq3 = np.polyfit(*D["it3"], 2); F["it3"]["R2cuad"] = bondad(D["it3"][1], np.polyval(cq3, D["it3"][0]))["R2"]
F["it3"]["cuad_a2"] = cq3[0]
le2 = lin(D["expo"][0], np.log(D["expo"][1])); F["expo"]["R2exp"] = bondad(D["expo"][1], np.exp(le2["a0"] + le2["a1"] * D["expo"][0]))["R2"]
F["expo"]["b_exp"] = le2["a1"]
cp = np.polyfit(*D["parab"], 2); F["parab"]["R2cuad"] = bondad(D["parab"][1], np.polyval(cp, D["parab"][0]))["R2"]; F["parab"]["a2"] = cp[0]
R["figs"] = F

if __name__ == "__main__":
    np.set_printoptions(precision=6, suppress=True)
    def show(k, v, ind=""):
        if isinstance(v, dict):
            print(f"{ind}[{k}]")
            for kk, vv in v.items():
                if kk in ("x", "y"): continue
                show(kk, vv, ind + "  ")
        elif isinstance(v, list):
            print(f"{ind}{k}:")
            for it in v: show("-", it, ind + "  ")
        else:
            print(f"{ind}{k} = {v}")
    for k, v in R.items(): show(k, v)
