#!/usr/bin/env python3
"""Verificación numérica de la Unidad 3 (Aproximación trigonométrica, Fourier, FFT, Nyquist, filtros).
Todo número que aparece en data/m3.js sale de acá. Se puede importar (gen_u3.py usa R y las funciones)
o correr directo para ver la salida:  python3 tools/gen/verif_u3.py
"""
import os, math, timeit
import numpy as np
from scipy.signal import find_peaks

FUENTES = "/home/claude/fuentes/tp"
R = {}          # resultados que usa el generador
PI = math.pi


def r(x, d=3):
    return float(round(float(x), d))


def picos(fr, amp, alto=0.15):
    """Picos locales de un espectro normalizado (máx = 1)."""
    a = amp / amp.max()
    p, _ = find_peaks(a, height=alto)
    return [(r(fr[i], 4), r(a[i], 2)) for i in p]


def alias(f, fs):
    """Frecuencia aparente (en [0, fs/2]) de una senoidal de f Hz muestreada a fs."""
    return abs(f - fs * round(f / fs))


# ---------------------------------------------------------------- 1. Ejemplo p12–13: s(t) = sin t + 1.2 cos 2t
R["ej_s"] = dict(A1=1.0, w1=1.0, f1=r(1 / (2 * PI), 4), A2=1.2, w2=2.0, f2=r(2 / (2 * PI), 4),
                 T1=r(2 * PI, 4), T2=r(PI, 4))
t = np.arange(0, 200 * PI, 0.05)
s = np.sin(t) + 1.2 * np.cos(2 * t)
S = np.abs(np.fft.rfft(s)) * 2 / len(s)
fr = np.fft.rfftfreq(len(s), 0.05)
R["ej_s"]["picos_fft"] = [(r(fr[i], 3), r(S[i], 2)) for i in find_peaks(S, height=0.5)[0]]

# ---------------------------------------------------------------- 2. Escalón (onda cuadrada) p5: coeficientes y sumas parciales
def sq(x):
    return np.sign(np.sin(x))

xs = np.linspace(-PI, PI, 200001)
bn = []
for n in range(1, 8):
    an_ = np.trapezoid(sq(xs) * np.cos(n * xs), xs) / PI
    bn_ = np.trapezoid(sq(xs) * np.sin(n * xs), xs) / PI
    bn.append((n, r(an_, 4), r(bn_, 4), r(4 / (n * PI), 4) if n % 2 else 0.0))
R["escalon_coef"] = bn
a0 = np.trapezoid(sq(xs), xs) / PI
R["escalon_a0"] = r(a0, 4)
xf = np.linspace(1e-4, PI - 1e-4, 200001)
sobre = []
for K in [1, 2, 3, 5, 10, 50, 200]:          # K términos impares
    ss = sum(4 / ((2 * m - 1) * PI) * np.sin((2 * m - 1) * xf) for m in range(1, K + 1))
    sobre.append((K, r(ss.max(), 4), r(ss.max() - 1, 4)))
R["escalon_gibbs"] = sobre

# ---------------------------------------------------------------- 3. Período de sin/cos (p6 dice T = Π)
R["periodo_sin"] = r(2 * PI, 4)
R["orto"] = r(np.trapezoid(np.sin(xs) * np.cos(xs), xs), 6)

# ---------------------------------------------------------------- 4. TDF a mano, N = 4, datos de la p18
x4 = np.array([1.2, 1.4, 1.3, 1.0])
N = 4
X4 = [sum(x4[n] * np.exp(-2j * PI * k * n / N) for n in range(N)) for k in range(N)]
R["tdf4_x"] = x4.tolist()
R["tdf4_X"] = [(r(z.real, 4), r(z.imag, 4)) for z in X4]
R["tdf4_ok"] = bool(np.allclose(X4, np.fft.fft(x4)))
R["tdf4_w"] = [(k, r(np.exp(-2j * PI * k / N).real, 3), r(np.exp(-2j * PI * k / N).imag, 3)) for k in range(4)]
R["tdf4_ops"] = dict(mult=N * N, sumas=N * (N - 1))
R["tdf4_ifft"] = [r(v, 4) for v in np.fft.ifft(np.fft.fft(x4)).real]

# ---------------------------------------------------------------- 5. Costos O(N²) vs O(N log2 N) (p19–20)
cost = []
for Nn in [1000, 10**6, 10**9]:
    n2 = Nn ** 2
    nl = Nn * math.log2(Nn)
    cost.append(dict(N=Nn, N2=n2, NlogN=nl, seg_N2=n2 * 1e-9, anios_N2=n2 * 1e-9 / (365.25 * 24 * 3600),
                     seg_NlogN=nl * 1e-9, cociente=n2 / nl))
R["costos"] = cost

# ---------------------------------------------------------------- 6. Cooley–Tukey: pares/impares, periodicidad, simetría
x8 = np.array([3., 1., 4., 1., 5., 9., 2., 6.])
N8 = 8
w = lambda N_, e: np.exp(-2j * PI * e / N_)
E = np.fft.fft(x8[0::2]); O = np.fft.fft(x8[1::2])
Xct = np.array([E[k % 4] + w(8, k) * O[k % 4] for k in range(8)])
R["ct_ok"] = bool(np.allclose(Xct, np.fft.fft(x8)))
R["ct_mitad"] = bool(np.allclose([Xct[k + 4] for k in range(4)], [E[k] - w(8, k) * O[k] for k in range(4)]))
R["ct_niveles"] = int(math.log2(N8))
R["ident_wNkN"] = r(abs(w(8, 3 * 8) - 1), 12)
R["periodicidad"] = bool(np.isclose(w(8, 3 * 5), w(8, 3 * (5 + 8))))
R["simetria"] = bool(np.isclose(w(8, 3 * (8 - 5)), np.conj(w(8, 3 * 5))))
R["w_2kn"] = bool(np.isclose(w(8, 2 * 3 * 1), w(4, 3 * 1)))


def fft_rec(x):
    """FFT radix-2 recursiva (decimación en el tiempo), para contar multiplicaciones complejas."""
    n = len(x)
    if n == 1:
        return x.astype(complex), 0
    e, me = fft_rec(x[0::2]); o, mo = fft_rec(x[1::2])
    tw = np.exp(-2j * PI * np.arange(n // 2) / n) * o
    return np.concatenate([e + tw, e - tw]), me + mo + n // 2


mults = []
for Nn in [8, 16, 1024]:
    xx = np.random.default_rng(0).standard_normal(Nn)
    Xr, m = fft_rec(xx)
    mults.append((Nn, Nn * Nn, m, int(Nn / 2 * math.log2(Nn)), bool(np.allclose(Xr, np.fft.fft(xx)))))
R["mults"] = mults

# ---------------------------------------------------------------- 7. tiempo_fft.ipynb (TP3 ej. 2)
import sympy
tiempos = []
for n in [166464, 166400, 2**16, 2**17, 1 + 2**16]:
    a = np.zeros(n)
    tt = min(timeit.repeat(lambda: np.fft.fft(a, len(a)), number=5, repeat=7)) / 5
    tiempos.append(dict(n=n, fact=sympy.factorint(n), ms=round(tt * 1e3, 2)))
R["tiempos"] = tiempos
R["numpy"] = np.__version__

# padding: siguiente potencia de 2
R["pad"] = [(n, 1 << (n - 1).bit_length()) for n in [1000, 1886, 39, 166464, 65537]]

# ---------------------------------------------------------------- 8. Nyquist
R["EV_P2"] = dict(fs=10, fN=10 / 2)
R["nyq_f2"] = dict(f=2, fs_min=2 * 2)

# muestreo.ipynb: f1 = 1, f2 = 2 Hz, desfase π/2, t en [0, 2]
mu = []
for fsm in [60, 10, 4, 2, 1]:
    tt = np.arange(0, 2 + 1 / fsm, 1 / fsm)
    fila = dict(fs=fsm, fN=fsm / 2, n=len(tt))
    for f in (1, 2):
        y = np.sin(2 * PI * f * tt + PI / 2)
        fila[f"f{f}_alias"] = r(alias(f, fsm), 3)
        fila[f"f{f}_muestras"] = sorted(set(np.round(y, 3).tolist()))[:6]
        fila[f"f{f}_por_periodo"] = r(fsm / f, 2)
    mu.append(fila)
R["muestreo"] = mu
# celda 2: suma de 1, 2, 5 Hz con fs = 60
fsm = 60
tt = np.arange(0, 2 + 1 / fsm, 1 / fsm)
y = sum(np.sin(2 * PI * f * tt + PI / 2) for f in (1, 2, 5))
Y = np.abs(np.fft.fft(y)); frq = np.fft.fftfreq(len(y), 1 / fsm); c = frq > 0
R["muestreo_c2"] = dict(n=len(tt), df=r(fsm / len(tt), 4), fN=fsm / 2, picos=picos(frq[c], Y[c], 0.5))

# Resolución de frecuencias (nyquist p4–6)
res = []
for (a_, b_, fsm) in [(0, 2 * PI, 10), (0, 2 * PI, 100), (0, 6 * PI, 6)]:
    tt = np.arange(a_, b_, 1.0 / fsm)
    n = tt.size
    y = sum(np.sin(2 * PI * f * tt) for f in (0.5, 1, 2))
    Y = np.abs(np.fft.fft(y)); frq = np.fft.fftfreq(n, 1.0 / fsm); c = frq > 0
    res.append(dict(b=r(b_, 3), fs=fsm, n=n, df=r(fsm / n, 3), fmax=r(frq[c].max(), 3), fN=fsm / 2,
                    picos=picos(frq[c], Y[c], 0.4)))
R["resolucion"] = res

# ET25-T3 ej. 4 / EV-P8: f = 10, 20, 100 con fs = 20
k = np.arange(40)
ev8 = []
for f in (10, 20, 100):
    sn = np.sin(2 * PI * f * k / 20); cs = np.cos(2 * PI * f * k / 20)
    ev8.append(dict(f=f, alias=alias(f, 20), rel=("<" if f < 10 else "=" if f == 10 else ">") + " fN",
                    sin_max=r(np.abs(sn).max(), 6), cos_vals=sorted(set(np.round(cs, 6).tolist()))))
R["EV_P8"] = ev8
# fase en el límite f = fs/2: amplitud aparente de sin(2π·10·t + φ) muestreada a 20 Hz
R["limite_fase"] = [(fi, r(np.abs(np.sin(2 * PI * 10 * k / 20 + fi * PI / 180)).max(), 3)) for fi in (0, 30, 45, 90)]
# con fs = 21 (apenas > 2·10): 10 Hz se recupera, 20 → alias 1
R["fs21"] = [(f, r(alias(f, 21), 3)) for f in (10, 20, 100)]

# ET25-T3 ej. 6: f = 10, 15, 100
R["ET_ej6"] = [(fs_, fs_ / 2, [(f, r(alias(f, fs_), 2)) for f in (10, 15, 100)]) for fs_ in (10, 20, 100, 200)]

# EV-P3 / tabla de alias genérica
R["alias_tabla"] = [(f, fs_, r(alias(f, fs_), 3)) for (f, fs_) in [(3, 4), (6, 10), (9, 10), (23, 10), (23, 20),
                                                                  (1.2, 2), (3, 2), (20, 2), (100, 20), (60, 50)]]

# ---------------------------------------------------------------- 9. TP3 ej. 3: muestra1 / muestra2
def cargar(nombre):
    d = np.loadtxt(os.path.join(FUENTES, nombre), delimiter=",")
    return d[:, 0], d[:, 1]

t1, y1 = cargar("muestra1.csv"); t2, y2 = cargar("muestra2.csv")


def senial_tp3(t):
    """Señal que reproduce ambos CSV (encontrada por mínimos cuadrados; error máx ~1e-14)."""
    return (5 * np.cos(2 * PI * 1.0 * t) + 5 * np.cos(2 * PI * 1.2 * t) + 10 * np.sin(2 * PI * 3.0 * t)
            + np.cos(2 * PI * 20 * t) - 2.5 * np.sin(2 * PI * 0.125 * t))


def espectro(t, y, normal=True):
    n = len(t); dt = t[1] - t[0]
    Y = np.abs(np.fft.fft(y)); frq = np.fft.fftfreq(n, dt); c = frq > 0
    a = Y[c] / (Y[c].max() if normal else 1)
    return frq[c], a


mm = {}
for nom, t_, y_ in (("m1", t1, y1), ("m2", t2, y2)):
    dt = t_[1] - t_[0]
    frq, a = espectro(t_, y_)
    # picos con buena resolución (zero padding sólo para ubicarlos)
    NP = 2 ** 18
    Yp = np.abs(np.fft.rfft(y_, NP)); fp = np.fft.rfftfreq(NP, dt)
    mm[nom] = dict(n=len(t_), dt=r(dt, 3), t0=r(t_[0], 3), tf=r(t_[-1], 3), fs=r(1 / dt, 3), fN=r(1 / (2 * dt), 3),
                   df=r(1 / (len(t_) * dt), 4), picos_bins=picos(frq, a, 0.2), picos_finos=picos(fp, Yp, 0.15),
                   err_modelo=float(np.abs(senial_tp3(t_) - y_).max()))
R["muestras"] = mm
R["tp3_comp"] = [(0.125, "−2,5·sin"), (1.0, "5·cos"), (1.2, "5·cos"), (3.0, "10·sin"), (20.0, "1·cos")]
R["tp3_alias_m1"] = [(f, r(alias(f, 2.0), 3)) for f, _ in R["tp3_comp"]]
kk = np.arange(len(t1))
R["m1_3Hz_max"] = float(np.abs(10 * np.sin(2 * PI * 3 * 0.5 * kk)).max())
R["m1_20Hz_vals"] = sorted(set(np.round(np.cos(2 * PI * 20 * 0.5 * kk), 6).tolist()))
R["m1_1Hz_vals"] = sorted(set(np.round(5 * np.cos(2 * PI * 1 * 0.5 * kk), 6).tolist()))

# ---------------------------------------------------------------- 10. EV-P5 recreado: tercera componente 23 Hz
def senial_p5(t):
    return np.sin(2 * PI * 0.7 * t) + np.sin(2 * PI * 1.25 * t) + np.sin(2 * PI * 23 * t)

p5 = {}
for nom, fs_ in (("m1", 10.0), ("m2", 60.0)):
    tt = np.arange(0, 20, 1 / fs_)
    frq, a = espectro(tt, senial_p5(tt))
    p5[nom] = dict(fs=fs_, fN=fs_ / 2, picos=picos(frq, a, 0.3))
R["EV_P5"] = p5
R["EV_P5_alias"] = [(fs_, r(alias(23, fs_), 2)) for fs_ in (8, 10, 20, 50, 60)]

# ---------------------------------------------------------------- 11. Señales filtradas (ET25-T3 ej. 5, EV-P6)
def filtrar_1d(y, dt, fc, tipo):
    Y = np.fft.fft(y); frq = np.fft.fftfreq(len(y), dt)
    m = (np.abs(frq) < fc) if tipo == "bajo" else (np.abs(frq) >= fc)
    return np.fft.ifft(Y * m).real


rng = np.random.default_rng(3)
tS = np.arange(0, 20, 0.01)
lento = 2 * np.sin(2 * PI * 0.15 * tS) + np.sin(2 * PI * 0.4 * tS)
ruido = 0.6 * rng.standard_normal(tS.size)
ej5_orig = lento + ruido
ej5_pb = filtrar_1d(ej5_orig, 0.01, 1.0, "bajo")
ej5_pa = filtrar_1d(ej5_orig, 0.01, 1.0, "alto")
R["ej5"] = dict(fc=1.0, err_pb=r(np.sqrt(np.mean((ej5_pb - lento) ** 2)), 3), std_ruido=r(ruido.std(), 3),
                err_pa=r(np.sqrt(np.mean((ej5_pa - lento) ** 2)), 3))
ev6_lento = 2 * np.sin(2 * PI * 0.12 * tS); ev6_rap = np.sin(2 * PI * 3 * tS)
ev6_orig = ev6_lento + ev6_rap
ev6_pa = filtrar_1d(ev6_orig, 0.01, 1.0, "alto")
ev6_pb = filtrar_1d(ev6_orig, 0.01, 1.0, "bajo")
R["EV_P6"] = dict(fc=1.0, amp_pa=r(np.abs(ev6_pa[200:-200]).max(), 2), err_pa=r(np.sqrt(np.mean((ev6_pa - ev6_rap) ** 2)), 4),
                  amp_pb=r(np.abs(ev6_pb[200:-200]).max(), 2))

# ---------------------------------------------------------------- 12. Imágenes (notebook imagen_con_filtro)
def filtro_pasa_bajo(D0, shape):
    rows, cols = shape
    base = np.zeros((rows, cols))
    for x in range(cols):
        for yy in range(rows):
            if math.hypot(yy - rows / 2, x - cols / 2) < D0:
                base[yy, x] = 1
    return base


def filtro_pasa_alto(D0, shape):
    return 1 - filtro_pasa_bajo(D0, shape)


def aplicar(img, filtro):
    c2 = np.fft.fft2(img); c3 = np.fft.fftshift(c2)
    c4 = np.fft.ifftshift(c3 * filtro); c5 = np.fft.ifft2(c4)
    return np.abs(c5), np.log(1 + np.abs(c3))


def tablero(n=128, celda=8, periodo=4):
    """Tablero sintético parecido al del examen: damero gris claro/blanco + diagonales oscuras."""
    i, j = np.indices((n, n)) // celda
    img = np.where((i + j) % 2 == 0, 0.78, 0.97)
    img = np.where((i - j) % periodo == 0, 0.30, img)
    return img


def figura_sintetica(n=128):
    yy, xx = np.indices((n, n))
    img = np.full((n, n), 0.25)
    img[(xx - 40) ** 2 + (yy - 45) ** 2 < 22 ** 2] = 0.85
    img[70:115, 70:118] = 0.6
    img[20:40, 80:120] = 0.95
    img += 0.08 * np.sin(2 * PI * xx / 64)
    return img


TB = tablero()
F = np.fft.fftshift(np.abs(np.fft.fft2(TB)))
Fc = F.copy(); Fc[64, 64] = 0
pk = np.argwhere(Fc > 0.25 * Fc.max())
R["tablero_picos"] = sorted({(int(a_ - 64), int(b_ - 64), r(math.hypot(a_ - 64, b_ - 64), 2)) for a_, b_ in pk},
                            key=lambda z: z[2])[:12]
D0_PB, D0_PA = 8, 24
tb_pb, tb_spec = aplicar(TB, filtro_pasa_bajo(D0_PB, TB.shape))
tb_pa, _ = aplicar(TB, filtro_pasa_alto(D0_PA, TB.shape))
ii, jj = np.indices(TB.shape)
diag = np.cos(2 * PI * (ii - jj) / 32)
R["tablero_pb_corr_diag"] = r(np.corrcoef((tb_pb - tb_pb.mean()).ravel(), diag.ravel())[0, 1], 3)
dist = np.minimum(np.minimum(ii % 8, 7 - ii % 8), np.minimum(jj % 8, 7 - jj % 8))   # distancia al borde de la casilla
R["tablero_pa_borde_vs_int"] = (r(tb_pa[dist == 0].mean(), 3), r(tb_pa[dist >= 2].mean(), 3))
tb_pa8, _ = aplicar(TB, filtro_pasa_alto(8, TB.shape))
R["tablero_pa8_borde_vs_int"] = (r(tb_pa8[dist == 0].mean(), 3), r(tb_pa8[dist >= 2].mean(), 3))
R["tablero_pa_media"] = r(tb_pa.mean(), 3)
R["tablero_pb_rango"] = (r(tb_pb.min(), 3), r(tb_pb.max(), 3))
R["D0"] = (D0_PB, D0_PA)
# zero padding de la imagen: la imagen reconstruida es la original más ceros (no la suaviza)
pad = np.abs(np.fft.ifft2(np.fft.fft2(TB, s=(256, 256))))
R["pad_img_ok"] = bool(np.allclose(pad[:128, :128], TB)) and bool(np.allclose(pad[128:, :], 0))

FS = figura_sintetica()
fs_pa, fs_spec = aplicar(FS, filtro_pasa_alto(10, FS.shape))
fs_pb, _ = aplicar(FS, filtro_pasa_bajo(10, FS.shape))
R["sint_pa_media"] = r(fs_pa.mean(), 3); R["sint_media"] = r(FS.mean(), 3)

# fft2 = fft por filas y después por columnas
A_ = np.random.default_rng(1).standard_normal((8, 6))
R["fft2_filas_cols"] = bool(np.allclose(np.fft.fft(np.fft.fft(A_, axis=1), axis=0), np.fft.fft2(A_)))
# energía en bajas frecuencias: imagen suave vs tablero
yy, xx = np.indices((128, 128))
suave = np.exp(-((xx - 64) ** 2 + (yy - 64) ** 2) / (2 * 30 ** 2))
def frac_baja(img, D0=8):
    P = np.abs(np.fft.fftshift(np.fft.fft2(img - img.mean()))) ** 2
    return r((P * filtro_pasa_bajo(D0, img.shape)).sum() / P.sum(), 3)
R["frac_baja"] = dict(suave=frac_baja(suave), tablero=frac_baja(TB))
# ubicación de la frecuencia cero antes y después de fftshift
R["dc_sin_shift"] = (0, 0); R["dc_con_shift"] = (64, 64)

# rayas verticales -> picos en el eje horizontal del espectro
rayas = 0.5 + 0.5 * np.sin(2 * PI * xx / 16)
Fr = np.fft.fftshift(np.abs(np.fft.fft2(rayas))); Fr[64, 64] = 0
R["rayas_picos"] = [(int(a_ - 64), int(b_ - 64)) for a_, b_ in np.argwhere(Fr > 0.5 * Fr.max())]


if __name__ == "__main__":
    import pprint
    for kname, v in R.items():
        print(f"== {kname}")
        pprint.pprint(v, width=130, compact=True)
