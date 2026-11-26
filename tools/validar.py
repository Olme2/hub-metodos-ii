#!/usr/bin/env python3
"""Valida un archivo de unidad data/mN.js contra el estándar del hub.
Lee data/config.js (preguntasPorUnidad, proporcionVerdaderas) si existe.

Uso: python3 tools/validar.py data/m2.js [data/m3.js ...]
Sale con código 1 si hay errores bloqueantes.
"""
import json, re, sys, statistics

PROHIBIDAS = [r"mir(a|á)r? (esto )?en el pdf", r"consult(a|á|ar) el pdf", r"ver el pdf", r"conviene mirar"]

def cargar_config(path_unidad):
    import os
    c = os.path.join(os.path.dirname(os.path.abspath(path_unidad)), "config.js")
    if not os.path.exists(c):
        return {}
    src = open(c, encoding="utf-8").read().strip()
    m = re.match(r"^HUB\.configurar\((.*)\);?\s*$", src, re.S)
    return json.loads(m.group(1)) if m else {}

def cargar(path):
    src = open(path, encoding="utf-8").read().strip()
    m = re.match(r"^ASO\.addUnit\((.*)\);?\s*$", src, re.S)
    if not m:
        raise SystemExit(f"{path}: el archivo debe ser exactamente ASO.addUnit({{...JSON...}});")
    return json.loads(m.group(1))

def plano(html):
    return re.sub(r"<[^>]+>", "", html or "")

def validar(path):
    errs, warns = [], []
    u = cargar(path)
    for k in ["id", "titulo", "corto", "fuente", "peso", "intro", "secciones", "conexiones",
              "fidelidad", "auditoria", "glosario", "quiz", "desarrollo"]:
        if k not in u:
            errs.append(f"falta la clave '{k}'")
    if errs:
        return errs, warns, u
    n = u["id"]
    sec_ids = [s["id"] for s in u["secciones"]]
    if len(set(sec_ids)) != len(sec_ids):
        errs.append("ids de sección duplicados")
    for s in u["secciones"]:
        if not s["id"].startswith(f"m{n}-s"):
            errs.append(f"id de sección mal formado: {s['id']}")
        if s.get("estrellas") not in (0, 1, 2, 3):
            errs.append(f"{s['id']}: estrellas debe ser 0..3")
        if len(plano(s.get("html", ""))) < 200:
            warns.append(f"{s['id']} '{s.get('titulo')}': menos de 200 caracteres de contenido")
    todo = json.dumps(u, ensure_ascii=False).lower()
    for p in PROHIBIDAS:
        if re.search(p, todo):
            errs.append(f"frase prohibida (derivar al PDF): /{p}/")
    for a in u["auditoria"]:
        if a.get("estado") == "cubierto" and a.get("seccion") not in sec_ids:
            errs.append(f"auditoría apunta a sección inexistente: {a}")
    for g in u["glosario"]:
        if g.get("seccion") and g["seccion"] not in sec_ids:
            errs.append(f"glosario '{g['termino']}' apunta a sección inexistente")
    for k in ("fuentes", "agregado"):
        if not u["fidelidad"].get(k):
            errs.append(f"fidelidad.{k} vacío")

    q = u["quiz"]
    cfg = cargar_config(path)
    objetivo = cfg.get("preguntasPorUnidad", 40)
    if len(q) != objetivo:
        errs.append(f"el quiz tiene {len(q)} preguntas (deben ser {objetivo})")
    ids = [x["id"] for x in q]
    if len(set(ids)) != len(ids):
        errs.append("ids de pregunta duplicados")
    largo_correcta = 0
    prop_correctas = []
    lc, li = [], []
    verdaderas = 0
    for x in q:
        if x.get("pide") not in ("verdaderas", "falsas", "una"):
            errs.append(f"{x['id']}: pide debe ser 'verdaderas', 'falsas' o 'una'")
            continue
        verdaderas += x["pide"] == "verdaderas"
        if x.get("seccion") not in sec_ids:
            errs.append(f"{x['id']}: seccion '{x.get('seccion')}' no existe")
        ops = x.get("opciones", [])
        if not 4 <= len(ops) <= 5:
            errs.append(f"{x['id']}: debe tener 4 o 5 opciones")
        corr = [o for o in ops if o["v"] == (x["pide"] != "falsas")]
        if x["pide"] == "una" and len(corr) != 1:
            errs.append(f"{x['id']}: pide 'una' exige exactamente una opción con v=true (hay {len(corr)})")
        if x["pide"] != "una" and not 1 <= len(corr) <= len(ops) - 1:
            errs.append(f"{x['id']}: debe haber entre 1 y {len(ops)-1} opciones correctas (hay {len(corr)})")
        if re.search(r"(opci[oó]n|letra)\s+\(?[a-e]\)?\b|(^|[\s(])[a-e]\)\s", x.get("explicacion", ""), re.I):
            errs.append(f"{x['id']}: la explicación nombra opciones por letra (se mezclan en cada intento)")
        if len(x.get("explicacion", "")) < 200:
            errs.append(f"{x['id']}: explicación de {len(x.get('explicacion',''))} caracteres (mínimo 200)")
        if ops:
            mas_larga = max(ops, key=lambda o: len(o["t"]))
            largo_correcta += mas_larga in corr
            prop_correctas.append(len(corr) / len(ops))
            lc += [len(o["t"]) for o in corr]
            li += [len(o["t"]) for o in ops if o not in corr]
    nq = max(len(q), 1)
    pct_larga = 100 * largo_correcta / nq
    esperado = 100 * statistics.mean(prop_correctas) if prop_correctas else 0
    ratio = (statistics.mean(lc) / statistics.mean(li)) if lc and li else 0
    pct_v = 100 * verdaderas / nq
    if abs(pct_larga - esperado) > 5:
        errs.append(f"balance: la más larga es correcta en {pct_larga:.0f}% (esperado {esperado:.0f}% ±5)")
    if not 0.9 <= ratio <= 1.1:
        errs.append(f"balance: largo medio correctas/incorrectas = {ratio:.2f} (debe estar entre 0,90 y 1,10)")
    meta_v = cfg.get("proporcionVerdaderas")  # p. ej. 67 si el patrón real pide 2/3 "verdaderas"
    if meta_v is not None and abs(pct_v - meta_v) > 10:
        warns.append(f"proporción 'verdaderas' = {pct_v:.0f}% (patrón real ≈ {meta_v}%)")
    for it in list(q) + list(u["desarrollo"]):
        f = it.get("fig")
        if f is None: continue
        if not (isinstance(f, str) and f.lstrip().startswith("<svg") and f.rstrip().endswith("</svg>")):
            errs.append(f"{it.get('id')}: fig debe ser un <svg>…</svg> completo")
        elif re.search(r"<script|on\w+=|<image[^>]+href=\"http", f, re.I):
            errs.append(f"{it.get('id')}: fig contiene script, eventos o recursos externos")
        elif len(f) > 150000:
            warns.append(f"{it.get('id')}: fig de {len(f)//1000} KB (conviene < 150 KB)")
    for d in u["desarrollo"]:
        if d.get("seccion") not in sec_ids:
            errs.append(f"desarrollo {d.get('id')}: sección inexistente")
        if len(plano(d.get("modelo", ""))) < 250:
            warns.append(f"desarrollo {d.get('id')}: respuesta modelo corta")
    metr = dict(preguntas=len(q), pct_larga=round(pct_larga), esperado=round(esperado),
                ratio=round(ratio, 2), pct_verdaderas=round(pct_v), secciones=len(sec_ids),
                glosario=len(u["glosario"]), desarrollo=len(u["desarrollo"]),
                practica=sum(1 for d in u["desarrollo"] if d.get("practica")), figs=sum(1 for it in list(q)+list(u["desarrollo"]) if it.get("fig")))
    u["_metricas"] = metr
    return errs, warns, u

if __name__ == "__main__":
    mal = False
    for p in sys.argv[1:]:
        e, w, u = validar(p)
        print(f"== {p}: {u.get('_metricas')}")
        for x in e: print("  ERROR", x)
        for x in w: print("  aviso", x)
        mal |= bool(e)
        if not e: print("  OK")
    sys.exit(1 if mal else 0)
