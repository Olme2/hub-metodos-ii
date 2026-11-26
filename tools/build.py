#!/usr/bin/env python3
"""Arma index.html (documento completo, para GitHub Pages) a partir de src/app.html."""
import os
root = os.path.join(os.path.dirname(__file__), "..")
import json, re
frag = open(os.path.join(root, "src", "app.html"), encoding="utf-8").read()
cfg_src = open(os.path.join(root, "data", "config.js"), encoding="utf-8").read().strip()
cfg = json.loads(re.match(r"^HUB\.configurar\((.*)\);?\s*$", cfg_src, re.S).group(1))
frag = frag.replace("{{TITULO}}", cfg["titulo"])
head_end = frag.index("</style>") + len("</style>")
head, body = frag[:head_end], frag[head_end:]
doc = f"""<!doctype html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<style>:root{{padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px);color-scheme:light}}[hidden]{{display:none!important}}img{{max-width:100%}}</style>
{head}
</head>
<body>{body}
</body>
</html>
"""
open(os.path.join(root, "index.html"), "w", encoding="utf-8").write(doc)
print("index.html listo")
