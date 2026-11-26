# Consola MN II — hub de estudio

Hub de estudio de Métodos Numéricos II (FACET-UNT) para el último parcial (1/12/2026): Ajuste de curvas (U2), Aproximación trigonométrica y Fourier (U3) y EDO (U4). Resúmenes con auditoría contra las diapositivas, quizzes con el formato de la cátedra, desarrollo (teoría y práctica), glosario y progreso.

## Abrir

Es un sitio estático: `index.html` + `data/`. Se publica en GitHub Pages (Settings → Pages → `main` / `(root)`), o se abre localmente con `python3 -m http.server` desde esta carpeta.

## Regenerar o ampliar

```
cd tools/gen
python3 verif_u4.py        # todos los números y órdenes, verificados (sympy/numpy)
python3 gen_u4.py          # escribe data/m4.js (incluye burden_u4.py)
cd ../..
python3 tools/validar.py data/m*.js
python3 tools/build.py
NODE_PATH=<carpeta con jsdom>/node_modules node tools/smoke.js .
```

Para mejorar una unidad se edita su generador y se compara la salida con la versión anterior antes de reemplazarla. No se reescribe desde cero.

## Fuentes

Diapositivas y TP de la cátedra 2026, exámenes 2025 y finales 2018/2020, Burden-Faires 10.ª ed. cap. 5 (ejercicios de práctica de EDO). Lo que no sale del material está marcado como "Agregado".
