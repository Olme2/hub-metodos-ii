# Especificación de una unidad — Métodos Numéricos II (último parcial 2026)


Vas a producir **un archivo**: `/home/claude/hub/data/mN.js` (N = número de unidad), con el resumen, el quiz, las preguntas de desarrollo y el glosario de esa unidad. Lo usa Valentino (Ing. Informática, FACET-UNT) para rendir **el último parcial (integral teórico escrito + evaluación práctica)** de **Métodos Numéricos II** el **martes 1/12/2026 al mediodía**. Todo en **español rioplatense (vos), denso, sin relleno**.

Antes de escribir, leé completo `/home/claude/calibracion/02-patron-evaluacion.md` (qué toma la cátedra y cómo arma los distractores). **No hay inventario previo: lo hacés vos** (ver "Índice oficial" abajo) y lo dejás en `/home/claude/calibracion/01-inventario-uN.md`.

## Fuente

- U2 Ajuste de curvas: `/home/claude/fuentes/ajuste.pdf` (23 diap.). TP: `/home/claude/fuentes/tp2.pdf`, datos `/home/claude/fuentes/tp/datos.csv`.
- U3 Aproximación trigonométrica: `/home/claude/fuentes/trig.pdf` (27 diap.) + `/home/claude/fuentes/nyquist.pdf` (10 diap., casi todo imagen). TP: `/home/claude/fuentes/tp3.pdf`, `/home/claude/fuentes/tp/muestra1.csv`, `muestra2.csv`, `/home/claude/fuentes/tp/notebooks-TP3.md`.
- U4 EDO: `/home/claude/fuentes/edo1.pdf` (17), `edo2.pdf` (30), `edo3.pdf` (18), `edo4.pdf` (18). TP: `/home/claude/fuentes/tp4.pdf`.
- Exámenes: `/home/claude/fuentes/parciales.pdf` (17 pág., imágenes) — ya transcriptos en el patrón; abrilos sólo si necesitás ver un gráfico.

### Índice oficial (inventario, lo hacés vos)
Casi **todas** las diapositivas tienen fórmulas o contenido en imagen (la extracción de texto pierde las fórmulas). Por eso:
1. Renderizá **todas** las diapositivas de tus PDFs, en hojas de contacto de 4 por imagen para ahorrar lecturas:
   `pdftoppm -r 60 -png <pdf> /home/claude/img/uN/<nombre>` y después `montage` (ImageMagick) o PIL para armar grillas 2×2; si una fórmula no se lee, rendereá esa diapositiva sola a `-r 120`.
2. Leé cada hoja con la herramienta Read. El índice oficial = el título de cada diapositiva (sin portada ni "FIN") + los objetivos si hay diapositiva de objetivos.
3. Guardá el inventario: número de diapositiva → título → qué contiene (fórmulas transcriptas), en `/home/claude/calibracion/01-inventario-uN.md`.

Texto ya extraído: `/home/claude/fuentes/txt/<archivo>.txt`

**Diapositivas con el contenido sólo en imagen:** la lista de tu unidad está en el inventario, marcada [img]. Renderizá **cada una** y leela con la herramienta Read:
```
pdftoppm -r 90 -png -f P -l P <pdf> /home/claude/img/mN-pP
```
Si no se lee, subí a `-r 140`. Transcribí su contenido al resumen. Si sigue sin leerse, anotalo en `fidelidad.huecos` y no lo completes de memoria.

## Reglas de fidelidad (no negociables)

1. **No inventar fórmulas, coeficientes de métodos, órdenes, condiciones ni números.** Fórmulas, notación, nombres de métodos, órdenes de error y condiciones salen **textuales** de las diapositivas (con la página). **Todo número** que aparezca (resultados de ejercicios, tablas de iteraciones, errores, coeficientes de un ajuste, frecuencias de un espectro, respuestas de una pregunta de cálculo, orden de un método deducido por Taylor) sale de una corrida de Python guardada en `/home/claude/hub/tools/gen/verif_uN.py`: se corre, se mira la salida y se escribe **lo que dio la corrida**. Si la cátedra usa una convención (aritmética de 3 dígitos, varianza dividiendo por n, f < fs/2 o f ≤ fs/2), se usa la de la cátedra y se dice cuál; si no la explicita, se dice que no está explicitada.
   Todo eso sale **textual** del material. Si la cátedra tiene un error tipográfico, se cita su forma y al lado va la correcta, dentro de un bloque `agregado`.
2. **Lo que no está en el material va dentro de `<div class="agregado">`.** Sólo se agrega cuando hace falta para responder algo que la cátedra toma o para entender una diapositiva.
3. **Cada tema del índice oficial (objetivos y títulos de diapositiva) tiene desarrollo propio** en alguna sección. Prohibido derivar al PDF: el validador lo rechaza.
4. **Las definiciones van textuales** de la diapositiva, en una caja `def` que indica la página.
5. Si un objetivo del índice no tiene diapositivas, va en `auditoria` con estado `hueco` y en `fidelidad.huecos`. No se rellena.

## Cómo se escribe

- Cada sección **arranca por el problema que el tema resuelve** (1–2 frases); después va el contenido. Nada de relleno.
- Peso según el patrón: ⭐⭐⭐ = aparece en la mayoría de los exámenes; ⭐⭐ = aparece alguna vez; ⭐ = nunca. Se cubre todo el índice, pero lo ⭐⭐⭐ va con más profundidad, ejemplos y una caja `trampa`.
- **Rótulo de las cajas `trampa`:** la etiqueta es "Así lo toma la cátedra" **sólo** si la caja cita un examen real (código del patrón, ej. P25-T1 MC7). Si es una predicción, la etiqueta es "Así lo tomaría la cátedra".
- Tablas comparativas donde haya cosas que se confunden entre sí.
- Secciones numeradas en el orden del PDF.

### HTML permitido dentro de `html`, `intro`, `conexiones` y `modelo`

Sólo estas piezas:
```html
<p>…</p>  <ul><li>…</li></ul>  <ol>…</ol>  <b>…</b>  <i>…</i>  <code>…</code>
<h4>Subtítulo dentro de la sección</h4>
<div class="def"><b>Término</b><p>Definición textual.</p><span class="src">U2 p8</span></div>
<div class="tbl"><table><thead><tr><th>…</th></tr></thead><tbody><tr><td>…</td></tr></tbody></table></div>
<pre class="cmd">Texto literal de la cátedra (comando, fórmula, código)</pre>
<div class="trampa"><b>Así lo toma la cátedra</b><p>Ejemplo real del parcial y por qué engaña.</p></div>
<div class="agregado"><b>Agregado</b><p>Lo que no está en el material.</p></div>
<div class="ej"><b>Ejemplo</b><p>…</p></div>
```
Fórmulas: en texto con Unicode (y′, yₙ₊₁, Σ, √, π, ≤, ·, h², O(h⁴)) o dentro de `<pre class="cmd">` para ecuaciones largas. Nada de LaTeX (no se renderiza).
Además se permite, en `html` de una sección, una figura propia: `<figure class="fig"><svg>…</svg></figure>` generada con `tools/gen/figs.py`.
Nada de estilos inline, scripts, `<img>`, `<br>` ni otras clases.

## Formato del archivo

Generalo con un **script Python guardado en `/home/claude/hub/tools/gen/`**, nunca en una carpeta temporal, usando `json.dumps(obj, ensure_ascii=False, indent=1)`. Mirá `ejemplos/gen_m2-aso2.py`. El archivo debe ser exactamente `ASO.addUnit(<json>);`

```json
{
 "id": 2,
 "titulo": "Título oficial de la unidad",
 "corto": "Nombre corto para el menú",
 "fuente": "archivo.pdf (32 diapositivas)",
 "peso": 3,
 "intro": "<p>Qué problema resuelve la unidad y qué es lo que más toman.</p>",
 "secciones": [
  {"id": "m2-s1", "titulo": "…", "estrellas": 3, "pags": "p3–7", "html": "…"}
 ],
 "conexiones": "<ul><li>Con la unidad 3: …</li></ul>",
 "fidelidad": {"fuentes": "…", "agregado": "…", "huecos": "…"},
 "auditoria": [{"tema": "Objetivo: …", "seccion": "m2-s4", "estado": "cubierto"}],
 "glosario": [{"termino": "…", "def": "…", "seccion": "m2-s7"}],
 "quiz": [
  {"id": "m2-q01", "seccion": "m2-s9", "tema": "…", "pide": "verdaderas",
   "sobre": "…", "opciones": [{"t": "Afirmación", "v": true}, {"t": "…", "v": false}],
   "explicacion": "≥200 caracteres…"}
 ],
 "desarrollo": [
  {"id": "m2-d01", "seccion": "m2-s9", "real": "P25-T1 D4 o null",
   "pregunta": "Enunciado estilo cátedra.", "modelo": "<p>Respuesta modelo completa.</p>"}
 ]
}
```

- `auditoria`: **una fila por cada objetivo oficial y por cada título de diapositiva distinto** (sin portada ni "FIN").
- `glosario`: 12–25 términos, los que un alumno confundiría.
- `desarrollo`: todas las preguntas reales de la unidad + variantes rotadas, y además ejercicios de **práctica** (ver abajo). Si el patrón tiene preguntas reales de tu unidad, van **todas** (con `real`), más variantes rotadas. La respuesta modelo cumple todo lo que el patrón dice que exige la cátedra: cada inciso (a, b, c, d) respondido por separado y rotulado; definición + justificación explícita ("porque…") en cada inciso que pide explicar; la ecuación cuando se pide; para clasificar un método: los 4 incisos (RK sí/no, orden, uno/multipaso, implícito/explícito) cada uno con su razón.
- **Pregunta ambigua del parcial** (admite dos lecturas, o el enunciado dice "si alguna"): la respuesta modelo presenta las dos lecturas y dice cuál escribir primero. Nunca una respuesta categórica que la cátedra podría corregir distinto.

## El quiz: 30 preguntas con el formato real

Formato de la cátedra según el patrón: **opción única con 4 opciones (estilo Moodle y evaluación TP2/TP3, con o sin gráfico) y bloques de verdadero/falso (4 afirmaciones, se marcan las verdaderas, todo o nada)**.

Tipos de pregunta que acepta el hub (campo `pide`):
- `"verdaderas"` / `"falsas"`: el hub muestra el `enunciado` de `config.js` con `{pide}` y `{sobre}`. Las correctas son las que coinciden con `pide`. Hay entre 1 y n−1 correctas y la corrección es todo o nada.
- `"una"`: opción única. El hub muestra `sobre` como pregunta y hay **exactamente una** opción con `v: true`.

Reglas:
- `v` es **si la afirmación es verdadera**.
- 4 o 5 opciones. Si el formato es de varias correctas, variá cuántas son correctas.
- Proporción de tipos igual a la del patrón: U2 y U3: ~70 % `"una"` y ~30 % `"verdaderas"`; U4: ~60 % `"una"` y ~40 % `"verdaderas"`. No uses `"falsas"`.
- Distribución de las preguntas por tema según las ⭐.
- **Distractores según el patrón**, con los mecanismos que documenta `02-patron`. Ninguno descartable de un vistazo y ninguno absurdo.
- No repitas la misma afirmación con la misma polaridad.
- **Las opciones se mezclan en cada intento:** la explicación **nunca** nombra letras ("la a", "opción b"). Cada afirmación se nombra por su contenido. El validador lo rechaza.
- **Explicación ≥200 caracteres** que diga por qué cada opción está bien o mal.
- **Balance de longitud (lo mide el validador):** el porcentaje de preguntas donde la opción más larga es correcta tiene que estar a ±5 puntos de la proporción media de opciones correctas. Con opción única de 4 opciones, eso es 20–30 %. Además, el largo medio de correctas / incorrectas tiene que quedar entre 0,90 y 1,10. Si no pasa, alargá distractores o acortá correctas y volvé a medir.

## Notas de patrón para tu unidad

Te llegan en el mensaje del coordinador.
(Qué preguntas reales hubo, textuales; qué temas son ⭐⭐⭐; distribución sugerida del quiz; errores tipográficos conocidos de la cátedra; diapositivas-imagen críticas.)

## Cierre obligatorio

1. `python3 /home/claude/hub/tools/validar.py /home/claude/hub/data/mN.js` **sin errores**.
2. Revisá la auditoría vos mismo: cada objetivo y cada título tiene que tener una sección con desarrollo real.
3. Respuesta final (para el coordinador, no para el alumno), breve:
   - métricas del validador;
   - diapositivas-imagen transcriptas y las que quedaron ilegibles;
   - agregados hechos;
   - dudas de fidelidad;
   - ruta del generador.


## Campos extra que soporta este hub (además de los del esquema)

- **`fig`** (opcional) en una pregunta de `quiz` o de `desarrollo`: un `<svg>…</svg>` completo generado con `/home/claude/hub/tools/gen/figs.py` (`from figs import nueva, svg`). Se muestra debajo del enunciado. Usalo para espectros, señales, nubes de puntos, imágenes sintéticas filtradas (con `imshow(..., rasterized=True)` y como mucho 128×128 píxeles), gráficos de error. < 150 KB cada uno. Nunca reproduzcas imágenes de terceros (las fotos del TP3, pósters, etc.): generá equivalentes propios.
- **`practica: true`** (opcional) en un ítem de `desarrollo`: lo marca como ejercicio práctico (sale en el filtro "Práctica"). La `pregunta` da todos los datos; el `modelo` muestra el procedimiento paso a paso (tabla con `<div class="tbl">` si hay iteraciones) y el resultado, todo salido de `verif_uN.py`. Si replica una pregunta real (evaluación TP o Moodle), `real` lleva el código; si es del TP, `real: null` y la pregunta empieza con "TP2 ej. 3:" etc.; si es una predicción, la respuesta modelo cierra con una caja `trampa` rotulada "Así lo tomaría la cátedra".
- En `quiz`, el campo `sobre` y `pregunta` son **texto plano** (se escapan): sin HTML. Fórmulas con Unicode.

## Verificación numérica obligatoria

- Script `/home/claude/hub/tools/gen/verif_uN.py`: calcula **todo** número que use la unidad y lo imprime. El generador `gen_uN.py` puede importarlo o copiar sus resultados, pero ningún número se escribe "a mano".
- Para el orden/consistencia de un método de EDO dado como fórmula: deducilo con sympy (desarrollo de Taylor del error de truncamiento local) en `verif_u4.py` y que coincida con lo que dice la cátedra.
- Corré el script, leé la salida, y si algo no coincide con lo que esperabas, investigá antes de escribir.
