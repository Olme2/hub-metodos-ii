// Prueba funcional del hub sin navegador.
// Uso:  npm install jsdom@24.1.3   (una vez, en cualquier carpeta fuera del repo)
//       node tools/smoke.js <carpeta-del-hub>
// Carga index.html, recorre cada resumen, hace un quiz completo de la primera unidad,
// abre glosario, progreso y desarrollo, y lista cualquier error de JavaScript.
const { JSDOM, VirtualConsole } = require("jsdom");
const path = require("path");
const site = path.resolve(process.argv[2] || ".");
const vc = new VirtualConsole(), errs = [];
vc.on("jsdomError", e => { if (!/Could not load link/.test(e.message)) errs.push(e.message); });
JSDOM.fromFile(path.join(site, "index.html"), { runScripts: "dangerously", resources: "usable", pretendToBeVisual: true, virtualConsole: vc }).then(dom => {
  const w = dom.window, d = w.document;
  w.Element.prototype.scrollIntoView = function () {}; w.scrollTo = () => {};
  const wait = ms => new Promise(r => setTimeout(r, ms));
  setTimeout(async () => {
    const out = [];
    try {
      const ids = w.HUB.cfg.unidades.map(u => u.id);
      out.push(`título: ${d.title} · unidades en config: ${ids.join(",")} · cargadas: ${Object.keys(w.HUB.units).join(",")}`);
      for (const id of ids) {
        w.location.hash = "m" + id; await wait(30);
        out.push(`  ${id}: ${d.querySelectorAll(".sec").length} secciones · ${d.querySelectorAll(".def").length} def · ${d.querySelectorAll(".agregado").length} agregados · auditoría ${d.querySelectorAll("details.fid tbody tr").length} filas`);
      }
      const first = ids.find(i => w.HUB.units[i]);
      w.location.hash = "m" + first; await wait(30); d.querySelector("#toQuiz").click(); await wait(30);
      const n = w.HUB.units[first].quiz.length;
      for (let k = 0; k < n; k++) {
        const inp = d.querySelector("input[data-k]"); inp.checked = true; inp.dispatchEvent(new w.Event("change"));
        d.querySelector("#chk").click(); await wait(5); d.querySelector("#next").click(); await wait(5);
      }
      out.push(`quiz unidad ${first}: resultado ${d.querySelector(".score b").textContent} · links a releer ${d.querySelectorAll(".flist a").length}`);
      w.location.hash = "glosario"; await wait(30); out.push(`glosario: ${d.querySelectorAll(".gl .t").length} términos`);
      w.location.hash = "progreso"; await wait(30); out.push(`progreso: ${d.querySelectorAll(".ptable tbody tr").length} filas`);
      w.location.hash = "quiz"; await wait(30); d.querySelector("#td").click(); await wait(30);
      out.push(`desarrollo: ${d.querySelectorAll("[data-dev]").length} preguntas`);
    } catch (e) { errs.push("PRUEBA: " + e.message); }
    console.log(out.join("\n"));
    console.log("ERRORES:", errs.length ? errs : "ninguno");
    process.exit(errs.length ? 1 : 0);
  }, 1500);
});
