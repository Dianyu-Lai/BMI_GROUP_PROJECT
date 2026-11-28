import { byId, el, getJSON, option, toolName } from "./ui.js";
import { renderResults } from "./results.js";
import { renderBenchmarks } from "./benchmarks.js";
import { initGeneration } from "./generate.js";

const form = byId("analysis-form");
const run = byId("run-analysis");
const upload = byId("fasta-file");
const status = byId("run-status");
const errorBox = byId("form-error");
const fileMeta = byId("file-meta");
const dropZone = byId("drop-zone");
const loadButton = byId("load-example");
const examplePicker = byId("example");
const aligner = byId("aligner");
const method = byId("method");
const preset = byId("preset");
const customPreset = preset.querySelector('[value="custom"]');
const customRvne = byId("custom-rvne");
const customContext = byId("custom-context");
const results = byId("results");
let input = null;
let config = null;
let busy = false;
let inputVersion = 0;
let releaseResults = null;
let benchmarksPromise = null;
let previousPreset = preset.value;
let customInitialized = false;

function error(message = "") {
  errorBox.hidden = !message;
  errorBox.textContent = message;
}

function updateForm() {
  const unaligned = form.elements.mode.value === "unaligned";
  const isExample = Boolean(input?.example);
  const isCustom = preset.value === "custom";
  const isScaw = method.value === "scaw-rvne";
  for (const radio of form.elements.mode) radio.disabled = busy || isExample;
  byId("aligner-field").hidden = !unaligned;
  const selectedTool = toolName(aligner.value);
  byId("aligner-status").textContent = config && !config.aligners[aligner.value]
    ? `${selectedTool} is unavailable on this machine.` : "";
  run.textContent = unaligned ? "Align & score" : "Score alignment";
  run.disabled = !input || busy;
  form.setAttribute("aria-busy", String(busy));
  for (const control of ["upload", "remove-file", "aligner", "method", "preset", "max-gap"]) {
    byId(control).disabled = busy;
  }
  examplePicker.disabled = loadButton.disabled = busy || !config;
  customRvne.hidden = !isCustom;
  customContext.hidden = !isCustom || !isScaw;
  customRvne.disabled = busy || customRvne.hidden;
  customContext.disabled = busy || customContext.hidden;
  const parameterValues = byId("parameter-values");
  parameterValues.hidden = isCustom || !config;
  if (parameterValues.hidden) return;

  const effective = config?.[preset.value === "frozen" ? "frozen" : "defaults"];
  const values = effective?.rvne ? [["RVNE · r / s1 / s2", `${effective.rvne.r} / ${effective.rvne.s1} / ${effective.rvne.s2}`]] : [];
  if (isScaw && effective?.context) {
    const context = effective.context;
    values.push(["Segmentation penalty multiplier", context.penalty_scale],
      ["Context window", "3 columns each side, clipped at segment boundaries"]);
  }
  parameterValues.replaceChildren(...values.map(([label, value]) =>
    el("div", {}, [el("dt", { text: label }), el("dd", { text: value })])));
}

function fillParameters(fields, defaults) {
  for (const field of fields.elements) {
    const value = defaults[field.name];
    field.value = Array.isArray(value) ? value.join(", ") : (value ?? "");
  }
}

function readParameters(fields) {
  return Object.fromEntries(Array.from(fields.elements, (field) => [field.name,
    field.value === "" ? null : field.valueAsNumber,
  ]));
}

function showInput() {
  fileMeta.hidden = !input;
  byId("drop-prompt").textContent = input ? "Replace FASTA file" : "Drag a protein FASTA here";
  byId("file-name").textContent = input?.filename ?? "";
  const count = input?.fasta.match(/^\s*>/gm)?.length ?? 0;
  byId("file-summary").textContent = input
    ? `${count} ${count === 1 ? "sequence" : "sequences"}${input.source ? ` · ${input.source}` : ""}`
    : "";
  updateForm();
}

async function loadFile(file) {
  if (!file || busy) return;
  const version = ++inputVersion;
  error();
  try {
    if (file.size > 1024 * 1024) throw new Error("Choose a FASTA file of up to 1 MiB.");
    const fasta = await file.text();
    if (version !== inputVersion) return;
    input = { filename: file.name, fasta };
    showInput();
  } catch (problem) {
    if (version === inputVersion) error(problem.message);
  }
}

byId("upload").addEventListener("click", () => upload.click());
upload.addEventListener("change", () => {
  loadFile(upload.files[0]);
  upload.value = "";
});
dropZone.addEventListener("dragover", (event) => {
  if (!event.dataTransfer.types.includes("Files")) return;
  event.preventDefault();
  dropZone.classList.toggle("dragging", !busy);
  event.dataTransfer.dropEffect = busy ? "none" : "copy";
});
dropZone.addEventListener("dragleave", () => dropZone.classList.remove("dragging"));
dropZone.addEventListener("drop", (event) => {
  event.preventDefault();
  dropZone.classList.remove("dragging");
  loadFile(event.dataTransfer.files[0]);
});
byId("remove-file").addEventListener("click", () => {
  inputVersion++;
  input = null;
  error();
  showInput();
  byId("upload").focus();
});
loadButton.addEventListener("click", async () => {
  const key = examplePicker.value;
  const version = ++inputVersion;
  error();
  loadButton.disabled = true;
  try {
    const loaded = await getJSON(`/api/examples/${key}`);
    if (version !== inputVersion || busy) return;
    input = { ...loaded, example: key };
    form.elements.mode.value = loaded.mode;
    showInput();
  } catch (problem) {
    if (version === inputVersion) error(problem.message);
  } finally {
    loadButton.disabled = busy;
  }
});
preset.addEventListener("change", () => {
  if (preset.value === "custom" && !customInitialized) {
    const initial = previousPreset === "frozen" ? config.frozen : config.defaults;
    fillParameters(customRvne, initial.rvne);
    fillParameters(customContext, initial.context);
    customInitialized = true;
  }
  previousPreset = preset.value;
});
byId("rvne-s2").addEventListener("input", (event) => {
  const field = event.currentTarget;
  field.setCustomValidity(field.valueAsNumber >= 1 ? "Enter an s2 value less than 1." : "");
});
form.addEventListener("invalid", (event) => {
  const advanced = byId("advanced");
  if (advanced.contains(event.target)) advanced.open = true;
}, true);
form.addEventListener("change", updateForm);

form.addEventListener("submit", async (event) => {
  event.preventDefault();
  if (!input || busy) return;
  const options = {
    ...(input.example ? { example: input.example } : { filename: input.filename, fasta: input.fasta }),
    mode: form.elements.mode.value,
    method: method.value,
    aligner: aligner.value,
    preset: preset.value,
    maxGap: byId("max-gap").valueAsNumber,
  };
  if (options.preset === "custom") {
    options.parameters = { rvne: readParameters(customRvne) };
    if (!customContext.hidden) options.parameters.context = readParameters(customContext);
  }
  error();
  busy = true;
  status.textContent = "Submitting…";
  updateForm();
  try {
    const { id } = await getJSON("/api/jobs", {
      method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(options),
    });
    let job;
    do {
      job = await getJSON(`/api/jobs/${id}`);
      if (job.state === "running") {
        status.textContent = `${job.stage}…`;
        await new Promise((resolve) => setTimeout(resolve, 250));
      }
    } while (job.state === "running");
    if (job.state === "error") throw new Error(job.error);
    results.hidden = false;
    releaseResults?.();
    releaseResults = renderResults(results, job.result, id, resetAnalysis);
    byId("setup").hidden = true;
    if (!byId("analyze-page").hidden) results.querySelector("h1").focus();
  } catch (problem) {
    results.hidden = true;
    const unsupportedPreset = options.preset === "custom"
      && problem.message === "Choose valid scoring settings and alignment tool.";
    error(unsupportedPreset
      ? "The running server does not support Custom settings. Restart it with python code/frontend/server.py, then reload this page."
      : problem.message);
  } finally {
    busy = false;
    status.textContent = "";
    updateForm();
  }
});

function resetAnalysis() {
  releaseResults?.();
  releaseResults = null;
  results.replaceChildren();
  results.hidden = true;
  byId("setup").hidden = false;
  input = null;
  error();
  showInput();
  byId("upload").focus();
  window.scrollTo(0, 0);
}

function navigate() {
  const page = ["#benchmarks", "#generate"].includes(location.hash) ? location.hash.slice(1) : "analyze";
  const isBenchmarks = page === "benchmarks";
  for (const name of ["analyze", "generate", "benchmarks"]) {
    byId(`${name}-page`).hidden = name !== page;
    const link = byId(`nav-${name}`);
    if (name === page) link.setAttribute("aria-current", "page");
    else link.removeAttribute("aria-current");
  }
  const bench = byId("benchmarks-page");
  if (isBenchmarks && !benchmarksPromise) {
    bench.replaceChildren(el("p", { role: "status", text: "Loading benchmark results…" }));
    benchmarksPromise = getJSON("/api/benchmarks").then((data) => renderBenchmarks(bench, data)).catch((problem) => {
      benchmarksPromise = null;
      const retry = el("button", { type: "button", className: "button", text: "Retry" });
      retry.addEventListener("click", navigate);
      bench.replaceChildren(el("p", { className: "notice error", role: "alert", text: problem.message }), retry);
    });
  }
}

getJSON("/api/config").then((data) => {
  config = data;
  customPreset.disabled = false;
  if (!config.frozen) {
    preset.value = "default";
    preset.options[0].disabled = true;
  }
  previousPreset = preset.value;
  if (!config.aligners.mafft && config.aligners.clustalw) aligner.value = "clustalw";
  if (config.examples) {
    examplePicker.replaceChildren(...Object.entries(config.examples).map(([key, label]) => option(key, label)));
  }
  updateForm();
}).catch((problem) => error(problem.message));
window.addEventListener("hashchange", navigate);
initGeneration((generated) => {
  if (busy) throw new Error("An analysis is running. Wait for it to finish before loading this MSA.");
  inputVersion++;
  releaseResults?.();
  releaseResults = null;
  results.replaceChildren();
  results.hidden = true;
  byId("setup").hidden = false;
  input = { filename: generated.filename, fasta: generated.fasta,
    source: `INDELible · ${generated.scenario} · seed ${generated.seed}` };
  form.elements.mode.value = "aligned";
  error();
  showInput();
  status.textContent = "Generated MSA loaded. Click Score alignment when ready.";
  location.hash = "#analyze";
  navigate();
  run.focus();
  window.scrollTo(0, 0);
});
updateForm();
navigate();
