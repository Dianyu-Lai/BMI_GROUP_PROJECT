import { byId, el, getJSON, option } from "./ui.js";

const notes = {
  conserved: "Uniform low mutation rate across the ancestor.",
  variable: "Uniform high mutation rate across the ancestor.",
  blocks: "Low, medium and high mutation-rate regions.",
  motifs: "Three conserved motifs in a variable background.",
  gradual: "A gradual rate pattern approximated by twelve constant-rate regions.",
  divergent: "A block pattern with divergence multiplied by 2.5 by the generator.",
  indels: "A block pattern with evolutionary indels. With both rates zero and no target supplied, the generator calibrates to 20% gap occupancy; calibration can take longer.",
  redundant: "A block pattern with a near-tip clade of related sequences. Blank redundancy uses the generator default of 0.7.",
  "conserved-variable-blocks": "Conserved flanks around a central variable block.",
  "single-motif": "A single central conserved motif in a variable background.",
  "three-rate-blocks": "Low, medium, high, medium and low mutation-rate regions.",
};

export function initGeneration(useMSA) {
  const form = byId("generation-form");
  const errorBox = byId("generation-error");
  const status = byId("generation-status");
  const resultPanel = byId("generation-result");
  let config;
  let busy = false;
  let generated = null;
  let scenario;
  const controls = new Map();

  function error(message = "") {
    errorBox.textContent = message;
    errorBox.hidden = !message;
  }

  function update() {
    for (const [key, { input, wrapper, spec }] of controls) {
      const applicable = !spec.scenarios || spec.scenarios.includes(scenario.value);
      wrapper.hidden = !applicable;
      input.disabled = busy || !applicable;
      if (key === "motif_length") input.max = Math.floor(controls.get("length").input.valueAsNumber / 4);
    }
    scenario.disabled = busy;
    byId("run-generation").disabled = busy;
    byId("use-generated").disabled = busy || !generated;
    form.setAttribute("aria-busy", String(busy));
    byId("scenario-note").textContent = notes[scenario.value] || "";
  }

  function validate() {
    for (const { input, spec } of controls.values()) {
      input.setCustomValidity("");
      if (input.disabled || input.tagName === "SELECT" || (spec.optional && input.value === "")) continue;
      const number = input.valueAsNumber;
      if (!Number.isFinite(number) || (spec.integer && !Number.isSafeInteger(number))) {
        input.setCustomValidity(`Enter ${spec.integer ? "an integer" : "a finite number"} for ${spec.label.toLowerCase()}.`);
      } else if (spec.exclusive && (number <= spec.min || number >= spec.max)) {
        input.setCustomValidity(`Enter a value strictly between ${spec.min} and ${spec.max}.`);
      }
    }
  }

  async function load() {
    error();
    byId("generation-retry").hidden = true;
    byId("generation-loading").hidden = false;
    try {
      config = await getJSON("/api/generation/config");
      scenario = el("select", { id: "generation-scenario", name: "scenario" },
        config.scenarios.map((name) => option(name, name)));
      scenario.value = config.scenario;
      byId("generation-basic").replaceChildren(el("div", { className: "field" }, [
        el("label", { for: scenario.id, text: "Simulation scenario" }), scenario,
      ]));
      byId("generation-advanced-fields").replaceChildren();
      controls.clear();
      for (const [key, spec] of Object.entries(config.fields)) {
        const id = `generation-${key}`;
        const input = spec.choices
          ? el("select", { id, name: key }, spec.choices.map((name) => option(name, name)))
          : el("input", { id, name: key, type: "number", min: spec.min, max: spec.max,
            step: spec.integer ? "1" : "any" });
        input.value = spec.default ?? "";
        input.required = !spec.optional;
        if (spec.optional) input.placeholder = "Generator default";
        const wrapper = el("div", { className: "field" }, [
          el("label", { for: id, text: `${spec.label}${spec.optional ? " (optional)" : ""}` }), input,
        ]);
        if (spec.note) {
          const noteId = `${id}-note`;
          input.setAttribute("aria-describedby", noteId);
          wrapper.append(el("p", { id: noteId, className: "field-note", text: spec.note }));
        }
        const basic = ["size", "length", "seed", "divergence"].includes(key);
        byId(basic ? "generation-basic" : "generation-advanced-fields").append(wrapper);
        controls.set(key, { input, wrapper, spec });
      }
      form.hidden = false;
      update();
      if (config.unavailable) error(config.unavailable);
    } catch (problem) {
      error(problem.message);
      byId("generation-retry").hidden = false;
    } finally {
      byId("generation-loading").hidden = true;
    }
  }

  byId("generation-retry").addEventListener("click", load);
  function changed() {
    update();
    validate();
    if (generated) status.textContent = "Parameters changed. Generate again to update the result below.";
  }
  form.addEventListener("input", changed);
  form.addEventListener("change", changed);
  form.addEventListener("invalid", (event) => {
    const advanced = byId("generation-advanced");
    if (advanced.contains(event.target)) advanced.open = true;
  }, true);
  form.addEventListener("submit", async (event) => {
    event.preventDefault();
    if (busy) return;
    validate();
    if (!form.reportValidity()) return;
    const parameters = { scenario: scenario.value };
    for (const [key, { input, spec }] of controls) {
      if (input.disabled || (spec.optional && input.value === "")) continue;
      parameters[key] = spec.choices ? input.value : input.valueAsNumber;
    }
    error();
    generated = null;
    resultPanel.hidden = true;
    busy = true;
    status.textContent = "Submitting generation…";
    update();
    try {
      const { id } = await getJSON("/api/generate", {
        method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(parameters),
      });
      let job;
      do {
        job = await getJSON(`/api/jobs/${id}`);
        if (job.state === "running") {
          status.textContent = `${job.stage}…`;
          await new Promise((resolve) => setTimeout(resolve, 250));
        }
      } while (job.state === "running");
      if (job.state !== "done") throw new Error(job.error || "Generation did not complete.");
      generated = job.result;
      const summary = [
        ["Sequences", generated.sequenceCount], ["Ancestor length", generated.ancestorLength],
        ["Aligned columns", generated.alignedColumns], ["Scenario", generated.scenario],
        ["Seed", generated.seed], ["Effective divergence", generated.divergence],
      ];
      byId("generation-summary").replaceChildren(...summary.map(([label, value]) =>
        el("div", {}, [el("dt", { text: label }), el("dd", { text: value })])));
      byId("generation-preview").textContent = generated.fasta.trim().split(/\r?\n/).slice(0, 6)
        .map((line) => line.startsWith(">") || line.length <= 120 ? line : `${line.slice(0, 120)}…`).join("\n");
      for (const [element, filename] of [["generation-fasta", "aligned.fasta"], ["generation-truth", "site_truth.tsv"]]) {
        byId(element).href = `/api/jobs/${id}/download/${filename}`;
        byId(element).setAttribute("download", filename);
      }
      resultPanel.hidden = false;
      status.textContent = "Generation complete.";
      if (!byId("generate-page").hidden) byId("generation-result-heading").focus();
    } catch (problem) {
      error(problem.message);
      status.textContent = "Generation failed.";
    } finally {
      busy = false;
      update();
    }
  });
  byId("use-generated").addEventListener("click", () => {
    if (!generated || busy) return;
    try {
      useMSA(generated);
    } catch (problem) {
      error(problem.message);
    }
  });
  load();
}
