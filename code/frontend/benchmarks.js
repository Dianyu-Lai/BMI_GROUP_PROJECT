import { el, linePath, option, percent, svg, valueRange } from "./ui.js";

const METHODS = {
  E0: "Original RVNE",
  E1: "Optimized RVNE",
  E2: "Uniform smoothing",
  E3: "Gaussian smoothing",
  E4: "Adaptive windows",
  E5: "Segmentation only",
  E6: "Default SCAW-RVNE",
  E7: "Retuned RVNE",
  E8: "Retuned SCAW-RVNE",
};
const METRICS = [
  ["spearman", "Spearman"],
  ["variable_recall", "Variable-site recall"],
  ["boundary_f1", "Boundary F1"],
  ["auroc", "AUROC"],
  ["auprc", "AUPRC"],
  ["scored_fraction", "Scored-column fraction"],
  ["medium_as_variable", "Medium called variable"],
  ["sensitivity", "Conserved-site recall"],
  ["seconds", "Runtime"],
];
const PATTERNS = {
  A: "Uniform conservation",
  B: "Central variable block",
  C: "Short conserved motif",
  D: "Three mutation rates",
};
const LABEL = Object.fromEntries(METRICS);
const methodLabel = (method) => METHODS[method] ? `${method} ${METHODS[method]}` : method;
const stat = (value) => value == null ? "N/A" : value.toFixed(4);
const fullStat = (value) => value == null ? { text: "N/A" } : { text: stat(value), title: String(value) };
const order = ["E8", "E7", "Shannon", "E5", "E6", "E4", "E3", "E2", "E1", "E0"];

function table(headers, rows) {
  const head = el("thead", {}, [el("tr", {}, headers.map((text) => el("th", { scope: "col", text })))]);
  const body = el("tbody", {}, rows.map(([name, ...values]) => el("tr", {}, [
    el("th", { scope: "row", text: name }),
    ...values.map((value) => el("td", typeof value === "object" ? value : { text: value })),
  ])));
  return el("div", { className: "table-scroll" }, [el("table", {}, [head, body])]);
}

function labeledSelect(label, options) {
  const select = el("select", {}, options.map(([value, text]) => option(value, text)));
  select.value = options[0]?.[0] ?? "";
  return [el("label", { className: "field" }, [el("span", { text: label }), select]), select];
}

function formatSource(path) {
  return el("p", { className: "source-note" }, ["Source ", el("code", { text: path })]);
}

function panel(title, ...children) {
  return el("section", { className: "panel bench-panel" }, [
    el("div", { className: "panel-heading" }, [el("h2", { text: title })]),
    ...children,
  ]);
}

function orderedMethods(methods) {
  const rank = (method) => order.includes(method) ? order.indexOf(method) : order.length;
  return [...methods].sort((a, b) => rank(a) - rank(b) || a.localeCompare(b));
}

function missingReason(pattern, metric) {
  if (["A", "conserved", "variable"].includes(pattern) && ["auroc", "auprc"].includes(metric)) {
    return "Each family has one truth class. The benchmark reports AUROC and AUPRC as undefined.";
  }
  if (["A", "conserved", "variable"].includes(pattern) && metric === "spearman") return "Truth is constant within each family.";
  if (metric === "variable_recall" && pattern === "A") return "Pattern A has no true variable sites.";
  if (metric === "boundary_f1") return "Boundary F1 requires predicted and true boundaries for evaluation.";
  return "The source contains no defined value for this method and metric.";
}

function barChart(rows, metric) {
  const signed = metric === "spearman";
  const domainLow = signed ? -1 : 0;
  const x = (value) => (value - domainLow) * 100 / (1 - domainLow);
  const bars = rows.map((row) => {
    const width = row.mean == null ? 0 : Math.abs(x(row.mean) - x(0));
    const start = row.mean == null ? 0 : Math.min(x(0), x(row.mean));
    const bar = el("div", { className: "bar-track" }, [
      ...(signed ? [el("i", { className: "bar-zero" })] : []),
      ...(row.mean == null ? [] : [el("i", { className: "bar", style: `left:${start}%;width:${width}%` })]),
    ]);
    return el("div", { className: "bar-row" }, [
      el("span", { text: methodLabel(row.method) }), bar, el("span", { className: "bar-value", text: stat(row.mean) }),
    ]);
  });
  const tickValues = signed ? ["−1", "−0.5", "0", "0.5", "1"] : ["0", "0.25", "0.5", "0.75", "1"];
  const axis = el("div", { className: "bar-row bar-axis" }, [
    el("span", { text: LABEL[metric] }),
    el("div", { className: "bar-ticks" }, tickValues.map((text) => el("span", { text }))),
    el("span"),
  ]);
  return el("div", { className: "bar-chart", role: "group", "aria-label": `${LABEL[metric]} by method` }, [...bars, axis]);
}

function motifPanel(data) {
  if (!data.motifs.length) return panel("Pattern C · Short motif recovery",
    el("p", { className: "bench-copy", text: "Motif-width summary unavailable." }));
  const methods = ["E2", "E3", "E7", "E8", "Shannon"];
  const widths = [3, 5, 9];
  const rows = methods.map((method) => [methodLabel(method), ...widths.map((width) => {
    const row = data.motifs.find((item) => item.method === method && item.width === width);
    return row ? `${percent(row.mean)} · n=${row.n}` : "N/A";
  })]);
  return panel("Pattern C · Short motif recovery",
    el("p", { className: "bench-copy", text: "Conserved-site recall at the frozen thresholds, by synthetic motif width. Uniform and Gaussian smoothing missed the 3- and 5-column motifs in this set." }),
    table(["Method", "3 columns", "5 columns", "9 columns"], rows), formatSource(data.motif_source));
}

function threeRatesPanel(dataset) {
  const methods = ["E7", "E8", "Shannon", "E5"];
  const metrics = ["medium_as_variable", "variable_recall", "spearman"];
  const rows = methods.map((method) => [methodLabel(method), ...metrics.map((metric) => {
    const row = dataset.rows.find((item) => item.pattern === "D" && item.method === method && item.metric === metric);
    if (!row || row.mean == null) return "N/A";
    return `${metric === "spearman" ? stat(row.mean) : percent(row.mean)} · n=${row.n}`;
  })]);
  return panel("Pattern D · Three mutation rates",
    el("p", { className: "bench-copy", text: "Binary recall treats medium- and high-rate sites as variable. The medium-only column isolates detection at the intermediate rate." }),
    table(["Method", "Medium called variable", "Variable-site recall", "Spearman"], rows), formatSource(dataset.source));
}

function runtimePanel(data) {
  if (!data.runtime.length) return panel("Runtime", el("p", { className: "bench-copy", text: "Recorded timings unavailable." }));
  return panel("Reference runtime",
    table(["Stage", "Sequences", "Columns", "Median seconds", "Repeats"],
      data.runtime.map((row) => [row.stage, row.sequences, row.columns, row.seconds.toFixed(6), row.repeats])),
    formatSource(data.runtime_source));
}

function casePanel(data) {
  if (!data?.sites?.length) return panel("Case study · UniProt globins",
    el("p", { className: "bench-copy", text: "The MAFFT and ClustalW comparison appears here when P02144 reference-coordinate outputs are available." }));
  const tracks = [["mafft", "MAFFT", "final-line"], ["clustalw", "ClustalW", "raw-line"]];
  const [low, high] = valueRange(data.sites.flatMap((row) => [row.mafft, row.clustalw]));
  const x = (index) => 24 + 740 * index / Math.max(1, data.sites.length - 1);
  const y = (value) => 141 - 118 * (value - low) / (high - low);
  const ticks = Array.from({ length: 4 }, (_, i) => low + i * (high - low) / 3);
  const chart = svg("svg", { class: "case-chart", viewBox: "0 0 800 166", role: "img",
    "aria-label": "MAFFT and ClustalW profiles by ungapped P02144 reference residue" }, [
    ...ticks.flatMap((value) => {
      const label = svg("text", { x: 29, y: y(value) - 5 });
      label.textContent = value.toFixed(2);
      return [svg("line", { x1: 24, x2: 764, y1: y(value), y2: y(value), class: "plot-grid" }), label];
    }),
    ...tracks.map(([key, , style]) => svg("path", { d: linePath(data.sites.map((row) => row[key]), x, y), class: style })),
  ]);
  return panel("Case study · UniProt globins",
    el("p", { className: "bench-copy", text: `Alignment sensitivity in P02144 coordinates · jointly scored ${data.joint}/${data.sites.length} residues (${percent(data.joint / data.sites.length)}).` }),
    chart,
    el("div", { className: "case-legend" }, [
      ...tracks.map(([, label, style]) => el("span", {}, [el("i", { className: `line-key ${style}` }), label])),
      el("span", { text: `P02144 reference residues 1–${data.sites.length}` }),
    ]),
    formatSource(data.source));
}

export function renderBenchmarks(target, data) {
  const [datasetField, dataset] = labeledSelect("Dataset", [
    ["supplemental", "Supplemental A–D"], ["reference", "Reference held-out"],
  ]);
  const [patternField, pattern] = labeledSelect("Pattern", []);
  const [metricField, metric] = labeledSelect("Metric", METRICS);
  const [methodField, method] = labeledSelect("Method", []);
  const comparison = el("div", { className: "bench-comparison" });
  const evidence = el("div", { className: "benchmark-evidence" });
  target.replaceChildren(
    el("div", { className: "page-heading" }, [
      el("h1", { text: "Benchmarks" }), el("p", { text: "Conservation scoring on synthetic protein alignments." }),
    ]),
    el("div", { className: "panel benchmark-filters" }, [datasetField, patternField, methodField, metricField]),
    comparison, evidence, casePanel(data.case),
  );

  function choosePatterns() {
    const groups = data.datasets[dataset.value].groups ?? {};
    const keys = dataset.value === "supplemental" ? Object.keys(PATTERNS) : Object.keys(groups).sort();
    pattern.replaceChildren(...keys.map((key) => option(key,
      PATTERNS[key] ? `Pattern ${key} · ${PATTERNS[key]}` : key[0].toUpperCase() + key.slice(1))));
    if (dataset.value === "supplemental") pattern.value = "D";
    chooseMethods();
  }

  function chooseMethods() {
    const current = method.value;
    const available = data.datasets[dataset.value].groups?.[pattern.value] ?? [];
    method.replaceChildren(option("", "All methods"), ...orderedMethods(available).map((key) => option(key, methodLabel(key))));
    if (available.includes(current)) method.value = current;
    refresh();
  }

  function refresh() {
    const current = data.datasets[dataset.value];
    const kind = metric.value;
    pattern.disabled = kind === "seconds";
    method.disabled = kind === "seconds";
    dataset.disabled = kind === "seconds";
    const isSupplemental = dataset.value === "supplemental";
    evidence.replaceChildren();
    if (kind === "seconds") {
      comparison.replaceChildren(runtimePanel(data));
      return;
    }
    if (!current.available) {
      comparison.replaceChildren(panel("Benchmark results", el("p", {
        className: "bench-copy", text: "The benchmark summary TSV is unavailable.",
      }), formatSource(current.source)));
      return;
    }
    const methods = orderedMethods(current.groups[pattern.value] ?? []).filter((key) => !method.value || method.value === key);
    const rows = methods.map((key) => current.rows.find((row) =>
      row.pattern === pattern.value && row.method === key && row.metric === kind
    ) ?? { method: key, mean: null, sd: null, n: 0 });
    const title = isSupplemental ? `Pattern ${pattern.value} · ${LABEL[kind]}` : `${LABEL[kind]} · ${pattern.value}`;
    const display = rows.length ? [
      barChart(rows, kind),
      ...(rows.some((row) => row.mean == null) ? [
        el("p", { className: "bench-copy", text: `N/A · ${missingReason(pattern.value, kind)}` }),
      ] : []),
      table(["Method", "Mean", "SD", "Families (n)"], rows.map((row) => [
        methodLabel(row.method), fullStat(row.mean), fullStat(row.sd), row.n,
      ])),
    ] : [el("p", { className: "bench-copy", text: "No methods are recorded for this pattern." })];
    comparison.replaceChildren(panel(title, ...display,
      el("p", { className: "bench-copy", text: !isSupplemental && kind === "variable_recall"
        ? "Variable-site recall is recorded as specificity. n counts families with a defined value, weighted equally."
        : "Equal weight per family. n counts families with a defined metric." }),
      formatSource(current.source)));
    if (isSupplemental && pattern.value === "C") evidence.append(motifPanel(data));
    if (isSupplemental && pattern.value === "D") evidence.append(threeRatesPanel(current));
  }

  dataset.addEventListener("change", choosePatterns);
  pattern.addEventListener("change", chooseMethods);
  method.addEventListener("change", refresh);
  metric.addEventListener("change", refresh);
  choosePatterns();
}
