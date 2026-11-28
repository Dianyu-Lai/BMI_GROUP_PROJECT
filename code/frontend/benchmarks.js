import { el, linePath, option, svg, valueRange } from "./ui.js";

const METHODS = {
  v2: "SCAW-RVNE",
  v2_segment_mean: "Segmentation only",
  fixed_uniform_7: "Uniform window (7 columns)",
  rvne: "Original RVNE",
  shannon: "Shannon",
  kl: "KL divergence",
  jsd: "Jensen–Shannon divergence",
  blosum62_sop: "BLOSUM62 sum of pairs",
  blosum62_cost: "BLOSUM62 substitution cost",
  adjusted_entropy: "Adjusted entropy",
  clustalw_consensus: "ClustalW consensus",
  al2co_variance: "AL2CO variance",
};
const METRICS = [
  ["auprc", "AUPRC"],
  ["spearman", "Spearman"],
  ["boundary_f1", "Boundary F1 (±2 columns)"],
  ["boundary_f1_exact", "Exact boundary F1"],
  ["auroc", "AUROC"],
  ["scored_fraction", "Scored-column fraction"],
];
const PATTERNS = {
  all: "All scenarios",
  conserved: "Uniform conservation",
  "conserved-variable-blocks": "Central variable block",
  "single-motif": "Short conserved motif",
  "three-rate-blocks": "Three mutation rates",
};
const LABEL = Object.fromEntries(METRICS);
const methodLabel = (method) => METHODS[method] ?? method;
const stat = (value) => value == null ? "N/A" : value.toFixed(4);
const fullStat = (value) => value == null ? { text: "N/A" } : { text: stat(value), title: String(value) };
const order = Object.keys(METHODS);

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

function panel(title, ...children) {
  return el("section", { className: "panel bench-panel" }, [
    el("div", { className: "panel-heading" }, [el("h2", { text: title })]),
    ...children,
  ]);
}

function datasetPanel(dataset) {
  const total = dataset.rows.find((row) =>
    row.pattern === "all" && row.method === "v2" && row.metric === "scored_fraction");
  const summary = panel("Dataset",
    table(["Cohort", "Families", "Sequences per family", "Aligned columns", "Composition"], [
      ["INDELible", total?.n ?? "N/A", 20, "120–145",
        "120 ancestral sites. Two uniformly conserved, two with a central variable block, four with short motifs (3 or 5 sites), and four with three rate tiers (two with indels)."],
    ]),
    el("p", { className: "bench-copy", text: "AUPRC is undefined for the two uniformly conserved families." }));
  summary.classList.add("bench-datasets");
  return summary;
}

function orderedMethods(methods) {
  const rank = (method) => order.includes(method) ? order.indexOf(method) : order.length;
  return [...methods].sort((a, b) => rank(a) - rank(b) || a.localeCompare(b));
}

function barChart(rows, metric, { titled = false, showCounts = false } = {}) {
  const signed = metric === "spearman";
  const domainLow = signed ? -1 : 0;
  const x = (value) => (value - domainLow) * 100 / (1 - domainLow);
  const bars = rows.map((row) => {
    const width = row.mean == null ? 0 : Math.abs(x(row.mean) - x(0));
    const start = row.mean == null ? 0 : Math.min(x(0), x(row.mean));
    const bar = el("div", { className: "bar-track" }, [
      ...(signed ? [el("i", { className: "bar-zero" })] : []),
      ...(row.mean == null ? [] : [el("i", {
        className: row.method === "v2" ? "bar" : "bar bar-comparator",
        style: `left:${start}%;width:${width}%`,
      })]),
    ]);
    const label = showCounts ? `${methodLabel(row.method)} (n=${row.n})` : methodLabel(row.method);
    return el("div", { className: "bar-row" }, [
      el("span", { text: label }), bar, el("span", { className: "bar-value", text: stat(row.mean) }),
    ]);
  });
  const tickValues = signed ? ["−1", "−0.5", "0", "0.5", "1"] : ["0", "0.25", "0.5", "0.75", "1"];
  const axis = el("div", { className: "bar-row bar-axis" }, [
    el("span", { text: titled ? "" : LABEL[metric] }),
    el("div", { className: "bar-ticks" }, tickValues.map((text) => el("span", { text }))),
    el("span"),
  ]);
  return el("div", { className: "bar-chart", role: "group", "aria-label": `${LABEL[metric]} by method` }, [...bars, axis]);
}

function methodComparison(data) {
  const baselineMethods = orderedMethods(new Set(data.baselines.map((row) => row.method)));
  if (!baselineMethods.length) return null;
  const current = data.datasets.fresh_indelible?.rows ?? [];
  const plots = ["auprc", "spearman"].map((metric) => {
    const scaw = current.find((row) => row.pattern === "all" && row.metric === metric && row.method === "v2");
    const rows = baselineMethods.map((method) =>
      data.baselines.find((row) => row.metric === metric && row.method === method) ?? { method, mean: null, n: 0 });
    return [metric, [scaw, ...rows]];
  });
  if (plots.some(([, rows]) => !rows[0])) return null;
  return panel("Method comparison",
    el("p", { className: "bench-copy", text: "Family means. n counts families with a defined score." }),
    el("div", { className: "bench-chart-grid" }, plots.map(([metric, rows]) => {
      const commonCount = rows.every((row) => row.n === rows[0].n);
      const title = commonCount ? `${LABEL[metric]} (n=${rows[0].n})` : LABEL[metric];
      return el("div", {}, [
        el("h3", { text: title }), barChart(rows, metric, { titled: true, showCounts: !commonCount }),
      ]);
    })));
}

function runtimePanel(data) {
  if (!data.runtime.length) return panel("Runtime", el("p", { className: "bench-copy", text: "Recorded timings unavailable." }));
  return panel("Context scoring runtime",
    table(["Method", "Columns", "Median seconds", "Repeats"],
      data.runtime.map((row) => [methodLabel(row.method), row.columns, row.seconds.toFixed(6), row.repeats])));
}

function casePanel(data) {
  if (!data?.sites?.length) return panel("Case study · UniProt globins",
    el("p", { className: "bench-copy", text: data?.error ?? "Run MAFFT and ClustalW to generate the P02144 profiles." }));
  const tracks = [["mafft", "MAFFT", "final-line"], ["clustalw", "ClustalW", "raw-line"]];
  const [low, high] = valueRange(data.sites.flatMap((row) => [row.mafft, row.clustalw]));
  const plotLeft = 52;
  const plotRight = 764;
  const x = (index) => plotLeft + (plotRight - plotLeft) * index / Math.max(1, data.sites.length - 1);
  const y = (value) => 141 - 118 * (value - low) / (high - low);
  const ticks = Array.from({ length: 4 }, (_, i) => low + i * (high - low) / 3);
  const chart = svg("svg", { class: "case-chart", viewBox: "0 0 800 166", role: "img",
    "aria-label": "MAFFT and ClustalW profiles by ungapped P02144 reference residue" }, [
    ...ticks.flatMap((value) => {
      const tickY = y(value);
      const label = svg("text", {
        x: plotLeft - 10, y: tickY, "text-anchor": "end", "dominant-baseline": "middle",
      });
      label.textContent = value.toFixed(2);
      return [svg("line", { x1: plotLeft, x2: plotRight, y1: tickY, y2: tickY, class: "plot-grid" }), label];
    }),
    ...tracks.map(([key, , style]) => svg("path", { d: linePath(data.sites.map((row) => row[key]), x, y), class: style })),
  ]);
  return panel("Case study · UniProt globins",
    el("p", { className: "bench-copy", text: `Jointly scored residues: ${data.joint}/${data.sites.length}${data.pearson == null ? "" : ` · Pearson r = ${data.pearson.toFixed(3)}`}` }),
    chart,
    el("div", { className: "case-legend" }, [
      ...tracks.map(([, label, style]) => el("span", {}, [el("i", { className: `line-key ${style}` }), label])),
      el("span", { text: `P02144 reference residues 1–${data.sites.length}` }),
    ]));
}

export function renderBenchmarks(target, data) {
  if (typeof data.source !== "string" || !Array.isArray(data.baselines)) {
    throw new Error("The local Python service is running an older version. Restart code/frontend/server.py, then retry.");
  }
  const dataset = data.datasets?.fresh_indelible;
  const heading = el("div", { className: "page-heading" }, [
    el("h1", { text: "SCAW-RVNE benchmarks" }),
  ]);
  if (!dataset) {
    target.replaceChildren(heading, panel("Benchmark results",
      el("p", { className: "bench-copy", text: "INDELible results are missing. Restore or regenerate the benchmark summary." })),
    casePanel(data.case));
    return;
  }

  const patterns = Object.keys(dataset.groups).sort();
  const [patternField, pattern] = labeledSelect("Scenario",
    patterns.map((key) => [key, PATTERNS[key] ?? key.replaceAll("-", " ")]));
  if (patterns.includes("all")) pattern.value = "all";
  const [metricField, metric] = labeledSelect("Metric", METRICS);
  const [methodField, method] = labeledSelect("Method", []);
  const comparison = el("div", { className: "bench-comparison" });
  const overview = methodComparison(data);
  target.replaceChildren(heading, datasetPanel(dataset), ...(overview ? [overview] : []),
    panel("Ablation experiments",
      el("div", { className: "benchmark-filters" }, [patternField, methodField, metricField]),
      comparison),
    casePanel(data.case), runtimePanel(data));

  function chooseMethods() {
    const current = method.value;
    const available = dataset.groups[pattern.value] ?? [];
    method.replaceChildren(option("", "All methods"), ...orderedMethods(available).map((key) => option(key, methodLabel(key))));
    if (available.includes(current)) method.value = current;
    refresh();
  }

  function refresh() {
    const kind = metric.value;
    const methods = orderedMethods(dataset.groups[pattern.value] ?? []).filter((key) => !method.value || method.value === key);
    const rows = methods.map((key) => dataset.rows.find((row) =>
      row.pattern === pattern.value && row.method === key && row.metric === kind
    ) ?? { method: key, mean: null, sd: null, n: 0 });
    const boundaryMetric = kind === "boundary_f1" || kind === "boundary_f1_exact";
    const missingBoundaries = boundaryMetric && methods.includes("fixed_uniform_7");
    comparison.replaceChildren(barChart(rows, kind),
      table(["Method", "Mean", "SD", "Families (n)"], rows.map((row) => [
        methodLabel(row.method), fullStat(row.mean), fullStat(row.sd), row.n,
      ])),
      ...(missingBoundaries ? [el("p", { className: "bench-copy", text: "The uniform window does not predict boundaries." })] : []));
  }

  pattern.addEventListener("change", chooseMethods);
  method.addEventListener("change", refresh);
  metric.addEventListener("change", refresh);
  chooseMethods();
}
