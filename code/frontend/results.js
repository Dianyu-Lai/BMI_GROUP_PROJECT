import { button, clamp, displayScore, el, linePath, methodName, percent, svg, toolName, valueRange } from "./ui.js";

const X = 1000;
const PLOT_HEIGHT = 154;
const ROW_HEIGHT = 27;
// Use the residue groups from the MView P1 palette for the SMART-style alignment.
const RESIDUE_GROUPS = [
  ["AGILMPV", "aa-hydrophobic", "Hydrophobic"],
  ["FHWY", "aa-aromatic", "Aromatic"],
  ["KR", "aa-basic", "Basic"],
  ["DE", "aa-acidic", "Acidic"],
  ["NQ", "aa-amides", "Amides"],
  ["ST", "aa-ser-thr", "Ser / Thr"],
  ["C", "aa-cysteine", "Cysteine"],
  ["-.", "aa-gap", "Gap"],
  ["BJOUXZ", "aa-other", "Other"],
];

const residueClasses = new Map(RESIDUE_GROUPS.flatMap(([letters, key]) => [...letters].map((letter) => [letter, key])));
const residueClass = (letter) => residueClasses.get(letter) ?? "aa-other";
const measure = (name, value) => el("div", {}, [el("dt", { text: name }), el("dd", { text: value })]);

function switcher(label, labels, select) {
  const group = el("div", { className: "switcher", role: "group", "aria-label": label });
  const switches = labels.map(([value, text], index) => {
    const control = button(text, { "aria-pressed": String(!index) });
    control.addEventListener("click", () => {
      switches.forEach((item) => item.setAttribute("aria-pressed", String(item === control)));
      select(value);
    });
    return control;
  });
  group.append(...switches);
  return group;
}

function toggle(text, checked, change) {
  const input = el("input", { type: "checkbox" });
  input.checked = checked;
  input.addEventListener("change", () => change(input.checked));
  return el("label", { className: "check-control" }, [input, text]);
}

function makeDownloadMenu(id, hasSegments) {
  const downloads = el("details", { className: "download-menu" }, [el("summary", { className: "button", text: "Download" })]);
  const files = [["alignment.fasta", "Aligned FASTA"], ["scores.tsv", "Site scores TSV"], ["settings.json", "Settings JSON"]];
  if (hasSegments) files.splice(2, 0, ["segments.tsv", "Segments TSV"]);
  const links = files.map(([file, label]) => {
    const link = el("a", { href: `/api/jobs/${id}/download/${file}`, download: file, text: label });
    link.addEventListener("click", () => { downloads.open = false; });
    return link;
  });
  downloads.append(el("div", { className: "download-links" }, links));
  return downloads;
}

function color(value, low, high) {
  const t = high === low ? 0.5 : clamp((value - low) / (high - low), 0, 1);
  return `rgb(${[230, 244, 248].map((start, i) => Math.round(start + t * ([18, 59, 86][i] - start))).join(",")})`;
}

function extent(values) {
  const scored = values.filter(Number.isFinite);
  return scored.length ? [Math.min(...scored), Math.max(...scored)] : null;
}

export function renderResults(target, data, id, resetAnalysis) {
  const length = data.sites.length;
  const isScaw = data.method === "scaw-rvne";
  const finalName = methodName(data.method);
  const minSpan = Math.min(8, length);
  let maxSpan = Math.min(64, length);
  let preferredSpan = 36;
  let span = Math.min(preferredSpan, length);
  let start = 0;
  let selected = 0;
  let activeSeries = "score";
  let view = "line";
  let showHeatTrack = false;
  let showSegments = data.segments.length > 0;
  let showTruth = data.truth.length > 0;
  let resizePending = false;
  let dragged = null;
  const isMotifExample = data.example === "pattern-c";
  const truthRegions = [];
  for (const site of data.truth) {
    const previous = truthRegions.at(-1);
    if (previous?.rate === site.rate) previous.end = site.position;
    else truthRegions.push({ start: site.position, end: site.position, rate: site.rate, conserved: site.conserved });
  }
  const series = {
    score: { label: finalName, style: "final-line", values: data.sites.map((row) => row.score) },
    raw_rvne: { label: "Raw RVNE", style: "raw-line", values: data.sites.map((row) => row.raw_rvne) },
  };
  for (const item of Object.values(series)) {
    item.bounds = extent(item.values);
    item.axis = valueRange(item.values);
  }
  const heading = el("h1", { text: "Conservation results", tabindex: "-1" });
  const newButton = button("New analysis");
  newButton.addEventListener("click", resetAnalysis);
  const summary = el("dl", { className: "result-summary" }, [
    measure("File", data.filename), measure("Sequences", data.ids.length),
    measure("Aligned columns", length), measure("Scored columns", `${data.scored}/${length} · ${percent(data.scored / length)}`),
    measure("Method", finalName), measure("Alignment", data.aligner
      ? `${toolName(data.aligner)}${data.alignerVersion ? ` ${data.alignerVersion}` : ""}` : "Pre-aligned FASTA"),
  ]);
  const presetName = { frozen: "Frozen benchmark settings", default: "Software defaults", custom: "Custom" };
  const metadata = [data.source, presetName[data.preset]].filter(Boolean);
  const rangeLabel = el("output", { className: "range-label", "aria-label": "Visible alignment columns" });
  const left = button("‹", { className: "button icon-button", "aria-label": "Pan left", title: "Pan left" });
  const right = button("›", { className: "button icon-button", "aria-label": "Pan right", title: "Pan right" });
  const zoom = el("input", { type: "range", min: minSpan, max: maxSpan, value: span, "aria-label": "Columns per view" });
  const overview = svg("svg", {
    class: "overview-plot", viewBox: `0 0 ${X} 44`, preserveAspectRatio: "none", tabindex: 0, role: "slider",
    "aria-valuemin": 1, "aria-label": "Alignment overview. Drag the visible window, or use arrow keys to pan.",
  });
  const overviewLabel = el("div", { className: "overview-title" }, [el("h2", { text: "Alignment overview" }), rangeLabel]);
  const overviewEnd = el("div", { className: "overview-ends" }, [el("span", { text: "1" }), el("span", { text: length })]);
  const map = el("div", { className: "overview-map" }, [overview, overviewEnd]);
  const rangeControls = el("div", { className: "range-controls" }, [
    left, right, el("label", { className: "zoom-label" }, [zoom, "Columns per view"]),
  ]);
  if (isMotifExample) {
    const motif = truthRegions.find((region) => region.conserved);
    const jump = button("Go to motif");
    jump.addEventListener("click", () => selectColumn(motif.start - 1, true));
    rangeControls.append(jump);
  }
  const overviewPanel = el("section", { className: "panel overview-panel" }, [overviewLabel, map, rangeControls]);

  const plot = svg("svg", {
    class: "line-plot", viewBox: `0 0 ${X} ${PLOT_HEIGHT}`, preserveAspectRatio: "none", tabindex: 0,
    role: "img", "aria-label": "Conservation line chart. Use the arrow keys to select a column.",
  });
  const yLabels = el("div", { className: "y-labels", "aria-hidden": "true" });
  const line = el("div", { className: "track line-track" }, [yLabels, plot]);
  const plotLabel = el("div", { className: "plot-label" });
  const scoreToggle = switcher("Score series", [["score", "Final score"], ["raw_rvne", "Raw RVNE"]], (next) => {
    activeSeries = next;
    refresh();
  });
  const viewToggle = switcher("Profile view", [["line", "Line"], ["heat", "Heatmap"]], (next) => {
    view = next;
    refresh();
  });
  const heat = el("div", { className: "heatmap" });
  const trackHeat = el("div", { className: "heatmap compact-heat" });
  const overlays = el("div", { className: "overlay-tracks" });
  const profile = el("section", { className: "panel profile-panel" }, [
    el("div", { className: "panel-heading profile-heading" }, [
      el("h2", { text: "Conservation profile" }), el("div", { className: "profile-controls" }, [viewToggle, scoreToggle]),
    ]),
    el("div", { className: "profile-info" }, [
      plotLabel, toggle("Score track above MSA", false, (value) => { showHeatTrack = value; refresh(); }),
    ]),
    line, heat, overlays,
  ]);
  const trackControls = el("div", { className: "layer-controls" });
  if (data.segments.length) trackControls.append(toggle("SCAW segments", showSegments, (value) => { showSegments = value; refresh(); }));
  if (data.truth.length) trackControls.append(toggle("Synthetic truth", showTruth, (value) => { showTruth = value; refresh(); }));
  if (trackControls.childElementCount) profile.append(trackControls);

  const ruler = el("div", { className: "msa-ruler track" });
  const plane = el("div", { className: "msa-plane", style: `height:${ROW_HEIGHT * data.ids.length}px` });
  const msa = el("div", {
    className: "msa-viewport", tabindex: 0, role: "region",
    "aria-label": "Multiple sequence alignment. Arrow keys select positions, Page Up and Page Down scroll sequences.",
  }, [plane]);
  const msaLegend = el("div", { className: "chemistry-legend", role: "group", "aria-label": "Residue classes" }, RESIDUE_GROUPS.map(
    ([letters, key, name]) => el("span", { className: "chemistry-key", title: `${name} (${[...letters].join(", ")})` },
      [el("i", { className: key, "aria-hidden": "true" }), name]),
  ));
  const msaPanel = el("section", { className: "panel msa-panel" }, [
    el("div", { className: "panel-heading msa-heading" }, [
      el("h2", { text: "Multiple sequence alignment" }), msaLegend,
    ]),
    trackHeat, ruler, msa,
  ]);
  const detail = el("aside", { className: "panel position-panel", "aria-label": "Selected position" });
  const positionLabel = el("p", { className: "position-number", "aria-live": "polite" });
  const positionValues = el("dl", { className: "position-values" });
  const countList = el("ul", { className: "residue-counts" });
  detail.append(
    el("div", { className: "panel-heading" }, [el("h2", { text: "Selected position" })]),
    el("div", { className: "position-body" }, [
      positionLabel, positionValues, el("h3", { text: "Residue counts" }), countList,
    ]),
  );
  const tracks = el("div", { className: "track-stack" }, [overviewPanel, profile, msaPanel]);
  target.replaceChildren(
    el("div", { className: "result-heading" }, [
      heading, el("div", { className: "actions" }, [newButton, makeDownloadMenu(id, isScaw)]),
    ]),
    summary, el("p", { className: "result-source", text: metadata.join(" · ") }),
    el("div", { className: "results-layout" }, [tracks, detail]),
  );

  const lastStart = () => length - span;
  const xAt = (index, count) => X * (index + 0.5) / count;

  function selectColumn(index, center = false) {
    selected = clamp(index, 0, length - 1);
    if (center) start = clamp(selected - Math.floor(span / 2), 0, lastStart());
    else if (selected < start) start = selected;
    else if (selected >= start + span) start = selected - span + 1;
    refresh();
  }

  function panTo(index) {
    start = clamp(index, 0, lastStart());
    refresh();
  }

  function keyboard(event) {
    if (event.key === "ArrowLeft") selectColumn(selected - 1);
    else if (event.key === "ArrowRight") selectColumn(selected + 1);
    else if (event.key === "Home") selectColumn(0);
    else if (event.key === "End") selectColumn(length - 1);
    else return;
    event.preventDefault();
  }

  function cellClick(event) {
    const cell = event.target.closest("[data-column]");
    if (cell) selectColumn(Number(cell.dataset.column));
  }

  function renderOverview() {
    const active = series[activeSeries];
    const [low, high] = active.axis;
    const path = linePath(active.values, (index) => xAt(index, length), (value) => 40 - 35 * (value - low) / (high - low));
    const leftEdge = start * X / length;
    const width = span * X / length;
    const children = [];
    if (showTruth) {
      for (const region of truthRegions.filter((item) => item.conserved)) {
        children.push(svg("rect", { x: (region.start - 1) * X / length, y: 0,
          width: (region.end - region.start + 1) * X / length, height: 44, class: "overview-motif" }));
      }
    }
    children.push(svg("path", { d: path, class: active.style }),
      svg("rect", { x: 0, y: 0, width: leftEdge, height: 44, class: "overview-mask" }),
      svg("rect", { x: leftEdge + width, y: 0, width: Math.max(0, X - leftEdge - width), height: 44, class: "overview-mask" }),
      svg("rect", { x: leftEdge, y: 1, width, height: 42, class: "overview-window" }));
    overview.replaceChildren(...children);
    overview.setAttribute("aria-valuemax", lastStart() + 1);
    overview.setAttribute("aria-valuenow", start + 1);
    overview.setAttribute("aria-valuetext", `Columns ${start + 1}–${start + span} of ${length}`);
    rangeLabel.textContent = `${start + 1}–${start + span} of ${length}`;
    left.disabled = start === 0;
    right.disabled = start === lastStart();
    zoom.value = span;
  }

  function renderLine() {
    const active = series[activeSeries];
    const [low, high] = active.axis;
    const y = (value) => PLOT_HEIGHT - 12 - (PLOT_HEIGHT - 24) * (value - low) / (high - low);
    const ticks = Array.from({ length: 4 }, (_, i) => low + i * (high - low) / 3);
    yLabels.replaceChildren(...ticks.map((value) => el("span", {
      text: value.toFixed(3), style: `top:${100 * y(value) / PLOT_HEIGHT}%`,
    })));
    const nodes = ticks.map((value) => svg("line", { x1: 0, x2: X, y1: y(value), y2: y(value), class: "plot-grid" }));
    if (showSegments) {
      for (const segment of data.segments) {
        const boundary = segment.start - 1;
        if (boundary <= start || boundary >= start + span) continue;
        const x = X * (boundary - start) / span;
        nodes.push(svg("line", { x1: x, x2: x, y1: 0, y2: PLOT_HEIGHT, class: "segment-boundary" }));
      }
    }
    const values = active.values.slice(start, start + span);
    const local = selected - start;
    const visible = local >= 0 && local < span;
    if (visible) {
      nodes.push(svg("rect", { x: X * local / span, y: 0, width: X / span, height: PLOT_HEIGHT, class: "plot-selection" }));
    }
    nodes.push(svg("path", { d: linePath(values, (i) => xAt(i, span), y), class: active.style }));
    values.forEach((value, i) => {
      if (value != null && (i === local || (values[i - 1] == null && values[i + 1] == null))) {
        nodes.push(svg("line", { x1: xAt(i, span), x2: xAt(i, span), y1: y(value) - 2, y2: y(value) + 2, class: "plot-point" }));
      }
    });
    plot.replaceChildren(...nodes);
  }

  function renderHeat(targetHeat, compact) {
    const keys = compact ? [activeSeries] : Object.keys(series);
    const focused = targetHeat.contains(document.activeElement) ? document.activeElement.dataset.series : null;
    targetHeat.replaceChildren(...keys.map((key) => {
      const active = series[key];
      const [low, high] = active.bounds ?? [0, 1];
      const cells = active.values.slice(start, start + span).map((value, index) => el("span", {
        className: `heat-cell${value == null ? " unscored" : ""}${start + index === selected ? " selected-cell" : ""}`,
        "data-column": start + index,
        title: `Column ${start + index + 1} · ${active.label} · ${displayScore(value)}`,
        ...(value == null ? {} : { style: `background:${color(value, low, high)}` }),
      }));
      const scoreRow = el("div", {
        className: "column-grid heat-row", tabindex: 0, role: "group",
        "aria-label": `${active.label} score track. Use arrow keys to select positions.`, "data-series": key,
      }, cells);
      scoreRow.addEventListener("keydown", keyboard);
      scoreRow.addEventListener("click", cellClick);
      const children = [el("div", { className: "track-label", text: active.label }), scoreRow];
      children.push(el("div", { className: "heat-scale" }, [
        el("span", { text: active.bounds ? low.toFixed(3) : "Unscored" }),
        ...(active.bounds ? [el("span", { className: "heat-gradient" }), el("span", { text: high.toFixed(3) })] : []),
      ]));
      return el("div", { className: "track heat-track" }, children);
    }), ...(compact ? [] : [
      el("div", { className: "unscored-key" }, [el("i", { className: "unscored" }), "Unscored"]),
    ]));
    if (focused) targetHeat.querySelector(`[data-series="${focused}"]`)?.focus({ preventScroll: true });
  }

  function renderOverlays() {
    const rows = [];
    if (showSegments) rows.push(["SCAW segments", data.segments, "segment"]);
    if (showTruth) rows.push(["Synthetic truth", truthRegions, "truth"]);
    overlays.replaceChildren(...rows.map(([label, regions, kind]) => {
      const blocks = [];
      for (const region of regions) {
        const from = Math.max(start + 1, region.start);
        const to = Math.min(start + span, region.end);
        if (from > to) continue;
        const name = kind === "segment" ? `Segment ${region.segment_id}`
          : (isMotifExample && region.conserved ? "Motif" : `Rate ${region.rate}`);
        blocks.push(el("span", {
          className: `overlay-block ${kind}${region.conserved ? " motif" : ""}`,
          style: `grid-column:${from - start}/span ${to - from + 1}`,
          text: name,
          title: `${name} · columns ${region.start}–${region.end}`,
        }));
      }
      return el("div", { className: "track overlay-track" }, [
        el("div", { className: "track-label", text: label }), el("div", { className: "column-grid" }, blocks),
      ]);
    }));
  }

  function renderRuler() {
    const positions = data.sites.slice(start, start + span).map((row, i) => {
      const label = (span <= 12 || i === 0 || row.position % 5 === 0) ? row.position : "";
      return el("span", { className: start + i === selected ? "selected-cell" : "", "data-column": start + i,
        title: `Alignment column ${row.position}`, text: label });
    });
    ruler.replaceChildren(el("div", { className: "track-label", text: "Alignment column" }),
      el("div", { className: "column-grid ruler-grid" }, positions));
  }

  function renderRows() {
    const first = Math.max(0, Math.floor(msa.scrollTop / ROW_HEIGHT) - 2);
    const last = Math.min(data.ids.length, first + Math.ceil(msa.clientHeight / ROW_HEIGHT) + 5);
    const rows = [];
    for (let row = first; row < last; row++) {
      const residues = [...data.sequences[row].slice(start, start + span)].map((letter, i) => el("span", {
        className: `residue ${residueClass(letter)}${start + i === selected ? " selected-cell" : ""}`,
        "data-column": start + i,
        title: `${data.ids[row]} · column ${start + i + 1} · ${letter}`,
        text: letter,
      }));
      rows.push(el("div", { className: "track msa-row", style: `top:${row * ROW_HEIGHT}px` }, [
        el("div", { className: "sequence-id", text: data.ids[row], title: data.ids[row] }),
        el("div", { className: "column-grid" }, residues),
      ]));
    }
    plane.replaceChildren(...rows);
  }

  function renderDetail() {
    const site = data.sites[selected];
    const truth = data.truth[selected];
    const counts = new Map();
    for (const sequence of data.sequences) {
      const residue = ".-".includes(sequence[selected]) ? "-" : sequence[selected];
      counts.set(residue, (counts.get(residue) ?? 0) + 1);
    }
    const values = [
      measure("Non-gap coverage", percent(site.nongap_coverage)),
      measure("Final score", displayScore(site.score)),
      measure("Raw RVNE", displayScore(site.raw_rvne)),
    ];
    if (site.window_size != null && isScaw) values.push(measure("Window width", `${site.window_size} columns`));
    if (site.segment_start != null) values.push(measure("Fitted segment", `${site.segment_start}–${site.segment_end}`));
    if (data.reference) values.push(measure(`${data.reference} residue`, data.referencePositions[selected] ?? "Gap"));
    if (truth) {
      const label = truth.conserved ? (isMotifExample ? "Conserved motif" : "Conserved")
        : (truth.rate < 1 ? "Intermediate" : "Variable");
      values.push(measure("True mutation rate", truth.rate), measure("Ground truth", label));
    }
    const residueCounts = [...counts.entries()].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))
      .map(([residue, count]) => el("li", { className: "count" }, [
        el("span", { className: residueClass(residue), text: residue }),
        el("span", { text: count, "aria-label": `${count} sequences` }),
      ]));
    const label = `Column ${site.position}`;
    if (positionLabel.textContent !== label) positionLabel.textContent = label;
    positionValues.replaceChildren(...values);
    countList.replaceChildren(...residueCounts);
  }

  function refresh() {
    tracks.style.setProperty("--columns", span);
    tracks.style.setProperty("--scrollbar-width", `${msa.offsetWidth - msa.clientWidth}px`);
    const active = series[activeSeries];
    plotLabel.replaceChildren(el("i", { className: `line-key ${active.style}` }),
      `${active.label} · ${active.bounds ? "higher is more conserved" : "no scored positions"}`);
    line.hidden = view !== "line";
    heat.hidden = view !== "heat";
    trackHeat.hidden = !showHeatTrack;
    renderOverview();
    if (view === "line") renderLine();
    if (view === "heat") renderHeat(heat, false);
    if (showHeatTrack) renderHeat(trackHeat, true);
    renderOverlays();
    renderRuler();
    renderRows();
    renderDetail();
  }

  function updateSpan() {
    const previous = span;
    span = clamp(preferredSpan, minSpan, maxSpan);
    const center = start + Math.floor(previous / 2);
    panTo(center - Math.floor(span / 2));
  }

  function overviewColumn(event) {
    const bounds = overview.getBoundingClientRect();
    return clamp(Math.floor(length * (event.clientX - bounds.left) / bounds.width), 0, length - 1);
  }

  overview.addEventListener("pointerdown", (event) => {
    if (event.button !== 0) return;
    const column = overviewColumn(event);
    dragged = column >= start && column < start + span ? column - start : Math.floor(span / 2);
    overview.setPointerCapture(event.pointerId);
    panTo(column - dragged);
  });
  overview.addEventListener("pointermove", (event) => { if (dragged != null) panTo(overviewColumn(event) - dragged); });
  overview.addEventListener("lostpointercapture", () => { dragged = null; });
  overview.addEventListener("keydown", (event) => {
    const next = { ArrowLeft: start - 1, ArrowRight: start + 1, Home: 0, End: lastStart() }[event.key];
    if (next == null) return;
    event.preventDefault();
    panTo(next);
  });
  left.addEventListener("click", () => panTo(start - Math.ceil(span / 3)));
  right.addEventListener("click", () => panTo(start + Math.ceil(span / 3)));
  zoom.addEventListener("input", () => {
    preferredSpan = Number(zoom.value);
    updateSpan();
  });
  plot.addEventListener("click", (event) => {
    const bounds = plot.getBoundingClientRect();
    const local = Math.floor((event.clientX - bounds.left) * span / bounds.width);
    selectColumn(start + clamp(local, 0, span - 1));
  });
  plot.addEventListener("keydown", keyboard);
  ruler.addEventListener("click", cellClick);
  msa.addEventListener("click", cellClick);
  msa.addEventListener("keydown", keyboard);
  msa.addEventListener("scroll", () => {
    if (resizePending) return;
    resizePending = true;
    requestAnimationFrame(() => {
      resizePending = false;
      renderRows();
    });
  });
  const resize = new ResizeObserver(() => {
    const available = tracks.clientWidth - parseFloat(getComputedStyle(tracks).getPropertyValue("--track-label")) - 16;
    if (available <= 0) return;
    const next = clamp(Math.floor(available / 16), minSpan, Math.min(64, length));
    if (maxSpan !== next) {
      maxSpan = next;
      zoom.max = maxSpan;
      updateSpan();
    } else renderRows();
  });
  resize.observe(tracks);
  refresh();
  return () => resize.disconnect();
}
