export const byId = (id) => document.getElementById(id);

export function el(tag, options = {}, children = []) {
  const element = document.createElement(tag);
  for (const [key, value] of Object.entries(options)) {
    if (key === "className") element.className = value;
    else if (key === "text") element.textContent = value;
    else element.setAttribute(key, value);
  }
  element.append(...children);
  return element;
}

export function svg(tag, attrs = {}, children = []) {
  const element = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [key, value] of Object.entries(attrs)) element.setAttribute(key, value);
  element.append(...children);
  return element;
}

export function option(value, text) {
  return el("option", { value, text });
}

export async function getJSON(url, options) {
  const response = await fetch(url, options);
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || "Could not reach the local analysis service.");
  return data;
}

export const clamp = (value, min, max) => Math.min(Math.max(min, value), max);
export const displayScore = (value) => value == null ? "Unscored" : value.toFixed(4);
export const percent = (value) => `${(100 * value).toFixed(1)}%`;
export const methodName = (id) => id === "rvne" ? "Original RVNE" : "SCAW-RVNE";
export const toolName = (id) => ({ mafft: "MAFFT", clustalw: "ClustalW" })[id];

export function button(text, attributes = {}) {
  return el("button", { type: "button", className: "button", text, ...attributes });
}

export function valueRange(values) {
  const valid = values.filter(Number.isFinite);
  if (!valid.length) return [0, 1];
  const low = Math.min(...valid);
  const high = Math.max(...valid);
  const padding = (high - low || Math.abs(high) || 1) * 0.06;
  return [low - padding, high + padding];
}

export function linePath(values, x, y) {
  let path = "";
  let next = "M";
  values.forEach((value, index) => {
    if (value == null) {
      next = "M";
      return;
    }
    path += `${next}${x(index).toFixed(3)} ${y(value).toFixed(3)}`;
    next = "L";
  });
  return path;
}
