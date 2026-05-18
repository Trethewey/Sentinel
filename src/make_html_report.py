"""Render the single-file HTML report from the post-processed per-sample table.

Reads:
  <PROJ>/results/per_sample_report.tsv
  <PROJ>/work/run_manifest.tsv         (for the optional `well` column)

Writes:
  <PROJ>/results/report.html

The output is a self-contained HTML document with the logo, an interactive
96-well plate view, a flagged-samples table, a per-sample contamination bar
plot, and a sortable all-samples table. Plotly is loaded from CDN.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd


PROJ = Path(os.environ.get("PROJ", "."))
RESULTS_DIR = PROJ / "results"
MANIFEST = PROJ / "work" / "run_manifest.tsv"
OUT = RESULTS_DIR / "report.html"


LOGO_SVG = """
<svg viewBox="200 20 250 270" width="170" height="74" role="img" aria-label="Sentinel">
  <path d="M274 38 L406 38 L406 156 Q406 198 340 218 Q274 198 274 156 Z" fill="#042c53"/>
  <rect x="315" y="42" width="50" height="8" rx="2" fill="#85b7eb"/>
  <rect x="287" y="51" width="62" height="8" rx="2" fill="#85b7eb"/>
  <rect x="273" y="60" width="54" height="8" rx="2" fill="#85b7eb"/>
  <rect x="258" y="69" width="60" height="8" rx="2" fill="#85b7eb"/>
  <rect x="256" y="78" width="58" height="8" rx="2" fill="#85b7eb"/>
  <rect x="265" y="87" width="52" height="8" rx="2" fill="#85b7eb"/>
  <rect x="274" y="96" width="64" height="8" rx="2" fill="#85b7eb"/>
  <rect x="298" y="105" width="56" height="8" rx="2" fill="#85b7eb"/>
  <rect x="320" y="114" width="58" height="8" rx="2" fill="#85b7eb"/>
  <rect x="337" y="123" width="66" height="8" rx="2" fill="#85b7eb"/>
  <rect x="359" y="132" width="54" height="8" rx="2" fill="#85b7eb"/>
  <rect x="363" y="141" width="62" height="8" rx="2" fill="#f09595"/>
  <rect x="365" y="150" width="58" height="8" rx="2" fill="#85b7eb"/>
  <rect x="354" y="159" width="60" height="8" rx="2" fill="#85b7eb"/>
  <rect x="341" y="168" width="52" height="8" rx="2" fill="#85b7eb"/>
  <rect x="314" y="177" width="62" height="8" rx="2" fill="#85b7eb"/>
  <rect x="300" y="186" width="46" height="8" rx="2" fill="#85b7eb"/>
  <path d="M274 38 L406 38 L406 156 Q406 198 340 218 Q274 198 274 156 Z" fill="none" stroke="#85b7eb" stroke-width="2"/>
  <line x1="425" y1="145" x2="433" y2="145" stroke="#f09595" stroke-width="1"/>
  <circle cx="439" cy="145" r="6" fill="#f09595"/>
  <text x="340" y="255" text-anchor="middle" font-family="Georgia, serif" font-size="34" font-weight="700" letter-spacing="10" fill="#b5d4f4">SENTINEL</text>
  <text x="340" y="278" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" letter-spacing="3" fill="#85b7eb">NGS CONTAMINATION DETECTION</text>
</svg>
"""


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Sentinel report</title>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js" charset="utf-8"></script>
<style>
:root {
  --bg: #0b1426; --bg-deep: #0f1b32; --card: #131f3a; --card-strong: #1a2944;
  --ink: #e8eef9; --muted: #7e8fae; --line: #233358; --line-soft: #1c2a47;
  --pass: #4ade80; --warn: #fbbf24; --fail: #f87171;
  --empty: #6b7589; --accent: #5aa9ff; --accent-soft: rgba(90,169,255,0.12);
}
* { box-sizing: border-box; }
body { margin: 0; font-family: ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif; color: var(--ink); background: var(--bg); line-height: 1.45; -webkit-font-smoothing: antialiased; }
header { background: linear-gradient(180deg, #0a1426 0%, #0e1932 100%); color: var(--ink); border-bottom: 1px solid var(--line); padding: 18px 32px; display: flex; align-items: center; gap: 14px; }
header .meta { margin-left: auto; color: var(--muted); font-size: 12px; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; text-align: right; }
main { padding: 24px 28px 64px; max-width: 1600px; margin: 0 auto; }
.split { display: grid; grid-template-columns: minmax(640px, 1fr) minmax(440px, 580px); gap: 8px; align-items: start; }
@media (max-width: 1280px) { .split { grid-template-columns: 1fr; } }
.split section { margin-bottom: 0; }
.split .table-scroll { max-height: 540px; }
section { background: var(--card); border: 1px solid var(--line); border-radius: 12px; padding: 18px 22px; margin-bottom: 14px; }
h2 { margin: 0 0 12px; font-size: 12px; font-weight: 600; color: var(--muted); letter-spacing: 1.2px; text-transform: uppercase; }
.toolbar { display: flex; gap: 10px; align-items: center; margin-bottom: 12px; flex-wrap: wrap; }
.toolbar label { font-size: 12px; color: var(--muted); }
.toolbar button { padding: 6px 14px; border-radius: 6px; border: 1px solid var(--line); background: var(--card-strong); color: var(--muted); cursor: pointer; font-size: 12px; font-weight: 500; letter-spacing: 0.3px; }
.toolbar button:hover { border-color: var(--accent); color: var(--ink); }
.toolbar button.active { background: var(--accent); color: #0a1426; border-color: var(--accent); font-weight: 600; }
.toggle-btn { min-width: 200px; text-align: center; transition: none; }
.toggle-btn.active { background: var(--accent); color: #0a1426; border-color: var(--accent); font-weight: 600; }
.toggle-btn:not(.active) { background: var(--card-strong); color: var(--muted); border-color: var(--line); }
.plate-wrap { display: flex; justify-content: center; padding: 18px 12px; background: var(--bg-deep); border-radius: 8px; border: 1px solid var(--line); align-items: center; gap: 18px; }
.inline-scale { display: flex; align-items: center; gap: 4px; font-size: 12px; }
.inline-scale .bar { display: inline-block; width: 110px; height: 10px; border-radius: 3px; background: linear-gradient(to right, rgb(165,214,167), rgb(255,183,77), rgb(229,57,53)); margin: 0 6px; }
svg.plate { width: 100%; max-width: 720px; height: auto; }
.row-label, .col-label { font: 700 16px/1 ui-sans-serif, system-ui, sans-serif; fill: #0a1426; text-anchor: middle; dominant-baseline: middle; }
.well-id { font: 700 10px/1 ui-sans-serif, system-ui, sans-serif; fill: #0a1426; text-anchor: middle; dominant-baseline: middle; pointer-events: none; }
.well-id.light { fill: #ffffff; }
.well-circle { stroke-width: 3; cursor: pointer; }
.well-circle.dim { opacity: 0.18; }
.well-circle.empty { fill: var(--empty); stroke: var(--line); stroke-width: 1.5; }
.well-circle:hover { stroke-width: 4; }
.well-circle.active { stroke-width: 5; filter: drop-shadow(0 0 8px rgba(90,169,255,0.8)); }
.donor-link { fill: none; stroke-linecap: round; opacity: 0.85; pointer-events: none; }
.swatch { width: 14px; height: 14px; border-radius: 50%; display: inline-block; border: 1px solid #ccc; }
.table-scroll { max-height: 540px; overflow-y: auto; border: 1px solid var(--line); border-radius: 8px; background: var(--card-strong); }
table.flagged { width: 100%; border-collapse: collapse; font-size: 12.5px; color: var(--ink); }
table.flagged th, table.flagged td { padding: 9px 12px; text-align: left; border-bottom: 1px solid var(--line-soft); }
table.flagged th { color: var(--muted); font-weight: 500; text-transform: uppercase; font-size: 10.5px; letter-spacing: 0.8px; position: sticky; top: 0; background: var(--card); z-index: 1; border-bottom: 1px solid var(--line); }
table.flagged td.num { text-align: right; font-variant-numeric: tabular-nums; }
table.flagged tr { cursor: pointer; }
table.flagged tr:nth-child(even) { background: rgba(255,255,255,0.02); }
table.flagged tr:hover { background: var(--accent-soft); }
table.flagged tr.row-selected { background: rgba(251,191,36,0.10); }
table.flagged tr.row-active { box-shadow: inset 3px 0 0 var(--accent); }
.verdict-pill { display: inline-block; padding: 2px 9px; border-radius: 999px; font-size: 10.5px; font-weight: 600; letter-spacing: 0.5px; }
.verdict-pill.pass { background: rgba(74,222,128,0.16); color: var(--pass); }
.verdict-pill.warn { background: rgba(251,191,36,0.18); color: var(--warn); }
.verdict-pill.fail { background: rgba(248,113,113,0.18); color: var(--fail); }
.toolbar button:disabled { opacity: 0.35; cursor: not-allowed; }
.filter-select, .filter-input { background: var(--card-strong); color: var(--ink); border: 1px solid var(--line); border-radius: 6px; padding: 5px 9px; font-size: 12px; }
.spark { display: inline-block; width: 90px; height: 16px; vertical-align: middle; }
.bar-bg { fill: var(--line-soft); }
.bar-fg.pass { fill: var(--pass); }
.bar-fg.warn { fill: var(--warn); }
.bar-fg.fail { fill: var(--fail); }
input[type="checkbox"] { accent-color: var(--accent); }
.tooltip { position: absolute; pointer-events: none; background: #1f2937; color: white; padding: 8px 12px; border-radius: 6px; font-size: 12px; line-height: 1.4; display: none; z-index: 100; box-shadow: 0 4px 12px rgba(0,0,0,0.15); max-width: 220px; }
::-webkit-scrollbar { width: 8px; }
::-webkit-scrollbar-track { background: var(--card-strong); }
::-webkit-scrollbar-thumb { background: var(--line); border-radius: 4px; }
</style>
</head>
<body>

<header>
__LOGO__
  <div class="meta">__META_HTML__</div>
</header>

<main>

<div class="split">
<section>
  <div class="toolbar" style="justify-content:space-between">
    <div style="display:flex;gap:10px;align-items:center">
      <button id="filter-btn" class="toggle-btn" data-state="all">Show flagged wells only</button>
      <button id="links-btn" class="toggle-btn" data-state="off">Show donor lines</button>
    </div>
    <div class="inline-scale" style="font-size:11px">
      <span style="color:var(--muted);margin-right:6px">contamination</span>
      <span style="color:var(--muted)">clean</span>
      <span class="bar"></span>
      <span style="color:var(--muted)">heavy</span>
    </div>
  </div>
  <div class="plate-wrap" style="margin-top:14px">
    <svg class="plate" id="plate" viewBox="0 0 880 580" aria-label="96-well plate">
      <rect x="20" y="20" width="840" height="540" rx="16" ry="16" fill="white" stroke="#cfd4dc" stroke-width="2"/>
      <g id="labels"></g>
      <g id="wells"></g>
      <g id="links-layer"></g>
    </svg>
    <div style="display:flex;flex-direction:column;gap:10px;font-size:12px;color:var(--ink)">
      <span style="display:inline-flex;align-items:center;gap:8px"><span class="swatch" style="border:3px solid var(--pass);background:white"></span> PASS</span>
      <span style="display:inline-flex;align-items:center;gap:8px"><span class="swatch" style="border:3px solid var(--warn);background:white"></span> WARN</span>
      <span style="display:inline-flex;align-items:center;gap:8px"><span class="swatch" style="border:3px solid var(--fail);background:white"></span> FAIL</span>
    </div>
  </div>
</section>

<section>
  <h2>Flagged Samples</h2>
  <div class="toolbar" style="justify-content:flex-end;margin-bottom:8px">
    <div style="display:flex;gap:8px;align-items:center">
      <button id="copy-sel" disabled>Copy sample IDs</button>
      <button id="csv-sel" disabled>Download as CSV</button>
    </div>
  </div>
  <div class="table-scroll">
    <table class="flagged" id="flagged-table">
      <thead><tr>
        <th style="width:34px"></th>
        <th>Well</th><th>Sample ID</th><th>Verdict</th><th>Type</th>
        <th class="num">Contamination</th><th class="num">Top-donor strength</th>
        <th>Inferred donor</th><th>Donor well</th>
      </tr></thead>
      <tbody id="flagged-rows"></tbody>
    </table>
  </div>
  <p style="margin:10px 0 0;font-size:12px;color:var(--muted)">
    Click a row to highlight that well on the plate. Tick the boxes to select samples for export.
  </p>
</section>
</div>

<section>
  <h2>Contamination level (every well, by plate order)</h2>
  <div id="bar-plot" style="height:300px"></div>
</section>

<section>
  <h2>All samples</h2>
  <div class="toolbar" style="gap:10px">
    <label>Verdict</label>
    <select id="all-verdict" class="filter-select">
      <option value="all">all</option><option value="FAIL">FAIL</option>
      <option value="WARN">WARN</option><option value="PASS">PASS</option>
    </select>
    <label>Type</label>
    <select id="all-type" class="filter-select">
      <option value="all">all</option>__TYPE_OPTIONS__
    </select>
    <label>Search</label>
    <input id="all-search" type="text" placeholder="well or sample id" class="filter-input">
  </div>
  <div class="table-scroll" style="max-height:520px">
    <table class="flagged" id="all-table">
      <thead><tr>
        <th>Well</th><th>Sample ID</th><th>Verdict</th><th>Type</th>
        <th class="num">Contamination</th><th class="num">Top-donor strength</th>
        <th>Top donor</th><th>Signal</th>
      </tr></thead>
      <tbody id="all-rows"></tbody>
    </table>
  </div>
</section>

</main>

<div class="tooltip" id="tt"></div>

<script>
const DATA = __DATA_JSON__;
const samples = DATA.samples;

const COL = {PASS:"#16a34a", WARN:"#d97706", FAIL:"#dc2626"};
const COL_DARK = {PASS:"#4ade80", WARN:"#fbbf24", FAIL:"#f87171"};
const ROWS = ["A","B","C","D","E","F","G","H"];
const COLS = [1,2,3,4,5,6,7,8,9,10,11,12];
const X0 = 110, Y0 = 90, DX = 60, DY = 58, R = 22;
const MAX_C = Math.max(0.01, ...samples.map(s => s.defl));
const PLATE_MAX = Math.max(0.05, MAX_C * 1.1);

function wellCenter(well){
  const r = well.charAt(0), c = parseInt(well.slice(1), 10);
  return {cx: X0 + (c - 1) * DX, cy: Y0 + (ROWS.indexOf(r)) * DY};
}
function gradient(t){
  t = Math.max(0, Math.min(1, t));
  const stops = [[0,[165,214,167]],[0.5,[255,183,77]],[1,[229,57,53]]];
  for (let i=1; i<stops.length; i++){
    if (t <= stops[i][0]){
      const lo = stops[i-1], hi = stops[i];
      const u = (t - lo[0]) / (hi[0] - lo[0]);
      return `rgb(${Math.round(lo[1][0]+(hi[1][0]-lo[1][0])*u)},${Math.round(lo[1][1]+(hi[1][1]-lo[1][1])*u)},${Math.round(lo[1][2]+(hi[1][2]-lo[1][2])*u)})`;
    }
  }
  return "rgb(127,0,0)";
}

const wellMap = {};
const sampleToWell = {};
samples.forEach(s => {
  if (s.well) { wellMap[s.well] = s; sampleToWell[s.id] = s.well; }
});

const labelsG = document.getElementById("labels");
COLS.forEach(c => labelsG.insertAdjacentHTML("beforeend", `<text class="col-label" x="${X0 + (c - 1) * DX}" y="${Y0 - 30}">${c}</text>`));
ROWS.forEach((r, i) => labelsG.insertAdjacentHTML("beforeend", `<text class="row-label" x="${X0 - 35}" y="${Y0 + i * DY}">${r}</text>`));

function render(){
  const wellsG = document.getElementById("wells");
  wellsG.innerHTML = "";
  ROWS.forEach((r, ri) => COLS.forEach(c => {
    const well = `${r}${c}`;
    const s = wellMap[well];
    const {cx, cy} = wellCenter(well);
    if (!s){
      wellsG.insertAdjacentHTML("beforeend", `<g data-well="${well}"><circle class="well-circle empty" cx="${cx}" cy="${cy}" r="${R}"></circle></g>`);
      return;
    }
    const t = s.defl / PLATE_MAX;
    const fill = gradient(t);
    const ring = COL[s.verdict];
    const shortId = (s.id || "").slice(-5);
    wellsG.insertAdjacentHTML("beforeend",
      `<g data-well="${well}">
        <circle class="well-circle" cx="${cx}" cy="${cy}" r="${R}" fill="${fill}" stroke="${ring}"></circle>
        <text class="well-id" x="${cx}" y="${cy + 1}">${shortId}</text>
      </g>`);
  }));
  attachHover();
  drawLinks();
}

let isolatedWell = null;

function applyFilter(state){
  isolatedWell = null;
  const keep = new Set();
  if (state === "flagged"){
    Object.values(wellMap).forEach(s => {
      if (s.verdict === "WARN" || s.verdict === "FAIL"){
        keep.add(s.well);
        const dw = sampleToWell[s.donor];
        if (dw) keep.add(dw);
      }
    });
  }
  document.querySelectorAll("g[data-well] .well-circle").forEach(c => {
    const w = c.parentNode.dataset.well;
    const dim = state === "flagged" && !keep.has(w);
    c.classList.toggle("dim", dim);
    c.classList.remove("active");
    const txt = c.parentNode.querySelector(".well-id");
    if (txt) txt.style.opacity = dim ? "0.15" : "1";
  });
  document.querySelectorAll("#flagged-rows tr.row-active").forEach(r => r.classList.remove("row-active"));
  drawLinks();
}

function drawLinks(){
  const layer = document.getElementById("links-layer");
  layer.innerHTML = "";
  const btn = document.getElementById("links-btn");
  if (btn && !btn.classList.contains("active")) return;
  Object.values(wellMap).forEach(s => {
    if (s.verdict === "PASS") return;
    const donorWell = sampleToWell[s.donor];
    if (!donorWell) return;
    const donor = wellCenter(donorWell), recipient = wellCenter(s.well);
    const dx = recipient.cx - donor.cx, dy = recipient.cy - donor.cy;
    const L = Math.hypot(dx, dy);
    if (L < 1) return;
    const ux = dx / L, uy = dy / L;
    const startX = donor.cx + ux * R, startY = donor.cy + uy * R;
    const endX = recipient.cx - ux * R, endY = recipient.cy - uy * R;
    const adj = isAdjacent(s.well, donorWell);
    layer.insertAdjacentHTML("beforeend",
      `<line class="donor-link" x1="${startX}" y1="${startY}" x2="${endX}" y2="${endY}" stroke="#e8eef9" stroke-width="${adj?2.5:1.6}" stroke-dasharray="${adj?'4 3':''}"/>`);
  });
}
function isAdjacent(w1, w2){
  const r1 = ROWS.indexOf(w1.charAt(0)), c1 = parseInt(w1.slice(1), 10);
  const r2 = ROWS.indexOf(w2.charAt(0)), c2 = parseInt(w2.slice(1), 10);
  return Math.max(Math.abs(r1 - r2), Math.abs(c1 - c2)) === 1;
}

const tt = document.getElementById("tt");
function attachHover(){
  document.querySelectorAll("g[data-well]").forEach(g => {
    g.addEventListener("mouseenter", () => {
      const w = g.dataset.well, s = wellMap[w];
      if (!s){ tt.innerHTML = `<div style="font-weight:600">${w}</div><div style="color:#9ca3af;font-size:11px">empty</div>`; }
      else {
        tt.innerHTML = `
          <div style="font-weight:600">${s.id} (${w})</div>
          <div style="color:#9ca3af;font-size:11px">${s.type} &middot; <b style="color:${COL[s.verdict]}">${s.verdict}</b></div>
          <div style="margin-top:4px">contamination: ${(s.defl*100).toFixed(2)}%</div>
          <div>top donor: ${s.donor || '-'}</div>`;
      }
      tt.style.display = "block";
    });
    g.addEventListener("mousemove", e => { tt.style.left = (e.pageX + 14) + "px"; tt.style.top = (e.pageY + 14) + "px"; });
    g.addEventListener("mouseleave", () => tt.style.display = "none");
  });
}

document.getElementById("filter-btn").addEventListener("click", function(){
  this.classList.toggle("active");
  applyFilter(this.classList.contains("active") ? "flagged" : "all");
});
document.getElementById("links-btn").addEventListener("click", function(){
  this.classList.toggle("active");
  drawLinks();
});

render();

const flagged = samples.filter(s => s.verdict === "WARN" || s.verdict === "FAIL")
  .sort((a, b) => (a.verdict === b.verdict ? b.defl - a.defl : (a.verdict === "FAIL" ? -1 : 1)));
const tbody = document.getElementById("flagged-rows");
flagged.forEach((s, i) => {
  const donorWell = sampleToWell[s.donor] || "-";
  tbody.insertAdjacentHTML("beforeend", `
    <tr data-well="${s.well || ''}" data-idx="${i}">
      <td><input type="checkbox" class="row-cb" data-id="${s.id}"></td>
      <td><b>${s.well || '-'}</b></td>
      <td><b>${s.id}</b></td>
      <td><span class="verdict-pill ${s.verdict.toLowerCase()}">${s.verdict}</span></td>
      <td>${s.type}</td>
      <td class="num">${(s.defl * 100).toFixed(2)}%</td>
      <td class="num">${(s.score * 100).toFixed(2)}%</td>
      <td>${s.donor || '-'}</td>
      <td>${donorWell}</td>
    </tr>`);
});
tbody.querySelectorAll("tr").forEach(tr => {
  tr.addEventListener("click", e => {
    if (e.target.tagName === "INPUT") return;
    const well = tr.dataset.well;
    if (!well) return;
    if (isolatedWell === well){
      const state = document.getElementById("filter-btn").classList.contains("active") ? "flagged" : "all";
      applyFilter(state);
      return;
    }
    isolateRow(tr, well);
  });
});
function isolateRow(tr, well){
  const s = wellMap[well];
  const donorWell = s ? sampleToWell[s.donor] : null;
  const keep = new Set([well]);
  if (donorWell) keep.add(donorWell);
  isolatedWell = well;
  document.querySelectorAll("g[data-well] .well-circle").forEach(c => {
    const w = c.parentNode.dataset.well;
    const dim = !keep.has(w);
    c.classList.toggle("dim", dim);
    c.classList.toggle("active", keep.has(w));
    const txt = c.parentNode.querySelector(".well-id");
    if (txt) txt.style.opacity = dim ? "0.10" : "1";
  });
  document.querySelectorAll("#flagged-rows tr.row-active").forEach(r => r.classList.remove("row-active"));
  tr.classList.add("row-active");
  const layer = document.getElementById("links-layer");
  layer.innerHTML = "";
  if (donorWell){
    const donor = wellCenter(donorWell), recipient = wellCenter(well);
    const dx = recipient.cx - donor.cx, dy = recipient.cy - donor.cy;
    const L = Math.hypot(dx, dy);
    if (L >= 1){
      const ux = dx / L, uy = dy / L;
      layer.insertAdjacentHTML("beforeend",
        `<line class="donor-link" x1="${donor.cx + ux * R}" y1="${donor.cy + uy * R}" x2="${recipient.cx - ux * R}" y2="${recipient.cy - uy * R}" stroke="#e8eef9" stroke-width="2.5"/>`);
    }
  }
}

const btnCopy = document.getElementById("copy-sel");
const btnCsv = document.getElementById("csv-sel");
function refreshSelection(){
  let n = 0;
  document.querySelectorAll(".row-cb").forEach(cb => {
    const tr = cb.closest("tr");
    if (cb.checked){ n++; tr.classList.add("row-selected"); } else tr.classList.remove("row-selected");
  });
  btnCopy.disabled = n === 0; btnCsv.disabled = n === 0;
}
document.querySelectorAll(".row-cb").forEach(cb => cb.addEventListener("change", refreshSelection));
btnCopy.addEventListener("click", () => {
  const ids = Array.from(document.querySelectorAll(".row-cb")).filter(c => c.checked).map(c => c.dataset.id).join("\n");
  navigator.clipboard.writeText(ids).then(() => { const o = btnCopy.textContent; btnCopy.textContent = "Copied"; setTimeout(() => btnCopy.textContent = o, 1200); });
});
btnCsv.addEventListener("click", () => {
  const rows = Array.from(document.querySelectorAll(".row-cb")).filter(c => c.checked).map(c => {
    const tr = c.closest("tr"), cells = tr.querySelectorAll("td");
    return [cells[1].textContent, cells[2].textContent, cells[3].textContent, cells[4].textContent,
            cells[5].textContent, cells[6].textContent, cells[7].textContent, cells[8].textContent];
  });
  const cols = ["well","sample_id","verdict","type","contamination","top_donor_strength","donor","donor_well"];
  const csv = [cols.join(",")].concat(rows.map(r => r.map(c => `"${c.trim()}"`).join(","))).join("\n");
  const blob = new Blob([csv], {type: "text/csv"});
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob); a.download = "sentinel_flagged_samples.csv"; a.click();
  URL.revokeObjectURL(a.href);
});

function wellSortKey(w){
  if (!w) return 1e9;
  const r = w.charAt(0), c = parseInt(w.slice(1), 10);
  return c * 100 + ROWS.indexOf(r);
}
const ordered = samples.slice().sort((a, b) => wellSortKey(a.well) - wellSortKey(b.well));
Plotly.newPlot("bar-plot", [{
  type: "bar",
  x: ordered.map(s => s.well ? `${s.well} (${s.id})` : s.id),
  y: ordered.map(s => s.defl * 100),
  marker: { color: ordered.map(s => COL_DARK[s.verdict]) },
  hovertemplate: "<b>%{x}</b><br>contamination: %{y:.2f}%<extra></extra>",
}], {
  margin: {t: 8, r: 12, b: 90, l: 50},
  xaxis: {tickangle: -55, automargin: true, fixedrange: true, gridcolor: "#233358", color: "#7e8fae", tickfont: {size: 10}},
  yaxis: {title: "contamination (%)", gridcolor: "#233358", zerolinecolor: "#233358", color: "#7e8fae"},
  paper_bgcolor: "#131f3a", plot_bgcolor: "#131f3a",
  font: {family: "ui-sans-serif", size: 11, color: "#e8eef9"},
  showlegend: false,
}, {displayModeBar: false, responsive: true});

function renderAllTable(){
  const fv = document.getElementById("all-verdict").value;
  const ft = document.getElementById("all-type").value;
  const q = document.getElementById("all-search").value.trim().toLowerCase();
  const body = document.getElementById("all-rows");
  body.innerHTML = ordered
    .filter(s => fv === "all" || s.verdict === fv)
    .filter(s => ft === "all" || s.type === ft)
    .filter(s => !q || (s.id||"").toLowerCase().includes(q) || (s.well||"").toLowerCase().includes(q))
    .map(s => {
      const w = Math.max(2, Math.round((s.defl / Math.max(MAX_C, 0.001)) * 86));
      const spark = `<svg class="spark" viewBox="0 0 90 16"><rect class="bar-bg" x="0" y="5" width="90" height="6" rx="2"/><rect class="bar-fg ${s.verdict.toLowerCase()}" x="0" y="5" width="${w}" height="6" rx="2"/></svg>`;
      return `<tr data-well="${s.well || ''}">
        <td><b>${s.well || '-'}</b></td>
        <td><b>${s.id}</b></td>
        <td><span class="verdict-pill ${s.verdict.toLowerCase()}">${s.verdict}</span></td>
        <td>${s.type}</td>
        <td class="num">${(s.defl * 100).toFixed(2)}%</td>
        <td class="num">${(s.score * 100).toFixed(2)}%</td>
        <td>${s.donor || '-'}</td>
        <td>${spark}</td>
      </tr>`;
    }).join("");
  body.querySelectorAll("tr").forEach(tr => {
    tr.addEventListener("click", () => {
      const well = tr.dataset.well;
      if (!well) return;
      if (isolatedWell === well){
        const state = document.getElementById("filter-btn").classList.contains("active") ? "flagged" : "all";
        applyFilter(state);
      } else { isolateRow(tr, well); }
    });
  });
}
document.getElementById("all-verdict").addEventListener("change", renderAllTable);
document.getElementById("all-type").addEventListener("change", renderAllTable);
document.getElementById("all-search").addEventListener("input", renderAllTable);
renderAllTable();
</script>
</body>
</html>
"""


def _default_well_assignment(n_samples: int) -> list:
    rows = list("ABCDEFGH")
    out = []
    for col in range(1, 13):
        for r in rows:
            if len(out) == n_samples:
                return out
            out.append(f"{r}{col}")
    return out


def build_samples(per_sample: pd.DataFrame) -> list:
    samples = []
    wells_col = per_sample["well"].fillna("").astype(str).tolist()
    if not any(w for w in wells_col):
        wells_col = _default_well_assignment(len(per_sample))
    for i, (sid, row) in enumerate(per_sample.iterrows()):
        defl = row.get("homalt_deflation")
        score = row.get("top_score_homalt")
        donor = row.get("top_donor_sample_id")
        samples.append({
            "id": str(sid),
            "well": (wells_col[i] if wells_col[i] else _default_well_assignment(len(per_sample))[i]),
            "type": str(row.get("sample_type") or "clinical"),
            "verdict": str(row.get("verdict") or "PASS"),
            "defl": float(defl) if defl is not None and not pd.isna(defl) else 0.0,
            "score": float(score) if score is not None and not pd.isna(score) else 0.0,
            "donor": "" if (donor is None or pd.isna(donor)) else str(donor),
        })
    return samples


def main():
    in_tsv = RESULTS_DIR / "per_sample_report.tsv"
    if not in_tsv.exists():
        raise FileNotFoundError(f"Expected {in_tsv}")

    per_sample = pd.read_csv(in_tsv, sep="\t").set_index("sample_id")
    if "well" not in per_sample.columns:
        per_sample["well"] = ""
    if "sample_type" not in per_sample.columns:
        per_sample["sample_type"] = "clinical"
    if "verdict" not in per_sample.columns:
        per_sample["verdict"] = "PASS"

    samples = build_samples(per_sample)

    types_present = sorted({s["type"] for s in samples})
    type_options = "".join(f'<option value="{t}">{t}</option>' for t in types_present)

    verdict_counts = per_sample["verdict"].value_counts().to_dict()
    meta_html = (
        f"{len(samples)} samples &middot; "
        f"<span style='color:#4ade80'>{verdict_counts.get('PASS', 0)} pass</span> &middot; "
        f"<span style='color:#fbbf24'>{verdict_counts.get('WARN', 0)} warn</span> &middot; "
        f"<span style='color:#f87171'>{verdict_counts.get('FAIL', 0)} fail</span>"
    )

    html = (HTML_TEMPLATE
            .replace("__LOGO__", LOGO_SVG.strip())
            .replace("__META_HTML__", meta_html)
            .replace("__TYPE_OPTIONS__", type_options)
            .replace("__DATA_JSON__", json.dumps({"samples": samples})))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(html, encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
