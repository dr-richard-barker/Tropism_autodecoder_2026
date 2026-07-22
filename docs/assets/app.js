/*
 * app.js — Tropism Autodecoder web tool (Phase 1, client-side).
 *
 * Pipeline: parse matrix -> map gene IDs -> per-sample gene ranks ->
 * singscore per tropism signature -> composed RdBu SVG figure + auto legend + exports.
 * Nothing leaves the browser.
 */

const TOOL_VERSION = "1.0.0";
const APP = {
  matrix: null,      // { genes:[], samples:[], values:[[]] }  values[gene][sample]
  scores: null,      // { samples:[], sigKeys:[], values:[[]], coverage:{}, calls:[] }
  fileName: null,
};

/* ------------------------------------------------------------------ *
 * Parsing
 * ------------------------------------------------------------------ */

function detectDelimiter(text) {
  const line = text.split(/\r?\n/).find((l) => l.trim().length) || "";
  const counts = { "\t": (line.match(/\t/g) || []).length,
                   ",": (line.match(/,/g) || []).length,
                   ";": (line.match(/;/g) || []).length };
  return Object.keys(counts).sort((a, b) => counts[b] - counts[a])[0] || ",";
}

function looksLikeGeneId(s) {
  if (!s) return false;
  const t = s.trim().toUpperCase().replace(/\.\d+$/, "");
  if (/^AT[1-5MC]G\d{5}$/.test(t)) return true;      // AGI locus
  if (Object.prototype.hasOwnProperty.call(SYMBOL_TO_AGI, t)) return true;
  return false;
}

// Parse a delimited expression table into { genes, samples, values }.
function parseMatrix(text) {
  const delim = detectDelimiter(text);
  const rows = text.split(/\r?\n/).filter((l) => l.trim().length)
                   .map((l) => l.split(delim).map((c) => c.trim().replace(/^"|"$/g, "")));
  if (rows.length < 2) throw new Error("File has too few rows to be an expression matrix.");

  const header = rows[0];
  const body = rows.slice(1);

  // Orientation A (default): first column = gene IDs, header[1..] = sample names.
  const firstColIds = body.map((r) => r[0]);
  const geneLikeInFirstCol = firstColIds.filter(looksLikeGeneId).length;
  // Orientation B: first row (after cell 0) = gene IDs, first column = sample names.
  const geneLikeInHeader = header.slice(1).filter(looksLikeGeneId).length;

  let genes, samples, values;
  if (geneLikeInHeader > geneLikeInFirstCol) {
    // samples x genes -> transpose to genes x samples
    genes = header.slice(1);
    samples = body.map((r) => r[0]);
    values = genes.map((_, gi) => samples.map((_, si) => toNum(body[si][gi + 1])));
  } else {
    genes = firstColIds;
    samples = header.slice(1);
    values = body.map((r) => samples.map((_, si) => toNum(r[si + 1])));
  }

  if (!samples.length) throw new Error("No sample columns detected.");
  if (!genes.length) throw new Error("No gene rows detected.");
  return { genes, samples, values };
}

function toNum(v) {
  const n = parseFloat(v);
  return Number.isFinite(n) ? n : NaN;
}

function normalizeGeneId(raw) {
  const t = String(raw || "").trim().toUpperCase().replace(/\.\d+$/, "");
  if (/^AT[1-5MC]G\d{5}$/.test(t)) return t;
  if (Object.prototype.hasOwnProperty.call(SYMBOL_TO_AGI, t)) return SYMBOL_TO_AGI[t];
  return t; // leave unknown as-is; it simply won't match a signature
}

/* ------------------------------------------------------------------ *
 * Scoring — singscore (rank-based single-sample enrichment)
 * ------------------------------------------------------------------ */

// Average ranks (1..N) for a numeric vector; NaN treated as lowest.
function averageRanks(vec) {
  const idx = vec.map((v, i) => [Number.isFinite(v) ? v : -Infinity, i])
                 .sort((a, b) => a[0] - b[0]);
  const ranks = new Array(vec.length);
  let i = 0;
  while (i < idx.length) {
    let j = i;
    while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++;
    const avg = (i + j) / 2 + 1; // 1-based average rank for the tie block
    for (let k = i; k <= j; k++) ranks[idx[k][1]] = avg;
    i = j + 1;
  }
  return ranks;
}

// singscore for one signature (up-regulated set) in one sample.
// Returns centered score in [-1, 1]; positive = coordinately high expression.
function singscore(ranks, memberRows, N) {
  if (!memberRows.length) return { score: NaN, matched: 0 };
  const meanRank = memberRows.reduce((s, r) => s + ranks[r], 0) / memberRows.length;
  const m = memberRows.length;
  const minMean = (m + 1) / 2;               // lowest possible mean rank
  const maxMean = N - (m - 1) / 2;            // highest possible mean rank
  const norm = (meanRank - minMean) / (maxMean - minMean); // 0..1
  return { score: 2 * norm - 1, matched: m };
}

function computeScores() {
  const { genes, samples, values } = APP.matrix;
  const N = genes.length;

  // Map each signature to the row indices present in the uploaded matrix.
  const geneRow = new Map();
  genes.forEach((g, i) => {
    const norm = normalizeGeneId(g);
    if (!geneRow.has(norm)) geneRow.set(norm, i);
  });

  const sigRows = TROPISM_SIGNATURES.map((sig) => {
    const rows = [], matchedGenes = [];
    sig.genes.forEach((agi) => {
      if (geneRow.has(agi)) { rows.push(geneRow.get(agi)); matchedGenes.push(agi); }
    });
    return { rows, matchedGenes };
  });

  // Per-sample ranks, then per-signature scores.
  const scoreMatrix = TROPISM_SIGNATURES.map(() => []);
  samples.forEach((_, si) => {
    const col = values.map((row) => row[si]);
    const ranks = averageRanks(col);
    TROPISM_SIGNATURES.forEach((sig, gi) => {
      scoreMatrix[gi].push(singscore(ranks, sigRows[gi].rows, N).score);
    });
  });

  const coverage = {};
  TROPISM_SIGNATURES.forEach((sig, gi) => {
    coverage[sig.key] = { matched: sigRows[gi].rows.length, total: sig.genes.length,
                          matchedGenes: sigRows[gi].matchedGenes };
  });

  // Dominant tropism call per sample.
  const calls = samples.map((s, si) => {
    const perSig = TROPISM_SIGNATURES.map((sig, gi) => ({ key: sig.key, name: sig.name,
                                                          v: scoreMatrix[gi][si] }))
                                     .filter((x) => Number.isFinite(x.v))
                                     .sort((a, b) => b.v - a.v);
    if (!perSig.length) return { sample: s, call: "n/a", gap: 0 };
    const gap = perSig.length > 1 ? perSig[0].v - perSig[1].v : perSig[0].v;
    return { sample: s, call: perSig[0].name, top: perSig[0].v, gap };
  });

  APP.scores = { samples, values: scoreMatrix, coverage, calls,
                 sigKeys: TROPISM_SIGNATURES.map((s) => s.key), genesTotal: N };
  return APP.scores;
}

/* ------------------------------------------------------------------ *
 * Color scale (RdBu diverging over [-1, 1])
 * ------------------------------------------------------------------ */

const RDBU_STOPS = [
  [-1.0, [103, 0, 31]],   [-0.6, [178, 24, 43]], [-0.3, [239, 138, 98]],
  [-0.1, [253, 219, 199]],[0.0, [247, 247, 247]],[0.1, [209, 229, 240]],
  [0.3, [103, 169, 207]], [0.6, [33, 102, 172]], [1.0, [5, 48, 97]],
];
// Note: RdBu convention here maps HIGH score -> BLUE, LOW -> RED per ColorBrewer RdBu.
// We flip so HIGH signature activity reads as warm RED (intuitive "on"): invert index.
function rdbuColor(v) {
  const x = Math.max(-1, Math.min(1, Number.isFinite(v) ? -v : 0)); // flip: high -> red
  for (let i = 1; i < RDBU_STOPS.length; i++) {
    if (x <= RDBU_STOPS[i][0]) {
      const [x0, c0] = RDBU_STOPS[i - 1], [x1, c1] = RDBU_STOPS[i];
      const t = (x - x0) / (x1 - x0 || 1);
      const c = c0.map((ch, k) => Math.round(ch + t * (c1[k] - ch)));
      return `rgb(${c[0]},${c[1]},${c[2]})`;
    }
  }
  return "rgb(5,48,97)";
}
function textOn(v) { return Math.abs(v) > 0.55 ? "#ffffff" : "#000000"; }

/* ------------------------------------------------------------------ *
 * Figure — one composed SVG (Panels A, B, C)
 * ------------------------------------------------------------------ */

function buildFigureSVG(figNumber) {
  const S = APP.scores;
  const sigs = TROPISM_SIGNATURES;
  const nS = S.samples.length, nSig = sigs.length;

  const W = 900;
  const padL = 150, padT = 90;
  const cell = 46, cellH = Math.max(20, Math.min(40, 320 / Math.max(nS, 1)));
  const heatW = nSig * cell, heatH = nS * cellH;
  const barTop = padT + heatH + 58;
  const barH = 140, barBase = barTop + barH;
  const H = barBase + 140;

  const esc = (s) => String(s).replace(/[&<>]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[c]));
  let svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" font-family="Helvetica, Arial, sans-serif" role="img" aria-label="Tropism signature decoding figure">`;
  svg += `<rect width="${W}" height="${H}" fill="#ffffff"/>`;
  svg += `<text x="24" y="34" font-size="20" font-weight="700" fill="#000">Figure ${figNumber}. Tropism signature decoding</text>`;
  svg += `<text x="24" y="56" font-size="12" fill="#000">Arabidopsis thaliana — ${nS} sample${nS === 1 ? "" : "s"}, ${S.genesTotal.toLocaleString()} genes. Rank-based single-sample enrichment (singscore).</text>`;

  // ---- Panel A: heatmap ----
  svg += `<text x="24" y="${padT - 18}" font-size="14" font-weight="700" fill="#000">A  Signature score heatmap</text>`;
  sigs.forEach((sig, gi) => {
    const x = padL + gi * cell + cell / 2;
    svg += `<text x="${x}" y="${padT - 4}" font-size="11" fill="#000" text-anchor="middle" transform="rotate(-25 ${x} ${padT - 4})">${esc(sig.name)}</text>`;
  });
  S.samples.forEach((s, si) => {
    const y = padT + si * cellH;
    svg += `<text x="${padL - 8}" y="${y + cellH / 2 + 4}" font-size="11" fill="#000" text-anchor="end">${esc(s.length > 20 ? s.slice(0, 19) + "…" : s)}</text>`;
    sigs.forEach((sig, gi) => {
      const v = S.values[gi][si];
      const x = padL + gi * cell;
      const fill = Number.isFinite(v) ? rdbuColor(v) : "#e6e6e6";
      svg += `<rect x="${x}" y="${y}" width="${cell - 2}" height="${cellH - 2}" fill="${fill}" stroke="#000" stroke-width="0.4"/>`;
      const label = Number.isFinite(v) ? v.toFixed(2) : "NA";
      svg += `<text x="${x + (cell - 2) / 2}" y="${y + cellH / 2 + 3.5}" font-size="9.5" fill="${Number.isFinite(v) ? textOn(v) : "#000"}" text-anchor="middle">${label}</text>`;
    });
  });

  // colorbar
  const cbX = padL + heatW + 40, cbY = padT, cbH = Math.min(heatH, 200), cbW = 16;
  for (let i = 0; i < cbH; i++) {
    const v = 1 - (2 * i) / cbH; // +1 top -> -1 bottom
    svg += `<rect x="${cbX}" y="${cbY + i}" width="${cbW}" height="1" fill="${rdbuColor(v)}"/>`;
  }
  svg += `<rect x="${cbX}" y="${cbY}" width="${cbW}" height="${cbH}" fill="none" stroke="#000" stroke-width="0.5"/>`;
  svg += `<text x="${cbX + cbW + 6}" y="${cbY + 8}" font-size="10" fill="#000">+1 high</text>`;
  svg += `<text x="${cbX + cbW + 6}" y="${cbY + cbH / 2 + 4}" font-size="10" fill="#000">0</text>`;
  svg += `<text x="${cbX + cbW + 6}" y="${cbY + cbH}" font-size="10" fill="#000">−1 low</text>`;

  // ---- Panel B: mean score bar chart ----
  svg += `<text x="24" y="${barTop - 18}" font-size="14" font-weight="700" fill="#000">B  Mean signature score across samples</text>`;
  const means = sigs.map((sig, gi) => {
    const vals = S.values[gi].filter(Number.isFinite);
    return vals.length ? vals.reduce((a, b) => a + b, 0) / vals.length : NaN;
  });
  const bw = 90, gap = 40, bx0 = padL - 30;
  const zeroY = barTop + barH / 2;
  svg += `<line x1="${bx0 - 10}" y1="${zeroY}" x2="${bx0 + nSig * (bw + gap)}" y2="${zeroY}" stroke="#000" stroke-width="0.7"/>`;
  svg += `<text x="${bx0 - 16}" y="${zeroY + 3}" font-size="9" fill="#000" text-anchor="end">0</text>`;
  sigs.forEach((sig, gi) => {
    const v = means[gi];
    const x = bx0 + gi * (bw + gap);
    if (Number.isFinite(v)) {
      const h = (Math.abs(v) * barH) / 2;
      const y = v >= 0 ? zeroY - h : zeroY;
      svg += `<rect x="${x}" y="${y}" width="${bw}" height="${h}" fill="${rdbuColor(v)}" stroke="#000" stroke-width="0.5"/>`;
      svg += `<text x="${x + bw / 2}" y="${v >= 0 ? y - 4 : y + h + 12}" font-size="11" fill="#000" text-anchor="middle">${v.toFixed(2)}</text>`;
    }
    svg += `<text x="${x + bw / 2}" y="${barBase + 16}" font-size="11" fill="#000" text-anchor="middle">${esc(sig.name)}</text>`;
  });

  // ---- Panel C: coverage strip ----
  const cy = barBase + 54;
  svg += `<text x="24" y="${cy - 12}" font-size="14" font-weight="700" fill="#000">C  Signature coverage (genes matched in your data)</text>`;
  sigs.forEach((sig, gi) => {
    const cov = S.coverage[sig.key];
    const x = 24 + gi * 220;
    const flag = sig.lowCoverage ? "  ⚠ low-coverage set" : "";
    svg += `<text x="${x}" y="${cy + 10}" font-size="11" fill="#000" font-weight="700">${esc(sig.name)}</text>`;
    svg += `<text x="${x}" y="${cy + 26}" font-size="10" fill="#000">${cov.matched}/${cov.total} genes${flag}</text>`;
  });

  svg += `<text x="24" y="${H - 12}" font-size="9" fill="#000">Tropism Autodecoder web tool v${TOOL_VERSION} · RdBu scale: red = high, white ≈ 0, blue = low · rank-based markers, not the auto-decoder.</text>`;
  svg += `</svg>`;
  return svg;
}

/* ------------------------------------------------------------------ *
 * Figure legend (auto-generated Markdown/plain text)
 * ------------------------------------------------------------------ */

function buildLegend(figNumber) {
  const S = APP.scores;
  const cov = S.coverage;
  const covStr = TROPISM_SIGNATURES.map((s) =>
    `${s.name} (${cov[s.key].matched}/${cov[s.key].total}: ${s.members})`).join("; ");
  const date = document.getElementById("run-date")?.textContent || "";
  const lowSets = TROPISM_SIGNATURES.filter((s) => s.lowCoverage).map((s) => s.name);

  return [
`**Figure ${figNumber}. Tropism signature decoding of ${S.samples.length} Arabidopsis thaliana transcriptome${S.samples.length === 1 ? "" : "s"}.**`,
`(**A**) Heatmap of rank-based single-sample enrichment scores (singscore; Foroutan et al., 2018) for four curated tropism signature gene sets, computed independently for each of the ${S.samples.length} sample column${S.samples.length === 1 ? "" : "s"} across ${S.genesTotal.toLocaleString()} genes. Scores are centered to [−1, +1]; positive (red) indicates coordinately high expression of the signature, negative (blue) low. (**B**) Mean signature score across all samples. (**C**) Number of signature genes detected in the uploaded matrix.`,
`Signatures: ${covStr}. Marker loci are canonical tropism regulators (AGI identifiers verified against TAIR/primary literature).`,
`Method: expression values were rank-transformed per sample (ties averaged); singscore rescales the mean rank of each signature's detected genes against the theoretical rank range for a set of that size. The method is scale-invariant, so raw counts, TPM/FPKM, or normalized values may be used without pre-processing.`,
lowSets.length ? `Note: ${lowSets.join(", ")} is a deliberately small, low-coverage signature (mirroring the limited hydrotropism transcriptomics available); interpret its score with caution.` : ``,
`Caveat: these are curated marker-enrichment scores, not the manuscript's variational auto-decoder / deconvolution output. They report relative tropism-signature activity within the uploaded samples. Color scale: ColorBrewer RdBu diverging.`,
`Generated by the Tropism Autodecoder web tool v${TOOL_VERSION}${date ? " on " + date : ""}.`,
  ].filter(Boolean).join("\n\n");
}

/* ------------------------------------------------------------------ *
 * Exports
 * ------------------------------------------------------------------ */

function download(filename, content, type) {
  const blob = new Blob([content], { type });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = filename; document.body.appendChild(a); a.click();
  a.remove(); URL.revokeObjectURL(url);
}

function scoresToCSV() {
  const S = APP.scores;
  const head = ["sample", ...TROPISM_SIGNATURES.map((s) => s.key), "dominant_call", "confidence_gap"];
  const lines = [head.join(",")];
  S.samples.forEach((s, si) => {
    const row = [csv(s)];
    TROPISM_SIGNATURES.forEach((sig, gi) => row.push(fmt(S.values[gi][si])));
    row.push(csv(S.calls[si].call), fmt(S.calls[si].gap));
    lines.push(row.join(","));
  });
  return lines.join("\n");
}
function csv(s) { return /[",\n]/.test(s) ? `"${String(s).replace(/"/g, '""')}"` : s; }
function fmt(v) { return Number.isFinite(v) ? v.toFixed(4) : "NA"; }

function svgToPng(svgText, scale, cb) {
  const img = new Image();
  const svg64 = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svgText)));
  img.onload = function () {
    const vb = svgText.match(/viewBox="0 0 (\d+) (\d+)"/);
    const w = vb ? +vb[1] : 900, h = vb ? +vb[2] : 700;
    const canvas = document.createElement("canvas");
    canvas.width = w * scale; canvas.height = h * scale;
    const ctx = canvas.getContext("2d");
    ctx.fillStyle = "#ffffff"; ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
    canvas.toBlob(cb, "image/png");
  };
  img.src = svg64;
}

/* ------------------------------------------------------------------ *
 * UI wiring
 * ------------------------------------------------------------------ */

function setStatus(msg, kind) {
  const el = document.getElementById("status");
  el.textContent = msg;
  el.className = "status " + (kind || "");
}

function runAnalysis(text, name) {
  try {
    APP.fileName = name;
    APP.matrix = parseMatrix(text);
    computeScores();
    renderResults();
    const totalMatched = Object.values(APP.scores.coverage).reduce((a, c) => a + c.matched, 0);
    if (totalMatched === 0) {
      setStatus(`Parsed ${APP.matrix.samples.length} samples × ${APP.matrix.genes.length} genes, but no signature genes matched. Check that gene IDs are AGI loci (e.g. AT1G70940) or recognised symbols.`, "warn");
    } else {
      setStatus(`Analysed ${APP.matrix.samples.length} sample(s) × ${APP.matrix.genes.length.toLocaleString()} genes · ${totalMatched} signature genes matched.`, "ok");
    }
  } catch (e) {
    setStatus("Could not analyse this file: " + e.message, "err");
  }
}

function currentMethod() {
  const el = document.querySelector('input[name="method"]:checked');
  return el ? el.value : "phase1";
}

function renderResults() {
  const figNumber = document.getElementById("fig-number").value || "1";
  let svg, legend, csv;
  if (currentMethod() === "phase2" && window.PHASE2 && PHASE2.ready) {
    const res = PHASE2.run(APP.matrix);          // real deconvolution + classifier
    svg = PHASE2.buildFigureSVG(res, figNumber);
    legend = PHASE2.buildLegend(res, figNumber);
    csv = PHASE2.scoresCSV(res);
  } else {
    svg = buildFigureSVG(figNumber);             // Phase 1 marker enrichment
    legend = buildLegend(figNumber);
    csv = scoresToCSV();
  }
  APP.active = { svg, legend, csv };
  document.getElementById("figure").innerHTML = svg;
  document.getElementById("legend").textContent = legend;
  document.getElementById("results").hidden = false;
  document.getElementById("results").scrollIntoView({ behavior: "smooth", block: "start" });
}

function init() {
  const dateEl = document.getElementById("run-date");
  if (dateEl) dateEl.textContent = new Date().toISOString().slice(0, 10);

  const fileInput = document.getElementById("file");
  const drop = document.getElementById("drop");

  fileInput.addEventListener("change", (e) => {
    const f = e.target.files[0];
    if (f) f.text().then((t) => runAnalysis(t, f.name));
  });

  ["dragover", "dragenter"].forEach((ev) => drop.addEventListener(ev, (e) => {
    e.preventDefault(); drop.classList.add("drag");
  }));
  ["dragleave", "drop"].forEach((ev) => drop.addEventListener(ev, (e) => {
    e.preventDefault(); drop.classList.remove("drag");
  }));
  drop.addEventListener("drop", (e) => {
    const f = e.dataTransfer.files[0];
    if (f) f.text().then((t) => runAnalysis(t, f.name));
  });

  document.getElementById("example-btn").addEventListener("click", () => {
    setStatus("Loading synthetic example…", "");
    fetch("assets/sample_data.csv").then((r) => r.text())
      .then((t) => runAnalysis(t, "sample_data.csv (synthetic demo)"))
      .catch(() => setStatus("Could not load the bundled example.", "err"));
  });

  document.getElementById("fig-number").addEventListener("input", () => {
    if (APP.scores) renderResults();
  });

  document.getElementById("dl-svg").addEventListener("click", () => {
    if (APP.active) download("tropism_figure.svg", APP.active.svg, "image/svg+xml");
  });
  document.getElementById("dl-png").addEventListener("click", () => {
    if (!APP.active) return;
    svgToPng(APP.active.svg, 2, (blob) => {
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a"); a.href = url; a.download = "tropism_figure.png";
      document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(url);
    });
  });
  document.getElementById("dl-csv").addEventListener("click", () => {
    if (APP.active) download("tropism_scores.csv", APP.active.csv, "text/csv");
  });
  document.getElementById("dl-legend").addEventListener("click", () => {
    if (APP.active) download("tropism_legend.md", APP.active.legend, "text/markdown");
  });

  // Re-render when the analysis method toggle changes.
  document.querySelectorAll('input[name="method"]').forEach((r) =>
    r.addEventListener("change", () => { if (APP.matrix) renderResults(); }));
}

if (typeof document !== "undefined") {
  document.addEventListener("DOMContentLoaded", init);
}

// Exposed for headless testing (Node). Harmless in the browser.
if (typeof module !== "undefined" && module.exports) {
  module.exports = { parseMatrix, computeScores, singscore, averageRanks,
                     normalizeGeneId, detectDelimiter, buildLegend, APP };
}
