/*
 * app.js — Tropism Autodecoder web tool (Phase 1, 2, and 3 client-side).
 *
 * Pipeline: parse matrix -> map gene IDs (Arabidopsis or cross-species translation) ->
 * per-sample gene ranks -> singscore / NNLS deconvolution / cross-species explorer ->
 * composed RdBu SVG figure + auto legend + exports.
 * Nothing leaves the browser.
 */

const TOOL_VERSION = "2.0.0";
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
  if (window.PHASE3 && PHASE3.orthomap) {
    const sp = PHASE3.selectedSpecies || "human";
    if (PHASE3.orthomap.mappings && PHASE3.orthomap.mappings[sp] && PHASE3.orthomap.mappings[sp][t]) {
      return true;
    }
  }
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

  if (geneLikeInHeader > geneLikeInFirstCol && geneLikeInHeader >= 2) {
    // Transposed layout (samples in rows, genes in columns).
    const genes = header.slice(1).map(normalizeGeneId);
    const samples = body.map((r) => r[0]);
    const values = genes.map((_, gi) =>
      samples.map((_, si) => {
        const v = parseFloat(body[si][gi + 1]);
        return Number.isFinite(v) ? v : 0;
      })
    );
    return { genes, samples, values };
  }

  // Standard layout (genes in rows, samples in columns).
  const samples = header.slice(1);
  if (!samples.length) throw new Error("Header has no sample columns.");
  const genes = [];
  const values = [];
  for (let r = 0; r < body.length; r++) {
    const row = body[r];
    if (row.length < 2) continue;
    const g = normalizeGeneId(row[0]);
    if (!g) continue;
    genes.push(g);
    const rowVals = [];
    for (let c = 1; c <= samples.length; c++) {
      const v = parseFloat(row[c]);
      rowVals.push(Number.isFinite(v) ? v : 0);
    }
    values.push(rowVals);
  }
  if (!genes.length) throw new Error("No recognisable gene IDs found in the first column.");
  return { genes, samples, values };
}

function normalizeGeneId(raw) {
  if (!raw) return "";
  let s = raw.trim().replace(/\.\d+$/, ""); // strip .1 transcript suffix
  const up = s.toUpperCase();
  if (SYMBOL_TO_AGI[up]) return SYMBOL_TO_AGI[up];
  if (SYMBOL_TO_AGI[s]) return SYMBOL_TO_AGI[s];
  return up;
}

/* ------------------------------------------------------------------ *
 * Scoring: ranks + singscore
 * ------------------------------------------------------------------ */

// Return 1-based ranks with average ties for an array of values.
function averageRanks(arr) {
  const n = arr.length;
  const indices = new Array(n);
  for (let i = 0; i < n; i++) indices[i] = i;
  indices.sort((a, b) => arr[a] - arr[b]);

  const ranks = new Float64Array(n);
  let i = 0;
  while (i < n) {
    let j = i;
    while (j + 1 < n && arr[indices[j + 1]] === arr[indices[i]]) j++;
    const avgRank = (i + 1 + j + 1) / 2.0;
    for (let k = i; k <= j; k++) ranks[indices[k]] = avgRank;
    i = j + 1;
  }
  return ranks;
}

// singscore (Foroutan et al., 2018): mean rank of target set rescaled against
// the theoretical min/max mean rank for a set of size k, centered to [-1, +1].
function singscore(sampleRanks, matchedIndices, totalGenes) {
  const k = matchedIndices.length;
  if (!k) return 0;
  let sum = 0;
  for (let i = 0; i < k; i++) sum += sampleRanks[matchedIndices[i]];
  const rMean = sum / k;
  const minMean = (k + 1) / 2.0;
  const maxMean = totalGenes - (k - 1) / 2.0;
  if (maxMean <= minMean) return 0;
  const norm01 = (rMean - minMean) / (maxMean - minMean);
  return (norm01 - 0.5) * 2.0;
}

function computeScores() {
  const { genes, samples, values } = APP.matrix;
  const N = genes.length;
  const S = samples.length;

  const geneToIdx = new Map();
  for (let i = 0; i < N; i++) {
    if (!geneToIdx.has(genes[i])) geneToIdx.set(genes[i], i);
  }

  const sigIndices = TROPISM_SIGNATURES.map((sig) => {
    const hits = [];
    sig.genes.forEach((g) => {
      if (geneToIdx.has(g)) hits.push(geneToIdx.get(g));
    });
    return hits;
  });

  const sampleCols = [];
  for (let si = 0; si < S; si++) {
    const col = new Float64Array(N);
    for (let gi = 0; gi < N; gi++) col[gi] = values[gi][si];
    sampleCols.push(col);
  }

  const scores = TROPISM_SIGNATURES.map(() => new Float64Array(S));
  for (let si = 0; si < S; si++) {
    const ranks = averageRanks(sampleCols[si]);
    for (let gi = 0; gi < TROPISM_SIGNATURES.length; gi++) {
      scores[gi][si] = singscore(ranks, sigIndices[gi], N);
    }
  }

  const coverage = {};
  TROPISM_SIGNATURES.forEach((sig, i) => {
    coverage[sig.key] = { matched: sigIndices[i].length, total: sig.genes.length };
  });

  const calls = [];
  for (let si = 0; si < S; si++) {
    let best = -Infinity, bestKey = "none", runnerUp = -Infinity;
    for (let gi = 0; gi < TROPISM_SIGNATURES.length; gi++) {
      const v = scores[gi][si];
      if (v > best) {
        runnerUp = best;
        best = v;
        bestKey = TROPISM_SIGNATURES[gi].name;
      } else if (v > runnerUp) {
        runnerUp = v;
      }
    }
    calls.push({ call: bestKey, score: best, gap: best - runnerUp });
  }

  APP.scores = {
    samples,
    sigKeys: TROPISM_SIGNATURES.map((s) => s.key),
    values: scores,
    coverage,
    calls,
  };
}

/* ------------------------------------------------------------------ *
 * Composed figure rendering (SVG, RdBu palette)
 * ------------------------------------------------------------------ */

function rdbuColor(val) {
  const v = Math.max(-1, Math.min(1, val));
  if (v === 0) return RDBU.white;
  if (v > 0) {
    if (v < 0.33) return RDBU.redPale;
    if (v < 0.66) return RDBU.redLight;
    if (v < 0.9) return RDBU.red;
    return RDBU.redDark;
  }
  const u = -v;
  if (u < 0.33) return RDBU.bluePale;
  if (u < 0.66) return RDBU.blueLight;
  if (u < 0.9) return RDBU.blue;
  return RDBU.blueDark;
}

function contrastFg(val) {
  const v = Math.abs(val);
  return v > 0.6 ? "#ffffff" : "#000000";
}

function buildFigureSVG(figNumber) {
  const { samples, values, coverage, calls } = APP.scores;
  const nSamp = samples.length;
  const nSig = TROPISM_SIGNATURES.length;

  const rowH = 26;
  const colW = 120;
  const labelW = 180;
  const pad = 20;

  const panelAW = labelW + nSig * colW + 40;
  const panelAH = 70 + nSamp * rowH;

  const barH = 22;
  const barW = 200;
  const panelBW = 360;
  const panelBH = 70 + nSig * 32;

  const rightPanelW = Math.max(panelBW, 340);
  const totalW = panelAW + rightPanelW + pad * 3;
  const totalH = Math.max(panelAH, panelBH + 120) + 90;

  let out = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${totalW} ${totalH}" width="100%" height="auto" style="background:#ffffff;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">`;
  out += `<rect width="${totalW}" height="${totalH}" fill="#ffffff"/>`;

  // Title
  out += `<text x="${pad}" y="28" font-size="16" font-weight="700" fill="#000000">Figure ${figNumber}. Arabidopsis tropism signature decoding across samples</text>`;
  out += `<text x="${pad}" y="46" font-size="12" fill="#555555">Single-sample singscore enrichment against literature-verified tropism marker sets (scale [-1, +1]). RdBu diverging palette.</text>`;

  // PANEL A: Heatmap
  const ax = pad, ay = 64;
  out += `<g transform="translate(${ax},${ay})">`;
  out += `<text x="0" y="0" font-size="13" font-weight="700" fill="#000000">(A) Single-sample tropism scores</text>`;

  // Column headers
  for (let gi = 0; gi < nSig; gi++) {
    const x = labelW + gi * colW + colW / 2;
    out += `<text x="${x}" y="22" font-size="12" font-weight="600" text-anchor="middle" fill="#000000">${TROPISM_SIGNATURES[gi].name}</text>`;
  }

  // Rows
  for (let si = 0; si < nSamp; si++) {
    const y = 30 + si * rowH;
    const name = samples[si].length > 22 ? samples[si].slice(0, 20) + "…" : samples[si];
    out += `<text x="${labelW - 10}" y="${y + rowH * 0.7}" font-size="11.5" text-anchor="end" fill="#000000">${name}</text>`;

    for (let gi = 0; gi < nSig; gi++) {
      const v = values[gi][si];
      const fill = rdbuColor(v);
      const fg = contrastFg(v);
      const cellX = labelW + gi * colW;
      out += `<rect x="${cellX}" y="${y}" width="${colW - 2}" height="${rowH - 2}" fill="${fill}" rx="2" stroke="#eeeeee" stroke-width="0.5"/>`;
      out += `<text x="${cellX + colW / 2}" y="${y + rowH * 0.68}" font-size="11" font-weight="600" text-anchor="middle" fill="${fg}">${v.toFixed(3)}</text>`;
    }
  }
  out += `</g>`;

  // PANEL B: Mean scores
  const bx = ax + panelAW + pad, by = ay;
  out += `<g transform="translate(${bx},${by})">`;
  out += `<text x="0" y="0" font-size="13" font-weight="700" fill="#000000">(B) Mean signature activity</text>`;
  const zeroX = 140 + barW / 2;

  for (let gi = 0; gi < nSig; gi++) {
    const y = 30 + gi * 32;
    let sum = 0;
    for (let si = 0; si < nSamp; si++) sum += values[gi][si];
    const mean = sum / nSamp;

    out += `<text x="130" y="${y + 15}" font-size="12" text-anchor="end" fill="#000000">${TROPISM_SIGNATURES[gi].name}</text>`;
    out += `<line x1="${140}" y1="${y + 11}" x2="${140 + barW}" y2="${y + 11}" stroke="#cccccc" stroke-width="1"/>`;
    out += `<line x1="${zeroX}" y1="${y}" x2="${zeroX}" y2="${y + barH}" stroke="#888888" stroke-width="1"/>`;

    const span = (mean * (barW / 2));
    const bxPos = mean >= 0 ? zeroX : zeroX + span;
    const bW = Math.max(1, Math.abs(span));
    const fill = rdbuColor(mean);
    out += `<rect x="${bxPos}" y="${y + 2}" width="${bW}" height="${barH - 4}" fill="${fill}" rx="2"/>`;
    out += `<text x="${mean >= 0 ? bxPos + bW + 6 : bxPos - 6}" y="${y + 15}" font-size="11" font-weight="600" text-anchor="${mean >= 0 ? "start" : "end"}" fill="#000000">${mean.toFixed(3)}</text>`;
  }

  // PANEL C: Dominant calls strip + coverage
  const cy = by + panelBH + 10;
  out += `<g transform="translate(0,${cy})">`;
  out += `<text x="0" y="0" font-size="13" font-weight="700" fill="#000000">(C) Calls &amp; marker coverage</text>`;

  let covY = 18;
  TROPISM_SIGNATURES.forEach((sig) => {
    const c = coverage[sig.key];
    const flag = c.matched < c.total ? ` (${c.matched}/${c.total} genes)` : ` (${c.total}/${c.total})`;
    const warn = c.matched < 3 ? " ⚠ low coverage" : "";
    out += `<text x="10" y="${covY}" font-size="11" fill="#444444">• ${sig.name}: ${flag}${warn}</text>`;
    covY += 16;
  });

  covY += 4;
  out += `<text x="10" y="${covY}" font-size="11" font-weight="600" fill="#000000">Dominant sample calls:</text>`;
  covY += 14;
  const maxCallsShown = Math.min(6, nSamp);
  for (let si = 0; si < maxCallsShown; si++) {
    const sName = samples[si].length > 14 ? samples[si].slice(0, 12) + "…" : samples[si];
    out += `<text x="18" y="${covY}" font-size="10.5" fill="#333333">${sName}: <tspan font-weight="600">${calls[si].call}</tspan> (gap ${calls[si].gap.toFixed(3)})</text>`;
    covY += 14;
  }
  if (nSamp > maxCallsShown) {
    out += `<text x="18" y="${covY}" font-size="10" fill="#666666">… and ${nSamp - maxCallsShown} more (see full CSV export)</text>`;
  }
  out += `</g>`;
  out += `</g>`;

  out += `</svg>`;
  return out;
}

/* ------------------------------------------------------------------ *
 * Auto-generated figure legend (Methods-style)
 * ------------------------------------------------------------------ */

function buildLegend(figNumber) {
  const { samples, coverage } = APP.scores;
  const dateStr = document.getElementById("run-date")?.textContent || new Date().toISOString().slice(0, 10);
  const nSamp = samples.length;
  const nGenes = APP.matrix.genes.length;

  const covSummary = TROPISM_SIGNATURES.map((sig) => {
    const c = coverage[sig.key];
    return `${sig.name} (${sig.members}; ${c.matched}/${c.total} genes matched)`;
  }).join("; ");

  return `Figure ${figNumber}. Deconvolution of Arabidopsis thaliana tropism signatures from bulk expression data.

(A) Sample-by-tropism signature enrichment heatmap computed from ${nSamp} sample(s) across ${nGenes.toLocaleString()} genes in the uploaded matrix. Within each sample, genes are rank-ordered by expression and scored using singscore (Foroutan et al., BMC Bioinformatics 19:404, 2018), centered to [-1, +1]. Positive values (red) denote coordinate enrichment of the signature gene set relative to the sample background; negative values (blue) denote low expression.

(B) Mean signature score per tropism across all samples, centered at zero.

(C) Signature marker coverage and dominant tropism calls. Signature definitions (AGI loci verified against TAIR/primary literature): ${covSummary}.

Generated by the Tropism Autodecoder web tool (v${TOOL_VERSION}) on ${dateStr}. Everything computed client-side; data remained in the browser. Note: scores reflect relative marker enrichment within each sample, not absolute flux.`;
}

function scoresToCSV() {
  const { samples, values, calls } = APP.scores;
  const lines = [];
  const header = ["sample", ...TROPISM_SIGNATURES.map((s) => s.key), "dominant_call", "confidence_gap"];
  lines.push(header.join(","));
  samples.forEach((s, si) => {
    const row = [csv(s)];
    TROPISM_SIGNATURES.forEach((sig, gi) => row.push(fmt(values[gi][si])));
    row.push(csv(calls[si].call), fmt(calls[si].gap));
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

function download(fileName, text, mime) {
  const blob = new Blob([text], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = fileName;
  document.body.appendChild(a); a.click(); a.remove();
  URL.revokeObjectURL(url);
}

/* ------------------------------------------------------------------ *
 * UI wiring
 * ------------------------------------------------------------------ */

function setStatus(msg, kind) {
  const el = document.getElementById("status");
  if (!el) return;
  el.textContent = msg;
  el.className = "status " + (kind || "");
}

function runAnalysis(text, name) {
  try {
    APP.fileName = name;
    let analysisText = text;
    let translationStats = null;

    if (window.PHASE3 && PHASE3.selectedSpecies && PHASE3.selectedSpecies !== "arabidopsis") {
      const trans = PHASE3.translate(text, PHASE3.selectedSpecies);
      if (trans.isCrossSpecies && trans.stats) {
        analysisText = trans.translatedText;
        translationStats = trans.stats;
      }
    }

    APP.matrix = parseMatrix(analysisText);
    computeScores();
    renderResults();

    const totalMatched = Object.values(APP.scores.coverage).reduce((a, c) => a + c.matched, 0);
    if (translationStats) {
      setStatus(`[${translationStats.fromSpecies.toUpperCase()}] Mapped ${translationStats.mappedCount}/${translationStats.totalGenes} genes (${translationStats.mappingRate}%) to ${translationStats.uniqueAgiMapped} Arabidopsis loci. Analysed ${APP.matrix.samples.length} sample(s).`, "ok");
    } else if (totalMatched === 0) {
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
  const method = currentMethod();
  let svg, legend, csv;

  if (method === "phase3" && window.PHASE3 && PHASE3.ready) {
    const p3Sec = document.getElementById("phase3");
    if (p3Sec) p3Sec.scrollIntoView({ behavior: "smooth", block: "start" });
    svg = buildFigureSVG(figNumber);
    legend = buildLegend(figNumber) + "\n\nNote: Phase 3 Cross-Species Functional Conservation is actively displayed in Section 05 below.";
    csv = scoresToCSV();
  } else if (method === "phase2" && window.PHASE2 && PHASE2.ready) {
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
  if (method !== "phase3") {
    document.getElementById("results").scrollIntoView({ behavior: "smooth", block: "start" });
  }
}

function init() {
  const dateEl = document.getElementById("run-date");
  if (dateEl) dateEl.textContent = new Date().toISOString().slice(0, 10);

  const fileInput = document.getElementById("file");
  const drop = document.getElementById("drop");

  if (fileInput) {
    fileInput.addEventListener("change", (e) => {
      const f = e.target.files[0];
      if (f) f.text().then((t) => runAnalysis(t, f.name));
    });
  }

  if (drop) {
    ["dragover", "dragenter"].forEach((ev) => drop.addEventListener(ev, (e) => {
      e.preventDefault(); drop.classList.add("drag");
    }));
    ["dragleave", "drop"].forEach((ev) => drop.addEventListener(ev, (e) => {
      e.preventDefault(); drop.classList.remove("drag");
    }));
    drop.addEventListener("drop", (e) => {
      e.preventDefault(); drop.classList.remove("drag");
      const f = e.dataTransfer.files[0];
      if (f) f.text().then((t) => runAnalysis(t, f.name));
    });
  }

  // Example button listener (Arabidopsis default)
  const exBtn = document.getElementById("example-btn");
  if (exBtn) {
    exBtn.addEventListener("click", () => {
      setStatus("Loading synthetic example…", "");
      fetch("assets/sample_data.csv").then((r) => r.text())
        .then((t) => runAnalysis(t, "sample_data.csv (synthetic demo)"))
        .catch(() => setStatus("Could not load the bundled example.", "err"));
    });
  }

  // Cross-species demo buttons
  document.querySelectorAll(".demo-species-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      const sp = btn.getAttribute("data-species");
      if (window.PHASE3 && PHASE3.demos && PHASE3.demos[sp]) {
        if (window.PHASE3) PHASE3.selectedSpecies = sp;
        const spRadio = document.querySelector(`input[name="input-species"][value="${sp}"]`);
        if (spRadio) spRadio.checked = true;
        setStatus(`Loading ${sp.toUpperCase()} spaceflight demo…`, "");
        runAnalysis(PHASE3.demos[sp], `${sp}_flight_demo.csv`);
      }
    });
  });

  // Species selector radios
  document.querySelectorAll('input[name="input-species"]').forEach((r) => {
    r.addEventListener("change", (e) => {
      if (window.PHASE3) PHASE3.selectedSpecies = e.target.value;
      const label = e.target.value.toUpperCase();
      setStatus(`Input species set to: ${label}. Upload a matrix or select a demo.`, "");
    });
  });

  const figInput = document.getElementById("fig-number");
  if (figInput) {
    figInput.addEventListener("input", () => {
      if (APP.scores) renderResults();
    });
  }

  const dlSvg = document.getElementById("dl-svg");
  if (dlSvg) dlSvg.addEventListener("click", () => {
    if (APP.active) download("tropism_figure.svg", APP.active.svg, "image/svg+xml");
  });

  const dlPng = document.getElementById("dl-png");
  if (dlPng) dlPng.addEventListener("click", () => {
    if (!APP.active) return;
    svgToPng(APP.active.svg, 2, (blob) => {
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a"); a.href = url; a.download = "tropism_figure.png";
      document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(url);
    });
  });

  const dlCsv = document.getElementById("dl-csv");
  if (dlCsv) dlCsv.addEventListener("click", () => {
    if (APP.active) download("tropism_scores.csv", APP.active.csv, "text/csv");
  });

  const dlLegend = document.getElementById("dl-legend");
  if (dlLegend) dlLegend.addEventListener("click", () => {
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
