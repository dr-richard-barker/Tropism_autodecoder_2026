/*
 * model.js — Phase 2 consumer for the Tropism Autodecoder web tool.
 *
 * Reproduces the manuscript pipeline's INFERENCE in the browser, using artifacts
 * exported by Code/export_web_artifacts.py into assets/model/:
 *   signature_matrix.bin (+ signature_index.json), stimulus_codes.json, classifier.json.
 *
 * Per uploaded sample:
 *   1. NNLS(signature_matrix, bulk) -> cell-type fractions          (as project_bulk.py)
 *   2. fractions · stimulus_codes  -> 32 stimulus activation scores  (as project_bulk.py)
 *   3. StandardScaler + elastic-net logistic on [183 fractions + 32 stim]
 *      -> P(Space Flight)                                            (as meta_classifier.py)
 *
 * Phase 2 stays disabled until a complete artifact bundle is present (manifest.phase2_ready).
 */
(function () {
  const PHASE2 = { ready: false, model: null, _cacheKey: null, _cache: null };
  window.PHASE2 = PHASE2;

  const MODEL_DIR = "assets/model";

  // ---- dense linear solve (Gaussian elimination, partial pivot) ---------- //
  function solveDense(A, b, k) {
    const M = Float64Array.from(A), y = Float64Array.from(b);
    for (let col = 0; col < k; col++) {
      let piv = col, mx = Math.abs(M[col * k + col]);
      for (let r = col + 1; r < k; r++) {
        const v = Math.abs(M[r * k + col]); if (v > mx) { mx = v; piv = r; }
      }
      if (mx < 1e-12) return null;
      if (piv !== col) {
        for (let c = 0; c < k; c++) { const t = M[col * k + c]; M[col * k + c] = M[piv * k + c]; M[piv * k + c] = t; }
        const t = y[col]; y[col] = y[piv]; y[piv] = t;
      }
      const d = M[col * k + col];
      for (let r = col + 1; r < k; r++) {
        const f = M[r * k + col] / d;
        if (f !== 0) { for (let c = col; c < k; c++) M[r * k + c] -= f * M[col * k + c]; y[r] -= f * y[col]; }
      }
    }
    const x = new Float64Array(k);
    for (let r = k - 1; r >= 0; r--) {
      let acc = y[r]; for (let c = r + 1; c < k; c++) acc -= M[r * k + c] * x[c];
      x[r] = acc / M[r * k + r];
    }
    return x;
  }

  // ---- NNLS (Lawson-Hanson, Gram form): min ||Sw - b||, w >= 0 ----------- //
  // AtA = S^T S (n x n), Atb = S^T b (n). Returns w (length n).
  function nnlsGram(AtA, Atb, n) {
    const x = new Float64Array(n), P = new Uint8Array(n), w = new Float64Array(n);
    const tol = 1e-9, outerMax = 3 * n;
    for (let i = 0; i < n; i++) w[i] = Atb[i];
    let iter = 0;
    while (iter++ < outerMax) {
      let j = -1, best = tol;
      for (let i = 0; i < n; i++) if (!P[i] && w[i] > best) { best = w[i]; j = i; }
      if (j < 0) break;
      P[j] = 1;
      let guard = 0;
      while (guard++ < outerMax) {
        const idx = []; for (let i = 0; i < n; i++) if (P[i]) idx.push(i);
        const k = idx.length;
        const sub = new Float64Array(k * k), rhs = new Float64Array(k);
        for (let a = 0; a < k; a++) {
          rhs[a] = Atb[idx[a]];
          for (let b2 = 0; b2 < k; b2++) sub[a * k + b2] = AtA[idx[a] * n + idx[b2]];
        }
        const s = solveDense(sub, rhs, k);
        if (!s) { P[j] = 0; x[j] = 0; break; }
        let minS = Infinity; for (let a = 0; a < k; a++) if (s[a] < minS) minS = s[a];
        if (minS > tol) { for (let a = 0; a < k; a++) x[idx[a]] = s[a]; break; }
        let alpha = Infinity;
        for (let a = 0; a < k; a++) if (s[a] <= tol) {
          const xi = x[idx[a]], denom = xi - s[a];
          if (denom > 0) { const r = xi / denom; if (r < alpha) alpha = r; }
        }
        if (!isFinite(alpha)) alpha = 0;
        for (let a = 0; a < k; a++) { const i = idx[a]; x[i] += alpha * (s[a] - x[i]); }
        for (let a = 0; a < k; a++) { const i = idx[a]; if (x[i] <= tol) { P[i] = 0; x[i] = 0; } }
      }
      for (let i = 0; i < n; i++) { let acc = Atb[i]; const row = i * n; for (let t = 0; t < n; t++) acc -= AtA[row + t] * x[t]; w[i] = acc; }
    }
    return x;
  }

  // ---------------------------- load artifacts ---------------------------- //
  PHASE2.load = async function () {
    try {
      const mres = await fetch(`${MODEL_DIR}/manifest.json`, { cache: "no-cache" });
      if (!mres.ok) return false;
      const manifest = await mres.json();
      if (!manifest.phase2_ready) { PHASE2.manifest = manifest; return false; }

      const [idx, clf, stim] = await Promise.all([
        fetch(`${MODEL_DIR}/signature_index.json`).then((r) => r.json()),
        fetch(`${MODEL_DIR}/classifier.json`).then((r) => r.json()),
        fetch(`${MODEL_DIR}/stimulus_codes.json`).then((r) => r.json()),
      ]);
      const buf = await fetch(`${MODEL_DIR}/signature_matrix.bin`).then((r) => r.arrayBuffer());
      const S = new Float32Array(buf);                 // row-major genes x clusters
      const [nGenes, nClusters] = idx.shape;
      if (S.length !== nGenes * nClusters) throw new Error("signature_matrix size mismatch");

      // signature gene -> row (normalized AGI)
      const geneRow = new Map();
      idx.genes.forEach((g, i) => { const k = norm(g); if (!geneRow.has(k)) geneRow.set(k, i); });

      PHASE2.model = { S, nGenes, nClusters, clusters: idx.clusters, geneRow, clf, stim, manifest };
      PHASE2.ready = true;
      return true;
    } catch (e) {
      console.warn("Phase 2 load failed:", e);
      return false;
    }
  };

  function norm(g) {
    // Prefer the app's normaliser; fall back to a local AGI cleanup.
    if (typeof normalizeGeneId === "function") return normalizeGeneId(g);
    return String(g || "").trim().toUpperCase().replace(/\.\d+$/, "");
  }

  // ------------------------------- scoring -------------------------------- //
  PHASE2.run = function (matrix) {
    if (PHASE2._cache && PHASE2._cacheKey === matrix) return PHASE2._cache;
    const M = PHASE2.model;
    const { samples, genes, values } = matrix;

    // rows shared between the uploaded matrix and the signature genes
    const pairs = []; // [sigRow, matrixRow]
    genes.forEach((g, mi) => { const r = M.geneRow.get(norm(g)); if (r !== undefined) pairs.push([r, mi]); });
    const nCommon = pairs.length;

    // decide log1p (bulk must be on the signature's log-normalized scale)
    let mx = 0;
    for (let si = 0; si < samples.length; si++) for (let p = 0; p < nCommon; p++) {
      const v = values[pairs[p][1]][si]; if (Number.isFinite(v) && v > mx) mx = v;
    }
    const doLog = mx > 20;

    // S_common (nCommon x nClusters) and its Gram matrix AtA (nClusters x nClusters)
    const nC = M.nClusters;
    const Sc = new Float64Array(nCommon * nC);
    for (let p = 0; p < nCommon; p++) {
      const srow = pairs[p][0] * nC;
      for (let c = 0; c < nC; c++) Sc[p * nC + c] = M.S[srow + c];
    }
    const AtA = new Float64Array(nC * nC);
    for (let p = 0; p < nCommon; p++) {
      const base = p * nC;
      for (let a = 0; a < nC; a++) {
        const va = Sc[base + a]; if (va === 0) continue;
        const rowA = a * nC;
        for (let b = a; b < nC; b++) AtA[rowA + b] += va * Sc[base + b];
      }
    }
    for (let a = 0; a < nC; a++) for (let b = a + 1; b < nC; b++) AtA[b * nC + a] = AtA[a * nC + b];
    let dmean = 0; for (let a = 0; a < nC; a++) dmean += AtA[a * nC + a];
    const ridge = (dmean / nC) * 1e-8 + 1e-12;
    for (let a = 0; a < nC; a++) AtA[a * nC + a] += ridge;

    // classifier feature layout
    const clf = M.clf, feat = clf.feature_names;
    const clusterIndex = new Map(); M.clusters.forEach((c, i) => clusterIndex.set(c, i));
    const stimDimOf = (name) => { const m = /^stim_dim_(\d+)$/.exec(name); return m ? +m[1] : -1; };
    const latent = M.stim[Object.keys(M.stim)[0]] ? M.stim[Object.keys(M.stim)[0]].length : 32;

    const out = { samples: [], fractions: [], stim: [], flightProb: [], nCommon,
                  doLog, clusters: M.clusters, auc: clf.cv_auc_mean, version: (M.manifest || {}).version };

    for (let si = 0; si < samples.length; si++) {
      const Atb = new Float64Array(nC);
      for (let p = 0; p < nCommon; p++) {
        let v = values[pairs[p][1]][si]; if (!Number.isFinite(v)) v = 0;
        if (doLog) v = Math.log1p(v);
        const base = p * nC;
        for (let c = 0; c < nC; c++) Atb[c] += Sc[base + c] * v;
      }
      const w = nnlsGram(AtA, Atb, nC);
      let sum = 0; for (let c = 0; c < nC; c++) sum += w[c];
      const frac = new Float64Array(nC); if (sum > 0) for (let c = 0; c < nC; c++) frac[c] = w[c] / sum;

      // stimulus activation = fractions · stimulus_codes (aligned by cluster name)
      const stim = new Float64Array(latent);
      for (let c = 0; c < nC; c++) {
        const code = M.stim[M.clusters[c]]; if (!code || frac[c] === 0) continue;
        for (let d = 0; d < latent; d++) stim[d] += frac[c] * code[d];
      }

      // feature vector in classifier order -> scale -> logistic
      let logit = clf.intercept;
      for (let k = 0; k < feat.length; k++) {
        const name = feat[k];
        let raw;
        const d = stimDimOf(name);
        if (d >= 0) raw = stim[d] || 0;
        else { const ci = clusterIndex.get(name); raw = ci === undefined ? 0 : frac[ci]; }
        const scaled = (raw - clf.scaler_mean[k]) / (clf.scaler_scale[k] || 1);
        logit += clf.coef[k] * scaled;
      }
      const prob = 1 / (1 + Math.exp(-logit)); // P(positive_class = classes_[1], "Space Flight")

      out.samples.push(samples[si]);
      out.fractions.push(frac); out.stim.push(stim); out.flightProb.push(prob);
    }

    PHASE2._cacheKey = matrix; PHASE2._cache = out;
    return out;
  };

  // ------------------------------ rendering ------------------------------- //
  function rdbu(v) { // v in [-1,1], high -> red (reuse app.js if available)
    if (typeof rdbuColor === "function") return rdbuColor(v);
    const x = Math.max(-1, Math.min(1, v));
    const r = x > 0 ? 178 : Math.round(247 + x * (33 - 247) * -1);
    return `rgb(${r},${Math.round(247 - Math.abs(x) * 200)},${x < 0 ? 172 : 43})`;
  }
  const esc = (s) => String(s).replace(/[&<>]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[c]));

  PHASE2.buildFigureSVG = function (res, figNumber) {
    const nS = res.samples.length;
    const W = 900, padL = 150;
    // top clusters by max fraction across samples
    const nC = res.clusters.length, maxFrac = new Float64Array(nC);
    res.fractions.forEach((f) => { for (let c = 0; c < nC; c++) if (f[c] > maxFrac[c]) maxFrac[c] = f[c]; });
    const topIdx = Array.from({ length: nC }, (_, i) => i).sort((a, b) => maxFrac[b] - maxFrac[a])
      .filter((i) => maxFrac[i] > 0).slice(0, 12);

    const rowH = Math.max(20, Math.min(38, 300 / Math.max(nS, 1)));
    const pAtop = 90, pAh = nS * rowH;
    const pBtop = pAtop + pAh + 70, cell = 46, pBh = nS * rowH;
    const stimTop = pBtop + pBh + 70, sc = 22, stimH = nS * rowH;
    const H = stimTop + stimH + 110;

    let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" font-family="Helvetica, Arial, sans-serif" role="img" aria-label="Auto-decoder deconvolution figure">`;
    s += `<rect width="${W}" height="${H}" fill="#fff"/>`;
    s += `<text x="24" y="34" font-size="20" font-weight="700" fill="#000">Figure ${esc(figNumber)}. Auto-decoder deconvolution &amp; spaceflight classification</text>`;
    s += `<text x="24" y="56" font-size="12" fill="#000">Arabidopsis thaliana — ${nS} sample${nS === 1 ? "" : "s"}, ${res.nCommon} signature genes matched${res.doLog ? " (log1p applied)" : ""}. NNLS deconvolution onto the 183-cluster atlas signature; elastic-net classifier.</text>`;

    // Panel A: Flight probability bars
    s += `<text x="24" y="${pAtop - 16}" font-size="14" font-weight="700" fill="#000">A  P(Space Flight) — elastic-net classifier${res.auc ? ` (nested-CV AUC ${res.auc.toFixed(3)})` : ""}</text>`;
    const barX = padL, barW = 480, halfX = barX + barW * 0.5;
    res.samples.forEach((name, si) => {
      const y = pAtop + si * rowH, p = res.flightProb[si];
      s += `<text x="${barX - 8}" y="${y + rowH / 2 + 4}" font-size="11" fill="#000" text-anchor="end">${esc(name.length > 20 ? name.slice(0, 19) + "…" : name)}</text>`;
      s += `<rect x="${barX}" y="${y + 3}" width="${barW}" height="${rowH - 8}" fill="#f2f2f2" stroke="#ccc" stroke-width="0.5"/>`;
      s += `<rect x="${barX}" y="${y + 3}" width="${barW * p}" height="${rowH - 8}" fill="${rdbu(-(2 * p - 1))}" stroke="#000" stroke-width="0.4"/>`;
      s += `<text x="${barX + barW + 8}" y="${y + rowH / 2 + 4}" font-size="11" fill="#000">${p.toFixed(3)}</text>`;
    });
    s += `<line x1="${halfX}" y1="${pAtop - 2}" x2="${halfX}" y2="${pAtop + pAh + 2}" stroke="#000" stroke-dasharray="3 3" stroke-width="0.7"/>`;
    s += `<text x="${halfX}" y="${pAtop + pAh + 14}" font-size="9" fill="#000" text-anchor="middle">0.5</text>`;

    // Panel B: top cell-type fractions heatmap
    s += `<text x="24" y="${pBtop - 16}" font-size="14" font-weight="700" fill="#000">B  Top cell-type fractions (NNLS)</text>`;
    topIdx.forEach((ci, k) => {
      const x = padL + k * cell + cell / 2;
      s += `<text x="${x}" y="${pBtop - 4}" font-size="9" fill="#000" text-anchor="middle" transform="rotate(-30 ${x} ${pBtop - 4})">${esc(res.clusters[ci])}</text>`;
    });
    res.samples.forEach((name, si) => {
      const y = pBtop + si * rowH;
      s += `<text x="${padL - 8}" y="${y + rowH / 2 + 4}" font-size="10" fill="#000" text-anchor="end">${esc(name.length > 18 ? name.slice(0, 17) + "…" : name)}</text>`;
      topIdx.forEach((ci, k) => {
        const v = res.fractions[si][ci], x = padL + k * cell;
        const g = Math.round(255 - Math.min(1, v / (maxFrac[ci] || 1)) * 210);
        s += `<rect x="${x}" y="${y}" width="${cell - 2}" height="${rowH - 2}" fill="rgb(${g},${g},255)" stroke="#000" stroke-width="0.3"/>`;
        s += `<text x="${x + (cell - 2) / 2}" y="${y + rowH / 2 + 3}" font-size="8.5" fill="${v > 0.25 ? "#fff" : "#000"}" text-anchor="middle">${(v * 100).toFixed(0)}</text>`;
      });
    });
    s += `<text x="${padL}" y="${pBtop + pBh + 14}" font-size="9" fill="#000">values = % of sample</text>`;

    // Panel C: stimulus activation heatmap (samples x latent), per-column z-scaled
    const latent = res.stim[0].length;
    s += `<text x="24" y="${stimTop - 16}" font-size="14" font-weight="700" fill="#000">C  Stimulus activation (${latent}-dim auto-decoder codes)</text>`;
    const colStats = [];
    for (let d = 0; d < latent; d++) {
      let m = 0; res.stim.forEach((v) => m += v[d]); m /= nS;
      let sd = 0; res.stim.forEach((v) => sd += (v[d] - m) ** 2); sd = Math.sqrt(sd / Math.max(1, nS)) || 1;
      colStats.push([m, sd]);
    }
    res.samples.forEach((name, si) => {
      const y = stimTop + si * rowH;
      s += `<text x="${padL - 8}" y="${y + rowH / 2 + 4}" font-size="10" fill="#000" text-anchor="end">${esc(name.length > 18 ? name.slice(0, 17) + "…" : name)}</text>`;
      for (let d = 0; d < latent; d++) {
        const z = Math.max(-1, Math.min(1, (res.stim[si][d] - colStats[d][0]) / colStats[d][1] / 2.5));
        s += `<rect x="${padL + d * sc}" y="${y}" width="${sc - 1}" height="${rowH - 2}" fill="${rdbu(z)}" stroke="#fff" stroke-width="0.2"/>`;
      }
    });
    s += `<text x="${padL}" y="${stimTop + stimH + 16}" font-size="9" fill="#000">stim_dim_0 … stim_dim_${latent - 1} · RdBu, per-column z-scaled (red high)</text>`;

    s += `<text x="24" y="${H - 12}" font-size="9" fill="#000">Phase 2 — reproduces the manuscript auto-decoder pipeline in-browser · model ${esc(res.version || "artifacts")}</text>`;
    s += `</svg>`;
    return s;
  };

  PHASE2.buildLegend = function (res, figNumber) {
    const pc = PHASE2.model.clf.positive_class || "Space Flight";
    return [
      `**Figure ${figNumber}. Auto-decoder deconvolution and spaceflight classification of ${res.samples.length} Arabidopsis thaliana sample(s).**`,
      `(**A**) Probability of ${pc} from the elastic-net logistic classifier${res.auc ? ` (nested 5-fold CV AUC = ${res.auc.toFixed(3)})` : ""}, computed on 183 cell-type fractions + ${res.stim[0].length} stimulus-activation scores. (**B**) Top cell-type fractions from non-negative least-squares (NNLS) deconvolution of each sample onto the 183-cluster Salk atlas signature matrix. (**C**) Stimulus activation = fractions weighted by the auto-decoder's per-cluster latent codes; shown per-column z-scaled (RdBu).`,
      `Method: this reproduces the manuscript pipeline's inference in the browser — NNLS onto the signature matrix (as in project_bulk.py), stimulus-code projection, then StandardScaler + elastic-net logistic regression (as in meta_classifier.py). ${res.nCommon} of the model's signature genes were detected in the uploaded matrix${res.doLog ? "; values were log1p-transformed to match the log-normalized reference" : ""}.`,
      `Note: for faithful probabilities, provide normalized bulk expression (e.g. DESeq2-normalized or CPM/TPM); NNLS fidelity depends on signature-gene coverage. Model artifacts: ${res.version || "exported bundle"}.`,
    ].join("\n\n");
  };

  PHASE2.scoresCSV = function (res) {
    const latent = res.stim[0].length, nC = res.clusters.length;
    const head = ["sample", "flight_probability", "dominant_cell_type", "dominant_fraction"];
    for (let d = 0; d < latent; d++) head.push(`stim_dim_${d}`);
    const lines = [head.join(",")];
    res.samples.forEach((name, si) => {
      const f = res.fractions[si]; let bi = 0; for (let c = 1; c < nC; c++) if (f[c] > f[bi]) bi = c;
      const row = [csv(name), res.flightProb[si].toFixed(4), csv(res.clusters[bi]), f[bi].toFixed(4)];
      for (let d = 0; d < latent; d++) row.push(res.stim[si][d].toFixed(4));
      lines.push(row.join(","));
    });
    return lines.join("\n");
  };
  function csv(s) { return /[",\n]/.test(s) ? `"${String(s).replace(/"/g, '""')}"` : s; }

  // ------------------------------ UI init --------------------------------- //
  function updateUI() {
    const radio = document.getElementById("method-phase2");
    const status = document.getElementById("method-status");
    if (PHASE2.ready) {
      if (radio) radio.disabled = false;
      if (status) status.textContent = "Full auto-decoder (Phase 2) available.";
    } else {
      if (radio) { radio.disabled = true; radio.checked = false; }
      const p1 = document.getElementById("method-phase1"); if (p1) p1.checked = true;
      if (status) status.textContent = "Full auto-decoder (Phase 2) unavailable — model artifacts not yet exported (run Code/export_web_artifacts.py).";
    }
  }

  if (typeof document !== "undefined") {
    document.addEventListener("DOMContentLoaded", () => {
      PHASE2.load().then(updateUI);
    });
  }

  if (typeof module !== "undefined" && module.exports) {
    module.exports = { nnlsGram, solveDense, PHASE2 };
  }
})();
