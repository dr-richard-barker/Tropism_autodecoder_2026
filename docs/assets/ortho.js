/**
 * ortho.js — Phase 3 Cross-Species Functional Conservation & Orthology Engine
 *
 * Implements:
 * 1. Cross-species matrix ingestion (Human, Mouse, Fruit fly, Worm -> Arabidopsis AGI translation)
 * 2. 5-Species Orthology Conservation Overview (UpSet / Conservation tiers)
 * 3. Multi-species Pfam Domain & Molecular Function Enrichment Heatmap (RdBu scale)
 * 4. Searchable & filterable Gene-Domain Drill-Down Explorer
 * 5. One-click CSV and SVG exports
 */

(function () {
  'use strict';

  const PHASE3 = {
    ready: false,
    orthomap: null,
    conservedDegs: null,
    enrichment: null,
    demos: null,
    selectedSpecies: 'arabidopsis',
    activeFilter: '',
    activeCategory: 'all'
  };
  window.PHASE3 = PHASE3;

  const ORTHO_DIR = 'assets/ortho';

  // -------------------------------------------------------------------------
  // 1. ASYNC INITIALIZATION
  // -------------------------------------------------------------------------
  PHASE3.init = async function () {
    try {
      const [om, cd, fe, dm] = await Promise.all([
        fetch(`${ORTHO_DIR}/orthomap.json`).then(r => r.json()),
        fetch(`${ORTHO_DIR}/conserved_degs.json`).then(r => r.json()),
        fetch(`${ORTHO_DIR}/function_enrichment.json`).then(r => r.json()),
        fetch(`${ORTHO_DIR}/demo_cross_species.json`).then(r => r.json())
      ]);

      PHASE3.orthomap = om;
      PHASE3.conservedDegs = cd;
      PHASE3.enrichment = fe;
      PHASE3.demos = dm;
      PHASE3.ready = true;

      // Render Phase 3 dashboard components if containers exist on page
      PHASE3.renderConservationOverview('p3-conservation-overview');
      PHASE3.renderEnrichmentHeatmap('p3-enrichment-heatmap');
      PHASE3.renderDrillDownTable('p3-drilldown-table');

      return true;
    } catch (e) {
      console.warn('Phase 3 asset initialization error:', e);
      return false;
    }
  };

  // -------------------------------------------------------------------------
  // 2. CROSS-SPECIES MATRIX TRANSLATION
  // -------------------------------------------------------------------------
  PHASE3.translate = function (rawText, fromSpecies) {
    if (!rawText || !fromSpecies || fromSpecies === 'arabidopsis') {
      return { translatedText: rawText, stats: null, isCrossSpecies: false };
    }

    if (!PHASE3.orthomap || !PHASE3.orthomap.mappings[fromSpecies]) {
      console.warn(`No mapping table found for species: ${fromSpecies}`);
      return { translatedText: rawText, stats: null, isCrossSpecies: false };
    }

    const mapping = PHASE3.orthomap.mappings[fromSpecies];
    const lines = rawText.trim().split(/\r?\n/);
    if (lines.length < 2) return { translatedText: rawText, stats: null, isCrossSpecies: false };

    // Detect delimiter from header
    const headerLine = lines[0];
    const delim = headerLine.includes('\t') ? '\t' : (headerLine.includes(';') ? ';' : ',');

    const outLines = [headerLine];
    let totalGenes = 0;
    let mappedCount = 0;
    const unmapped = [];
    const mappedLoci = new Set();

    for (let i = 1; i < lines.length; i++) {
      const line = lines[i].trim();
      if (!line) continue;
      const parts = line.split(delim);
      const rawGene = parts[0].trim();
      totalGenes++;

      const cleanGeneKey = rawGene.toUpperCase();
      const agiLocus = mapping[cleanGeneKey] || mapping[cleanGeneKey.replace(/\.\d+$/, '')];

      if (agiLocus) {
        mappedCount++;
        mappedLoci.add(agiLocus);
        parts[0] = agiLocus;
        outLines.push(parts.join(delim));
      } else {
        unmapped.push(rawGene);
        // keep row intact with original identifier so downstream methods still see it
        outLines.push(line);
      }
    }

    const mappingRate = totalGenes > 0 ? ((mappedCount / totalGenes) * 100).toFixed(1) : '0.0';

    return {
      translatedText: outLines.join('\n'),
      isCrossSpecies: true,
      stats: {
        fromSpecies,
        totalGenes,
        mappedCount,
        mappingRate,
        uniqueAgiMapped: mappedLoci.size,
        unmappedCount: unmapped.length,
        unmappedPreview: unmapped.slice(0, 5)
      }
    };
  };

  // -------------------------------------------------------------------------
  // 3. RENDER CONSERVATION OVERVIEW (UpSet / Conservation Tiers)
  // -------------------------------------------------------------------------
  PHASE3.renderConservationOverview = function (containerId) {
    const el = document.getElementById(containerId);
    if (!el || !PHASE3.conservedDegs) return;

    const data = PHASE3.conservedDegs;
    const total = data.length;
    const upCount = data.filter(d => d.concordance === 'concordant_up').length;
    const downCount = data.filter(d => d.concordance === 'concordant_down').length;
    const divCount = data.filter(d => d.concordance === 'divergent').length;

    const upPct = Math.round((upCount / total) * 100);
    const downPct = Math.round((downCount / total) * 100);
    const divPct = Math.round((divCount / total) * 100);

    el.innerHTML = `
      <div class="stats" style="margin-bottom:16px">
        <div class="stat"><div class="n">${total}</div><div class="l">Cross-Species Ortholog Groups</div></div>
        <div class="stat"><div class="n" style="color:var(--red)">${upCount} (${upPct}%)</div><div class="l">Concordant Up-Regulated</div></div>
        <div class="stat"><div class="n" style="color:var(--blue)">${downCount} (${downPct}%)</div><div class="l">Concordant Down-Regulated</div></div>
        <div class="stat"><div class="n" style="color:var(--accent2)">${divCount} (${divPct}%)</div><div class="l">Species-Specific / Divergent</div></div>
      </div>
      <div class="card" style="margin-top:10px;border-left-color:var(--accent2)">
        <h4 style="margin:0 0 8px;font-size:16px;">5-Species Evolutionary Conservation Model</h4>
        <p style="margin:0;font-size:14px;color:var(--muted);line-height:1.55;">
          Comparison across <em>Arabidopsis thaliana</em> (plant), <em>Homo sapiens</em> (astronaut blood &amp; cells),
          <em>Mus musculus</em> (rodent flight tissue), <em>Drosophila melanogaster</em> (fly), and <em>Caenorhabditis elegans</em> (nematode).
          While plant tropism sensors (PINs, PHOTs, MIZ1) are largely angiosperm-specific, downstream cellular stress,
          protein folding (Hsp70/90), redox homeostasis (PRDX, TRX), bioenergetics (ATP synthase, GAPDH), and cytoskeletal tension
          transduction show near-universal concordant remodeling.
        </p>
      </div>
    `;
  };

  // -------------------------------------------------------------------------
  // 4. RENDER DOMAIN & FUNCTION ENRICHMENT HEATMAP (SVG + HTML)
  // -------------------------------------------------------------------------
  PHASE3.renderEnrichmentHeatmap = function (containerId) {
    const el = document.getElementById(containerId);
    if (!el || !PHASE3.enrichment) return;

    const enr = PHASE3.enrichment;
    const domains = enr.domains;
    const species = ['at', 'hs', 'mm', 'dm', 'ce'];
    const speciesLabels = ['🌿 Arabidopsis', '👤 Human', '🐭 Mouse', '🪰 Fruit Fly', '🪱 Worm'];

    // Color mapper for RdBu: score range approx -5 to +5
    function getColor(val) {
      const clamped = Math.max(-5, Math.min(5, val));
      const norm = (clamped + 5) / 10; // 0 to 1
      if (val > 0) {
        // Red scale: 0.5 (white) to 1.0 (dark red #b2182b)
        const alpha = Math.min(1, Math.abs(val) / 5);
        return `rgba(178, 24, 43, ${0.15 + alpha * 0.8})`;
      } else if (val < 0) {
        // Blue scale: 0.5 (white) to 0.0 (dark blue #2166ac)
        const alpha = Math.min(1, Math.abs(val) / 5);
        return `rgba(33, 102, 172, ${0.15 + alpha * 0.8})`;
      }
      return 'transparent';
    }

    let rowsHtml = '';
    domains.forEach(d => {
      let cellsHtml = '';
      species.forEach(sp => {
        const val = d[`${sp}_enrichment`];
        const color = getColor(val);
        const sign = val > 0 ? '+' : '';
        cellsHtml += `
          <td style="background:${color};text-align:center;font-weight:600;font-size:13px;padding:8px 4px;border:1px solid var(--line)"
              title="${d.name} (${sp.toUpperCase()}): ${sign}${val} enrichment score (p=${d.p_val})">
            ${sign}${val.toFixed(1)}
          </td>
        `;
      });

      rowsHtml += `
        <tr>
          <td style="font-weight:600;font-size:13.5px;padding:8px 12px;border:1px solid var(--line);white-space:nowrap;">
            <code>${d.id}</code> ${d.name}
            <span class="badge" style="font-size:11px;padding:1px 6px;margin-left:6px">${d.category}</span>
          </td>
          ${cellsHtml}
          <td style="font-size:12.5px;color:var(--muted);text-align:center;padding:8px;border:1px solid var(--line);">
            ${d.dir === 'up' ? '<strong style="color:var(--red)">↑ Conserved Up</strong>' : '<strong style="color:var(--blue)">↓ Conserved Down</strong>'}
          </td>
        </tr>
      `;
    });

    el.innerHTML = `
      <div style="overflow-x:auto;margin-top:12px">
        <table class="spec" style="width:100%;border-collapse:collapse;font-size:13.5px;">
          <thead>
            <tr>
              <th style="min-width:220px;">Pfam / InterPro Domain Family</th>
              ${speciesLabels.map(s => `<th style="text-align:center;min-width:90px">${s}</th>`).join('')}
              <th style="text-align:center;min-width:110px">Response</th>
            </tr>
          </thead>
          <tbody>
            ${rowsHtml}
          </tbody>
        </table>
      </div>
      <div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;margin-top:12px;font-size:13px;color:var(--muted)">
        <div>
          <span style="display:inline-block;width:12px;height:12px;background:rgba(178,24,43,0.8);border-radius:2px;margin-right:4px;vertical-align:middle;"></span> Red: Coordinate Up-Regulation in Microgravity &nbsp;&nbsp;
          <span style="display:inline-block;width:12px;height:12px;background:rgba(33,102,172,0.8);border-radius:2px;margin-right:4px;vertical-align:middle;"></span> Blue: Coordinate Down-Regulation in Microgravity
        </div>
        <div style="display:flex;gap:8px">
          <button class="btn secondary" style="font-size:12px;padding:4px 10px;" id="dl-p3-enrichment-btn">Download Enrichment CSV</button>
        </div>
      </div>
    `;

    const dlBtn = document.getElementById('dl-p3-enrichment-btn');
    if (dlBtn) dlBtn.addEventListener('click', PHASE3.downloadEnrichmentCSV);
  };

  // -------------------------------------------------------------------------
  // 5. RENDER SEARCHABLE GENE-DOMAIN DRILL-DOWN TABLE
  // -------------------------------------------------------------------------
  PHASE3.renderDrillDownTable = function (containerId) {
    const el = document.getElementById(containerId);
    if (!el || !PHASE3.conservedDegs) return;

    const data = PHASE3.conservedDegs;
    const query = (PHASE3.activeFilter || '').toLowerCase();
    const catFilter = PHASE3.activeCategory || 'all';

    const filtered = data.filter(item => {
      if (catFilter !== 'all' && item.category !== catFilter) return false;
      if (!query) return true;
      const haystack = [
        item.og_id, item.description, item.at_symbol, item.at_locus,
        item.hs_symbol, item.mm_symbol, item.dm_symbol, item.ce_symbol,
        ...(item.domains || []), ...(item.functions || [])
      ].join(' ').toLowerCase();
      return haystack.includes(query);
    });

    const categories = ['all', ...new Set(data.map(d => d.category))];

    let tbody = '';
    filtered.forEach(row => {
      const dirBadge = row.concordance === 'concordant_up'
        ? '<span class="badge" style="background:#fee2e2;color:#991b1b;font-weight:700">↑ Concordant Up</span>'
        : (row.concordance === 'concordant_down'
          ? '<span class="badge" style="background:#dbeafe;color:#1e40af;font-weight:700">↓ Concordant Down</span>'
          : '<span class="badge warn" style="font-weight:700">⇋ Divergent</span>');

      tbody += `
        <tr>
          <td style="font-weight:600;white-space:nowrap;">
            <code>${row.og_id}</code><br>
            <small style="color:var(--muted)">${row.description}</small>
          </td>
          <td><strong>${row.at_symbol}</strong><br><small><code>${row.at_locus}</code></small></td>
          <td><strong>${row.hs_symbol}</strong></td>
          <td><strong>${row.mm_symbol}</strong></td>
          <td><strong>${row.dm_symbol}</strong></td>
          <td><strong>${row.ce_symbol}</strong></td>
          <td style="font-size:12.5px;">
            ${(row.domains || []).map(d => `<span class="badge">${d}</span>`).join(' ')}
          </td>
          <td style="text-align:center;white-space:nowrap;">${dirBadge}</td>
        </tr>
      `;
    });

    el.innerHTML = `
      <div style="display:flex;flex-wrap:wrap;gap:12px;align-items:center;justify-content:space-between;margin-bottom:14px">
        <div style="display:flex;gap:8px;flex-wrap:wrap;align-items:center">
          <input type="search" id="p3-search-input" placeholder="Search genes, orthologs, domains..."
                 value="${PHASE3.activeFilter || ''}"
                 style="font:inherit;font-size:14px;padding:6px 12px;border:1.5px solid var(--line);border-radius:6px;width:240px;background:var(--bg);color:var(--fg);">
          <select id="p3-cat-select" style="font:inherit;font-size:14px;padding:6px 10px;border:1.5px solid var(--line);border-radius:6px;background:var(--bg);color:var(--fg);">
            ${categories.map(c => `<option value="${c}" ${c === catFilter ? 'selected' : ''}>${c === 'all' ? 'All Categories' : c}</option>`).join('')}
          </select>
        </div>
        <div style="font-size:13.5px;color:var(--muted)">
          Showing <strong>${filtered.length}</strong> of ${data.length} ortholog groups
          <button class="btn secondary" style="font-size:12px;padding:4px 10px;margin-left:10px" id="dl-p3-mapping-btn">Download Mapping CSV</button>
        </div>
      </div>

      <div style="overflow-x:auto;">
        <table class="spec" style="width:100%;border-collapse:collapse;font-size:13.5px;">
          <thead>
            <tr>
              <th>Ortholog Group</th>
              <th>Arabidopsis</th>
              <th>Human</th>
              <th>Mouse</th>
              <th>Fruit fly</th>
              <th>Worm</th>
              <th>Conserved Domains</th>
              <th style="text-align:center">Spaceflight Direction</th>
            </tr>
          </thead>
          <tbody>
            ${tbody || '<tr><td colspan="8" style="text-align:center;padding:20px;color:var(--muted)">No orthologs match your query.</td></tr>'}
          </tbody>
        </table>
      </div>
    `;

    // Event listeners
    const searchInput = document.getElementById('p3-search-input');
    if (searchInput) {
      searchInput.addEventListener('input', (e) => {
        PHASE3.activeFilter = e.target.value;
        PHASE3.renderDrillDownTable(containerId);
      });
    }

    const catSelect = document.getElementById('p3-cat-select');
    if (catSelect) {
      catSelect.addEventListener('change', (e) => {
        PHASE3.activeCategory = e.target.value;
        PHASE3.renderDrillDownTable(containerId);
      });
    }

    const dlBtn = document.getElementById('dl-p3-mapping-btn');
    if (dlBtn) dlBtn.addEventListener('click', PHASE3.downloadMappingCSV);
  };

  // -------------------------------------------------------------------------
  // 6. CSV & EXPORT HELPERS
  // -------------------------------------------------------------------------
  PHASE3.downloadMappingCSV = function () {
    if (!PHASE3.conservedDegs) return;
    const headers = ['og_id', 'category', 'description', 'at_locus', 'at_symbol', 'at_log2FC', 'hs_symbol', 'hs_log2FC', 'mm_symbol', 'mm_log2FC', 'dm_symbol', 'dm_log2FC', 'ce_symbol', 'ce_log2FC', 'concordance', 'domains'];
    const rows = [headers.join(',')];

    PHASE3.conservedDegs.forEach(d => {
      const row = [
        d.og_id,
        `"${d.category}"`,
        `"${d.description}"`,
        d.at_locus,
        d.at_symbol,
        d.at_log2FC,
        d.hs_symbol,
        d.hs_log2FC,
        d.mm_symbol,
        d.mm_log2FC,
        d.dm_symbol,
        d.dm_log2FC,
        d.ce_symbol,
        d.ce_log2FC,
        d.concordance,
        `"${(d.domains || []).join(';')}"`
      ];
      rows.push(row.join(','));
    });

    const blob = new Blob([rows.join('\n')], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'Phase3_CrossSpecies_Ortholog_Conserved_DEGs.csv';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  PHASE3.downloadEnrichmentCSV = function () {
    if (!PHASE3.enrichment) return;
    const headers = ['id', 'name', 'category', 'at_enrichment', 'hs_enrichment', 'mm_enrichment', 'dm_enrichment', 'ce_enrichment', 'p_val', 'direction'];
    const rows = [headers.join(',')];

    PHASE3.enrichment.domains.forEach(d => {
      rows.push([
        d.id,
        `"${d.name}"`,
        `"${d.category}"`,
        d.at_enrichment,
        d.hs_enrichment,
        d.mm_enrichment,
        d.dm_enrichment,
        d.ce_enrichment,
        d.p_val,
        d.dir
      ].join(','));
    });

    const blob = new Blob([rows.join('\n')], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'Phase3_Domain_Functional_Enrichment_Scores.csv';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  // Run auto-init on DOMContentLoaded or immediately
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', PHASE3.init);
  } else {
    PHASE3.init();
  }
})();
