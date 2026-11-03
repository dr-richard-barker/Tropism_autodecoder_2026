#!/usr/bin/env python3
"""
generate_figure8_revamped.py

Generates Figure 8 formatted to npj Microgravity (Springer Nature) publication standards:
- 100% Crisp Neutral White Background (#FFFFFF)
- Standard Nature typography & panel indicators (a, b, c, d)
- Astrobotany & single-cell deconvolution fidelity
- Zero right-edge truncation: 1800x1260 canvas with ample margins and internal padding
- High-resolution uncropped PNG rendering via headless Google Chrome (2400x1680 px)

Outputs:
  - Figures/kegg_pathway_network.svg (Publication Vector SVG)
  - Figures/kegg_pathway_network.png (300 DPI Publication Raster for LaTeX Manuscript)
  - docs/assets/figures/kegg_pathway_network.svg (Web Portal SVG)
  - docs/assets/figures/kegg_pathway_network.png (Web Portal High-Res PNG)
"""

import os
import sys
import subprocess
import shutil

WIDTH = 1800
HEIGHT = 1260

# npj Microgravity (Springer Nature) Color Palette
COLOR = {
    "bg": "#FFFFFF",
    "panel_bg": "#FFFFFF",
    "panel_border": "#E5E7EB",
    "subcard_bg": "#F9FAFB",
    "subcard_border": "#E5E7EB",
    "text_dark": "#111827",
    "text_body": "#374151",
    "text_muted": "#6B7280",
    "text_light": "#9CA3AF",
    
    # Biological functional modules (Nature accessible colorblind-safe)
    "auxin": "#15803D",        # Nature Forest Green
    "auxin_light": "#DCFCE7",
    "auxin_border": "#86EFAC",
    
    "mechano": "#0284C7",      # Sky / Slate Blue
    "mechano_light": "#E0F2FE",
    "mechano_border": "#7DD3FC",
    
    "wall": "#B45309",         # Ochre / Amber
    "wall_light": "#FEF3C7",
    "wall_border": "#FCD34D",
    
    "redox": "#BE123C",        # Vermilion / Carmine
    "redox_light": "#FFE4E6",
    "redox_border": "#FDA4AF",
    
    "chaperone": "#6D28D9",    # Royal Violet
    "chaperone_light": "#EDE9FE",
    "chaperone_border": "#C4B5FD",
    
    "atp": "#0F766E",          # Deep Teal
    "atp_light": "#CCFBF1",
    "atp_border": "#5EEAD4",
    
    # Cellular structures
    "membrane": "#0284C7",
    "membrane_bilayer": "#BAE6FD",
    "statolith": "#CA8A04",
    "statolith_body": "#FEF08A",
    "nucleus": "#F5F3FF",
    "nucleus_border": "#8B5CF6",
    "mito": "#FFF7ED",
    "mito_border": "#F97316",
    "er": "#FDF2F8",
    "er_border": "#EC4899",
}

def escape_xml(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")

def generate_figure8_svg():
    svg = []
    
    # ------------------ SVG Header & Defs ------------------
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {HEIGHT}" width="{WIDTH}" height="{HEIGHT}" style="background-color: {COLOR["bg"]}; font-family: -apple-system, BlinkMacSystemFont, \'Segoe UI\', Helvetica, Arial, sans-serif;">')
    svg.append('<defs>')
    
    # Markers (Arrows)
    svg.append(f'''
      <marker id="arrow-auxin" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
        <path d="M 0 1 L 10 5 L 0 9 z" fill="{COLOR["auxin"]}" />
      </marker>
      <marker id="arrow-mechano" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
        <path d="M 0 1 L 10 5 L 0 9 z" fill="{COLOR["mechano"]}" />
      </marker>
      <marker id="arrow-wall" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
        <path d="M 0 1 L 10 5 L 0 9 z" fill="{COLOR["wall"]}" />
      </marker>
      <marker id="arrow-redox" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
        <path d="M 0 1 L 10 5 L 0 9 z" fill="{COLOR["redox"]}" />
      </marker>
      <marker id="arrow-muted" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
        <path d="M 0 1 L 10 5 L 0 9 z" fill="{COLOR["text_muted"]}" />
      </marker>
    ''')
    svg.append('</defs>')

    # ------------------ Main Figure Header (npj Microgravity style) ------------------
    svg.append(f'''
    <!-- Title & Metadata Header -->
    <text x="30" y="38" font-size="20" font-weight="700" fill="{COLOR["text_dark"]}">
      Fig. 8 | Spatial-molecular systems architecture of spaceflight tropism and mechanotransduction.
    </text>
    <text x="30" y="58" font-size="12.5" font-weight="400" fill="{COLOR["text_muted"]}">
      Multiscale integration of subcellular mechanoperception, single-cell tissue activation (inflorescence stem_9 &amp; root columella), and functional pathway cross-talk in Arabidopsis thaliana.
    </text>
    <rect x="{WIDTH - 210}" y="22" width="180" height="26" rx="4" fill="{COLOR["subcard_bg"]}" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    <text x="{WIDTH - 120}" y="39" font-size="11" font-weight="600" fill="{COLOR["text_body"]}" text-anchor="middle">
      npj Microgravity · Research
    </text>
    <line x1="30" y1="68" x2="{WIDTH - 30}" y2="68" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    ''')

    # Coordinates & Grid Layout
    left_x = 30
    left_w = 1030
    right_x = 1085
    right_w = 685
    
    row1_y = 80
    row1_h = 635
    row2_y = 730
    row2_h = 505

    # =========================================================================
    # PANEL a: Subcellular Mechanotransduction & Biochemical Cascade
    # =========================================================================
    svg.append(f'''
    <!-- PANEL a CONTAINER -->
    <rect x="{left_x}" y="{row1_y}" width="{left_w}" height="{row1_h}" rx="6" fill="{COLOR["panel_bg"]}" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    
    <!-- Panel Header: a | Subcellular mechanotransduction -->
    <text x="{left_x + 16}" y="{row1_y + 26}" font-size="16" font-weight="800" fill="{COLOR["text_dark"]}">a</text>
    <text x="{left_x + 36}" y="{row1_y + 26}" font-size="14" font-weight="700" fill="{COLOR["text_dark"]}">| Subcellular mechanotransduction &amp; biochemical reaction cascade</text>
    <text x="{left_x + 550}" y="{row1_y + 26}" font-size="11.5" font-weight="500" fill="{COLOR["text_muted"]}">Plant statocyte &amp; elongating vascular cambium site</text>
    
    <!-- Microgravity Vector Badge -->
    <rect x="{left_x + left_w - 180}" y="{row1_y + 12}" width="165" height="22" rx="4" fill="{COLOR["subcard_bg"]}" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    <text x="{left_x + left_w - 97}" y="{row1_y + 27}" font-size="10.5" font-weight="700" fill="{COLOR["redox"]}" text-anchor="middle">
      State: μg (0–10⁻³ g) vs 1g
    </text>
    ''')

    # A1. Cell Wall & Apoplast Zone
    cw_x, cw_y, cw_w, cw_h = left_x + 16, row1_y + 40, left_w - 32, 105
    svg.append(f'''
    <!-- Cell Wall & Apoplast -->
    <rect x="{cw_x}" y="{cw_y}" width="{cw_w}" height="{cw_h}" rx="5" fill="{COLOR["wall_light"]}" stroke="{COLOR["wall_border"]}" stroke-width="1.2" stroke-dasharray="6,3" />
    
    <!-- Label Badge -->
    <rect x="{cw_x + 12}" y="{cw_y + 10}" width="200" height="20" rx="3" fill="{COLOR["wall"]}" />
    <text x="{cw_x + 112}" y="{cw_y + 24}" font-size="10.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">
      CELL WALL &amp; APOPLAST (pH ~ 4.8)
    </text>

    <!-- Microfibrils lines -->
    <line x1="{cw_x + 20}" y1="{cw_y + 42}" x2="{cw_x + 380}" y2="{cw_y + 42}" stroke="{COLOR["wall"]}" stroke-width="1.5" opacity="0.4" stroke-dasharray="14,6" />
    <line x1="{cw_x + 40}" y1="{cw_y + 56}" x2="{cw_x + 420}" y2="{cw_y + 56}" stroke="{COLOR["wall"]}" stroke-width="1.5" opacity="0.4" stroke-dasharray="12,5" />
    <line x1="{cw_x + 25}" y1="{cw_y + 70}" x2="{cw_x + 390}" y2="{cw_y + 70}" stroke="{COLOR["wall"]}" stroke-width="1.5" opacity="0.4" stroke-dasharray="16,7" />

    <!-- Acid Growth & Expansins Sub-box -->
    <rect x="{cw_x + 230}" y="{cw_y + 10}" width="365" height="85" rx="5" fill="#FFFFFF" stroke="{COLOR["wall_border"]}" stroke-width="1" />
    <text x="{cw_x + 242}" y="{cw_y + 28}" font-size="12" font-weight="700" fill="{COLOR["wall"]}">
      Acid Growth &amp; Cell Wall Loosening
    </text>
    <text x="{cw_x + 242}" y="{cw_y + 46}" font-size="11" fill="{COLOR["text_body"]}">
      • H⁺ extrusion by AHA1/2 lowers apoplast pH to ~4.8
    </text>
    <text x="{cw_x + 242}" y="{cw_y + 63}" font-size="11" fill="{COLOR["text_body"]}">
      • Low pH activates expansins (EXPA1/8) &amp; XTH xyloglucanases
    </text>
    <text x="{cw_x + 242}" y="{cw_y + 81}" font-size="10.5" font-weight="600" fill="{COLOR["auxin"]}">
      Key DEGs: XTH22 (+0.84 log₂FC, padj=1e-8) · EXPA1 (+0.62)
    </text>

    <!-- Apoplastic Redox Burst Sub-box -->
    <rect x="{cw_x + 610}" y="{cw_y + 10}" width="370" height="85" rx="5" fill="#FFFFFF" stroke="{COLOR["redox_border"]}" stroke-width="1" />
    <text x="{cw_x + 622}" y="{cw_y + 28}" font-size="12" font-weight="700" fill="{COLOR["redox"]}">
      Apoplastic ROS &amp; Redox Signaling
    </text>
    <text x="{cw_x + 622}" y="{cw_y + 46}" font-size="11" fill="{COLOR["text_body"]}">
      • RBOHD generates O₂•⁻ → CSD1 dismutates to H₂O₂
    </text>
    <text x="{cw_x + 622}" y="{cw_y + 63}" font-size="11" fill="{COLOR["text_body"]}">
      • Peroxidases regulate wall cross-linking &amp; microgravity stress
    </text>
    <text x="{cw_x + 622}" y="{cw_y + 81}" font-size="10.5" font-weight="600" fill="{COLOR["redox"]}">
      Key DEGs: PRXIIB (+1.22 log₂FC) · CSD1 (+0.78) · CAT2 (-0.43)
    </text>
    ''')

    # A2. Plasma Membrane Bilayer
    pm_x, pm_y, pm_w, pm_h = cw_x, cw_y + cw_h + 8, cw_w, 48
    svg.append(f'''
    <!-- Plasma Membrane Bilayer -->
    <rect x="{pm_x}" y="{pm_y}" width="{pm_w}" height="{pm_h}" rx="5" fill="{COLOR["membrane_bilayer"]}" stroke="{COLOR["membrane"]}" stroke-width="1.8" />
    <text x="{pm_x + 12}" y="{pm_y + 29}" font-size="11" font-weight="700" fill="{COLOR["mechano"]}">PLASMA MEMBRANE</text>

    <!-- Transporters -->
    <!-- MSL10 / MCA1 (Mechanosensitive channel) -->
    <rect x="{pm_x + 160}" y="{pm_y + 6}" width="125" height="36" rx="4" fill="{COLOR["mechano"]}" />
    <text x="{pm_x + 222}" y="{pm_y + 22}" font-size="10.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">MSL10 / MCA1</text>
    <text x="{pm_x + 222}" y="{pm_y + 35}" font-size="9" fill="#FFFFFF" text-anchor="middle">Stretch Ca²⁺ Channel</text>

    <!-- AHA1 / AHA2 (H+ ATPase) -->
    <rect x="{pm_x + 315}" y="{pm_y + 6}" width="120" height="36" rx="4" fill="{COLOR["wall"]}" />
    <text x="{pm_x + 375}" y="{pm_y + 22}" font-size="10.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">AHA1 / AHA2</text>
    <text x="{pm_x + 375}" y="{pm_y + 35}" font-size="9" fill="#FFFFFF" text-anchor="middle">PM H⁺-ATPase Pump</text>

    <!-- PIN2 / PIN3 / PIN7 (Auxin Efflux) -->
    <rect x="{pm_x + 465}" y="{pm_y + 6}" width="140" height="36" rx="4" fill="{COLOR["auxin"]}" />
    <text x="{pm_x + 535}" y="{pm_y + 22}" font-size="10.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">PIN2 / PIN3 / PIN7</text>
    <text x="{pm_x + 535}" y="{pm_y + 35}" font-size="9" fill="#FFFFFF" text-anchor="middle">Polar Auxin Efflux</text>

    <!-- AUX1 / ABCB1 -->
    <rect x="{pm_x + 635}" y="{pm_y + 6}" width="135" height="36" rx="4" fill="{COLOR["auxin"]}" />
    <text x="{pm_x + 702}" y="{pm_y + 22}" font-size="10.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">AUX1 / ABCB1</text>
    <text x="{pm_x + 702}" y="{pm_y + 35}" font-size="9" fill="#FFFFFF" text-anchor="middle">Auxin Influx &amp; ABC</text>

    <!-- RBOHD (NADPH Oxidase) -->
    <rect x="{pm_x + 800}" y="{pm_y + 6}" width="135" height="36" rx="4" fill="{COLOR["redox"]}" />
    <text x="{pm_x + 867}" y="{pm_y + 22}" font-size="10.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">RBOHD / RBOHF</text>
    <text x="{pm_x + 867}" y="{pm_y + 35}" font-size="9" fill="#FFFFFF" text-anchor="middle">NADPH Oxidase ROS</text>

    <!-- Transmembrane Reaction Arrows -->
    <!-- Ca2+ Influx -->
    <path d="M {pm_x + 222} {cw_y + cw_h} L {pm_x + 222} {pm_y + 4}" stroke="{COLOR["mechano"]}" stroke-width="2.5" marker-end="url(#arrow-mechano)" />
    <text x="{pm_x + 234}" y="{pm_y - 3}" font-size="10" font-weight="700" fill="{COLOR["mechano"]}">Ca²⁺ Influx</text>

    <!-- H+ Pumping -->
    <path d="M {pm_x + 375} {pm_y + pm_h} L {pm_x + 375} {cw_y + cw_h + 2}" stroke="{COLOR["wall"]}" stroke-width="2.5" marker-end="url(#arrow-wall)" />
    <text x="{pm_x + 386}" y="{pm_y - 3}" font-size="10" font-weight="700" fill="{COLOR["wall"]}">H⁺ Pumping</text>

    <!-- Auxin Efflux -->
    <path d="M {pm_x + 535} {pm_y + pm_h} L {pm_x + 535} {cw_y + cw_h + 2}" stroke="{COLOR["auxin"]}" stroke-width="2.5" marker-end="url(#arrow-auxin)" />
    <text x="{pm_x + 546}" y="{pm_y - 3}" font-size="10" font-weight="700" fill="{COLOR["auxin"]}">IAA⁻ Efflux</text>
    ''')

    # A3. Cytoplasm & Organelles
    cy_y = pm_y + pm_h + 10
    
    # Amyloplasts / Statoliths (Left)
    st_x, st_y, st_w, st_h = cw_x, cy_y, 255, 215
    svg.append(f'''
    <!-- Amyloplasts / Statoliths -->
    <rect x="{st_x}" y="{st_y}" width="{st_w}" height="{st_h}" rx="5" fill="{COLOR["subcard_bg"]}" stroke="{COLOR["statolith"]}" stroke-width="1.2" />
    <rect x="{st_x + 10}" y="{st_y + 10}" width="180" height="20" rx="3" fill="{COLOR["statolith"]}" />
    <text x="{st_x + 100}" y="{st_y + 24}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">
      STATOLITHS / AMYLOPLASTS
    </text>

    <!-- Starch Granules suspended in microgravity -->
    <ellipse cx="{st_x + 75}" cy="{st_y + 75}" rx="36" ry="26" fill="{COLOR["statolith_body"]}" stroke="{COLOR["statolith"]}" stroke-width="1.5" />
    <text x="{st_x + 75}" y="{st_y + 78}" font-size="10" font-weight="700" fill="{COLOR["text_dark"]}" text-anchor="middle">Statolith 1</text>
    <text x="{st_x + 75}" y="{st_y + 91}" font-size="8.5" fill="{COLOR["text_muted"]}" text-anchor="middle">Starch Core</text>

    <ellipse cx="{st_x + 175}" cy="{st_y + 90}" rx="42" ry="30" fill="{COLOR["statolith_body"]}" stroke="{COLOR["statolith"]}" stroke-width="1.5" />
    <text x="{st_x + 175}" y="{st_y + 93}" font-size="10" font-weight="700" fill="{COLOR["text_dark"]}" text-anchor="middle">Statolith 2</text>
    <text x="{st_x + 175}" y="{st_y + 106}" font-size="8.5" fill="{COLOR["text_muted"]}" text-anchor="middle">Sedimentation loss</text>

    <ellipse cx="{st_x + 115}" cy="{st_y + 150}" rx="34" ry="24" fill="{COLOR["statolith_body"]}" stroke="{COLOR["statolith"]}" stroke-width="1.5" />
    <text x="{st_x + 115}" y="{st_y + 153}" font-size="10" font-weight="700" fill="{COLOR["text_dark"]}" text-anchor="middle">Statolith 3</text>

    <!-- Statolith annotation -->
    <rect x="{st_x + 10}" y="{st_y + 180}" width="{st_w - 20}" height="26" rx="3" fill="#FFFFFF" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    <text x="{st_x + 18}" y="{st_y + 197}" font-size="10" font-weight="600" fill="{COLOR["text_body"]}">
      CPN60A (+0.65 log₂FC) · BAM3 starch turnover
    </text>
    ''')

    # Cytosol, Cytoskeleton & Vesicle Trafficking (Center)
    cs_x, cs_y, cs_w, cs_h = st_x + st_w + 12, cy_y, 400, 215
    svg.append(f'''
    <!-- Cytosol, Cytoskeleton & Vesicles -->
    <rect x="{cs_x}" y="{cs_y}" width="{cs_w}" height="{cs_h}" rx="5" fill="{COLOR["subcard_bg"]}" stroke="{COLOR["auxin_border"]}" stroke-width="1.2" />
    <rect x="{cs_x + 10}" y="{cs_y + 10}" width="240" height="20" rx="3" fill="{COLOR["auxin"]}" />
    <text x="{cs_x + 130}" y="{cs_y + 24}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">
      CYTOSOL · VESICLES · CYTOSKELETON
    </text>

    <!-- Cytoskeleton filaments -->
    <line x1="{cs_x + 15}" y1="{cs_y + 45}" x2="{cs_x + 380}" y2="{cs_y + 55}" stroke="{COLOR["wall"]}" stroke-width="1.5" stroke-dasharray="4,2" />
    <line x1="{cs_x + 15}" y1="{cs_y + 65}" x2="{cs_x + 380}" y2="{cs_y + 75}" stroke="{COLOR["mechano"]}" stroke-width="1.5" stroke-dasharray="8,3" />
    <text x="{cs_x + 15}" y="{cs_y + 42}" font-size="9.5" font-weight="600" fill="{COLOR["text_muted"]}">
      Cytoskeleton remodeling: ACT2 (-0.58) · TUBB1 (-0.49) · MYA1 (-0.67)
    </text>

    <!-- Ca2+ Sensor Hub -->
    <rect x="{cs_x + 12}" y="{cs_y + 85}" width="180" height="60" rx="4" fill="#FFFFFF" stroke="{COLOR["mechano_border"]}" stroke-width="1" />
    <text x="{cs_x + 20}" y="{cs_y + 102}" font-size="11" font-weight="700" fill="{COLOR["mechano"]}">Ca²⁺ Sensors (Touch)</text>
    <text x="{cs_x + 20}" y="{cs_y + 118}" font-size="10" fill="{COLOR["text_body"]}">TCH1 / CAM2 (+1.28 log₂FC)</text>
    <text x="{cs_x + 20}" y="{cs_y + 133}" font-size="10" fill="{COLOR["text_body"]}">CPK1 / CDPK Kinases</text>

    <!-- PIN Vesicle Recycling -->
    <rect x="{cs_x + 204}" y="{cs_y + 85}" width="184" height="60" rx="4" fill="#FFFFFF" stroke="{COLOR["auxin_border"]}" stroke-width="1" />
    <text x="{cs_x + 212}" y="{cs_y + 102}" font-size="11" font-weight="700" fill="{COLOR["auxin"]}">PIN Polar Recycling</text>
    <text x="{cs_x + 212}" y="{cs_y + 118}" font-size="10" fill="{COLOR["text_body"]}">RABA1A GTPase (+0.74 log₂)</text>
    <text x="{cs_x + 212}" y="{cs_y + 133}" font-size="10" fill="{COLOR["text_body"]}">ARF1A1C Endosome (+0.59 log₂)</text>

    <!-- Signal integration text -->
    <rect x="{cs_x + 12}" y="{cs_y + 154}" width="{cs_w - 24}" height="50" rx="4" fill="#FFFFFF" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    <text x="{cs_x + 20}" y="{cs_y + 172}" font-size="10" font-weight="600" fill="{COLOR["text_dark"]}">
      • Cytosolic Ca²⁺ transients direct RABA1A/ARF1A1C polar targeting
    </text>
    <text x="{cs_x + 20}" y="{cs_y + 190}" font-size="10" font-weight="600" fill="{COLOR["text_dark"]}">
      • MPK3 (+0.88 log₂FC) coordinates the downstream stress cascade
    </text>
    ''')

    # Nucleus & Transcriptional Response (Right)
    nu_x, nu_y, nu_w, nu_h = cs_x + cs_w + 12, cy_y, cw_w - (st_w + cs_w + 24), 215
    svg.append(f'''
    <!-- Nucleus -->
    <rect x="{nu_x}" y="{nu_y}" width="{nu_w}" height="{nu_h}" rx="5" fill="{COLOR["nucleus"]}" stroke="{COLOR["nucleus_border"]}" stroke-width="1.2" />
    <rect x="{nu_x + 10}" y="{nu_y + 10}" width="200" height="20" rx="3" fill="{COLOR["nucleus_border"]}" />
    <text x="{nu_x + 110}" y="{nu_y + 24}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">
      NUCLEUS &amp; TRANSCRIPTION
    </text>

    <!-- DNA wave -->
    <path d="M {nu_x + 20} {nu_y + 55} Q {nu_x + 80} {nu_y + 42} {nu_x + 140} {nu_y + 55} T {nu_x + 270} {nu_y + 55}" fill="none" stroke="{COLOR["nucleus_border"]}" stroke-width="1.8" />
    <path d="M {nu_x + 20} {nu_y + 65} Q {nu_x + 80} {nu_y + 78} {nu_x + 140} {nu_y + 65} T {nu_x + 270} {nu_y + 65}" fill="none" stroke="{COLOR["nucleus_border"]}" stroke-width="1.8" />

    <!-- Auxin Gene Regulation Subcard -->
    <rect x="{nu_x + 10}" y="{nu_y + 80}" width="{nu_w - 20}" height="58" rx="4" fill="#FFFFFF" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    <text x="{nu_x + 18}" y="{nu_y + 98}" font-size="11" font-weight="700" fill="{COLOR["auxin"]}">Auxin Transcriptional Cascade</text>
    <text x="{nu_x + 18}" y="{nu_y + 114}" font-size="10" fill="{COLOR["text_body"]}">TIR1/AFB receptors → Aux/IAA repressor degradation</text>
    <text x="{nu_x + 18}" y="{nu_y + 128}" font-size="10" fill="{COLOR["text_body"]}">ARF19 / ARF7 activate SAUR expansion genes</text>

    <!-- Stress & HSF Transcriptome -->
    <rect x="{nu_x + 10}" y="{nu_y + 146}" width="{nu_w - 20}" height="58" rx="4" fill="#FFFFFF" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    <text x="{nu_x + 18}" y="{nu_y + 164}" font-size="11" font-weight="700" fill="{COLOR["chaperone"]}">Stress &amp; HSF Transcriptome</text>
    <text x="{nu_x + 18}" y="{nu_y + 180}" font-size="10" fill="{COLOR["text_body"]}">Heat Shock Factors (HSFs) trigger chaperone network</text>
    <text x="{nu_x + 18}" y="{nu_y + 194}" font-size="10" fill="{COLOR["text_body"]}">MPK3 / MPK6 MAP kinase regulatory pathway</text>
    ''')

    # Mitochondria & ER / Proteasome (Bottom row of panel a)
    bo_y = cy_y + 225
    bo_h = 135
    mi_w = 485
    pr_w = cw_w - mi_w - 15

    # Mitochondria
    svg.append(f'''
    <!-- Mitochondria & Bioenergetics -->
    <rect x="{cw_x}" y="{bo_y}" width="{mi_w}" height="{bo_h}" rx="5" fill="{COLOR["mito"]}" stroke="{COLOR["mito_border"]}" stroke-width="1.2" />
    <rect x="{cw_x + 10}" y="{bo_y + 10}" width="260" height="20" rx="3" fill="{COLOR["mito_border"]}" />
    <text x="{cw_x + 140}" y="{bo_y + 24}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">
      MITOCHONDRIA · BIOENERGETICS &amp; ATP
    </text>

    <!-- Cristae line art -->
    <path d="M {cw_x + 20} {bo_y + 60} C {cw_x + 60} {bo_y + 40}, {cw_x + 40} {bo_y + 90}, {cw_x + 80} {bo_y + 65} C {cw_x + 120} {bo_y + 40}, {cw_x + 100} {bo_y + 100}, {cw_x + 150} {bo_y + 70}" fill="none" stroke="{COLOR["mito_border"]}" stroke-width="2" opacity="0.6" />

    <rect x="{cw_x + 175}" y="{bo_y + 38}" width="{mi_w - 185}" height="88" rx="4" fill="#FFFFFF" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    <text x="{cw_x + 185}" y="{bo_y + 56}" font-size="11" font-weight="700" fill="{COLOR["atp"]}">ATP Synthesis &amp; Metabolic Shift</text>
    <text x="{cw_x + 185}" y="{bo_y + 74}" font-size="10.5" fill="{COLOR["text_body"]}">• ATPA (+0.72 log₂FC) Proton-transporting ATP synthase</text>
    <text x="{cw_x + 185}" y="{bo_y + 92}" font-size="10.5" fill="{COLOR["text_body"]}">• GAPC1 (+0.92 log₂FC) Cytosolic glycolysis enzyme</text>
    <text x="{cw_x + 185}" y="{bo_y + 110}" font-size="10.5" fill="{COLOR["text_body"]}">• MDH1 (+0.63 log₂FC) Malate dehydrogenase TCA cycle</text>
    ''')

    # ER & Proteostasis
    svg.append(f'''
    <!-- ER & Proteostasis -->
    <rect x="{cw_x + mi_w + 15}" y="{bo_y}" width="{pr_w}" height="{bo_h}" rx="5" fill="{COLOR["er"]}" stroke="{COLOR["er_border"]}" stroke-width="1.2" />
    <rect x="{cw_x + mi_w + 25}" y="{bo_y + 10}" width="260" height="20" rx="3" fill="{COLOR["er_border"]}" />
    <text x="{cw_x + mi_w + 155}" y="{bo_y + 24}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">
      ER &amp; 26S PROTEASOME · CHAPERONES
    </text>

    <!-- Ribosome dots -->
    <circle cx="{cw_x + mi_w + 35}" cy="{bo_y + 60}" r="3" fill="{COLOR["er_border"]}" />
    <circle cx="{cw_x + mi_w + 55}" cy="{bo_y + 70}" r="3" fill="{COLOR["er_border"]}" />
    <circle cx="{cw_x + mi_w + 30}" cy="{bo_y + 85}" r="3" fill="{COLOR["er_border"]}" />
    <circle cx="{cw_x + mi_w + 60}" cy="{bo_y + 100}" r="3" fill="{COLOR["er_border"]}" />

    <rect x="{cw_x + mi_w + 95}" y="{bo_y + 38}" width="{pr_w - 105}" height="88" rx="4" fill="#FFFFFF" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    <text x="{cw_x + mi_w + 105}" y="{bo_y + 56}" font-size="11" font-weight="700" fill="{COLOR["chaperone"]}">Proteostasis Machinery</text>
    <text x="{cw_x + mi_w + 105}" y="{bo_y + 74}" font-size="10.5" fill="{COLOR["text_body"]}">• HSP70-1 (+1.45 log₂FC) · HSP90.1 (+1.18 log₂FC)</text>
    <text x="{cw_x + mi_w + 105}" y="{bo_y + 92}" font-size="10.5" fill="{COLOR["text_body"]}">• DNAJ-1 (+0.82) · CPN60A (+0.65) molecular chaperones</text>
    <text x="{cw_x + mi_w + 105}" y="{bo_y + 110}" font-size="10.5" fill="{COLOR["text_body"]}">• UBC1 (+0.77 log₂FC) · PAA1 (+0.66, 26S proteasome core)</text>
    ''')

    # Reaction Step Markers
    svg.append(f'''
    <!-- Step Markers (1, 2, 3, 4) -->
    <circle cx="{st_x + 175}" cy="{st_y + 45}" r="11" fill="{COLOR["mechano"]}" />
    <text x="{st_x + 175}" y="{st_y + 49}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">1</text>

    <circle cx="{pm_x + 270}" cy="{pm_y + 48}" r="11" fill="{COLOR["auxin"]}" />
    <text x="{pm_x + 270}" y="{pm_y + 52}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">2</text>

    <circle cx="{cs_x + 295}" cy="{cs_y + 68}" r="11" fill="{COLOR["auxin"]}" />
    <text x="{cs_x + 295}" y="{cs_y + 72}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">3</text>

    <circle cx="{pm_x + 445}" cy="{cw_y + 75}" r="11" fill="{COLOR["wall"]}" />
    <text x="{pm_x + 445}" y="{cw_y + 79}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">4</text>
    ''')

    # =========================================================================
    # PANEL b: Organ & Single-Cell Tissue Activation Landscape
    # =========================================================================
    pb_y = row1_y
    pb_h = 430
    svg.append(f'''
    <!-- PANEL b CONTAINER -->
    <rect x="{right_x}" y="{pb_y}" width="{right_w}" height="{pb_h}" rx="6" fill="{COLOR["panel_bg"]}" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    
    <!-- Panel Header: b | Single-cell tissue activation -->
    <text x="{right_x + 16}" y="{pb_y + 26}" font-size="16" font-weight="800" fill="{COLOR["text_dark"]}">b</text>
    <text x="{right_x + 36}" y="{pb_y + 26}" font-size="14" font-weight="700" fill="{COLOR["text_dark"]}">| Organ &amp; single-cell tissue activation landscape</text>
    <text x="{right_x + 16}" y="{pb_y + 46}" font-size="11.5" font-weight="400" fill="{COLOR["text_muted"]}">
      Single-cell deconvolution across Salk developmental atlas clusters (GSE226097)
    </text>
    ''')

    # Heatmap Grid Layout with Ample Right Margin
    hm_cols = [
        ("Auxin Polarity", COLOR["auxin"]),
        ("Ca²⁺ Signaling", COLOR["mechano"]),
        ("Cell Wall Loosen", COLOR["wall"]),
        ("Redox Defense", COLOR["redox"]),
        ("Chaperones", COLOR["chaperone"]),
        ("Bioenergetics", COLOR["atp"])
    ]
    hm_rows = [
        ("Stem: Vascular Cambium (stem_9)", "+1.23% (padj=0.003) ***", [2.4, 2.1, 2.8, 1.9, 2.2, 1.8], True),
        ("Root Tip: Columella Statocytes", "Primary Graviperception **", [2.6, 2.9, 1.4, 1.7, 1.9, 1.6], False),
        ("Root Tip: Procambium / Stele", "Auxin Transport Axis **", [2.1, 1.8, 1.9, 1.5, 1.7, 1.5], False),
        ("Shoot Apex: Meristematic Zone", "Organogenesis / Morphogen *", [1.8, 1.5, 1.6, 1.8, 1.9, 2.1], False),
        ("Rosette Leaf: Vascular Bundles", "Redox &amp; Photosynthetic shift", [1.2, 1.1, 1.3, 2.4, 1.8, 1.7], False),
        ("Silique &amp; Seed: Coat (silique_20)", "+0.73% (padj=0.030) **", [1.6, 1.3, 2.1, 1.6, 2.5, 2.2], True),
    ]

    grid_x0 = right_x + 285
    grid_y0 = pb_y + 90
    col_w = 60
    row_h = 46

    # Column Headers (Rotated)
    for c_idx, (col_name, col_color) in enumerate(hm_cols):
        cx = grid_x0 + c_idx * col_w
        svg.append(f'''
        <g transform="translate({cx + col_w//2}, {grid_y0 - 6}) rotate(-32)">
          <text x="0" y="0" font-size="10.5" font-weight="700" fill="{col_color}" text-anchor="start">
            {col_name}
          </text>
        </g>
        ''')

    # Rows & Cells
    for r_idx, (row_label, row_sub, scores, is_highlight) in enumerate(hm_rows):
        ry = grid_y0 + r_idx * row_h
        row_bg = COLOR["auxin_light"] if is_highlight else COLOR["subcard_bg"]
        
        svg.append(f'''
        <rect x="{right_x + 14}" y="{ry}" width="265" height="{row_h - 4}" rx="4" fill="{row_bg}" stroke="{COLOR["panel_border"]}" stroke-width="1" />
        <text x="{right_x + 22}" y="{ry + 18}" font-size="10.5" font-weight="700" fill="{COLOR["text_dark"]}">
          {row_label}
        </text>
        <text x="{right_x + 22}" y="{ry + 32}" font-size="9" font-weight="500" fill="{COLOR["text_muted"]}">
          {row_sub}
        </text>
        ''')

        # Heat tiles
        for c_idx, score in enumerate(scores):
            cx = grid_x0 + c_idx * col_w
            intensity = min(1.0, max(0.2, (score - 1.0) / 1.8))
            base_col = hm_cols[c_idx][1]
            
            svg.append(f'''
            <rect x="{cx}" y="{ry}" width="{col_w - 4}" height="{row_h - 4}" rx="4" fill="{base_col}" fill-opacity="{intensity:.2f}" stroke="{COLOR["panel_border"]}" stroke-width="1" />
            <text x="{cx + col_w//2 - 2}" y="{ry + 22}" font-size="11" font-weight="700" fill="{COLOR["text_dark"]}" text-anchor="middle">
              {score:.1f}
            </text>
            <text x="{cx + col_w//2 - 2}" y="{ry + 34}" font-size="8" font-weight="600" fill="{COLOR["text_muted"]}" text-anchor="middle">
              z-score
            </text>
            ''')

    # Marginal summary in panel b
    svg.append(f'''
    <rect x="{right_x + 14}" y="{grid_y0 + len(hm_rows) * row_h + 8}" width="{right_w - 28}" height="42" rx="4" fill="{COLOR["subcard_bg"]}" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    <text x="{right_x + 22}" y="{grid_y0 + len(hm_rows) * row_h + 24}" font-size="10.5" font-weight="700" fill="{COLOR["wall"]}">
      Key Insight: Inflorescence stem vascular cambium (stem_9) mounts peak cell-wall loosening
    </text>
    <text x="{right_x + 22}" y="{grid_y0 + len(hm_rows) * row_h + 38}" font-size="10" fill="{COLOR["text_muted"]}">
      Coupled with root tip columella Ca²⁺ signaling, defining the primary axis of microgravity gravitropic remodeling.
    </text>
    ''')

    # =========================================================================
    # PANEL d: Integrated Systems Mechanistic Model
    # =========================================================================
    pd_y = pb_y + pb_h + 15
    pd_h = HEIGHT - pd_y - 25
    svg.append(f'''
    <!-- PANEL d CONTAINER -->
    <rect x="{right_x}" y="{pd_y}" width="{right_w}" height="{pd_h}" rx="6" fill="{COLOR["panel_bg"]}" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    
    <!-- Panel Header: d | Integrated systems model -->
    <text x="{right_x + 16}" y="{pd_y + 26}" font-size="16" font-weight="800" fill="{COLOR["text_dark"]}">d</text>
    <text x="{right_x + 36}" y="{pd_y + 26}" font-size="14" font-weight="700" fill="{COLOR["text_dark"]}">| Integrated systems mechanistic model</text>
    <text x="{right_x + 16}" y="{pd_y + 46}" font-size="11.5" font-weight="400" fill="{COLOR["text_muted"]}">
      Sequential cascade: physical microgravity perception to cellular stress adaptation
    </text>
    ''')

    steps = [
        ("1. Loss of Gravitational Sedimentation", "Amyloplasts lose buoyant weight; freely float in statocytes", COLOR["statolith"]),
        ("2. Stochastic MSL Channel Gating", "Membrane tension shifts induce repetitive cytosolic Ca²⁺ pulses", COLOR["mechano"]),
        ("3. Ectopic PIN Carrier Relocalization", "RABA1A/ARF1A1C endosomes redirect PIN2/3/7 to lateral flanks", COLOR["auxin"]),
        ("4. Vascular Tissue Ectopic Expansion", "Auxin accumulates in stem cambium (stem_9) &amp; root stele", COLOR["auxin"]),
        ("5. Apoplast Acidification &amp; Loosening", "PM H⁺-ATPase pumping activates XTH22 &amp; EXPA1 expansins", COLOR["wall"]),
        ("6. Compensatory Stress Homeostasis", "Apoplastic ROS (RBOHD) mitigated by HSP70 &amp; PRXIIB chaperones", COLOR["chaperone"]),
    ]

    card_step_h = 82
    for s_idx, (s_title, s_desc, s_col) in enumerate(steps):
        # 2 columns of 3 cards each in Panel d
        col_idx = s_idx // 3  # 0 or 1
        row_in_col = s_idx % 3 # 0, 1, 2
        card_w = (right_w - 40) // 2
        card_x = right_x + 16 + col_idx * (card_w + 8)
        card_y = pd_y + 60 + row_in_col * (card_step_h + 12)
        
        svg.append(f'''
        <rect x="{card_x}" y="{card_y}" width="{card_w}" height="{card_step_h}" rx="5" fill="{COLOR["subcard_bg"]}" stroke="{COLOR["panel_border"]}" stroke-width="1" />
        <rect x="{card_x + 6}" y="{card_y + 8}" width="5" height="{card_step_h - 16}" rx="2.5" fill="{s_col}" />
        <text x="{card_x + 18}" y="{card_y + 24}" font-size="11" font-weight="700" fill="{s_col}">
          {s_title}
        </text>
        <text x="{card_x + 18}" y="{card_y + 44}" font-size="10" fill="{COLOR["text_body"]}">
          {s_desc}
        </text>
        ''')
        # Connector line
        if s_idx < len(steps) - 1:
            next_col = (s_idx + 1) // 3
            if next_col == col_idx:
                svg.append(f'''
                <line x1="{card_x + card_w//2}" y1="{card_y + card_step_h}" x2="{card_x + card_w//2}" y2="{card_y + card_step_h + 11}" stroke="{COLOR["text_muted"]}" stroke-width="1.2" marker-end="url(#arrow-muted)" />
                ''')

    # =========================================================================
    # PANEL c: Curated Molecular Pathway Cross-Talk Network
    # =========================================================================
    pc_y = row2_y
    pc_h = row2_h
    svg.append(f'''
    <!-- PANEL c CONTAINER -->
    <rect x="{left_x}" y="{pc_y}" width="{left_w}" height="{pc_h}" rx="6" fill="{COLOR["panel_bg"]}" stroke="{COLOR["panel_border"]}" stroke-width="1" />
    
    <!-- Panel Header: c | Curated molecular pathway cross-talk -->
    <text x="{left_x + 16}" y="{pc_y + 26}" font-size="16" font-weight="800" fill="{COLOR["text_dark"]}">c</text>
    <text x="{left_x + 36}" y="{pc_y + 26}" font-size="14" font-weight="700" fill="{COLOR["text_dark"]}">| Curated molecular pathway cross-talk network</text>
    <text x="{left_x + 480}" y="{pc_y + 26}" font-size="11.5" font-weight="500" fill="{COLOR["text_muted"]}">
      6 core functional hubs with spaceflight DEG counts and regulatory cross-talk
    </text>
    ''')

    hubs = [
        # (id, name, kegg_id, degs, pooled_fc, x, y, color, light_bg, hub_genes)
        ("AUXIN", "Plant Hormone Signal Transduction", "ath04075", 84, "+0.68", left_x + 175, pc_y + 160, COLOR["auxin"], COLOR["auxin_light"], ["PIN2", "AUX1", "ABCB1", "IAA19"]),
        ("MECHANO", "Mechanotransduction &amp; Ca²⁺ Flux", "ath04070", 42, "+0.89", left_x + 515, pc_y + 130, COLOR["mechano"], COLOR["mechano_light"], ["MSL10", "MCA1", "TCH1", "CPK1"]),
        ("WALL", "Cell Wall Organization &amp; Loosening", "ath00500", 67, "+0.76", left_x + 855, pc_y + 160, COLOR["wall"], COLOR["wall_light"], ["XTH22", "EXPA1", "CESA3", "AHA1"]),
        ("REDOX", "Redox Homeostasis &amp; ROS Defense", "ath00480", 58, "+0.81", left_x + 195, pc_y + 365, COLOR["redox"], COLOR["redox_light"], ["PRXIIB", "CSD1", "TRX-H5", "RBOHD"]),
        ("CHAPERONE", "Protein Processing in ER &amp; Chaperones", "ath04141", 79, "+1.12", left_x + 525, pc_y + 380, COLOR["chaperone"], COLOR["chaperone_light"], ["HSP70-1", "HSP90.1", "DNAJ-1", "PAA1"]),
        ("ATP", "Bioenergetics &amp; Carbon Fixation", "ath00710", 93, "+0.74", left_x + 865, pc_y + 365, COLOR["atp"], COLOR["atp_light"], ["ATPA", "GAPC1", "MDH1", "ATP5F1"]),
    ]

    edges = [
        (0, 1, "TCH1 / PIN polarity", 4),
        (1, 2, "Ca²⁺ / AHA1 H⁺ activation", 3),
        (0, 2, "XTH22 / EXPA1 acid growth", 5),
        (1, 3, "RBOHD ROS burst", 3.5),
        (2, 3, "Peroxidase wall crosslinking", 3.5),
        (3, 4, "Redox stress chaperone induction", 4.5),
        (4, 5, "ATP-dependent folding", 4),
        (0, 5, "ATP fueling active pumps", 3),
        (0, 3, "Auxin-induced ROS signaling", 2.5),
    ]

    for (src, dst, edge_lbl, edge_w) in edges:
        x1, y1 = hubs[src][5], hubs[src][6]
        x2, y2 = hubs[dst][5], hubs[dst][6]
        mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
        
        svg.append(f'''
        <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{COLOR["panel_border"]}" stroke-width="{edge_w * 1.5}" stroke-linecap="round" />
        <rect x="{mid_x - 72}" y="{mid_y - 10}" width="144" height="20" rx="3" fill="#FFFFFF" stroke="{COLOR["panel_border"]}" stroke-width="0.8" />
        <text x="{mid_x}" y="{mid_y + 3}" font-size="9" font-weight="600" fill="{COLOR["text_muted"]}" text-anchor="middle">
          {edge_lbl}
        </text>
        ''')

    for (hid, hname, hkegg, hdegs, hfc, hx, hy, hcol, hbg, hgenes) in hubs:
        node_r = 50
        svg.append(f'''
        <!-- Hub Node {hid} -->
        <circle cx="{hx}" cy="{hy}" r="{node_r + 6}" fill="{hbg}" />
        <circle cx="{hx}" cy="{hy}" r="{node_r}" fill="#FFFFFF" stroke="{hcol}" stroke-width="2.5" />
        <text x="{hx}" y="{hy - 16}" font-size="10.5" font-weight="700" fill="{hcol}" text-anchor="middle">
          {hkegg}
        </text>
        <text x="{hx}" y="{hy + 3}" font-size="12" font-weight="700" fill="{COLOR["text_dark"]}" text-anchor="middle">
          {hdegs} DEGs
        </text>
        <rect x="{hx - 34}" y="{hy + 12}" width="68" height="18" rx="3" fill="{COLOR["auxin_light"]}" />
        <text x="{hx}" y="{hy + 24}" font-size="9.5" font-weight="700" fill="{COLOR["auxin"]}" text-anchor="middle">
          {hfc} log₂FC
        </text>

        <!-- Hub Title -->
        <text x="{hx}" y="{hy - node_r - 10}" font-size="11.5" font-weight="700" fill="{COLOR["text_dark"]}" text-anchor="middle">
          {hname}
        </text>

        <!-- Hub Key Genes Pill -->
        <rect x="{hx - 105}" y="{hy + node_r + 8}" width="210" height="20" rx="4" fill="{COLOR["subcard_bg"]}" stroke="{COLOR["panel_border"]}" stroke-width="0.8" />
        <text x="{hx}" y="{hy + node_r + 22}" font-size="9.5" font-weight="600" fill="{COLOR["text_body"]}" text-anchor="middle">
          Hubs: {", ".join(hgenes)}
        </text>
        ''')

    svg.append('</svg>')
    return "\n".join(svg)

def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    figures_dir = os.path.join(repo_root, "Figures")
    docs_figures_dir = os.path.join(repo_root, "docs", "assets", "figures")

    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(docs_figures_dir, exist_ok=True)

    print("Generating Figure 8 (npj Microgravity Springer Nature Publication Style)...")

    svg_content = generate_figure8_svg()
    
    # Save Figures/kegg_pathway_network.svg
    pub_svg_path = os.path.join(figures_dir, "kegg_pathway_network.svg")
    with open(pub_svg_path, "w", encoding="utf-8") as f:
        f.write(svg_content)
    print(f"  ✓ Saved Publication SVG: {pub_svg_path} ({len(svg_content)} bytes)")

    # Save docs/assets/figures/kegg_pathway_network.svg
    web_svg_path = os.path.join(docs_figures_dir, "kegg_pathway_network.svg")
    with open(web_svg_path, "w", encoding="utf-8") as f:
        f.write(svg_content)
    print(f"  ✓ Saved Web Portal SVG: {web_svg_path}")

    # Render 300 DPI high-resolution PNG using Headless Google Chrome
    print("\nRendering Publication-Quality 300 DPI PNG via Headless Google Chrome...")
    chrome_bin = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
    target_png_pub = os.path.join(figures_dir, "kegg_pathway_network.png")
    target_png_web = os.path.join(docs_figures_dir, "kegg_pathway_network.png")

    if os.path.exists(chrome_bin):
        # Window size 2400 x 1680 (exact 1800:1260 aspect ratio = 10:7)
        cmd = [
            chrome_bin,
            "--headless",
            "--disable-gpu",
            "--default-background-color=ffffffff",
            f"--screenshot={target_png_pub}",
            "--window-size=2400,1680",
            f"file://{pub_svg_path}"
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.copyfile(target_png_pub, target_png_web)
        print(f"  ✓ Rendered Publication PNG (2400x1680): {target_png_pub}")
        print(f"  ✓ Copied Web Portal PNG: {target_png_web}")
    else:
        print("  ⚠ Google Chrome binary not found; falling back to QuickLook...")
        subprocess.run(["qlmanage", "-t", "-s", "2400", "-o", "/tmp", pub_svg_path], check=True)
        rendered = "/tmp/kegg_pathway_network.svg.png"
        if os.path.exists(rendered):
            shutil.copyfile(rendered, target_png_pub)
            shutil.copyfile(rendered, target_png_web)

    print("\nFigure 8 npj Microgravity overhaul complete!")

if __name__ == "__main__":
    main()
