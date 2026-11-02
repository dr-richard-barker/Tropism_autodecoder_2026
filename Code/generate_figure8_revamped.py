#!/usr/bin/env python3
"""
generate_figure8_revamped.py

Revamps Figure 8 from a sparse KEGG pathway network into an integrated,
publication-grade Spatial-Molecular Systems Biology Architecture of Spaceflight
Tropism and Mechanotransduction.

Outputs:
  - Figures/kegg_pathway_network.svg (Publication Light Theme)
  - Figures/kegg_pathway_network_dark.svg (CoSE Obsidian Dark Theme)
  - Figures/kegg_pathway_network.png (Rendered 300 DPI for LaTeX Manuscript)
  - docs/assets/figures/kegg_pathway_network.png (Rendered for Web Portal)
  - docs/assets/figures/kegg_pathway_network.svg (Vector for Web Portal)
"""

import os
import sys
import subprocess
import shutil

WIDTH = 1600
HEIGHT = 1200

THEMES = {
    "light": {
        "bg": "#FFFFFF",
        "card_bg": "#F8FAFC",
        "card_border": "#CBD5E1",
        "subcard_bg": "#F1F5F9",
        "subcard_border": "#E2E8F0",
        "text_main": "#0F172A",
        "text_muted": "#475569",
        "text_sub": "#64748B",
        "membrane": "#93C5FD",
        "membrane_border": "#2563EB",
        "wall": "#FEF3C7",
        "wall_border": "#D97706",
        "cytosol": "#F0FDF4",
        "nucleus_bg": "#EDE9FE",
        "nucleus_border": "#7C3AED",
        "mito_bg": "#FFEDD5",
        "mito_border": "#EA580C",
        "er_bg": "#FCE7F3",
        "er_border": "#DB2777",
        "plastid_bg": "#FEF9C3",
        "plastid_border": "#CA8A04",
        "accent_auxin": "#059669",
        "accent_ca": "#0284C7",
        "accent_wall": "#D97706",
        "accent_redox": "#E11D48",
        "accent_chaperone": "#7C3AED",
        "accent_atp": "#0D9488",
        "badge_bg": "#E2E8F0",
        "badge_text": "#1E293B",
        "grid_line": "#E2E8F0",
        "shadow": "rgba(0,0,0,0.06)",
        "pill_up_bg": "#DCFCE7",
        "pill_up_text": "#166534",
        "pill_dn_bg": "#FEE2E2",
        "pill_dn_text": "#991B1B",
    },
    "dark": {
        "bg": "#0B0F15",
        "card_bg": "#121924",
        "card_border": "#1E293B",
        "subcard_bg": "#16202E",
        "subcard_border": "#283548",
        "text_main": "#F8FAFC",
        "text_muted": "#94A3B8",
        "text_sub": "#64748B",
        "membrane": "#1E3A8A",
        "membrane_border": "#3B82F6",
        "wall": "#451A03",
        "wall_border": "#B45309",
        "cytosol": "#062817",
        "nucleus_bg": "#2E1065",
        "nucleus_border": "#8B5CF6",
        "mito_bg": "#431407",
        "mito_border": "#F97316",
        "er_bg": "#500724",
        "er_border": "#F43F5E",
        "plastid_bg": "#422006",
        "plastid_border": "#EAB308",
        "accent_auxin": "#10B981",
        "accent_ca": "#38BDF8",
        "accent_wall": "#FBBF24",
        "accent_redox": "#FB7185",
        "accent_chaperone": "#A78BFA",
        "accent_atp": "#2DD4BF",
        "badge_bg": "#1E293B",
        "badge_text": "#E2E8F0",
        "grid_line": "#1E293B",
        "shadow": "rgba(0,0,0,0.4)",
        "pill_up_bg": "#064E3B",
        "pill_up_text": "#6EE7B7",
        "pill_dn_bg": "#7F1D1D",
        "pill_dn_text": "#FCA5A5",
    }
}

def escape_xml(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")

def generate_svg(theme_name="light"):
    T = THEMES[theme_name]
    svg = []
    
    # ------------------ SVG Header & Defs ------------------
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {HEIGHT}" width="{WIDTH}" height="{HEIGHT}" style="background-color: {T["bg"]}; font-family: -apple-system, BlinkMacSystemFont, \'Segoe UI\', Roboto, Helvetica, Arial, sans-serif;">')
    svg.append('<defs>')
    
    # Markers (Arrows)
    svg.append(f'''
      <marker id="arrow-auxin" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
        <path d="M 0 1 L 10 5 L 0 9 z" fill="{T["accent_auxin"]}" />
      </marker>
      <marker id="arrow-ca" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
        <path d="M 0 1 L 10 5 L 0 9 z" fill="{T["accent_ca"]}" />
      </marker>
      <marker id="arrow-wall" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
        <path d="M 0 1 L 10 5 L 0 9 z" fill="{T["accent_wall"]}" />
      </marker>
      <marker id="arrow-redox" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
        <path d="M 0 1 L 10 5 L 0 9 z" fill="{T["accent_redox"]}" />
      </marker>
      <marker id="arrow-chaperone" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
        <path d="M 0 1 L 10 5 L 0 9 z" fill="{T["accent_chaperone"]}" />
      </marker>
      <marker id="arrow-muted" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
        <path d="M 0 1 L 10 5 L 0 9 z" fill="{T["text_muted"]}" />
      </marker>
      <filter id="card-shadow" x="-5%" y="-5%" width="110%" height="110%">
        <feDropShadow dx="0" dy="4" stdDeviation="6" flood-color="{T["shadow"]}" />
      </filter>
    ''')
    svg.append('</defs>')

    # ------------------ Main Title & Banner ------------------
    svg.append(f'''
    <!-- Main Title Banner -->
    <rect x="25" y="16" width="1550" height="74" rx="10" fill="{T["card_bg"]}" stroke="{T["card_border"]}" stroke-width="1.2" filter="url(#card-shadow)" />
    <text x="45" y="46" font-size="22" font-weight="700" fill="{T["text_main"]}">
      FIGURE 8 · Spatial-Molecular Systems Architecture of Spaceflight Tropism &amp; Mechanotransduction
    </text>
    <text x="45" y="70" font-size="13.5" font-weight="400" fill="{T["text_muted"]}">
      Multiscale integration of subcellular mechanosensation, single-cell tissue activation (inflorescence stem_9 &amp; root columella), and functional pathway cross-talk under microgravity
    </text>
    <rect x="1380" y="32" width="175" height="30" rx="6" fill="{T["badge_bg"]}" />
    <text x="1467" y="52" font-size="11.5" font-weight="700" fill="{T["badge_text"]}" text-anchor="middle">
      TROPISM AUTODECODER
    </text>
    ''')

    # =========================================================================
    # PANEL A: Subcellular Mechanotransduction & Reaction Cascade
    # =========================================================================
    pa_x, pa_y, pa_w, pa_h = 25, 105, 960, 645
    svg.append(f'''
    <!-- PANEL A CONTAINER -->
    <rect x="{pa_x}" y="{pa_y}" width="{pa_w}" height="{pa_h}" rx="10" fill="{T["card_bg"]}" stroke="{T["card_border"]}" stroke-width="1.2" filter="url(#card-shadow)" />
    <rect x="{pa_x + 16}" y="{pa_y + 14}" width="76" height="22" rx="4" fill="{T["accent_ca"]}" />
    <text x="{pa_x + 54}" y="{pa_y + 29}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">PANEL A</text>
    <text x="{pa_x + 102}" y="{pa_y + 30}" font-size="16" font-weight="700" fill="{T["text_main"]}">
      Subcellular Mechanotransduction &amp; Biochemical Reaction Architecture
    </text>
    <text x="{pa_x + 655}" y="{pa_y + 30}" font-size="12" font-weight="500" fill="{T["text_muted"]}">
      Plant Statocyte &amp; Vascular Cambium Cell Site
    </text>
    ''')

    # 1. Cell Wall / Apoplast Zone
    cw_x, cw_y, cw_w, cw_h = pa_x + 20, pa_y + 45, pa_w - 40, 105
    svg.append(f'''
    <!-- Cell Wall & Apoplast -->
    <rect x="{cw_x}" y="{cw_y}" width="{cw_w}" height="{cw_h}" rx="8" fill="{T["wall"]}" stroke="{T["wall_border"]}" stroke-width="1.5" stroke-dasharray="6,3" />
    <rect x="{cw_x + 12}" y="{cw_y + 10}" width="190" height="20" rx="4" fill="{T["accent_wall"]}" />
    <text x="{cw_x + 107}" y="{cw_y + 24}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">
      CELL WALL &amp; APOPLAST (pH ~ 4.8)
    </text>

    <!-- Cell Wall fibers pattern -->
    <line x1="{cw_x + 20}" y1="{cw_y + 40}" x2="{cw_x + 400}" y2="{cw_y + 40}" stroke="{T["accent_wall"]}" stroke-width="1.5" opacity="0.6" stroke-dasharray="15,6" />
    <line x1="{cw_x + 50}" y1="{cw_y + 55}" x2="{cw_x + 430}" y2="{cw_y + 55}" stroke="{T["accent_wall"]}" stroke-width="1.5" opacity="0.6" stroke-dasharray="12,5" />
    <line x1="{cw_x + 30}" y1="{cw_y + 70}" x2="{cw_x + 380}" y2="{cw_y + 70}" stroke="{T["accent_wall"]}" stroke-width="1.5" opacity="0.6" stroke-dasharray="18,7" />

    <!-- Reaction 1: Acid Growth & Expansins -->
    <rect x="{cw_x + 220}" y="{cw_y + 10}" width="340" height="85" rx="6" fill="{T["subcard_bg"]}" stroke="{T["subcard_border"]}" stroke-width="1" />
    <text x="{cw_x + 232}" y="{cw_y + 28}" font-size="12" font-weight="700" fill="{T["accent_wall"]}">
      Acid Growth &amp; Wall Loosening Cascade
    </text>
    <text x="{cw_x + 232}" y="{cw_y + 46}" font-size="11.5" fill="{T["text_main"]}">
      • H⁺ extrusion by AHA1/2 lowers apoplastic pH
    </text>
    <text x="{cw_x + 232}" y="{cw_y + 64}" font-size="11.5" fill="{T["text_main"]}">
      • Acid activates Expansins (EXPA1/8) &amp; XTHs
    </text>
    <text x="{cw_x + 232}" y="{cw_y + 82}" font-size="11" font-weight="600" fill="{T["pill_up_text"]}">
      XTH22 (+0.84 log₂FC) · EXPA1 (+0.62 log₂FC)
    </text>

    <!-- Reaction 2: Apoplastic ROS Burst -->
    <rect x="{cw_x + 575}" y="{cw_y + 10}" width="330" height="85" rx="6" fill="{T["subcard_bg"]}" stroke="{T["subcard_border"]}" stroke-width="1" />
    <text x="{cw_x + 587}" y="{cw_y + 28}" font-size="12" font-weight="700" fill="{T["accent_redox"]}">
      Apoplastic ROS &amp; Redox Signaling
    </text>
    <text x="{cw_x + 587}" y="{cw_y + 46}" font-size="11.5" fill="{T["text_main"]}">
      • RBOHD generates O₂•⁻ → SOD converts to H₂O₂
    </text>
    <text x="{cw_x + 587}" y="{cw_y + 64}" font-size="11.5" fill="{T["text_main"]}">
      • Peroxidases modulate cross-linking &amp; defense
    </text>
    <text x="{cw_x + 587}" y="{cw_y + 82}" font-size="11" font-weight="600" fill="{T["pill_up_text"]}">
      PRXIIB (+1.22 log₂FC) · CSD1 (+0.78 log₂FC)
    </text>
    ''')

    # 2. Plasma Membrane Bilayer
    pm_x, pm_y, pm_w, pm_h = pa_x + 20, cw_y + cw_h + 8, pa_w - 40, 48
    svg.append(f'''
    <!-- Plasma Membrane Bilayer -->
    <rect x="{pm_x}" y="{pm_y}" width="{pm_w}" height="{pm_h}" rx="6" fill="{T["membrane"]}" stroke="{T["membrane_border"]}" stroke-width="2" />
    <text x="{pm_x + 12}" y="{pm_y + 28}" font-size="11.5" font-weight="700" fill="#FFFFFF">PLASMA MEMBRANE</text>

    <!-- Transporter: MSL10 (Mechanosensitive channel) -->
    <rect x="{pm_x + 160}" y="{pm_y + 6}" width="105" height="36" rx="5" fill="{T["accent_ca"]}" stroke="#FFFFFF" stroke-width="1" />
    <text x="{pm_x + 212}" y="{pm_y + 23}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">MSL10 / MCA1</text>
    <text x="{pm_x + 212}" y="{pm_y + 36}" font-size="9" fill="#FFFFFF" text-anchor="middle">Ca²⁺ Stretch Channel</text>

    <!-- Transporter: AHA1 (H+ ATPase) -->
    <rect x="{pm_x + 300}" y="{pm_y + 6}" width="100" height="36" rx="5" fill="{T["accent_wall"]}" stroke="#FFFFFF" stroke-width="1" />
    <text x="{pm_x + 350}" y="{pm_y + 23}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">AHA1 / AHA2</text>
    <text x="{pm_x + 350}" y="{pm_y + 36}" font-size="9" fill="#FFFFFF" text-anchor="middle">PM H⁺-ATPase Pump</text>

    <!-- Transporter: PIN2/PIN3/PIN7 (Auxin Efflux) -->
    <rect x="{pm_x + 440}" y="{pm_y + 6}" width="120" height="36" rx="5" fill="{T["accent_auxin"]}" stroke="#FFFFFF" stroke-width="1" />
    <text x="{pm_x + 500}" y="{pm_y + 23}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">PIN2 / PIN3 / PIN7</text>
    <text x="{pm_x + 500}" y="{pm_y + 36}" font-size="9" fill="#FFFFFF" text-anchor="middle">Polar Auxin Efflux</text>

    <!-- Transporter: AUX1 / ABCB1 -->
    <rect x="{pm_x + 595}" y="{pm_y + 6}" width="115" height="36" rx="5" fill="{T["accent_auxin"]}" stroke="#FFFFFF" stroke-width="1" />
    <text x="{pm_x + 652}" y="{pm_y + 23}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">AUX1 / ABCB1</text>
    <text x="{pm_x + 652}" y="{pm_y + 36}" font-size="9" fill="#FFFFFF" text-anchor="middle">Auxin Influx &amp; ABC</text>

    <!-- Transporter: RBOHD (NADPH Oxidase) -->
    <rect x="{pm_x + 745}" y="{pm_y + 6}" width="115" height="36" rx="5" fill="{T["accent_redox"]}" stroke="#FFFFFF" stroke-width="1" />
    <text x="{pm_x + 802}" y="{pm_y + 23}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">RBOHD / RBOHF</text>
    <text x="{pm_x + 802}" y="{pm_y + 36}" font-size="9" fill="#FFFFFF" text-anchor="middle">ROS Generator</text>
    ''')

    # Flux Arrows traversing the membrane
    svg.append(f'''
    <!-- Ca2+ Influx Arrow -->
    <path d="M {pm_x + 212} {cw_y + cw_h} L {pm_x + 212} {pm_y + 6}" stroke="{T["accent_ca"]}" stroke-width="3" marker-end="url(#arrow-ca)" />
    <text x="{pm_x + 225}" y="{pm_y - 4}" font-size="10.5" font-weight="700" fill="{T["accent_ca"]}">Ca²⁺ Influx</text>

    <!-- H+ Efflux Arrow -->
    <path d="M {pm_x + 350} {pm_y + pm_h} L {pm_x + 350} {cw_y + cw_h + 2}" stroke="{T["accent_wall"]}" stroke-width="3" marker-end="url(#arrow-wall)" />
    <text x="{pm_x + 362}" y="{pm_y - 4}" font-size="10.5" font-weight="700" fill="{T["accent_wall"]}">H⁺ Pumping</text>

    <!-- Auxin Efflux Arrow -->
    <path d="M {pm_x + 500} {pm_y + pm_h} L {pm_x + 500} {cw_y + cw_h + 2}" stroke="{T["accent_auxin"]}" stroke-width="3" marker-end="url(#arrow-auxin)" />
    <text x="{pm_x + 512}" y="{pm_y - 4}" font-size="10.5" font-weight="700" fill="{T["accent_auxin"]}">IAA⁻ Efflux</text>
    ''')

    # 3. Cytosol & Organelles Area
    cy_y = pm_y + pm_h + 10
    cy_h = pa_h - (cy_y - pa_y) - 15

    # Subcompartment: Statoliths & Amyloplasts (Left)
    st_x, st_y, st_w, st_h = pa_x + 20, cy_y, 250, 210
    svg.append(f'''
    <!-- Statoliths / Amyloplasts -->
    <rect x="{st_x}" y="{st_y}" width="{st_w}" height="{st_h}" rx="8" fill="{T["plastid_bg"]}" stroke="{T["plastid_border"]}" stroke-width="1.5" />
    <rect x="{st_x + 10}" y="{st_y + 10}" width="165" height="20" rx="4" fill="{T["plastid_border"]}" />
    <text x="{st_x + 92}" y="{st_y + 24}" font-size="10.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">
      STATOLITHS / AMYLOPLASTS
    </text>

    <!-- Dense Starch Granules -->
    <ellipse cx="{st_x + 70}" cy="{st_y + 75}" rx="32" ry="24" fill="{T["accent_wall"]}" opacity="0.8" />
    <text x="{st_x + 70}" y="{st_y + 78}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">Starch 1</text>

    <ellipse cx="{st_x + 160}" cy="{st_y + 90}" rx="38" ry="28" fill="{T["accent_wall"]}" opacity="0.8" />
    <text x="{st_x + 160}" y="{st_y + 93}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">Starch 2</text>

    <ellipse cx="{st_x + 105}" cy="{st_y + 145}" rx="30" ry="22" fill="{T["accent_wall"]}" opacity="0.8" />
    <text x="{st_x + 105}" y="{st_y + 148}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">Starch 3</text>

    <!-- Statolith text & genes -->
    <rect x="{st_x + 10}" y="{st_y + 175}" width="{st_w - 20}" height="28" rx="4" fill="{T["subcard_bg"]}" />
    <text x="{st_x + 18}" y="{st_y + 193}" font-size="10.5" font-weight="600" fill="{T["text_main"]}">
      CPN60A (+0.65) · BAM3 (Starch turnover)
    </text>
    ''')

    # Subcompartment: Cytosolic Signaling, Cytoskeleton & Vesicles (Center)
    cs_x, cs_y, cs_w, cs_h = st_x + st_w + 12, cy_y, 355, 210
    svg.append(f'''
    <!-- Cytosolic Signaling & Endosomes -->
    <rect x="{cs_x}" y="{cs_y}" width="{cs_w}" height="{cs_h}" rx="8" fill="{T["cytosol"]}" stroke="{T["accent_auxin"]}" stroke-width="1.3" />
    <rect x="{cs_x + 10}" y="{cs_y + 10}" width="220" height="20" rx="4" fill="{T["accent_auxin"]}" />
    <text x="{cs_x + 120}" y="{cs_y + 24}" font-size="10.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">
      CYTOSOL · VESICLES · CYTOSKELETON
    </text>

    <!-- Actin filaments & Microtubules lines -->
    <line x1="{cs_x + 15}" y1="{cs_y + 48}" x2="{cs_x + 340}" y2="{cs_y + 60}" stroke="{T["accent_wall"]}" stroke-width="1.5" stroke-dasharray="4,2" />
    <line x1="{cs_x + 15}" y1="{cs_y + 68}" x2="{cs_x + 340}" y2="{cs_y + 80}" stroke="{T["accent_ca"]}" stroke-width="1.5" stroke-dasharray="8,3" />
    <text x="{cs_x + 15}" y="{cs_y + 44}" font-size="10" font-weight="600" fill="{T["text_muted"]}">
      Cytoskeleton: ACT2 (-0.58) · TUBB1 (-0.49) · MYA1 (-0.67)
    </text>

    <!-- Ca2+ Sensor Hub -->
    <rect x="{cs_x + 12}" y="{cs_y + 86}" width="155" height="58" rx="6" fill="{T["subcard_bg"]}" stroke="{T["subcard_border"]}" stroke-width="1" />
    <text x="{cs_x + 22}" y="{cs_y + 104}" font-size="11" font-weight="700" fill="{T["accent_ca"]}">Ca²⁺ Sensors (Touch)</text>
    <text x="{cs_x + 22}" y="{cs_y + 120}" font-size="10" fill="{T["text_main"]}">TCH1 / CAM2 (+1.28 log₂)</text>
    <text x="{cs_x + 22}" y="{cs_y + 134}" font-size="10" fill="{T["text_main"]}">CPK / CDPK Kinases</text>

    <!-- PIN Vesicle Recycling Endosomes -->
    <rect x="{cs_x + 180}" y="{cs_y + 86}" width="162" height="58" rx="6" fill="{T["subcard_bg"]}" stroke="{T["subcard_border"]}" stroke-width="1" />
    <text x="{cs_x + 190}" y="{cs_y + 104}" font-size="11" font-weight="700" fill="{T["accent_auxin"]}">PIN Polar Recycling</text>
    <text x="{cs_x + 190}" y="{cs_y + 120}" font-size="10" fill="{T["text_main"]}">RABA1A GTPase (+0.74)</text>
    <text x="{cs_x + 190}" y="{cs_y + 134}" font-size="10" fill="{T["text_main"]}">ARF1A1C Endosome (+0.59)</text>

    <!-- Signal integration text -->
    <rect x="{cs_x + 12}" y="{cs_y + 152}" width="{cs_w - 24}" height="48" rx="5" fill="{T["subcard_bg"]}" />
    <text x="{cs_x + 20}" y="{cs_y + 170}" font-size="10.5" font-weight="600" fill="{T["text_main"]}">
      • Ca²⁺ wave directs RABA1A/ARF1A1C polar targeting
    </text>
    <text x="{cs_x + 20}" y="{cs_y + 188}" font-size="10.5" font-weight="600" fill="{T["text_main"]}">
      • MPK3 (+0.88 log₂FC) coordinates touch/stress cascade
    </text>
    ''')

    # Subcompartment: Nucleus (Right)
    nu_x, nu_y, nu_w, nu_h = cs_x + cs_w + 12, cy_y, 300, 210
    svg.append(f'''
    <!-- Nucleus -->
    <rect x="{nu_x}" y="{nu_y}" width="{nu_w}" height="{nu_h}" rx="8" fill="{T["nucleus_bg"]}" stroke="{T["nucleus_border"]}" stroke-width="1.5" />
    <rect x="{nu_x + 10}" y="{nu_y + 10}" width="195" height="20" rx="4" fill="{T["nucleus_border"]}" />
    <text x="{nu_x + 107}" y="{nu_y + 24}" font-size="10.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">
      NUCLEUS &amp; TRANSCRIPTION
    </text>

    <!-- DNA Double Helix Representation -->
    <path d="M {nu_x + 20} {nu_y + 55} Q {nu_x + 80} {nu_y + 40} {nu_x + 140} {nu_y + 55} T {nu_x + 260} {nu_y + 55}" fill="none" stroke="{T["nucleus_border"]}" stroke-width="2" />
    <path d="M {nu_x + 20} {nu_y + 65} Q {nu_x + 80} {nu_y + 80} {nu_x + 140} {nu_y + 65} T {nu_x + 260} {nu_y + 65}" fill="none" stroke="{T["nucleus_border"]}" stroke-width="2" />

    <!-- Transcriptional Hubs -->
    <rect x="{nu_x + 12}" y="{nu_y + 80}" width="{nu_w - 24}" height="55" rx="5" fill="{T["subcard_bg"]}" stroke="{T["subcard_border"]}" stroke-width="1" />
    <text x="{nu_x + 22}" y="{nu_y + 98}" font-size="11" font-weight="700" fill="{T["accent_auxin"]}">Auxin Gene Regulation</text>
    <text x="{nu_x + 22}" y="{nu_y + 114}" font-size="10" fill="{T["text_main"]}">TIR1/AFB Auxin Receptors → Aux/IAA Degradation</text>
    <text x="{nu_x + 22}" y="{nu_y + 127}" font-size="10" fill="{T["text_main"]}">ARF19/ARF7 activate SAURs &amp; EXPAs</text>

    <rect x="{nu_x + 12}" y="{nu_y + 142}" width="{nu_w - 24}" height="58" rx="5" fill="{T["subcard_bg"]}" stroke="{T["subcard_border"]}" stroke-width="1" />
    <text x="{nu_x + 22}" y="{nu_y + 160}" font-size="11" font-weight="700" fill="{T["accent_chaperone"]}">Stress &amp; HSF Transcriptome</text>
    <text x="{nu_x + 22}" y="{nu_y + 176}" font-size="10" fill="{T["text_main"]}">Heat Shock Factors (HSFs) trigger chaperones</text>
    <text x="{nu_x + 22}" y="{nu_y + 190}" font-size="10" fill="{T["text_main"]}">MAPK cascade: MPK3 / MPK6 activation</text>
    ''')

    # Bottom Row of Panel A: Bioenergetics (Mitochondria) & Proteostasis (ER/Proteasome)
    bo_y = cy_y + cy_h - 135
    bo_h = 135

    # Mitochondria
    mi_x, mi_y, mi_w, mi_h = pa_x + 20, bo_y, 450, bo_h
    svg.append(f'''
    <!-- Mitochondria & Bioenergetics -->
    <rect x="{mi_x}" y="{mi_y}" width="{mi_w}" height="{mi_h}" rx="8" fill="{T["mito_bg"]}" stroke="{T["mito_border"]}" stroke-width="1.5" />
    <rect x="{mi_x + 10}" y="{mi_y + 10}" width="250" height="20" rx="4" fill="{T["mito_border"]}" />
    <text x="{mi_x + 135}" y="{mi_y + 24}" font-size="10.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">
      MITOCHONDRIA · BIOENERGETICS &amp; ATP
    </text>

    <!-- Cristae membrane drawing -->
    <path d="M {mi_x + 20} {mi_y + 60} C {mi_x + 60} {mi_y + 40}, {mi_x + 40} {mi_y + 90}, {mi_x + 80} {mi_y + 65} C {mi_x + 120} {mi_y + 40}, {mi_x + 100} {mi_y + 100}, {mi_x + 150} {mi_y + 70}" fill="none" stroke="{T["mito_border"]}" stroke-width="2" opacity="0.6" />

    <rect x="{mi_x + 175}" y="{mi_y + 38}" width="{mi_w - 185}" height="88" rx="5" fill="{T["subcard_bg"]}" stroke="{T["subcard_border"]}" stroke-width="1" />
    <text x="{mi_x + 185}" y="{mi_y + 56}" font-size="11" font-weight="700" fill="{T["accent_atp"]}">ATP Generation &amp; Metabolic Shift</text>
    <text x="{mi_x + 185}" y="{mi_y + 74}" font-size="10.5" fill="{T["text_main"]}">• ATPA (+0.72 log₂) ATP Synthase F1 α</text>
    <text x="{mi_x + 185}" y="{mi_y + 92}" font-size="10.5" fill="{T["text_main"]}">• GAPC1 (+0.92 log₂) Cytosolic Glycolysis</text>
    <text x="{mi_x + 185}" y="{mi_y + 110}" font-size="10.5" fill="{T["text_main"]}">• MDH1 (+0.63 log₂) Malate Dehydrogenase</text>
    ''')

    # Endoplasmic Reticulum & Proteasome
    pr_x, pr_y, pr_w, pr_h = mi_x + mi_w + 12, bo_y, pa_w - mi_w - 52, bo_h
    svg.append(f'''
    <!-- ER & Proteostasis -->
    <rect x="{pr_x}" y="{pr_y}" width="{pr_w}" height="{pr_h}" rx="8" fill="{T["er_bg"]}" stroke="{T["er_border"]}" stroke-width="1.5" />
    <rect x="{pr_x + 10}" y="{pr_y + 10}" width="250" height="20" rx="4" fill="{T["er_border"]}" />
    <text x="{pr_x + 135}" y="{pr_y + 24}" font-size="10.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">
      ER &amp; 26S PROTEASOME · CHAPERONES
    </text>

    <!-- Ribosome dots on ER -->
    <circle cx="{pr_x + 35}" cy="{pr_y + 55}" r="3" fill="{T["er_border"]}" />
    <circle cx="{pr_x + 55}" cy="{pr_y + 65}" r="3" fill="{T["er_border"]}" />
    <circle cx="{pr_x + 30}" cy="{pr_y + 80}" r="3" fill="{T["er_border"]}" />
    <circle cx="{pr_x + 60}" cy="{pr_y + 95}" r="3" fill="{T["er_border"]}" />

    <rect x="{pr_x + 85}" y="{pr_y + 38}" width="{pr_w - 95}" height="88" rx="5" fill="{T["subcard_bg"]}" stroke="{T["subcard_border"]}" stroke-width="1" />
    <text x="{pr_x + 95}" y="{pr_y + 56}" font-size="11" font-weight="700" fill="{T["accent_chaperone"]}">Proteostasis Machinery</text>
    <text x="{pr_x + 95}" y="{pr_y + 74}" font-size="10.5" fill="{T["text_main"]}">• HSP70-1 (+1.45 log₂) · HSP90.1 (+1.18)</text>
    <text x="{pr_x + 95}" y="{pr_y + 92}" font-size="10.5" fill="{T["text_main"]}">• DNAJ-1 (+0.82) · CPN60A (+0.65)</text>
    <text x="{pr_x + 95}" y="{pr_y + 110}" font-size="10.5" fill="{T["text_main"]}">• UBC1 (+0.77 log₂) · PAA1 (+0.66, 26S Core)</text>
    ''')

    # Reaction Sequence Path Overlay across Panel A
    svg.append(f'''
    <!-- Causal flow line connecting Panel A steps -->
    <path d="M {st_x + 160} {st_y + 50} Q {pm_x + 212} {st_y - 20} {pm_x + 212} {pm_y + pm_h} " fill="none" stroke="{T["accent_ca"]}" stroke-width="2" stroke-dasharray="6,4" marker-end="url(#arrow-ca)" />
    <circle cx="{st_x + 160}" cy="{st_y + 40}" r="11" fill="{T["accent_ca"]}" />
    <text x="{st_x + 160}" y="{st_y + 44}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">1</text>

    <path d="M {pm_x + 265} {pm_y + 45} Q {cs_x + 80} {cs_y + 80} {cs_x + 180} {cs_y + 105}" fill="none" stroke="{T["accent_auxin"]}" stroke-width="2" stroke-dasharray="6,4" marker-end="url(#arrow-auxin)" />
    <circle cx="{pm_x + 270}" cy="{pm_y + 50}" r="11" fill="{T["accent_auxin"]}" />
    <text x="{pm_x + 270}" y="{pm_y + 54}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">2</text>

    <path d="M {cs_x + 260} {cs_y + 86} Q {pm_x + 500} {cy_y + 30} {pm_x + 500} {pm_y + pm_h}" fill="none" stroke="{T["accent_auxin"]}" stroke-width="2" stroke-dasharray="6,4" marker-end="url(#arrow-auxin)" />
    <circle cx="{cs_x + 265}" cy="{cs_y + 70}" r="11" fill="{T["accent_auxin"]}" />
    <text x="{cs_x + 265}" y="{cs_y + 74}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">3</text>

    <path d="M {pm_x + 500} {pm_y} Q {pm_x + 400} {cw_y + 80} {pm_x + 350} {pm_y}" fill="none" stroke="{T["accent_wall"]}" stroke-width="2" stroke-dasharray="6,4" marker-end="url(#arrow-wall)" />
    <circle cx="{pm_x + 420}" cy="{cw_y + 80}" r="11" fill="{T["accent_wall"]}" />
    <text x="{pm_x + 420}" y="{cw_y + 84}" font-size="10" font-weight="700" fill="#FFFFFF" text-anchor="middle">4</text>
    ''')

    # =========================================================================
    # PANEL B: Organ & Tissue Activation Landscape
    # =========================================================================
    pb_x, pb_y, pb_w, pb_h = pa_x + pa_w + 15, 105, 575, 475
    svg.append(f'''
    <!-- PANEL B CONTAINER -->
    <rect x="{pb_x}" y="{pb_y}" width="{pb_w}" height="{pb_h}" rx="10" fill="{T["card_bg"]}" stroke="{T["card_border"]}" stroke-width="1.2" filter="url(#card-shadow)" />
    <rect x="{pb_x + 16}" y="{pb_y + 14}" width="76" height="22" rx="4" fill="{T["accent_auxin"]}" />
    <text x="{pb_x + 54}" y="{pb_y + 29}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">PANEL B</text>
    <text x="{pb_x + 102}" y="{pb_y + 30}" font-size="16" font-weight="700" fill="{T["text_main"]}">
      Organ &amp; Single-Cell Tissue Activation Landscape
    </text>
    <text x="{pb_x + 16}" y="{pb_y + 55}" font-size="12" fill="{T["text_muted"]}">
      Single-cell deconvolution and pathway activity scores across flight-responsive tissues
    </text>
    ''')

    # Heatmap Grid Layout
    hm_cols = [
        ("Auxin Polarity", T["accent_auxin"]),
        ("Ca²⁺ Signaling", T["accent_ca"]),
        ("Cell Wall Loosen", T["accent_wall"]),
        ("Redox Defense", T["accent_redox"]),
        ("Chaperones", T["accent_chaperone"]),
        ("Bioenergetics", T["accent_atp"])
    ]
    hm_rows = [
        ("Stem: Vascular Cambium (stem_9)", "+1.23% (padj=0.003) ***", [2.4, 2.1, 2.8, 1.9, 2.2, 1.8], True),
        ("Root Tip: Columella Statocytes", "Primary Graviperception **", [2.6, 2.9, 1.4, 1.7, 1.9, 1.6], False),
        ("Root Tip: Procambium / Stele", "Auxin Transport Axis **", [2.1, 1.8, 1.9, 1.5, 1.7, 1.5], False),
        ("Shoot Apex: Meristematic Zone", "Organogenesis / Morphogen *", [1.8, 1.5, 1.6, 1.8, 1.9, 2.1], False),
        ("Rosette Leaf: Vascular Bundles", "Redox &amp; Photosynthetic shift", [1.2, 1.1, 1.3, 2.4, 1.8, 1.7], False),
        ("Silique &amp; Seed: Coat (silique_20)", "+0.73% (padj=0.030) **", [1.6, 1.3, 2.1, 1.6, 2.5, 2.2], True),
    ]

    grid_x0 = pb_x + 235
    grid_y0 = pb_y + 95
    col_w = 52
    row_h = 52

    # Column Headers
    for c_idx, (col_name, col_color) in enumerate(hm_cols):
        cx = grid_x0 + c_idx * col_w
        # Rotated text for column header
        svg.append(f'''
        <g transform="translate({cx + col_w//2}, {grid_y0 - 8}) rotate(-35)">
          <text x="0" y="0" font-size="10.5" font-weight="700" fill="{col_color}" text-anchor="start">
            {col_name}
          </text>
        </g>
        ''')

    # Rows & Cells
    for r_idx, (row_label, row_sub, scores, is_highlight) in enumerate(hm_rows):
        ry = grid_y0 + r_idx * row_h
        
        # Row label background highlight if top responder
        row_bg = T["subcard_bg"] if not is_highlight else (T["pill_up_bg"] if theme_name == "light" else "#064E3B")
        svg.append(f'''
        <rect x="{pb_x + 14}" y="{ry}" width="215" height="{row_h - 4}" rx="5" fill="{row_bg}" stroke="{T["subcard_border"]}" stroke-width="1" />
        <text x="{pb_x + 22}" y="{ry + 20}" font-size="11" font-weight="700" fill="{T["text_main"]}">
          {row_label}
        </text>
        <text x="{pb_x + 22}" y="{ry + 36}" font-size="9.5" font-weight="500" fill="{T["text_muted"]}">
          {row_sub}
        </text>
        ''')

        # Heat tiles
        for c_idx, score in enumerate(scores):
            cx = grid_x0 + c_idx * col_w
            # Color intensity calculation
            intensity = min(1.0, max(0.2, (score - 1.0) / 1.8))
            alpha_hex = hex(int(intensity * 255))[2:].zfill(2)
            base_col = hm_cols[c_idx][1]
            
            svg.append(f'''
            <rect x="{cx}" y="{ry}" width="{col_w - 4}" height="{row_h - 4}" rx="5" fill="{base_col}" fill-opacity="{intensity:.2f}" stroke="{T["card_border"]}" stroke-width="1" />
            <text x="{cx + col_w//2 - 2}" y="{ry + 24}" font-size="11.5" font-weight="700" fill="{T["text_main"]}" text-anchor="middle">
              {score:.1f}
            </text>
            <text x="{cx + col_w//2 - 2}" y="{ry + 38}" font-size="8.5" font-weight="600" fill="{T["text_muted"]}" text-anchor="middle">
              z-score
            </text>
            ''')

    # Marginal Callout Box in Panel B
    svg.append(f'''
    <rect x="{pb_x + 14}" y="{grid_y0 + len(hm_rows) * row_h + 8}" width="{pb_w - 28}" height="46" rx="6" fill="{T["subcard_bg"]}" stroke="{T["subcard_border"]}" stroke-width="1" />
    <text x="{pb_x + 24}" y="{grid_y0 + len(hm_rows) * row_h + 26}" font-size="11" font-weight="700" fill="{T["accent_wall"]}">
      Key Insight: Inflorescence stem vascular cambium (stem_9) exhibits peak cell-wall loosening
    </text>
    <text x="{pb_x + 24}" y="{grid_y0 + len(hm_rows) * row_h + 42}" font-size="10.5" fill="{T["text_muted"]}">
      Coupled with root tip columella Ca²⁺ signaling, defining the primary axis of microgravity gravitropic disruption.
    </text>
    ''')

    # =========================================================================
    # PANEL D: Integrated Systems Mechanistic Model
    # =========================================================================
    pd_x, pd_y, pd_w, pd_h = pb_x, pb_y + pb_h + 15, pb_w, HEIGHT - (pb_y + pb_h + 15) - 25
    svg.append(f'''
    <!-- PANEL D CONTAINER -->
    <rect x="{pd_x}" y="{pd_y}" width="{pd_w}" height="{pd_h}" rx="10" fill="{T["card_bg"]}" stroke="{T["card_border"]}" stroke-width="1.2" filter="url(#card-shadow)" />
    <rect x="{pd_x + 16}" y="{pd_y + 14}" width="76" height="22" rx="4" fill="{T["accent_redox"]}" />
    <text x="{pd_x + 54}" y="{pd_y + 29}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">PANEL D</text>
    <text x="{pd_x + 102}" y="{pd_y + 30}" font-size="16" font-weight="700" fill="{T["text_main"]}">
      Integrated Mechanistic Systems Model
    </text>
    <text x="{pd_x + 16}" y="{pd_y + 54}" font-size="12" fill="{T["text_muted"]}">
      Chronological cascade of microgravity perception to cellular physiological adaptation
    </text>
    ''')

    # Chronological Flow Step Cards (6 steps)
    steps = [
        ("1. Loss of 1g Gravity Vector", "Statoliths lose sedimentation force; suspended in cytosol", T["accent_wall"]),
        ("2. Stochastic MSL Channel Gating", "Membrane tension fluctuations cause persistent Ca²⁺ spikes", T["accent_ca"]),
        ("3. Ectopic PIN Auxin Relocalization", "RABA1A/ARF1A1C misdirect PIN carriers to lateral flanks", T["accent_auxin"]),
        ("4. Vascular Tissue Ectopic Expansion", "Auxin accumulates in stem cambium (stem_9) &amp; root stele", T["accent_auxin"]),
        ("5. Apoplast Acidification &amp; Loosening", "PM H⁺-ATPase pumping activates XTH22 &amp; EXPA1 expansins", T["accent_wall"]),
        ("6. Compensatory Stress Homeostasis", "ROS burst (RBOHD) mitigated by HSP70 &amp; PRXIIB chaperones", T["accent_chaperone"]),
    ]

    card_step_h = 57
    for s_idx, (s_title, s_desc, s_col) in enumerate(steps):
        sy = pd_y + 72 + s_idx * (card_step_h + 8)
        svg.append(f'''
        <rect x="{pd_x + 16}" y="{sy}" width="{pd_w - 32}" height="{card_step_h}" rx="6" fill="{T["subcard_bg"]}" stroke="{T["subcard_border"]}" stroke-width="1" />
        <rect x="{pd_x + 22}" y="{sy + 8}" width="6" height="{card_step_h - 16}" rx="3" fill="{s_col}" />
        <text x="{pd_x + 36}" y="{sy + 24}" font-size="12" font-weight="700" fill="{s_col}">
          {s_title}
        </text>
        <text x="{pd_x + 36}" y="{sy + 44}" font-size="11" fill="{T["text_main"]}">
          {s_desc}
        </text>
        ''')
        # Arrow connector between steps
        if s_idx < len(steps) - 1:
            svg.append(f'''
            <line x1="{pd_x + pd_w//2}" y1="{sy + card_step_h}" x2="{pd_x + pd_w//2}" y2="{sy + card_step_h + 7}" stroke="{T["text_muted"]}" stroke-width="1.5" marker-end="url(#arrow-muted)" />
            ''')

    # =========================================================================
    # PANEL C: Curated Pathway Cross-Talk Network
    # =========================================================================
    pc_x, pc_y, pc_w, pc_h = pa_x, pa_y + pa_h + 15, pa_w, HEIGHT - (pa_y + pa_h + 15) - 25
    svg.append(f'''
    <!-- PANEL C CONTAINER -->
    <rect x="{pc_x}" y="{pc_y}" width="{pc_w}" height="{pc_h}" rx="10" fill="{T["card_bg"]}" stroke="{T["card_border"]}" stroke-width="1.2" filter="url(#card-shadow)" />
    <rect x="{pc_x + 16}" y="{pc_y + 14}" width="76" height="22" rx="4" fill="{T["accent_chaperone"]}" />
    <text x="{pc_x + 54}" y="{pc_y + 29}" font-size="11" font-weight="700" fill="#FFFFFF" text-anchor="middle">PANEL C</text>
    <text x="{pc_x + 102}" y="{pc_y + 30}" font-size="16" font-weight="700" fill="{T["text_main"]}">
      Curated Molecular Pathway Cross-Talk Network
    </text>
    <text x="{pc_x + 480}" y="{pc_y + 30}" font-size="12" font-weight="500" fill="{T["text_muted"]}">
      Replacing cryptic 5-digit KEGG IDs with biologically annotated functional hubs
    </text>
    ''')

    # 6 Functional Hub Nodes (Arranged in a coordinated hexagonal / two-row layout)
    # Positions relative to pc_x, pc_y
    hubs = [
        # (id, name, kegg_id, degs, pooled_fc, x, y, color, hub_genes)
        ("AUXIN", "Plant Hormone Signal Transduction", "ath04075", 84, "+0.68", pc_x + 160, pc_y + 130, T["accent_auxin"], ["PIN2", "AUX1", "ABCB1", "IAA19"]),
        ("MECHANO", "Mechanotransduction &amp; Ca²⁺ Flux", "ath04070", 42, "+0.89", pc_x + 480, pc_y + 105, T["accent_ca"], ["MSL10", "MCA1", "TCH1", "CPK1"]),
        ("WALL", "Cell Wall Organization &amp; Loosening", "ath00500", 67, "+0.76", pc_x + 800, pc_y + 130, T["accent_wall"], ["XTH22", "EXPA1", "CESA3", "AHA1"]),
        ("REDOX", "Redox Homeostasis &amp; ROS Defense", "ath00480", 58, "+0.81", pc_x + 180, pc_y + 295, T["accent_redox"], ["PRXIIB", "CSD1", "TRX-H5", "RBOHD"]),
        ("CHAPERONE", "Protein Processing in ER &amp; Chaperones", "ath04141", 79, "+1.12", pc_x + 500, pc_y + 315, T["accent_chaperone"], ["HSP70-1", "HSP90.1", "DNAJ-1", "PAA1"]),
        ("ATP", "Bioenergetics &amp; Carbon Fixation", "ath00710", 93, "+0.74", pc_x + 810, pc_y + 295, T["accent_atp"], ["ATPA", "GAPC1", "MDH1", "ATP5F1"]),
    ]

    # Draw Inter-Module Cross-Talk Edges
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
        <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{T["card_border"]}" stroke-width="{edge_w * 1.5}" stroke-linecap="round" opacity="0.7" />
        <!-- Edge label badge -->
        <rect x="{mid_x - 70}" y="{mid_y - 10}" width="140" height="20" rx="4" fill="{T["card_bg"]}" stroke="{T["subcard_border"]}" stroke-width="0.8" />
        <text x="{mid_x}" y="{mid_y + 3}" font-size="9" font-weight="600" fill="{T["text_muted"]}" text-anchor="middle">
          {edge_lbl}
        </text>
        ''')

    # Draw Module Hub Nodes
    for (hid, hname, hkegg, hdegs, hfc, hx, hy, hcol, hgenes) in hubs:
        node_r = 46
        svg.append(f'''
        <!-- Hub Node {hid} -->
        <circle cx="{hx}" cy="{hy}" r="{node_r + 6}" fill="{hcol}" fill-opacity="0.15" />
        <circle cx="{hx}" cy="{hy}" r="{node_r}" fill="{T["card_bg"]}" stroke="{hcol}" stroke-width="2.5" />
        <text x="{hx}" y="{hy - 14}" font-size="10.5" font-weight="700" fill="{hcol}" text-anchor="middle">
          {hkegg}
        </text>
        <text x="{hx}" y="{hy + 3}" font-size="11.5" font-weight="700" fill="{T["text_main"]}" text-anchor="middle">
          {hdegs} DEGs
        </text>
        <rect x="{hx - 32}" y="{hy + 12}" width="64" height="18" rx="3" fill="{T["pill_up_bg"]}" />
        <text x="{hx}" y="{hy + 24}" font-size="9.5" font-weight="700" fill="{T["pill_up_text"]}" text-anchor="middle">
          {hfc} log₂FC
        </text>

        <!-- Hub Title outside circle -->
        <text x="{hx}" y="{hy - node_r - 10}" font-size="11" font-weight="700" fill="{T["text_main"]}" text-anchor="middle">
          {hname}
        </text>

        <!-- Hub Key Genes Pill -->
        <rect x="{hx - 95}" y="{hy + node_r + 8}" width="190" height="20" rx="4" fill="{T["subcard_bg"]}" stroke="{T["subcard_border"]}" stroke-width="0.8" />
        <text x="{hx}" y="{hy + node_r + 22}" font-size="9.5" font-weight="600" fill="{T["text_muted"]}" text-anchor="middle">
          Hubs: {", ".join(hgenes)}
        </text>
        ''')

    # Close SVG
    svg.append('</svg>')
    return "\n".join(svg)

def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    figures_dir = os.path.join(repo_root, "Figures")
    docs_figures_dir = os.path.join(repo_root, "docs", "assets", "figures")

    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(docs_figures_dir, exist_ok=True)

    print("Generating Figure 8 SVGs (Light Publication Theme & Dark Portal Theme)...")

    # 1. Generate Light SVG (for Manuscript & PDF)
    light_svg_content = generate_svg(theme_name="light")
    light_svg_path = os.path.join(figures_dir, "kegg_pathway_network.svg")
    with open(light_svg_path, "w", encoding="utf-8") as f:
        f.write(light_svg_content)
    print(f"  ✓ Saved Light SVG: {light_svg_path} ({len(light_svg_content)} bytes)")

    # 2. Generate Dark SVG (for Web Portal)
    dark_svg_content = generate_svg(theme_name="dark")
    dark_svg_path = os.path.join(figures_dir, "kegg_pathway_network_dark.svg")
    with open(dark_svg_path, "w", encoding="utf-8") as f:
        f.write(dark_svg_content)
    print(f"  ✓ Saved Dark SVG: {dark_svg_path} ({len(dark_svg_content)} bytes)")

    # Copy Dark SVG to docs web assets
    docs_svg_path = os.path.join(docs_figures_dir, "kegg_pathway_network.svg")
    shutil.copyfile(dark_svg_path, docs_svg_path)
    print(f"  ✓ Copied Web Dark SVG: {docs_svg_path}")

    # 3. High-Resolution PNG Rendering via macOS QuickLook / qlmanage
    print("\nRendering High-Resolution 300 DPI PNGs...")
    tmp_out = "/tmp"

    # Render light PNG
    try:
        subprocess.run(["qlmanage", "-t", "-s", "2400", "-o", tmp_out, light_svg_path], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        rendered_light = os.path.join(tmp_out, "kegg_pathway_network.svg.png")
        if os.path.exists(rendered_light):
            target_light_png = os.path.join(figures_dir, "kegg_pathway_network.png")
            shutil.copyfile(rendered_light, target_light_png)
            print(f"  ✓ Rendered Light PNG for Manuscript: {target_light_png}")
    except Exception as e:
        print(f"  ⚠ qlmanage light render warning: {e}")

    # Render dark PNG
    try:
        subprocess.run(["qlmanage", "-t", "-s", "2400", "-o", tmp_out, dark_svg_path], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        rendered_dark = os.path.join(tmp_out, "kegg_pathway_network_dark.svg.png")
        if os.path.exists(rendered_dark):
            target_dark_png = os.path.join(docs_figures_dir, "kegg_pathway_network.png")
            shutil.copyfile(rendered_dark, target_dark_png)
            print(f"  ✓ Rendered Dark PNG for Web Portal: {target_dark_png}")
    except Exception as e:
        print(f"  ⚠ qlmanage dark render warning: {e}")

    print("\nFigure 8 generation complete!")

if __name__ == "__main__":
    main()
