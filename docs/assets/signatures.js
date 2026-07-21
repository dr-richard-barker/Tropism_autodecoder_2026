/*
 * signatures.js — curated tropism signature gene sets for the Tropism Autodecoder web tool.
 *
 * Phase 1 signatures. Every AGI locus below was verified against TAIR / primary
 * literature before inclusion — no invented identifiers. To swap in signatures
 * exported from the manuscript pipeline, replace TROPISM_SIGNATURES with your own
 * { name, color, genes:[AGI...], note } objects; the app code needs no changes.
 *
 * `genes` are AGI locus IDs (canonical key). `SYMBOL_TO_AGI` lets uploaded matrices
 * that use gene symbols resolve to the same loci.
 */

const TROPISM_SIGNATURES = [
  {
    key: "gravitropism",
    name: "Gravitropism",
    subtitle: "gravity sensing / auxin transport",
    genes: ["AT5G57090", "AT1G70940", "AT1G23080", "AT2G38120", "AT5G14090"],
    members: "PIN2, PIN3, PIN7, AUX1, LAZY1",
  },
  {
    key: "phototropism",
    name: "Phototropism",
    subtitle: "blue-light directed growth",
    genes: ["AT3G45780", "AT5G58140", "AT5G64330", "AT2G30520", "AT2G02950"],
    members: "PHOT1, PHOT2, NPH3, RPT2, PKS1",
  },
  {
    key: "thigmotropism",
    name: "Thigmotropism",
    subtitle: "touch / mechanostimulation",
    genes: ["AT5G37780", "AT5G37770", "AT2G41100", "AT5G57560", "AT5G12080"],
    members: "TCH1, TCH2, TCH3, TCH4, MSL10",
  },
  {
    key: "hydrotropism",
    name: "Hydrotropism",
    subtitle: "water-seeking (low-coverage set)",
    genes: ["AT2G41660", "AT1G13980"],
    members: "MIZ1, GNOM/MIZ2",
    lowCoverage: true,
  },
];

// Symbol -> AGI aliases so symbol-labelled matrices still match.
const SYMBOL_TO_AGI = {
  // gravitropism
  PIN2: "AT5G57090", EIR1: "AT5G57090", AGR1: "AT5G57090",
  PIN3: "AT1G70940",
  PIN7: "AT1G23080",
  AUX1: "AT2G38120",
  LAZY1: "AT5G14090", LZY1: "AT5G14090",
  // phototropism
  PHOT1: "AT3G45780", NPH1: "AT3G45780",
  PHOT2: "AT5G58140", NPL1: "AT5G58140",
  NPH3: "AT5G64330",
  RPT2: "AT2G30520",
  PKS1: "AT2G02950",
  // thigmotropism
  TCH1: "AT5G37780", CAM2: "AT5G37780",
  TCH2: "AT5G37770", CML24: "AT5G37770",
  TCH3: "AT2G41100", CML12: "AT2G41100",
  TCH4: "AT5G57560", XTH22: "AT5G57560",
  MSL10: "AT5G12080",
  // hydrotropism
  MIZ1: "AT2G41660",
  GNOM: "AT1G13980", MIZ2: "AT1G13980", EMB30: "AT1G13980",
};

// RdBu diverging endpoints (ColorBrewer). Used for the accessible red-white-blue scale.
const RDBU = {
  redDark: "#67001f", red: "#b2182b", redLight: "#ef8a62", redPale: "#fddbc7",
  white: "#f7f7f7",
  bluePale: "#d1e5f0", blueLight: "#67a9cf", blue: "#2166ac", blueDark: "#053061",
};

if (typeof module !== "undefined") {
  module.exports = { TROPISM_SIGNATURES, SYMBOL_TO_AGI, RDBU };
}
