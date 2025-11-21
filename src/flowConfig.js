// src/flowConfig.js

// 公共样式：让每个节点大一点、字体大一点
const nodeStyle = {
  width: 260,
  height: 70,
  fontSize: 20,
  fontWeight: "bold",
  display: "flex",
  justifyContent: "center",
  alignItems: "center",
  border: "1px solid #555",
  borderRadius: 6,
  background: "#2d2d2d",    // 节点暗色
  color: "#eaeaea",         // 白色字体
};

export const nodes = [
  { id: "chip_pre", position: { x: -400, y: -450 }, data: { label: "pre_process" }, style: nodeStyle },
  { id: "astro",    position: { x: -400, y: -300 }, data: { label: "proc_astrometry" }, style: nodeStyle },
  { id: "source",   position: { x: -400, y: -150 }, data: { label: "proc_source" }, style: nodeStyle },
  { id: "fourier",  position: { x: -400, y: 0 },    data: { label: "proc_Fourier_T" }, style: nodeStyle },
  { id: "psf",      position: { x: -400, y: 150 },  data: { label: "PSF_reconstruction" }, style: nodeStyle },
  { id: "shear",    position: { x: -400, y: 300 },  data: { label: "gen_shear_cat" }, style: nodeStyle },
  { id: "combine",  position: { x: -400, y: 450 },  data: { label: "combine_expo_catalog" }, style: nodeStyle },
];

export const mainTitleMap = {
  chip_pre: "pre_process",
  astro: "proc_astrometry",
  source: "proc_source",
  fourier: "proc_Fourier_T",
  psf: "PSF_reconstruction",
  shear: "gen_shear_cat",
  combine: "combine_expo_catalog",
};

export const edges = [
  { id: "e1", source: "chip_pre", target: "astro" },
  { id: "e2", source: "astro",    target: "source" },
  { id: "e3", source: "source",   target: "fourier" },
  { id: "e4", source: "fourier",  target: "psf" },
  { id: "e5", source: "psf",      target: "shear" },
  { id: "e6", source: "shear",    target: "combine" },
];

const subNodeStyle = {
  width: 260,
  height: 55,
  fontSize: 18,
  fontWeight: "bold",
  display: "flex",
  justifyContent: "center",
  alignItems: "center",
  background: "#2d2d2d",
  color: "#eaeaea",
  border: "1px solid #666",
  borderRadius: 6,
};

export const subFlowMap = {
  chip_pre: {
    nodes: [
      { id: "Set_Background",        position: { x: 0, y: -150 }, data:{ label: "Set_Background" },        style: subNodeStyle },
      { id: "Star_match",            position: { x: 0, y: -60 },  data:{ label: "Star_match" },            style: subNodeStyle },
      { id: "Locate_defects",        position: { x: 0, y: 30 },   data:{ label: "Locate_defects" },        style: subNodeStyle },
      { id: "Merge_defects",         position: { x: 0, y: 120 },  data:{ label: "Merge_defects" },         style: subNodeStyle },
      { id: "Gen_norm",              position: { x: 0, y: 210 },  data:{ label: "Gen_norm.fits" },         style: subNodeStyle },
    ],
    edges: [
      { id: "p1", source: "Set_Background",  target: "Star_match" },
      { id: "p2", source: "Star_match",      target: "Locate_defects" },
      { id: "p3", source: "Locate_defects",  target: "Merge_defects" },
      { id: "p4", source: "Merge_defects",   target: "Gen_norm" },
    ],
  },

  astro: {
    nodes: [
      { id: "get_astrometry",  position: { x: 0, y: -80 }, data:{ label: "get_astrometry" },  style: subNodeStyle },
      { id: "update_norm",     position: { x: 0, y: 20 },  data:{ label: "update_norm" },     style: subNodeStyle },
    ],
    edges: [
      { id: "a1", source: "get_astrometry", target: "update_norm" },
    ],
  },

  source: {
    nodes: [
      { id: "get_weight",               position: { x: 0, y: -150 }, data:{ label: "get_weight" },               style: subNodeStyle },
      { id: "get_expo_catalog",         position: { x: 0, y: -60 },  data:{ label: "get_expo_catalog" },         style: subNodeStyle },
      { id: "gen_source_ext_catalog",   position: { x: 0, y: 30 },   data:{ label: "gen_source_ext_catalog" },   style: subNodeStyle },
      { id: "gen_star_candidate_direct",position: { x: 0, y: 120 },  data:{ label: "gen_star_candidate_direct" },style: subNodeStyle },
    ],
    edges: [
      { id: "s1", source: "get_weight",             target: "get_expo_catalog" },
      { id: "s2", source: "get_expo_catalog",       target: "gen_source_ext_catalog" },
      { id: "s3", source: "gen_source_ext_catalog", target: "gen_star_candidate_direct" },
    ],
  },

  fourier: {
    nodes: [
      { id: "get_gal_power",  position: { x: 0, y: -60 }, data:{ label: "get_gal_power" },  style: subNodeStyle },
      { id: "get_star_power", position: { x: 0, y: 40 },  data:{ label: "get_star_power" }, style: subNodeStyle },
    ],
    edges: [
      { id: "f1", source: "get_gal_power", target: "get_star_power" },
    ],
  },

  psf: {
    nodes: [
      { id: "read_in_candidates", position: { x: 0, y: -180 }, data:{ label: "read_in_candidates" }, style: subNodeStyle },
      { id: "star_selection",     position: { x: 0, y: -90 },  data:{ label: "star_selection" },     style: subNodeStyle },
      { id: "plot_star_expo",     position: { x: 0, y: 0 },    data:{ label: "plot_star_expo" },     style: subNodeStyle },
      { id: "plot_stars",         position: { x: 0, y: 90 },   data:{ label: "plot_stars" },         style: subNodeStyle },
      { id: "make_PSF_local_fit", position: { x: 0, y: 180 },  data:{ label: "make_PSF_local_fit" }, style: subNodeStyle },
    ],
    edges: [
      { id: "ps1", source: "read_in_candidates", target: "star_selection" },
      { id: "ps2", source: "star_selection",     target: "plot_star_expo" },
      { id: "ps3", source: "plot_star_expo",     target: "plot_stars" },
      { id: "ps4", source: "plot_stars",         target: "make_PSF_local_fit" },
    ],
  },

  shear: {
    nodes: [
      { id: "Load_PSF",            position: { x: 0, y: -150 }, data:{ label: "Load_PSF" },            style: subNodeStyle },
      { id: "read_astrometry_para",position: { x: 0, y: -60 },  data:{ label: "read_astrometry_para" },style: subNodeStyle },
      { id: "read_source_info",    position: { x: 0, y: 30 },   data:{ label: "read_source_info" },    style: subNodeStyle },
      { id: "field_distortion_PU", position: { x: 0, y: 120 },  data:{ label: "field_distortion_PU" }, style: subNodeStyle },
      { id: "get_shear",           position: { x: 0, y: 210 },  data:{ label: "get_shear" },           style: subNodeStyle },
    ],
    edges: [
      { id: "sh1", source: "Load_PSF",            target: "read_astrometry_para" },
      { id: "sh2", source: "read_astrometry_para",target: "read_source_info" },
      { id: "sh3", source: "read_source_info",    target: "field_distortion_PU" },
      { id: "sh4", source: "field_distortion_PU", target: "get_shear" },
    ],
  },

  combine: {
    nodes: [
      { id: "read_sheardat",        position: { x: 0, y: -60 }, data: { label: "read_sheardat" },        style: subNodeStyle },
      { id: "Filtering_Combining",  position: { x: 0, y: 40 },  data: { label: "Filtering_Combining" },  style: subNodeStyle },
    ],
    edges: [
      { id: "c1", source: "read_sheardat", target: "Filtering_Combining" },
    ],
  },
};
