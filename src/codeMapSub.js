const codeMapSub = {

  Set_Background: { start: 79, end: 105 },
  Star_match: { start: 108, end: 119 },
  Locate_defects: { start: 120, end: 121 },
  Merge_defects: { start: 122, end: 123 },
  Gen_norm: { start: 127, end: 172 },

  get_astrometry: { start: 14, end: 21 },
  update_norm: { start: 22, end: 47 },

  get_weight: { start: 35, end: 77 },
  get_expo_catalog: { start: 80, end: 81 },
  gen_source_ext_catalog: { start: 99, end: 110 },
  gen_star_candidate_direct:{ start: 112, end: 113 },

  get_gal_power: { start: 31, end: 107 },
  get_star_power: { start: 112, end: 173 },

  read_in_candidates: { start: 34, end: 165 },
  star_selection: { start: 167, end: 337 },
  plot_star_expo: { start: 339, end: 404 },
  plot_stars: { start: 406, end: 490 },
  make_PSF_local_fit: { start: 492, end: 585 },
  
  Load_PSF: { start: 47, end: 61 },
  read_astrometry_para: { start: 72, end: 73 },
    read_source_info: { start: 78, end: 108 },
    field_distortion_PU: { start: 140, end: 147 },
  get_shear:{ start: 148, end: 158 },

    read_sheardat: { start: 17, end: 58 },
  Filtering_Combining: { start: 73, end: 99 },
};

export default codeMapSub;
