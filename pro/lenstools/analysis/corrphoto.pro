PRO corrphoto, struct, clrs, e1corr, e2corr, corr

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: corrphoto, struct, clrs, e1corr, e2corr, corr, dilution=dilution'
      return
  ENDIF 

  getphotosmear, struct, clrs, corr

  e1corr = struct.m_e1[clrs] - corr*struct.m_e1_psf[clrs]
  e2corr = struct.m_e2[clrs] - corr*struct.m_e2_psf[clrs]

END 
