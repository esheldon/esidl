PRO getphotosmear, struct, clrs, corr

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: getphotocorr, struct, clrs, corr'
      return
  ENDIF 

  corr = struct.m_rr_cc_psf[clrs]/struct.m_rr_cc[clrs]*(4./struct.m_cr4_psf[clrs] -1.)/(4./struct.m_cr4[clrs]-1.)

END 
