PRO photoz_select_struct, struct, useind

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: photoz_select_struct, struct, useind'
      return
  ENDIF 

  tzmin = 0.02                  ;0.02 problems at low z
  tzmax = 0.8                   ;because a default is around 1
  qualmin = 0
  qualmax = 12
  zerrmin = 0.01
  zerrmax = 0.4

  useind=where(struct.photoz_z GT tzmin AND $
               struct.photoz_z LT tzmax AND $
               struct.photoz_z NE 0.0 AND $
               struct.objc_prob_psf LT (1.-!hardprobcut) AND $
               struct.photoz_zerr GT zerrmin AND $
               struct.photoz_zerr LT zerrmax AND $
               struct.photoz_quality GT qualmin AND $
               struct.photoz_quality LT qualmax, nuse)


END 
