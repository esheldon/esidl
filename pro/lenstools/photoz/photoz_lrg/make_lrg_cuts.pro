PRO make_lrg_cuts, scat, lrgflags, nostar=nostar

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: make_lrg_cuts, scat, lrgflags, /nostar'
      return
  ENDIF 

  ;; This is a modified version of Ryan's make_lrg_flags

  n_cat = n_elements(scat)
  lrgflags = lonarr(n_cat)

  base = 2L

  OBJC_GAL = base^0
  LRG_GAL = base^1
  LRG = base^2
  LRG_HI_DEEP = base^3
  LRG_LO_DEEP = base^4
  LRG_HI_MED = base^5
  LRG_LO_MED = base^6
  LRG_SPECTRO = base^7

  gmr = scat.grmodel
  rmi = scat.rimodel

  ;; this should be cmodel actually
  rpetro = scat.rpetro
  ipetro = scat.ipetro

  C = 0.7
  cparr = C*(gmr) + (1.0-C)*4.0*(rmi - 0.177)
  cperp = rmi - gmr/4.0 - 0.177
  dperp = rmi - gmr/8.0

  cparr_min = 1.6
  dperp_min = 0.2
  dperp_max = 1.0
  gmr_cut = 1.0

  IF keyword_set(nostar) THEN BEGIN 
      stellar_max = 15.0/11.0*dperp + (2.35 - 15.0/11.0)
      stellar_min = 9.0/7.0*dperp + (2.05 - 9.0/7.0)
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; This selects LRG's according to scranton et al. ISW
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  w_lrg =     where(cparr GT cparr_min AND $
                    gmr GT gmr_cut AND $
                    ipetro LT 21.0 AND $
                    dperp GE dperp_min AND $
                    dperp LT dperp_max, n_lrg)

  w_lo_med =  where(cparr GT cparr_min AND $
                    gmr GT gmr_cut AND $
                    dperp GT dperp_min AND $
                    dperp LT dperp_min + 0.2, n_lo_med)

  w_hi_med =  where(cparr GT cparr_min AND $
                    gmr GT gmr_cut AND $
                    dperp GT dperp_min + 0.2 AND $
                    dperp LT dperp_min + 0.4, n_hi_med)

  w_lo_deep = where(cparr GT cparr_min AND $
                    gmr GT gmr_cut AND $
                    dperp GT dperp_min + 0.4 AND $
                    dperp LT dperp_min + 0.6, n_lo_deep)

  w_hi_deep = where(cparr GT cparr_min AND $
                    gmr GT gmr_cut AND $
                    dperp GT dperp_min + 0.6 AND $
                    dperp LT dperp_max, n_hi_deep)


  w_spectro = where(cperp LT 0.2 AND $
                    rpetro LT (13.1 + cparr/0.3), n_spectro)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; extra cuts to remove stars from the sample
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(nostar) THEN BEGIN 
      IF n_lrg NE 0 THEN BEGIN 
          w2_lrg = where(cparr[w_lrg] LT stellar_min OR $
                         cparr[w_lrg] GT stellar_max, n_lrg)
          IF n_lrg NE 0 THEN w_lrg=w_lrg[w2_lrg]
      ENDIF 

      IF n_lo_med NE 0 THEN BEGIN 
          w2_lo_med = where(cparr[w_lo_med] LT stellar_min OR $
                            cparr[w_lo_med] GT stellar_max, n_lo_med)
          IF n_lo_med NE 0 THEN w_lo_med=w_lo_med[w2_lo_med]
      ENDIF 

      IF n_hi_med NE 0 THEN BEGIN 
          w2_hi_med = where(cparr[w_hi_med] LT stellar_min OR $
                            cparr[w_hi_med] GT stellar_max, n_hi_med)
          IF n_hi_med NE 0 THEN w_hi_med=w_hi_med[w2_hi_med]
      ENDIF 

      IF n_lo_deep NE 0 THEN BEGIN 
          w2_lo_deep = where(cparr[w_lo_deep] LT stellar_min OR $
                             cparr[w_lo_deep] GT stellar_max, n_lo_deep)
          IF n_lo_deep NE 0 THEN w_lo_deep=w_lo_deep[w2_lo_deep]
      ENDIF 
      
      IF n_hi_deep NE 0 THEN BEGIN 
          w2_hi_deep = where(cparr[w_hi_deep] LT stellar_min OR $
                             cparr[w_hi_deep] GT stellar_max, n_hi_deep)
          IF n_hi_deep NE 0 THEN w_hi_deep=w_hi_deep[w2_hi_deep]
      ENDIF 

  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set the bits
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_lrg GT 0 THEN BEGIN
      lrgflags[w_lrg] = lrgflags[w_lrg] + LRG
  ENDIF

  IF n_lo_med GT 0 THEN BEGIN
      lrgflags[w_lo_med] = lrgflags[w_lo_med] + LRG_LO_MED
  ENDIF

  IF n_hi_med GT 0 THEN BEGIN
      lrgflags[w_hi_med] = lrgflags[w_hi_med] + LRG_HI_MED
  ENDIF

  IF n_lo_deep GT 0 THEN BEGIN
      lrgflags[w_lo_deep] = lrgflags[w_lo_deep] + LRG_LO_DEEP
  ENDIF

  IF n_hi_deep GT 0 THEN BEGIN
      lrgflags[w_hi_deep] = lrgflags[w_hi_deep] + LRG_HI_DEEP
  ENDIF

  IF n_spectro GT 0 THEN BEGIN
      lrgflags[w_spectro] = lrgflags[w_spectro] + LRG_SPECTRO
  ENDIF

END 
