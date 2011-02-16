PRO make_stripe_mask, stripe, regress=regress, tsgals=tsgals, outdir=outdir, indir=indir, plotscat=plotscat

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; THIS ROUTINE IS OBSOLETE. USE MAKE_MULTISTRIPE_MASK INSTEAD
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: make_stripe_mask, stripe, regress=regress, tsgals=tsgals, outdir=outdir, indir=indir, plotscat=plotscat'
      return
  ENDIF 

  sdssidl_setup,/silent
  setup_mystuff
  IF n_elements(outdir) EQ 0 THEN $
    outdir = sdssidl_config('SHAPECORR_DIR') + 'masks/'
  IF keyword_set(regress) THEN BEGIN 
      outfile = outdir + 'mask-spectra-regress-stripe'+ntostr(long(stripe))+'.sav'
  ENDIF ELSE IF keyword_set(tsgals) THEN BEGIN
      outfile = outdir + 'mask-tsgal-stripe'+ntostr(long(stripe))+'.sav'
  ENDIF ELSE BEGIN  
      outfile = outdir + 'mask-spectra-stripe'+ntostr(long(stripe))+'.sav'
  ENDELSE 

  IF keyword_set(tsgals) THEN BEGIN 
      get_tsgals, stripe, lcat, indir=indir
  ENDIF ELSE BEGIN 
      get_spectra_lcat, stripe, lcat, indir=indir, /spec
      IF stripe EQ 82 THEN BEGIN 
          make_tsflag_struct, tsf
          tsf.galaxy='Y'
          tsflag_select, lcat, tsf, keep
          lcat=lcat[keep]
      ENDIF 
  ENDELSE 

  IF tag_exist(lcat,'clambda') THEN BEGIN 
      lambda = lcat.clambda
      eta = lcat.ceta
  ENDIF ELSE BEGIN 
      eq2survey, lcat.ra, lcat.dec, lambda, eta
  ENDELSE 

  IF keyword_set(plotscat) THEN BEGIN 
      get_scat, stripe, [1,2,3], scat
      ox = scat.clambda
      oy = scat.ceta
  ENDIF 

  make_mask, lambda, eta, mask, ox=ox, oy=oy

  apply_mask, mask, lambda, bad, good

  oplot, lambda[good], eta[good], color=!green, psym=8

  print,'saving mask file: ',outfile
  save, mask, filename=outfile, /verbose

  setzero, lcat

END
