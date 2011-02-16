PRO make_multistripe_mask, stripes, outdir=outdir, indir=indir, plotscat=plotscat, scat=scat, lcat=lcat, noprimary_bound_cut=noprimary_bound_cut

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: make_multistripe_mask, stripes, outdir=outdir, indir=indir, plotscat=plotscat, scat=scat, lcat=lcat, /noprimary_bound_cut'
      return
  ENDIF 

  IF n_elements(stripes) EQ 1 THEN BEGIN 
      IF ((stripes[0] EQ 10) OR (stripes[0] EQ 76) OR $
          (stripes[0] EQ 82) OR (stripes[0] EQ 86)) THEN BEGIN 
          noprimary_bound_cut = 1
      ENDIF
  ENDIF 

  stripe_string = stripearr2string(stripes)

  sdssidl_setup,/silent
  setup_mystuff
  IF n_elements(outdir) EQ 0 THEN $
    outdir = sdssidl_config('SHAPECORR_DIR') + 'masks/'
  outfile = outdir + 'mask-spectra-stripe'+stripe_string+'.sav'
  print,'makeing mask file: ',outfile

  IF n_elements(lcat) EQ 0 THEN BEGIN 
      get_spectra_lcat, stripes, lcat, indir=indir, /spec
  ENDIF 

  make_tsflag_struct, tsf
  tsf.galaxy='Y'
  tsflag_select, lcat, tsf, keep
  lcat=lcat[keep]

  IF tag_exist(lcat,'clambda') THEN BEGIN 
      lambda = lcat.clambda
      eta = lcat.ceta
  ENDIF ELSE BEGIN 
      eq2survey, lcat.ra, lcat.dec, lambda, eta
  ENDELSE 

  IF keyword_set(plotscat) THEN BEGIN 
      IF n_elements(scat) EQ 0 THEN BEGIN 
          get_scat, stripes, [1,2,3], scat
      ENDIF 
      ox = scat.clambda
      oy = scat.ceta
  ENDIF 

  ;; get primary bounds
  primary_bound_multi, stripes, bound_arr

  make_mask_sdss2d, lambda, eta, bound_arr, mask, ox=ox, oy=oy, $
                    noprimary_bound_cut=noprimary_bound_cut

  apply_mask2d, mask, lambda, eta, bad, good

  oplot, lambda[good], eta[good], color=!green, psym=8
  print,ntostr(n_elements(bad))+' thrown out by mask'

  print,'saving mask file: ',outfile
  save, mask, filename=outfile, /verbose

END
