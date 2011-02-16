PRO read_stripe_mask, stripe, mask, regress=regress, tsgals=tsgals, indir=indir

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: read_stripe_mask, stripe, mask, regress=regress, tsgals=tsgals, indir=indir'
      return
  ENDIF 

  stripe_string = stripearr2string(stripe)

  sdssidl_setup,/silent
  ;setup_mystuff
  IF n_elements(indir) EQ 0 THEN $
    indir = sdssidl_config('SHAPECORR_DIR') + 'masks/'
  IF NOT keyword_set(regress) AND NOT keyword_set(tsgals) THEN BEGIN 
      infile = indir + 'mask-spectra-stripe'+stripe_string+'.sav'
  ENDIF ELSE IF keyword_set(tsgals) THEN BEGIN
      infile = indir + 'mask-tsgal-stripe'+stripe_string+'.sav'
  ENDIF ELSE BEGIN 
      infile = indir + 'mask-spectra-regress-stripe'+stripe_string+'.sav'
  ENDELSE 

  print
  print,'Reading mask file: ',infile
  IF fexist(infile) THEN restore, infile $
  ELSE message,'File does not exist: '+infile

END 
