PRO vagc_inwindow, lss=lss, letter=letter, post=post, sample=sample

  tt = systime(1)

  sdss_shapecorr_dir = sdssidl_config('shapecorr_dir')

  maskDir = sdss_shapecorr_dir + 'masks/'
  maskFile = maskDir + 'sphpoly_mask_sample14.dat'

  catname = vagc_catname(lss=lss, letter=letter, post=post, sample=sample)

  lssDir = sdssidl_config('lss_dir') + $
    sample + "/" + letter + "/" + post + "/"
  lssFile = 'post_catalog.'+sample +letter+post+'.fits'
  lssWindowFile = repstr(lssFile, '.fits', '_inwindow.fit')

  lssFile = lssDir + lssFile
  lssWindowFile = maskDir + lssWindowFile

  print
  print,'lssFile: ',lssFile
  print,'Output lssWindowFile: ',lssWindowFile

  print
  lss = mrdfits(lssFile, 1)
  nlss = n_elements(lss)

  ss = create_struct('completeness', 0.0, $
                     'poly_id', 0L, $
                     'poly_area', 0d)
  window = replicate(ss, nlss)

  apply_sphpoly_mask, lss.ra, lss.dec, masked, unmasked, $
    maskFile=maskFile, $
    completeness=completeness, poly_id=poly_id, poly_area=poly_area

  window.completeness = completeness
  window.poly_id = poly_id
  window.poly_area = poly_area

  print
  print,'Writing window file: ',lssWindowFile
  mwrfits, window, lssWindowFile, /create
  
return


  vagc_getstripes, stripes, nstripe, $
    lss=lss, letter=letter, post=post, sample=sample
  catname = vagc_catname(lss=lss, letter=letter, post=post, sample=sample)
  maskdir = sdss_shapecorr_dir + 'masks/'
  
  
  ss = create_struct('completeness', 0.0, $
                     'poly_id', 0L, $
                     'poly_area', 0.0)

GOTO, jump
  FOR i=0L, nstripe-1 DO BEGIN 
      stripe = stripes[i]
      sstr = stripearr2string(stripe)

      outfile = maskdir + 'stripe'+sstr+'_'+catname+'_inwindow.fit'
      print,'Outfile = ',outfile

      get_spectra_lcat, stripe, lcat, $
        lss=lss, letter=letter, post=post, sample=sample

      apply_sphpoly_mask, lcat.ra, lcat.dec, masked, unmasked, $
        maskfile=maskfile, $
        completeness=completeness, poly_id=poly_id, poly_area=poly_area

      nlcat = n_elements(lcat)
      outstruct = replicate(ss, nlcat)
      outstruct.completeness = completeness
      outstruct.poly_id = poly_id
      outstruct.poly_area = poly_area
      
      help,lcat,unmasked
      mwrfits, outstruct, outfile, /create

  ENDFOR 
jump:
  get_spectra_lcat, stripes, lcat, $
    lss=lss, letter=letter, post=post, sample=sample
  
  apply_sphpoly_mask, lcat.ra, lcat.dec, masked, unmasked, $
    maskfile=maskfile, $
    completeness=completeness, poly_id=poly_id, poly_area=poly_area
  
  nlcat = n_elements(lcat)
  outstruct = replicate(ss, nlcat)
  outstruct.completeness = completeness
  outstruct.poly_id = poly_id
  outstruct.poly_area = poly_area

  sstr = stripearr2string(stripes)
  
  outfile = maskdir + 'stripe'+sstr+'_'+catname+'_inwindow.fit'
  print,'Outfile = ',outfile

  help,lcat,unmasked
  mwrfits, outstruct, outfile, /create

  ptime,systime(1)-tt
END 
