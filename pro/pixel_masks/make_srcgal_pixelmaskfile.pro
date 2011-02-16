PRO make_srcgal_pixelmaskfile, rlrgMask=rlrgMask

  ;; currently no photoz for 42,43 so no point here
  ;; Also, 82 is masked currently.

  IF keyword_set(rlrgMask) THEN BEGIN 

      maskFile = $
        sdssidl_config('SHAPECORR_DIR')+$
        'masks/pixel_stripe_princeton.mask_simple'
      soFile = '/net/cheops2/home/esheldon/ccode/pixel_masks/apply_maskIDL.so'
      stripes = [9,10,11,12,13,14,15,27,28,29,30,31,32,33,34,35,36,37,76,82,86]

  ENDIF ELSE BEGIN 
      stripes = [9,10,11,12,13,14,15,27,28,29,30,31,32,33,34,35,36,37,76,86]
  ENDELSE 
  nst = n_elements(stripes)

  columns = ['clambda','ceta']
  ss = create_struct('index', 0L, $
                     'maskFlags', 0)

  clr = [1,2,3]
  hirata=1

  FOR i=0L, nst-1 DO BEGIN 

      print,'------------------------------------------------'
      print,'Stripe = '+ntostr(stripes[i])

      srcgal_pixelmaskfile_name, stripes[i], clr, outfile, $
        hirata=hirata, rlrgMask=rlrgMask
      
      print
      print,'Will write to file: ',outfile

      get_scat, stripes[i], clr, scat, $
        columns=columns, count=count, hirata=hirata

      apply_pixel_mask, $
        scat.clambda, scat.ceta, pbad, pgood, maskFlags, /combined,$
        soFile=soFile, maskFile=maskFile

      npass = n_elements(pgood)
      print,$
        ntostr(npass)+'/'+ntostr(count)+' sources unmasked'

      outstruct = replicate(ss, count)
      outstruct.index = lindgen(count)
      outstruct.maskFlags = maskFlags

      print
      print,'Writing to file: ',outfile
      mwrfits, outstruct, outfile, /create

  ENDFOR 

END 
