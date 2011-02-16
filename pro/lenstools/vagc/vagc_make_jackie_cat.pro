PRO vagc_make_jackie_cat

  sdss_shapecorr_dir = sdssidl_config('shapecorr_dir')
  outfile = sdss_shapecorr_dir + 'combined/jackie/jackiecat_kcorr_rmag.st'

;  GOTO, jump

  stripes = [9,10,11,12,13,14,15,$
             27,28,29,30,31,32,33,34,35,36,37, $
             76,86]

  columns=['ra','dec','absPetroMag','z','completeness','plate','kcorrpetro','petroflux']
  get_spectra_lcat, stripes, lcat, /lss, count=nTot,$
    columns=columns
    
  lcat = lcat[0:99]
  
  eq2csurvey, lcat.ra, lcat.dec, clambda, ceta
  w=where(lcat.completeness GT 0, ngood)
  print
  print,'Threw out '+ntostr((ntot-ngood))+' / '+ntostr(ntot)+' in completeness'
  
  apply_pixel_mask,clambda[w],ceta[w], m, um, /nomissf
  ngood2 = n_elements(um)

  print
  print,'Threw out '+ntostr((nGood-nGood2))+' in pixel mask'

  lcat = lcat[w[um]]

  outst = create_struct('ra', 0d, $
                        'dec', 0d, $
                        'abs_petro_mag', fltarr(5), $
                        'rpetro_mag', 0.0, $
                        'kcorr_petro', fltarr(5), $
                        'z', 0.0, $
                        'completeness', 0.0, $
                        'plate', 0L)

  ;; don't want all 8 mags
  outst = replicate(outst, n_elements(lcat))
  outst.ra = lcat.ra
  outst.dec = lcat.dec

  outst.abs_petro_mag[0] = lcat.absPetroMag[0]
  outst.abs_petro_mag[1] = lcat.absPetroMag[1]
  outst.abs_petro_mag[2] = lcat.absPetroMag[2]
  outst.abs_petro_mag[3] = lcat.absPetroMag[3]
  outst.abs_petro_mag[4] = lcat.absPetroMag[4]

  outst.kcorr_petro[0] = lcat.kcorrpetro[0]
  outst.kcorr_petro[1] = lcat.kcorrpetro[1]
  outst.kcorr_petro[2] = lcat.kcorrpetro[2]
  outst.kcorr_petro[3] = lcat.kcorrpetro[3]
  outst.kcorr_petro[4] = lcat.kcorrpetro[4]

  outst.rpetro_mag = 22.5 - 2.5*alog10(lcat.petroFlux[2])

  outst.z = lcat.z
  outst.completeness = lcat.completeness
  outst.plate = lcat.plate
  print
  print,'Writing file: ',outFile
  write_idlstruct, outst, outFile, /ascii
return
jump:

  lssRandDir = $
    sdssidl_config('lss_dir') + $
    sdssidl_config('lss_vers') + "/safe/random/"

  nrand = 20
  randNums = lindgen(nrand)

  FOR i=0L, nrand-1 DO BEGIN 

      randNum = randNums[i]

      randFile = $
        lssRandDir + "random-"+ntostr(randnum)+".sample14safe.fits"
      outFile = sdss_shapecorr_dir + 'combined/jackie/jackiecat_rand'+strn(randNum,padchar='0',len=2)+'.st'

      columns=['ra','dec','fgot']
      print
      print,'Reading file: ',randFile
      r=mrdfits(randFile,1,columns=columns)

      eq2csurvey, r.ra, r.dec, clambda, ceta
      w=where(r.fgot GT 0)
      apply_pixel_mask,clambda[w],ceta[w], m, um, /nomissf

      r = r[w[um]]

      outst = create_struct('ra', 0d, $
                            'dec', 0d, $
                            'completeness', 0.0)
      outst = replicate(outst, n_elements(r))
      
      outst.ra = r.ra
      outst.dec = r.dec
      outst.completeness = r.fgot



  ENDFOR 


END 
