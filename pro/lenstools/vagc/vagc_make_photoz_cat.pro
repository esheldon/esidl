PRO vagc_make_photoz_cat, lrg=lrg

  stripes = [9,10,11,12,13,14,15,16,$
             26,27,28,29,30,31,32,33,34,35,36,37, $
             42,43,44,$
             76,82,86]

  columns=['ra', 'dec', $
           'z','z_err', 'modelFlux', 'modelFlux_Ivar','extinction']

  IF keyword_set(lrg) THEN BEGIN 
      outfile = $
        sdssidl_config('shapecorr_dir') + 'combined/photozcat_lrg_input.st'
      get_spectra_lcat, stripes, lcat, /lrg, count=nTot,$
        columns=columns
  ENDIF ELSE BEGIN 
      outfile = $
        sdssidl_config('shapecorr_dir') + 'combined/photozcat_input2.st'
      get_spectra_lcat, stripes, lcat, /lss, count=nTot,$
        columns=columns
  ENDELSE 

  wgood = where(lcat.modelFlux_Ivar[0] GT 0.0 AND $
                lcat.modelFlux_Ivar[1] GT 0.0 AND $
                lcat.modelFlux_Ivar[2] GT 0.0 AND $
                lcat.modelFlux_Ivar[3] GT 0.0 AND $
                lcat.modelFlux_Ivar[4] GT 0.0 AND $
                $
                lcat.modelFlux[0] GT 0.0 AND $
                lcat.modelFlux[1] GT 0.0 AND $
                lcat.modelFlux[2] GT 0.0 AND $
                lcat.modelFlux[3] GT 0.0 AND $
                lcat.modelFlux[4] GT 0.0, nlcat)

  lcat = lcat[wgood]

  print
  print,ntostr(nlcat)+'/'+ntostr(nTot)+'  passed cuts'

  outst = create_struct('ra',0d, $
                        'dec',0d, $
                        'z', 0.0, $
                        'z_err', 0.0, $
                        'modelmag', fltarr(5),$
                        'modelmag_err', fltarr(5) )

  outst = replicate(outst, nlcat)

  outst.ra = lcat.ra
  outst.dec = lcat.dec

  outst.z = lcat.z
  outst.z_err = lcat.z_err
 
  outst.modelmag[0] = $
    22.5 - 2.5*alog10(lcat.modelFlux[0]) - lcat.extinction[0]
  outst.modelmag[1] = $
    22.5 - 2.5*alog10(lcat.modelFlux[1]) - lcat.extinction[1]
  outst.modelmag[2] = $
    22.5 - 2.5*alog10(lcat.modelFlux[2]) - lcat.extinction[2]
  outst.modelmag[3] = $
    22.5 - 2.5*alog10(lcat.modelFlux[3]) - lcat.extinction[3]
  outst.modelmag[4] = $
    22.5 - 2.5*alog10(lcat.modelFlux[4]) - lcat.extinction[4]

  sigma = 1./sqrt(lcat.modelFlux_Ivar[0])
  outst.modelmag_err[0] = 2.5/alog(10.0)*sigma/lcat.modelFlux[0]

  sigma = 1./sqrt(lcat.modelFlux_Ivar[1])
  outst.modelmag_err[1] = 2.5/alog(10.0)*sigma/lcat.modelFlux[1]

  sigma = 1./sqrt(lcat.modelFlux_Ivar[2])
  outst.modelmag_err[2] = 2.5/alog(10.0)*sigma/lcat.modelFlux[2]

  sigma = 1./sqrt(lcat.modelFlux_Ivar[3])
  outst.modelmag_err[3] = 2.5/alog(10.0)*sigma/lcat.modelFlux[3]

  sigma = 1./sqrt(lcat.modelFlux_Ivar[4])
  outst.modelmag_err[4] = 2.5/alog(10.0)*sigma/lcat.modelFlux[4]

  print
  print,'Writing file: ',outFile
  write_idlstruct, outst, outFile, /ascii

END 
