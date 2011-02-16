PRO match_spec2tsgal, stripe, notsame=notsame, tsgaldir=tsgaldir, specdir=specdir

  ;; This program is for making the w(theta) catalog. It will
  ;; contain all tsgals, with the redshift info for the galaxies
  ;; that have spectra

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax:  match_spec2tsgal, stripe, notsame=notsame'
      print
      print,'/notsame means use close_match_radec instead of photo_match'
      return
  ENDIF 

  get_tsgals, stripe, tsgals, hdr1, name=tsname, indir=tsgaldir
  FXHCLEAN, hdr1
  get_spectra_lcat, stripe, specgals, indir=specdir

  help,tsgals,specgals

  ep = 1.0                      ;arcsec
  ep = 1.0/3600.                ;degrees

  ;; select "main" galaxies
  main=where( (specgals.primtarget AND 2L^6) NE 0,nmain)

  IF keyword_set(notsame) THEN BEGIN 
      close_match_radec, specgals[main].ra,specgals[main].dec,tsgals.ra,tsgals.dec,$
        spmatch,tsmatch,ep,1
  ENDIF ELSE BEGIN 
      photo_match,specgals[main].run,specgals[main].rerun,specgals[main].camcol,$
        specgals[main].field,specgals[main].id,$
        tsgals.run,tsgals.rerun,tsgals.camcol,$
        tsgals.field,tsgals.id,$
        spmatch, tsmatch
  ENDELSE 

  IF spmatch[0] EQ -1 THEN message,'No matches found!'

  spmatch=main[spmatch]

  nmatch=n_elements(spmatch)
  print
  print,'Found '+ntostr(nmatch)+' matches'
  print

  ;; copy in spectra stuff
  tsgals[tsmatch].plateid = specgals[spmatch].plateid
  tsgals[tsmatch].fiberid = specgals[spmatch].fiberid
  tsgals[tsmatch].z1d = specgals[spmatch].z1d
  tsgals[tsmatch].z1d_error = specgals[spmatch].z1d_error
  tsgals[tsmatch].kcorr = specgals[spmatch].kcorr
  tsgals[tsmatch].absmag = specgals[spmatch].absmag
  tsgals[tsmatch].lum = specgals[spmatch].lum

  outname=repstr(tsname,'tsgal','tsgal_spec')
  print
  print,'Output file: ',outname
  mwrfits, tsgals, outname, hdr1, /create

  setzero,tsgals,specgals

END 
