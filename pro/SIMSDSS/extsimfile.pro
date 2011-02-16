
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME  extract_sim
;
; KEYWORDS:
;     number:  testimage_number.fit
;     fwhm:    fwhm to send to sextractor for star-gal separation
;     phdir:   directory to store temporary files for sextractor
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO extsimfile, number=number, fwhm=fwhm, phdir=phdir


;;;; run sextractor on some simulated frames

  IF (NOT keyword_set(fwhm)) THEN fwhm = 1.5
  dir = '/sdss4/data1/esheldon/'
  IF n_elements(number) EQ 0 THEN BEGIN 
    e=exist(dir+'testimage_*',count)
    begjj = 0
    endjj = count-1
  ENDIF ELSE BEGIN
    count = n_elements(number)
    begjj = number[0]-1
    endjj = number[count-1]-1
  ENDELSE

  nbad = 0
  FOR jj = begjj, endjj DO BEGIN
    num = strtrim(string(jj+1),2)
    imname = dir + 'testimage_'+num+'.fit'
    print
    print,'Reading From ',imname
    fits_info, imname, /silent, n_ext=nfields
    ;;;  fits_info doesn't count first one
    nfields=nfields+1
    print,'Number of Fields in File: ',nfields
    FOR i=0,nfields-1 DO BEGIN 
      im=0
      im = mrdfits(imname,i,hdr)
      s=size(im)
      sx=s[1]
      sy=s[2]
      IF (sx EQ 2048 AND sy EQ 1489) THEN BEGIN
        sdss_extract, im, tmpcat, fwhm=fwhm
        IF (i EQ 0) THEN BEGIN
          cat = tmpcat 
        ENDIF ELSE BEGIN
          concat_structs,cat,tmpcat,tmp
          cat = tmp
        ENDELSE
      ENDIF ELSE BEGIN
        print,'ERROR:  BAD IMAGE'
        nbad = nbad+1
      ENDELSE
    ENDFOR
    mwrfits, cat, dir+'testcat_'+num+'.fit',/create
  ENDFOR

  print, 'Number of bad images:  ',nbad



return
END







