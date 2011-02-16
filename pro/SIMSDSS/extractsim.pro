
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME  extractsim
;
; 
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO extractsim, nfields, cat, aratio=aratio, posangle=posangle, fwhm=fwhm, $
                tmpdir=tmpdir,dir=dir, catname=catname,fit=fit,psfile=psfile

  IF n_params() EQ 0 THEN BEGIN
      print,'-Syntax: extractsim, nfields, cat, aratio=aratio, posangle=posangle, fwhm=fwhm,tmpdir=tmpdir, dir=dir,fit=fit,catname=catname,psfile=psfile'
      return
  ENDIF
  
  s=systime(1)

  IF NOT keyword_set(tmpdir) THEN tmpdir = '/tmp/'
  IF NOT keyword_set(psfile) THEN psfile=0

  ngal=1000
  nstar=400
  IF NOT keyword_set(fwhm) THEN fwhm = 1.5

  ;; These defaults give e1=.2 e2=0 45 (deg gives e1=0 e2=.2)
  IF NOT keyword_set(aratio) THEN aratio = .816497
  IF NOT keyword_set(posangle) THEN posangle = 0.0 

;;;; run sextractor on some simulated frames

  FOR i=0,nfields-1 DO BEGIN 
      print,'Doing field ',strtrim(string(i+1),2),' Of ',$
            strtrim(string(nfields),2)
      im=0
      sdssframe, ngal, nstar, fwhm, aratio, posangle, im, gals, stars

      sdss_extract, im, tmpcat, fwhm=fwhm, tmpdir=tmpdir

      tmpcat.sfield = i+1
      matchsim, tmpcat, gals, stars, /plt,psfile=psfile
      IF (i EQ 0) THEN BEGIN
          cat = tmpcat 
      ENDIF ELSE BEGIN
          concat_structs,cat,tmpcat,tmp
          cat = tmp
      ENDELSE

  ENDFOR

  IF keyword_set(fit) THEN BEGIN
      IF NOT keyword_set(dir) THEN dir='/sdss4/data1/esheldon/'
      IF NOT keyword_set(catname) THEN BEGIN
          s=systime(1)
          catname = 'testcat_'+strmid(strtrim(string(s),2), 7, 2)+'.fit'
      ENDIF 
      mwrfits, cat, dir+catname,/create
  ENDIF

  print,strtrim(string( (systime(1)-s)/60.0 ), 2),' Total Minutes'

  return
END








