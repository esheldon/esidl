PRO extractfpc, nfields, phdir=phdir

  IF n_params() EQ 0 THEN BEGIN 
    print,'-Syntax: extractfpc, number, dir=dir'
    return
  ENDIF

  ;; phdir is where temporary files are kept
  IF (NOT keyword_set(phdir)) THEN phdir='/sdss4/data1/esheldon/'
  dir = '/sdss4/data1/run259/'
  outdir = '/sdss4/data1/esheldon/'

  ;;; Seeing for run 259 ~ 1.1 or so
  fwhm = 1.1
  color = 'r'
  camcol = '3'
  runstr = run2string(259)
  fieldzero = 210

  FOR i=0,nfields - 1 DO BEGIN
    fieldstr = field2string(fieldzero + i)
    fname = dir+'fpC-'+runstr+'-'+color+camcol+'-'+fieldstr+'.fit'
    print
    print,'Processing File: ',fname
    im=0
    im = mrdfits(fname, 0, hdr,/silent)
    bzero = fxpar(hdr,'BZERO')
    bscale = fxpar(hdr,'BSCALE')
    im = im*bscale + bzero
    sdss_extract, im, tmpcat, fwhm=fwhm, dir=phdir

    IF (i EQ 0) THEN BEGIN
      cat = tmpcat 
    ENDIF ELSE BEGIN
      concat_structs,cat,tmpcat,tmp
      cat = tmp
    ENDELSE
  ENDFOR
;  outname = outdir+'fpC-'+runstr+'-'+color+camcol+'-sexcat.fit'
  outname = outdir + 'tmp-fpc.fit'
  mwrfits, cat, outname,/create


return
END
