PRO get_tsgals, stripe, lcat, hdr1, hdr0, indir=indir, count=count, name=name, spec=spec

  np=n_params()
  IF np LT 1 THEN BEGIN 
      print,'-Syntax: get_tsgals, stripe, lcat, hdr1, hdr0, indir=indir, count=count, name=name, spec=spec'
      return
  ENDIF 

  colors = ['u','g','r','i','z']
  stripe_string = ntostr(long(stripe))

  tc = '.fit'

  IF n_elements(indir) EQ 0 THEN BEGIN
      defsysv, '!SDSS_SHAPECORR_DIR', exists=exists
      IF NOT exists THEN sdssidl_setup, /silent
      combdir = !SDSS_SHAPECORR_DIR + 'combined/'
  ENDIF ELSE BEGIN
      combdir = indir
  ENDELSE 
   
  IF keyword_set(spec) THEN BEGIN 
      name = combdir + 'stripe'+stripe_string+'_tsgal_spec'+tc
  ENDIF ELSE BEGIN 
      name = combdir + 'stripe'+stripe_string+'_tsgal'+tc
  ENDELSE 

  IF NOT fexist(name) THEN BEGIN 
      message,'File not found: '+name
  ENDIF 

  print,'Reading File: ',name

  IF np EQ 4 THEN hdr0=headfits(name)
  lcat = mrdfits(name, 1, hdr1)
  count=n_elements(lcat)

  return
END 
