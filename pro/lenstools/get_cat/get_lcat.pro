PRO get_lcat, stripe, clr, lcat, hdr1, hdr0, indir=indir, typecut=typecut, name=name, count=count, magcut=magcut

  np=n_params()
  IF np LT 2 THEN BEGIN 
      print,'-Syntax: get_lcat, stripe, clr, lcat, hdr1, hdr0, indir=indir, typecut=typecut, name=name, count=count, magcut=magcut'
      return
  ENDIF 

  colors = ['u','g','r','i','z']
  stripe_string = stripe2string(long(stripe))

  IF n_elements(indir) EQ 0 THEN BEGIN
      combdir = sdssidl_config('SHAPECORR_DIR') + 'combined/'
  ENDIF ELSE BEGIN
      combdir = indir
  ENDELSE 
  
  IF n_elements(magcut) EQ 0 THEN BEGIN 
      mcut = 18.
      mcutstr = ntostr( rnd(mcut,1), 4)
      fcutstr = ''
  ENDIF ELSE BEGIN 
      mcutstr = ntostr( rnd(magcut,1), 4)
      fcutstr = '_magcut'+mcutstr
  ENDELSE 

  IF NOT keyword_set(typecut) THEN tc = fcutstr+'.fit' ELSE tc = fcutstr+'_typecut.fit'

  name = combdir + 'stripe'+stripe_string+'_lensgal_'+colors[clr]+tc
    
  IF NOT fexist(name) THEN BEGIN 
      message,'File not found: '+name
  ENDIF 

  print
  print,'Reading File: ',name
  print

  IF np EQ 5 THEN hdr0=headfits(name)
  lcat = mrdfits(name, 1, hdr1)
  count=n_elements(lcat)

  return
END 
