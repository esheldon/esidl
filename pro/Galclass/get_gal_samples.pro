PRO get_gal_samples,l,s,e,b,f

  IF n_params() LT 3 THEN BEGIN
      print,'get_gal_sample,pstruct,s,e,b'
      print,'    returns indices of spirals and ellipticals'
      return
  END

  make_clflag_struct,cl
  cl.failed='Y'
  clflag_select,l,cl,f
  make_clflag_struct,cl
  cl.failed='N'
  cl.ellip='Y'
  clflag_select,l,cl,e
  cl.ellip='D'
  cl.ellip_likely='Y'
  clflag_select,l,cl,e1
  cl.ellip_likely='D'
  cl.spiral_likely='Y'
  clflag_select,l,cl,s1
  cl.spiral_likely='D'
  cl.spiral='Y'
  clflag_select,l,cl,s

  IF n_elements(e1) GT 1 THEN BEGIN 
      e=[e,e1]
      e=e[sort(e)]
  ENDIF
  IF n_elements(s1) GT 1 THEN BEGIN 
      s=[s,s1]
      s=s[sort(s)]
  ENDIF

  make_clflag_struct,cl
  cl.brg='Y'
  clflag_select,l,cl,b

  return
END 
