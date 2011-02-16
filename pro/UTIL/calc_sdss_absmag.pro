FUNCTION calc_sdss_absmag, kmag, z, lumsolar=lumsolar, omegamat=omegamat

  IF n_params() LT 2 THEN BEGIN
      print,'-Syntax: absmag=calc_sdss_absmag(kmag, z, absmag, omegamat=, lumsolar=) '
      print,'kmag should already be k-corrected'
      print,'lumsolar in units of 10^10 solar lum'
      print
      message,'Halting'
  ENDIF 

  nobj = n_elements(z)

  dang = angdist_lambda(z, dlum=dlum, omegamat=omegamat) ;Mpc
  dlum = dlum*1.e6              ;pc

  dlum2 = fltarr(5, nobj)

  dlum2[0,*]=dlum & dlum2[1,*]=dlum & dlum2[2,*]=dlum 
  dlum2[3,*]=dlum & dlum2[4,*]=dlum
  
  absmag = kmag - 5.*alog10(dlum2/10.0)

  wbad = where(absmag LT -30 OR absmag GT -5, nbad)
  IF nbad NE 0 THEN absmag[wbad] = -9999.0

  IF arg_present(lumsolar) THEN BEGIN 
      lumsolar = calc_sdss_lumsolar(absmag)
  ENDIF 

  return,absmag

END 
