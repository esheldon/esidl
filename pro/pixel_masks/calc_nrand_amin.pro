PRO calc_nrand_amin, z, radius, AreaMin, Pmax, nrand

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: calc_nrand_amin, z, radius, Amin, Pmax, nrand'
      return
  ENDIF 

  nz = n_elements(z)
  nrand = lonarr(nz)

  DL = angdist_lambda(z)        ; Mpc
  angle = radius/DL             ; Radians
  angle = angle*180./!pi        ; degrees
  Area = !pi*angle^2

  Pmiss = 1. - AreaMin/Area

  wbad=where(Pmiss LT 0, nbad, comp=wgood, ncomp=ngood)

  IF nbad NE 0 THEN nrand[wbad] = 0
  IF ngood NE 0 THEN nrand[wgood] = alog10(Pmax)/alog10(Pmiss[wgood])

END 
