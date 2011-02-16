PRO sismodel, sigvrange, modelx, model, sigv

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: cutmodel_rass, sigvrange, modelx, model, sigv'
      return
  ENDIF  

  rmax = 3500.                  ;kpc
  rmin = 10.

  nx = 1000.
  modelx = arrscl( findgen(nx), rmin, rmax )

  nsig = 500L
  sigv = arrscl( findgen(nsig), sigvrange[0], sigvrange[1] )

  model = fltarr( nsig, nx )

  FOR isig=0L, nsig-1 DO BEGIN 

      ;; Msolar/pc^2
      model[isig, *] = sigmasis( sigv[isig], modelx, /core )

  ENDFOR 

  return
END 
