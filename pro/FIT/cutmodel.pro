PRO cutmodel, sigvrange, cutrange, modelx, model, sigv, cutoff

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: cutmodel_rass, sigvrange, cutrange, modelx, model, sigv, cutoff'
      return
  ENDIF  

  rmax = 3500.
  rmin = 10.

  nx = 1000.
  
  modelx = arrscl( findgen(nx), rmin, rmax )

  nsig = 100L
  ncut = 100L

  sigv = arrscl( findgen(nsig), sigvrange[0], sigvrange[1] )
  cutoff = arrscl( findgen(ncut), cutrange[0], cutrange[1] )

  model = fltarr(ncut, nsig, nx )

  FOR isig=0L, nsig-1 DO BEGIN 
      FOR icut=0L, ncut-1 DO BEGIN 

          ;; Msolar/pc^2
          model[icut, isig, *] = sigdiffsis_trunc( sigv[isig], cutoff[icut], modelx, /core )

      ENDFOR 
  ENDFOR 

  return
END 
