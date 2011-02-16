PRO powmodel, normrange, powrange, xmin, xmax, modelx, model, norm, power, $
              nnorm=nnorm, npower=npower, nx=nx

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: normrange, powrange, xmin, xmax, modelx, model, norm, power, nnorm=nnorm, npower=npower'
      return
  ENDIF  

;  xmax = 3.500                  ;Mpc
;  xmin = .010
  
  IF n_elements(nnorm) EQ 0 THEN nnorm = 50L
  IF n_elements(npower) EQ 0 THEN npower = 50L
  IF n_elements(nx) EQ 0 THEN nx=1000

  modelx = arrscl( findgen(nx), xmin, xmax )
  norm = arrscl( findgen(nnorm), normrange[0], normrange[1] )
  power = arrscl( findgen(npower), powrange[0], powrange[1] )

  model = fltarr(npower, nnorm, nx )

  FOR inorm=0L, nnorm-1 DO BEGIN 
      FOR ipow=0L, npower-1 DO BEGIN 

          ;; Msolar/pc^2
          model[ipow, inorm, *] = norm[inorm]*(modelx)^power[ipow]

      ENDFOR 
  ENDFOR 

  return
END 
