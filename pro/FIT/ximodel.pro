PRO ximodel, gamrange, ngam, r0range, nr0, radrange, nrad, $
             gam, r0, rad, xi

  IF n_params() LT 6 THEN BEGIN  
      print,'-Syntax: ximodel, gamrange, ngam, r0range, nr0, radrange, nrad, $'
      print,'            r0, gam, rad, xi'
  ENDIF 

  ;; parameters
  ;; r0range in Mpc
  gam = arrscl( findgen(ngam), gamrange[0], gamrange[1] )
  r0 = arrscl( findgen(nr0), r0range[0], r0range[1] )

  ;; radius in Mpc: generate in log space
  logradrange = alog10(radrange)
  lograd = arrscl( findgen(nrad), logradrange[0], logradrange[1] )
  rad = 10.^lograd
  
  ;; the model
  xi = fltarr( ngam, nr0, nrad )

  FOR igam = 0L, ngam-1 DO BEGIN 
      tgam = gam[igam]
      
      FOR ir0=0L, nr0-1 DO BEGIN 

          tr0 = r0[ir0]
          
          xi[igam, ir0, *] = (tr0/rad)^tgam

      ENDFOR
 
  ENDFOR 

END 
