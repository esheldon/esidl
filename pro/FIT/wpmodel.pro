PRO wpmodel, gamrange, ngam, r0range, nr0, rprange, nrp, omegam, $
             gam, r0, rp, wp

  ;; IMPORTANT: this is designed to fit the the density contrast. Thus,
  ;; a factor of 1/[ (2.-alpha)/alpha ]*rhobar is applied to the model to 
  ;; convert the model to a density contrast

  IF n_params() LT 7 THEN BEGIN  
      print,'-Syntax: wpmodel, gamrange, ngam, r0range, nr0, rprange, nrp, omega, $'
      print,'            r0, gam, rp, wp'
  ENDIF 

  ;; parameters
  ;; r0range in Mpc
  gam = arrscl( findgen(ngam), gamrange[0], gamrange[1] )
  r0 = arrscl( findgen(nr0), r0range[0], r0range[1] )

  ;; radius in Mpc: generate in log space
  logrprange = alog10(rprange)
  logrp = arrscl( findgen(nrp), logrprange[0], logrprange[1] )
  rp = 10.^logrp
  wp = fltarr( ngam, nr0, nrp )

  rhobar = omegam*!rhocrit

  FOR igam = 0L, ngam-1 DO BEGIN 
      tgam = gam[igam]
      alpha = tgam - 1.

      fgam = gamma(0.5)*gamma(0.5*(tgam-1.))/gamma(0.5*tgam)
      falpha = (2. - alpha)/alpha

      omg = 1. - tgam
      tmpfunc = fgam/falpha*rp^omg
      FOR ir0=0L, nr0-1 DO BEGIN 

          tr0 = r0[ir0]
          
          wp[igam, ir0, *] = (tr0^tgam)*tmpfunc

      ENDFOR
 
  ENDFOR 

  ;; convert to density Msolar/pc^2
  wp = wp*(rhobar/1.e12)

END 
