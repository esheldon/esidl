PRO test_deltasig_invert, dops=dops, adderror=adderror, nrealize=nrealize
  
  IF n_elements(nrealize) EQ 0 THEN nrealize = 1

  IF keyword_set(dops) THEN BEGIN
      IF keyword_set(adderror) THEN BEGIN 
          name='~/plots/test_invert_deltasig_adderror_N1.ps'
          WHILE fexist(name) DO name = newname(name)
      ENDIF ELSE BEGIN 
          name='~/plots/test_invert_deltasig.ps'
      ENDELSE 
      print,name
  ENDIF 
  
 

  ;; get radii from iro's points
  dir = '~/lensout/iro/'
  xifile = dir + 'gm3d_18.0_to_21.6_sdssw_z0.1.dat'
  dfile = dir +  'gm2d_18.0_to_21.6_sdssw_z0.1.dat'
  readcol, xifile, rr, xi, xierrp, junk, xierrjack, xierrmax 
  readcol, dfile, meanr, dsig, sigm, sig, $
           dsigerrp, junk, dsigjack, dsigerrmax

; $\alpha = 0.76 \pm 0.05$ and $A = (4.3 \pm 0.4) h $M$_{\sun}$pc$^{-2}$

  r0 = 6.0                      ; Mpc
  gam = 1.8
  
  omegam = 0.27
  rhobar = omegam*!rhocrit ; in msun/Mpc^3

  G = gamma(0.5)*gamma( 0.5*(gam-1.0) )/gamma(0.5*gam)
  A = (gam-1.0)/(3.0-gam)*rhobar*r0^gam*G ; in Msun/Mpc^2
  A = A/1.e12                   ; In Msun/pc^2
  alpha = gam-1.0

  print
  print,'A = '+ntostr(A)+' Mpc'
  
  deltasig = A*meanr^(-alpha)

  forprint,dsigerrmax/deltasig
  print
  perc_err = mean(dsigerrmax/deltasig)


  deltasig_err = perc_err*deltasig

  ;; Now add error to deltasig

  nrad = n_elements(meanr)

  seed = long(systime(1))  
  IF keyword_set(adderror) THEN BEGIN 

      xi = fltarr(nrad-1)
      xierr = xi
      xitmp = fltarr(nrealize, nrad-1)
      FOR i=0L, nrealize-1 DO BEGIN 
 
          add_errors = randomu(seed, nrad, /normal)*deltasig_err 

          deltasig_use = deltasig + add_errors

          ;; create outliers
          ;;  deltasig[3] = deltasig[3] + add_errors[3]
          ;;  deltasig[9] = deltasig[9] + add_errors[9]
          covariance = diagonal_array(add_errors^2)

          ;; dave's code
          print,'Calling delta_sigma2der_sigma'
          
          delta_sigma2der_sigma, meanr, deltasig_use, dersig, covariance, $
                                 covder, derr
          print,'Calling der_sigma2xi'
          rrange = [min(meanr), max(meanr)]
          der_sigma2xi, meanr, dersig, r3, txi, corrfac, covder, covxi, $
                        xi_interior=xi_interior, rrange=rrange


          xitmp[i, *] = txi

      ENDFOR 

      FOR i=0L, nrad-2 DO BEGIN 
          mom = moment(xitmp[*,i])
          xi[i] = mom[0]
          xierr[i] = sqrt(mom[1]/nrealize)
      ENDFOR 

  ENDIF ELSE BEGIN 
      covariance = diagonal_array(replicate(1.0, nrad))
      deltasig_use = deltasig
      ;; dave's code
      print,'Calling delta_sigma2der_sigma'
      
      delta_sigma2der_sigma, meanr, deltasig, dersig, covariance, $
                             covder, derr
      print,'Calling der_sigma2xi'
      rrange = [min(meanr), max(meanr)]
      der_sigma2xi, meanr, dersig, r3, xi, corrfac, covder, covxi, $
                    xi_interior=xi_interior, rrange=rrange

      nr3 = n_elements(r3)
      xierr = fltarr(nr3)
      FOR i=0L, nr3-1 DO xierr[i] = sqrt(covxi[i,i])

  ENDELSE 


  !p.multi=0



;  aploterror, 1.0, r3, xi, xierr, psym=8, /xlog, /ylog, $
;              yrange = [0.1, 1.e5], ystyle=1+2


  IF keyword_set(dops) THEN begplot,name=name,/color
 
  setup_mystuff

  ;; plots
  IF !d.name EQ 'PS' THEN oclr = !blue ELSE oclr = !green

  
  model = (r0/r3)^gam
  ratio = xi/model
  ratio = ratio[1:nrad-2]
  ratio_r = r3[1:nrad-2]
  IF NOT keyword_set(adderror) THEN BEGIN 
      plottwo, r3, xi, r3, ratio, $
               frac1=0.75, /xlog, /topylog, toppsym=8, botpsym=8, $
               xoplot = r3, yoplot = model, oplotcolor=oclr, $
               botyrange = [0.97, 1.01], botystyle=1, topaspect=1, $
               topytitle = !csym.xi+'!Dgm!N(r)', xtitle=!xigmxtitle
  ENDIF ELSE BEGIN  
      ratioerr = ratio*xierr/xi
      plottwoerr, r3, xi, xierr, r3, ratio, ratioerr, $
                  frac1=0.75, /xlog, /topylog, toppsym=8, botpsym=8, $
                  xoplot = r3, yoplot = model, oplotcolor=oclr, $
                  botystyle=2, topaspect=1, $
                  topytitle = !csym.xi+'!Dgm!N(r)', xtitle=!xigmxtitle
      legend,'Nrealize = '+ntostr(nrealize),/bottom,/left,box=0
  ENDELSE 
  oplot, [0.01, 1000], [1,1], color=oclr

;  aplot, 1.0, r3, xi, psym=8, /xlog, /ylog, $
;              yrange = [0.1, 1.e5], ystyle=1+2

;  oplot, r3, (r0/r3)^gam, color=oclr

  IF keyword_set(dops) THEN endplot

return
END 
