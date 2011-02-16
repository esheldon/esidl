
PRO tdl_fitgauss, tmax, tsdev, denscont, like, AA

  nsig = 1.0
  w=where(denscont GE (tmax-nsig*tsdev) AND denscont LE (tmax+nsig*tsdev))
  
  ;;guess = [1.0, tmax, tsdev, 0.01]
  guess = [1.0, tmax, tsdev]
  nterms = n_elements(guess)
  tyfit = gaussfit(denscont[w], like[w], AA,$
                   nterms=nterms, estimates = guess)

END 


PRO tdl_fitquad, denscont, loglike, meanguess, sdevguess, tmean, tsdev

  varguess = sdevguess^2
  Aguess = fltarr(3)
  Aguess[2] = 1./(2.*varguess)
  Aguess[1] = -meanguess/varguess
  ;; this guess is correct if we do loglike = loglike-max(loglike)
  ;; before calling this program
  Aguess[0] = meanguess^2/(2.*varguess)

  ;; svdfit defaults to polynomial

  A = svdfit(denscont, -loglike, 3, A=Aguess,yfit=yfit)

  tvar = 1./(2.*A[2])
  tsdev = sqrt(tvar)
  tmean = -A[1]*tvar

;  print,meanguess,tmean
;  print,sdevguess,tsdev
;  print,Aguess[0], A[0]

;  plot, denscont,-loglike
;  oplot, denscont, yfit, color=!red

;key=get_kbrd(1)

END 

PRO tdl_getmean_sdev, denscont, like, tmean, tsdev, domax=domax

  npts = 100

  tnorm = qgauss(like, denscont, npts)
  IF keyword_set(domax) THEN BEGIN 
      junk = max(like, wmax)
      tmean = denscont[wmax]
  ENDIF ELSE BEGIN 
      tmean = qgauss(like*denscont, denscont, npts)/tnorm
  ENDELSE 
  tvar = qgauss((denscont-tmean)^2*like, denscont,npts)/tnorm
  tsdev = sqrt(tvar)
  
END 

PRO test_denscont_likelihood, zL, Nrand, $
                              prior=prior, trueprior=trueprior,$
                              truezs=truezs, justprior=justprior, $
                              doplot=doplot, zeroprior=zeroprior, $
                              meansigcritinv=meansigcritinv, sdsszs=sdsszs

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: test_denscont_likelihood, zL, Nrand, denscont, denscont_likelihood, prior=prior, trueprior=trueprior, truezs=truezs, justprior=justprior, sig=sig, med=med, mode=mode, doplot=doplot'
      return 
  ENDIF 

  myusersym, 'fill_circle'
  time=systime(1)

  IF keyword_set(truezs) THEN adds = '_true' ELSE adds='_est'

  IF keyword_set(doplot) THEN BEGIN 
      outdir = '/net/cheops1/data0/esheldon/denscont_likelihood/'
      outfile = outdir + $
        'test_denscont_like_zL'+ntostr(zl,4,/round)+'_Nrand'+ntostr(nrand)+$
        adds+'.ps'
      begplot,name=outfile,/color
  ENDIF 

  !p.multi=[0,0,2]
  !p.charsize = 1

  ;;;;;;;;;;;;;;;;;;;
  ;; Units
  ;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; error on ellipticites: this dominates error
  ;; on measurement
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  etanerr = replicate(0.0, nrand)
  shapenoise = 0.32

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The lens model
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  offset = 2.

  ndenscont = 200  
  deltasig = 10. + offset

  ;; range for single measurements
  ;; for 0.32
  mindenscont1 = -5000.0
  maxdenscont1 =  5000.0
      
  mindenscont = 0.0
  maxdenscont = 20.0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; sources drawn from broad gaussian
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  min_redshift = -1000.

  IF keyword_set(sdsszs) THEN BEGIN 
      dir = '/net/cheops1/data0/corrected/sigmacrit/'
      t=mrdfits(dir + 'nzstruct_stripe10.fit',1)

      trueprior_pofzs = t.rnz
      trueprior_zs = t.z

      tnorm = qgauss(trueprior_pofzs, trueprior_zs, 100)
      trueprior_pofzs = trueprior_pofzs/tnorm

  ENDIF ELSE BEGIN 
      
      trueprior_meanzs = 0.5
      trueprior_sigzs = 0.15
      min_pzs = trueprior_meanzs - 3.5*trueprior_sigzs
      max_pzs = trueprior_meanzs + 3.5*trueprior_sigzs
      ntrueprior = 200
      trueprior_zs = arrscl( findgen(ntrueprior), min_pzs, max_pzs )
      
      IF keyword_set(zeroprior) THEN BEGIN 
          min_redshift = 0.0
          trueprior_pofzs = gaussprob(trueprior_zs, $
                                      trueprior_meanzs, trueprior_sigzs,$
                                      minx = min_redshift)
      ENDIF ELSE BEGIN 
          min_redshift = -1000.
          trueprior_pofzs = gaussprob(trueprior_zs, $
                                      trueprior_meanzs, trueprior_sigzs)
      ENDELSE 
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; we have this distribution of errors in each estimated redshift
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sigzs = 0.05

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;
  ;; Now generate true redshifts, errors
  ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  genrand, trueprior_pofzs, trueprior_zs, nrand, true_zsvals

  ;; errors will scale with redshift in some way (this helps
  ;; us avoid crossing redshift zero)
  
  sigzsvals = replicate(sigzs, nrand)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; estimated redshifts: add random error
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  zsvals = true_zsvals + randomu(seed, nrand,/normal)*sigzsvals

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;
  ;; The shear of each object, based on true redshift
  ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  siginv = sigmacritinv(zl, true_zsvals)
  shear = deltasig*siginv
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;
  ;; add random errors to ellipticity measurements
  ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; our tangential ellipticities.

  evals = shear*2.0 + $
    randomu(seed, nrand,/normal)*etanerr + $
    randomu(seed, nrand,/normal)*shapenoise 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; values of the density contrast to evaluate at
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  denscont1 = arrscl(dindgen(ndenscont), mindenscont1, maxdenscont1)
  denscont = arrscl(dindgen(ndenscont), mindenscont, maxdenscont)
  loglike = dblarr(ndenscont)
  tmp_loglike = dblarr(ndenscont)
  tmp_like = dblarr(ndenscont)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; For zs integral
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  npts = 200

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plot meanzs plus distribution of errors
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF !d.name EQ 'X' THEN tclr = !green ELSE tclr = !blue

  IF !d.name EQ 'X' THEN window,2
  plothist, true_zsvals, bin=0.01,/norm
  ;plothist, zsvals, bin=0.01, /norm, /overplot, color=!red
  IF n_elements(input_true_zsvals) EQ 0 THEN BEGIN 
      oplot,trueprior_zs,trueprior_pofzs,color=tclr
  ENDIF 

  oplot, [zL, zL], [0.,10000.], color=!red

  IF keyword_set(zeroprior) THEN BEGIN 
      photoz_dist, zsvals, sigzsvals, 200, prior_zs, prior_pofzs, $
                   minz=min_redshift
  ENDIF ELSE BEGIN  
      photoz_dist, zsvals, sigzsvals, 200, prior_zs, prior_pofzs
  ENDELSE 

  ;;photoz_dist, true_zsvals, sigzsvals, 200, prior_zs, prior_pofzs

  oplot, prior_zs, prior_pofzs, color=!magenta

  legend,'Zl = '+ntostr(zL,4),/clear
  legend,['Input distribution',$
          'Randomly sampled: N = '+ntostr(nrand),$
          'Combined P(z) with noise'], line=[0,0,0],$
         colors = [tclr, !p.color, !magenta],/right,charsize=0.7,/clear

  IF keyword_set(justprior) THEN print,'Just sending prior'

  IF keyword_set(trueprior) THEN BEGIN
      print,'Sending true prior'
      sendprior_pofzs = trueprior_pofzs
      sendprior_zs = trueprior_pofzs
      sendprior_pofzs = sendprior_pofzs/qgauss(sendprior_pofzs, sendprior_zs, 1000)
  ENDIF ELSE IF keyword_set(prior) THEN BEGIN 
      print,'Sending estimated prior'
      sendprior_pofzs = prior_pofzs
      sendprior_zs = prior_zs
      sendprior_pofzs = sendprior_pofzs/qgauss(sendprior_pofzs, sendprior_zs, 1000)
  ENDIF 

  IF keyword_set(truezs) THEN BEGIN 
      legmes = 'True'
      print,'Sending true redshift dist.'
      send_zsvals = true_zsvals
  ENDIF ELSE BEGIN 
      legmes = 'Est'
      print,'Sending estimated redshift dist.'
      send_zsvals = zsvals
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; calculate the likelihoods
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  modnum = 100
  
  colors = arrscl( findgen(nrand), 255, 65280 )

  IF !d.name EQ 'X' THEN BEGIN
      window,0
      pos1 = [0.1,0.6,0.95,0.95]
      pos2 = [0.1,0.1,0.95,0.5]
  ENDIF 

  nperlens = 10
  nlens = nrand/nperlens
  meanvals = fltarr(nlens)
  sigvals  = meanvals
  
  mean_denscont = dblarr(nrand)
  err_denscont = dblarr(nrand)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; number of sigma to integrate over in zs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nsig = 4.5

  i = 0
  FOR j=0L, nlens-1 DO BEGIN 

      tmp_loglike[*] = 0d
      FOR k=0L, nperlens-1 DO BEGIN 
          
 

          IF keyword_set(justprior) THEN BEGIN 
              minzs = min(sendprior_zs, max=maxzs) > min_redshift
          ENDIF ELSE BEGIN 
              minzs = zsvals[i] - nsig*sigzsvals[i] > min_redshift; > (zL+0.001)
              maxzs = zsvals[i] + nsig*sigzsvals[i]
          ENDELSE 
          ;;fac = sqrt(2.)
          fac = 1.0
          IF keyword_set(trueprior) OR keyword_set(prior) THEN BEGIN 
              denscont_likelihood, zL, minzs, maxzs, npts, $
                                   send_zsvals[i], sigzsvals[i]*fac, $
                                   evals[i], etanerr[i], shapenoise, $
                                   denscont1, tl, mtl, $
                                   prior_zs = sendprior_zs, $
                                   prior_pofzs=sendprior_pofzs, $
                                   justprior=justprior, $
                                   meansigcritinv=meansigcritinv
          ENDIF ELSE BEGIN 
              
              denscont_likelihood, zL, minzs, maxzs, npts, $
                                   send_zsvals[i], sigzsvals[i]*fac, $
                                   evals[i], etanerr[i], shapenoise, $
                                   denscont1, tl, mtl, $
                                   meansigcritinv=meansigcritinv
              
          ENDELSE 

          tmp_loglike = tmp_loglike + alog(tl)

          mean_denscont[i] = tmdenscont
          err_denscont[i] = terrdenscont

         i = i+1L
      ENDFOR 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; renormalize
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      tmp_loglike = tmp_loglike - max(tmp_loglike)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; interplate to finer grid for overall likelihood
      ;; function
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      interp_tmp_loglike = interpol(tmp_loglike, denscont1, $
                                    denscont, /quadratic)
      loglike = loglike + interp_tmp_loglike

      ;; linear for plotting
      tmp_like = (exp(1d))^(tmp_loglike)
      tlike = (exp(1d))^(interp_tmp_loglike)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; get mean, sdev: these should be sufficient statistics
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      tdl_getmean_sdev, denscont1, tmp_like, tmean, tsdev
      meanvals[j] = tmean
      sigvals[j] = tsdev

      ;;tdl_fitquad, denscont1, tmp_loglike, meanvals[j], sigvals[j], $
      ;;             tmean2, tsdev2

      

      IF ((j MOD 10) EQ 0) THEN BEGIN
          print,'Nlens = ',ntostr(j)+'/'+ntostr(nlens-1)

          IF !d.name EQ 'X' THEN BEGIN 
              ;;like = (exp(1d))^(loglike-max(loglike))
              like = tmp_like
              ;;ml = max(like)
              IF j EQ 0 THEN BEGIN
                  plot, denscont1, tmp_like, yrange=[0,1.2], pos=pos1
                  oplot, denscont1, tmp_like, color=colors[i]

                  oplot, denscont, tlike, psym=8, symsize=0.5,color=colors[i]
;                  tmeanval = interpol(tmp_like, denscont1, tmean)
;                  oplot,[tmean],[tmeanval],color=!magenta, psym=8

                  plot, denscont, (exp(1d))^(loglike-max(loglike)),$
                        pos=pos2, /noerase,yrange=[0,1.2]

                  oplot, denscont, (exp(1d))^(loglike-max(loglike)),$
                         color=colors[i]
;stop
              ENDIF ELSE BEGIN 
                  plot, denscont1, tmp_like, color=colors[i], pos=pos1, $
                        /noerase, xstyle=4,ystyle=4, yrange=[0,1.2]

                  oplot, denscont, tlike, psym=8, symsize=0.5,color=colors[i]
;                  tmeanval = interpol(like, denscont1, tmean)
;                  oplot,[tmean],[tmeanval],color=!magenta, psym=8

                  plot, denscont, (exp(1d))^(loglike-max(loglike)),$
                        color=colors[i], pos=pos2, /noerase, xstyle=4,ystyle=4
;stop
              ENDELSE 
          ENDIF 
          
      ENDIF 

  ENDFOR 

  IF !d.name EQ 'X' THEN BEGIN
      wset,2
      !p.multi=[1,0,2]
  ENDIF 
  ;;key=get_kbrd(1)

  loglike = loglike - max(loglike)
  denscont_likelihood = (exp(1d))^( loglike )

  myusersym, 'fill_circle'
  plot, denscont, denscont_likelihood,/ynozero, $
        yrange = [0,1.2], psym=8;, xrange=[deltasig-5., deltasig+5.]
  oplot, denscont, denscont_likelihood
  oplot, [deltasig, deltasig], [0, 10000], color=!red

  norm = qgauss(denscont_likelihood, denscont, 100)
  denscont_likelihood = denscont_likelihood/norm

  meand = qgauss(denscont*denscont_likelihood, denscont, 100)
  vard  = qgauss((denscont-meand)^2*denscont_likelihood, denscont, 100)
  sdevd = sqrt(vard)

  meand = total(denscont*denscont_likelihood)/total(denscont_likelihood)
  oplot,[meand,meand],[0,10],color=tclr

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; using mean values of each lens likelihood function
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  weights = 1./sigvals^2
  wtot = total(weights)
  wmean = total(weights*meanvals)/wtot
  werr = sqrt(1./wtot)

  wyfit = gaussprob(denscont, wmean, werr)
  wyfit = wyfit/max(wyfit)
  oplot, denscont, wyfit, color=!cyan

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Using mean denscont from each source
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  w=where(err_denscont GT 0.0 AND finite(err_denscont))
  weights = 1./err_denscont[w]^2
  wtot = total(weights)
  wmean_s = total(weights*mean_denscont[w])/wtot
  werr_s = sqrt(1./wtot)

  wyfit = gaussprob(denscont, wmean_s, werr_s)
  wyfit = wyfit/max(wyfit)
  oplot, denscont, wyfit, color=!green

  print,'true: '+ntostr(deltasig)
  print,'mean: '+ntostr(meand)+' '+!plusminus+' '+ntostr(sdevd)
  print,'wmean: '+ntostr(wmean)+' '+!plusminus+' '+ntostr(werr)
  diff = deltasig-meand
  diffnsig = abs(diff/sdevd)
  print,'diff: '+ntostr(diff)+' (nsig='+ntostr(diffnsig)+')'


  legend,['Input Value: '+ntostr(deltasig,5,/round),$
          'Mean: '+ntostr(meand,5,/round)+!csym.plusminus+$
          ntostr(sdevd,5,/round)],$
         line=[0,0],$
         color=[!red, tclr],/clear

  legend,legmes,/right,/clear

  ptime,systime(1)-time

  !p.multi=0

  IF keyword_set(doplot) THEN endplot


END 
