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



PRO test_denscont_likelihood_inputz, zLens, true_zsvals, sigzsvals, $
                                     etanerr, input_evals=input_evals,$
                                     doplot=doplot, $
                                     meansigcritinv=meansigcritinv, $
                                     fracatlens=fracatlens

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: '
      return 
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; lens redshifts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nrand = n_elements(true_zsvals)
  nzlens = n_elements(zlens)

  IF nzlens EQ nrand THEN BEGIN 
      zL = zLens
  ENDIF ELSE IF nzlens EQ 1 THEN BEGIN 
      zL = replicate(zLens, nrand)
  ENDIF ELSE BEGIN 
      message,'Size of zLens is incorrect'
  ENDELSE 

  myusersym, 'fill_circle'

  
  time=systime(1)

  IF keyword_set(truezs) THEN adds = '_true' ELSE adds='_est'

  outdir = '/net/cheops1/data0/esheldon/denscont_likelihood/'

  IF keyword_set(meansigcritinv) THEN meanstr = '_mean' ELSE meanstr=''
  IF n_elements(fracatlens) NE 0 THEN BEGIN
      atlens=fracatlens
      clstr='_clust' 
  ENDIF ELSE BEGIN
      atlens=0.0
      clstr=''
  ENDELSE 
  
  IF n_elements(input_evals) NE 0 THEN estr='_etrue' ELSE estr=''

  IF nzlens EQ 1 THEN BEGIN 
      zlstr = ntostr(zlens,4,/round)
      zlstr = repstr(zlstr, '.', '_')
      fitfile = outdir + $
        'test_denscont_like'+meanstr+clstr+estr+'_zL'+zlstr+'_nrand'+ntostr(nrand)+$
        adds+'_N1.fit'
  ENDIF ELSE BEGIN 
      fitfile = outdir + $
        'test_denscont_like'+meanstr+clstr+estr+'_manyzL_nrand'+ntostr(nrand)+$
        adds+'_N1.fit'
  ENDELSE 
      psfile = repstr(fitfile, '.fit', '.ps')


  WHILE fexist(fitfile) OR fexist(psfile) DO BEGIN
      fitfile = newname(fitfile)
      psfile  = newname(psfile)
  ENDWHILE 

  print,'Fits File: ',fitfile

  IF keyword_set(doplot) THEN BEGIN
      print,'PS file: ',psfile
      begplot,name=psfile,/color
  ENDIF 

  !p.multi=[0,0,2]
  !p.charsize = 1


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The lens model
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ndenscont = 100

  deltasig = 10.

  mindenscont1 = -5000.0
  maxdenscont1 =  5000.0

  mindenscont = 0.0
  maxdenscont = 20.0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; estimated redshifts: add random error
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  zsvals = true_zsvals + randomu(seed, nrand,/normal)*sigzsvals

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;
  ;; The shear of each object, based on true redshift
  ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  siginv = sigmacritinv( (true_zsvals > zL), zL, $
                         /use_lambda, /pcmsolar) >0.0
  shear = deltasig*siginv

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;
  ;; add random errors to ellipticity measurements
  ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  shapenoise = 0.32

  ;; our tangential ellipticities.
  
  IF n_elements(input_evals) NE 0 THEN BEGIN 
      IF n_elements(input_evals) NE nrand THEN message,'input_evals must be size nrand'
      ;; these should already contain shape noise
      evals = shear*2.0 + input_evals; + $
;        randomu(seed, nrand,/normal)*etanerr
      mindenscont = -10.0
      maxdenscont = 30.0
  ENDIF ELSE BEGIN 
      evals = shear*2.0 + $
        randomu(seed, nrand,/normal)*etanerr + $
        randomu(seed, nrand,/normal)*shapenoise 
  ENDELSE 

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

  npts = 100

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plot meanzs plus distribution of errors
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF !d.name EQ 'X' THEN tclr = !green ELSE tclr = !blue

  IF !d.name EQ 'X' THEN window,2

  photoz_dist, zsvals, sigzsvals, 200, prior_zs, prior_pofzs

  xrange = fltarr(2)
  xrange[0] = min(prior_zs, max=ttmm) > (-0.5)
  xrange[1] = ttmm < 1.5
  plothist, true_zsvals, xhist, yhist, bin=0.01,/norm, xrange=xrange

  IF nzlens EQ 1 THEN BEGIN
      oplot, [zLens, zLens], [0.,10000.], color=!red
      legend,'Zl = '+ntostr(zLens,4),/clear
  ENDIF ELSE BEGIN 
      plothist, zL, bin=0.01, color=!red, /overplot, peak=max(yhist)
  ENDELSE 

  oplot, prior_zs, prior_pofzs, color=!magenta

  
  legend,['Input distribution',$
          'Randomly sampled: N = '+ntostr(nrand),$
          'Combined P(z) with noise'], line=[0,0,0],$
         colors = [tclr, !p.color, !magenta],/right,charsize=0.7,/clear


  legmes = 'Est'
  print,'Sending estimated redshift dist.'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; calculate the likelihoods
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
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
  
  meanetan = fltarr(nlens)
  meandeltasig = meanetan
  meanscritinv = meanetan

  mean_denscont = dblarr(nrand)
  err_denscont = dblarr(nrand)

  tmp_mdenscont = dblarr(nperlens)
  tmp_errdenscont = dblarr(nperlens)

  scritinvsum = 0d

  nsig = 4.5

  i = 0
  FOR j=0L, nlens-1 DO BEGIN 

      tmp_loglike[*] = 0d
      FOR k=0L, nperlens-1 DO BEGIN 
          

          minzs = zsvals[i] - nsig*sigzsvals[i] 
          maxzs = zsvals[i] + nsig*sigzsvals[i]


          denscont_likelihood, zL[i], minzs, maxzs, npts, $
                               zsvals[i], sigzsvals[i], $
                               evals[i], etanerr[i], shapenoise, $
                               denscont1, tl, mtl, $
                               tmdenscont, terrdenscont, $
                               prior_zs = sendprior_zs, $
                               prior_pofzs=sendprior_pofzs, $
                               justprior=justprior, $
                               meansigcritinv=meansigcritinv

          tmp_loglike = tmp_loglike + alog(tl)
          
          tmp_mdenscont[k] = tmdenscont
          tmp_errdenscont[k] = terrdenscont
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
      ;; get mean, sdev of true likelihood: these should be sufficient 
      ;; statistics
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;      tdl_getmean_sdev, denscont1, tmp_like, tmean, tsdev
;      meanvals[j] = tmean
;      sigvals[j] = tsdev

      ;; just mean/sdev over individual sources
      w=where(tmp_errdenscont GT 0 AND finite(tmp_errdenscont),nw)
      IF nw EQ 0 THEN BEGIN 
          meanvals[j] = 0.0
          sigvals[j] = 1./0.
      ENDIF ELSE BEGIN 
          weights = 1./tmp_errdenscont[w]^2
          wtot = total(weights)
          meanvals[j] = total(weights*tmp_mdenscont[w])/wtot
          sigvals[j] = sqrt(1./wtot)
      ENDELSE 

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
  print,'wmean_s: '+ntostr(wmean_s)+' '+!plusminus+' '+ntostr(werr_s)
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

  !p.multi=0

  IF nzlens EQ 1 THEN sendzL = zLens ELSE sendzL = -1.
  ss = create_struct('input_denscont', deltasig, $
                     'zlens', sendzL, $
                     'fracatlens', atlens, $
                     'nrand', nrand, $
                     'denscont', denscont, $
                     'denscont_likelihood', denscont_likelihood, $
                     'mean', meand, $
                     'err', sdevd, $
                     'wmean', wmean, $
                     'werr', werr, $
                     'wmean_s', wmean, $
                     'werr_s', werr)

  print,'Output Fits File: ',fitfile
  mwrfits, ss, fitfile, /create

END 
