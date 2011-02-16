PRO calc_r0_gamma_wgm, combstruct, $
                       gamma, r0, $
                       gammalow, gammahigh, r0low, r0high, chisq_surf, $
                       psfile_front=psfile_front, $
                       dolegend=dolegend, $
                       wuse=wuse, $
                       rebin=rebin, $
                       r0range=r0range, $
                       gamrange=gamrange, $
                       nr0=nr0, ngam=ngam, yfit=yfit, replace=replace, $
                       minchisq=minchisq, nocov=nocov, $
                       input_gamma=input_gamma, input_r0=input_r0, $
                       _extra=_extra

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: calc_r0_gamma_wgm, '
      return
  ENDIF 

  ;; will use r0,err from here as guess
  denscont2wgm, combstruct, wgm, wgmerr, r0, r0errl, r0errh, r0low, r0high,$
                /silent

  gamfac = 4.5
  r0fac = 10.0

;  gamfac = 6.5
;  r0fac = 12.0

  IF keyword_set(rebin) THEN nocov=1

  IF n_elements(r0range) EQ 0 THEN BEGIN 
      maxr0err = max([r0errl, r0errh])
      r0range = r0 + r0fac*[-maxr0err, maxr0err]
  ENDIF 
  IF n_elements(gamrange) EQ 0 THEN BEGIN 
      mgam = 1. + abs(combstruct.power)
      range2error, combstruct.powlow[0], combstruct.power, combstruct.powhigh[0],$
                   perrl, perrh
      gamrange = mgam + gamfac*[-perrl, perrh]
  ENDIF 

  IF n_elements(ngam) EQ 0 THEN ngam = 400L
  IF n_elements(nr0) EQ 0 THEN nr0 = 400L

  IF n_elements(wuse) EQ 0 THEN BEGIN 
      rprange = [0.8*min(combstruct.meanr/1000.), $
                 1.2*max(combstruct.meanr/1000.)]
  ENDIF ELSE BEGIN 
      rprange = [0.8*min(combstruct.meanr[wuse]/1000.), $
                 1.2*max(combstruct.meanr[wuse]/1000.)]
  ENDELSE 

  print,'Creating wp model'
  nrp = 200
  wpmodel, gamrange, ngam, r0range, nr0, rprange, nrp, !omegam, $
           gamvals, r0vals, rp, wp

  IF tag_exist(combstruct, 'covariance') AND NOT keyword_set(nocov) THEN BEGIN 
      print,'Using covariance'
      errsend = combstruct.covariance
  ENDIF ELSE BEGIN 
      print,'Not using covariance'
      IF keyword_set(rebin) THEN errsend = combstruct.sigmaerr_rebin $
      ELSE errsend = combstruct.sigmaerr
  ENDELSE 
  print,'Fitting for r0, gamma'

  IF keyword_set(rebin) THEN BEGIN 
      sigma = combstruct.sigma_rebin
      meanr = combstruct.meanr_rebin/1000.
  ENDIF ELSE BEGIN 
      sigma = combstruct.sigma
      meanr = combstruct.meanr/1000.
  ENDELSE 

  IF n_elements(psfile_front) NE 0 THEN BEGIN 
      xtit=!tsym.gamma
      ytit='r!D0!N [h!U'+!tsym.minus+'1!N Mpc]'
      names = [!tsym.gamma, 'r!D0!N']
  ENDIF ELSE BEGIN 
      xtit=!csym.gamma
      ytit='r!D0!N [h!U'+!csym.minus+'1!N Mpc]'
      names = [!csym.gamma, 'r!D0!N']
  ENDELSE 
  
  ;; use %%BoundingBox: 8 85 450 520
  chisq_conf,meanr,sigma,errsend,rp,wp,$
             gamvals,r0vals,$
             chisq_surf,gamma,r0,$
             gammalow, gammahigh, $
             r0low,r0high, $
             yfit=yfit, $
             errlow1=gerrlow,errhigh1=gerrhigh,$
             errlow2=r0errlow,errhigh2=r0errhigh, $
             xtit=xtit,$
             ytit=ytit,$
             dolegend=dolegend,names=names,nkeep=[4,4], $
             likelihood = tlike, xstyle=1+2, ystyle=1+2, /center, $
             psfile_front=psfile_front, wuse=wuse, $
             minchisq=minchisq, degfree=degfree, $
             input_best1=input_gamma, input_best2=input_r0, $
             _extra=_extra


  IF keyword_set(replace) THEN BEGIN 

      IF tag_exist(combstruct,'wgm_chisq') THEN combstruct.wgm_chisq=minchisq $
      ELSE combstruct = create_struct(combstruct,'wgm_chisq',minchisq)

      IF tag_exist(combstruct,'wgm_degfree') THEN combstruct.wgm_degfree=degfree $
      ELSE combstruct = create_struct(combstruct,'wgm_degfree',degfree)

      IF tag_exist(combstruct,'wgm_gamvals') THEN combstruct.wgm_gamvals=gamvals $
      ELSE combstruct = create_struct(combstruct,'wgm_gamvals',gamvals)

      IF tag_exist(combstruct,'wgm_r0vals') THEN combstruct.wgm_r0vals=r0vals $
      ELSE combstruct = create_struct(combstruct,'wgm_r0vals',r0vals)

      IF tag_exist(combstruct,'wgm_chisq_surf') THEN combstruct.wgm_chisq_surf=chisq_surf $
      ELSE combstruct = create_struct(combstruct,'wgm_chisq_surf',chisq_surf)


;      combstruct = create_struct(combstruct, $
;                                 'wgm_chisq', minchisq, $
;                                 'wgm_degfree', degfree, $
;                                 'wgm_gamvals', gamvals, $
;                                 'wgm_r0vals', r0vals, $
;                                 'wgm_chisq_surf', chisq_surf)

      combstruct.r0 = r0
      combstruct.r0low = r0low
      combstruct.r0high = r0high
      combstruct.gamma = gamma
      combstruct.gammalow = gammalow
      combstruct.gammahigh = gammahigh
      combstruct.yfit[*] = 0
      combstruct.yfit[wuse] = yfit

  ENDIF 

END 
