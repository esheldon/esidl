PRO fit_ximodel, struct, $
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

  gamguess = 1.8
  r0guess = 6.0
  Nguess = r0guess^gamguess

  aguess = [Nguess,gamguess]

  meanr = struct.r3
  xi = struct.xi
  xierr = struct.xierr
  IF keyword_set(nocov) THEN cov = xierr ELSE cov = struct.covxi

  IF n_elements(wuse) EQ 0 THEN wuse = lindgen(n_elements(meanr))

  wuse2 = where(struct.xi[wuse] GT 0)
  fitpower, struct.r3[wuse[wuse2]], struct.xi[wuse[wuse2]], struct.xierr[wuse[wuse2]], aguess, $
            yfit, aout, sigma,$
            /silent

  gamfac = 6
  r0fac =  15.0

  fgam = abs(aout[1])
  gamrange = fgam + gamfac*[-sigma[1], sigma[1]]

  Norm = aout[0]
;  Normmax = Norm + r0fac*sigma[0] 
;  Normmin = Norm - r0fac*sigma[1] > 1

  fr0 = Norm^(1./fgam)
;  fr0_min = Normmin
  sigma_r0 = (1/fgam)*sigma[0]*Norm^(1/fgam - 1)

  r0range = fr0 + r0fac*[-sigma_r0, sigma_r0]

  print
  print,'First gam: ',ntostr(fgam)+' '+!plusminus+' '+ntostr(sigma[1])
  print,'First r0:  ',ntostr(fr0)+' '+!plusminus+' '+ntostr(sigma_r0)

  IF n_elements(nr0) EQ 0 THEN nr0 = 400
  IF n_elements(ngam) EQ 0 THEN ngam = 400

  radrange = [0.9*min(struct.r3), 1.1*max(struct.r3)]
  nrad = 200

  ximodel, gamrange, ngam, r0range, nr0, radrange, nrad, $
           gamvals, r0vals, radmodel, ximodel

  
  ;;help,meanr,xi,cov,radmodel,ximodel

  IF n_elements(psfile_front) NE 0 THEN BEGIN 
      xtit=!tsym.gamma
      ytit='r!D0!N [h!U'+!tsym.minus+'1!N Mpc]'
      names = [!tsym.gamma, 'r!D0!N']
  ENDIF ELSE BEGIN 
      xtit=!csym.gamma
      ytit='r!D0!N [h!U'+!csym.minus+'1!N Mpc]'
      names = [!csym.gamma, 'r!D0!N']
  ENDELSE 

  meanr = struct.r3
  chisq_conf,meanr,xi,cov,radmodel,ximodel,$
             gamvals,r0vals,$
             chisq_surf, gamma, r0,$
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

  print
  print,'Final gam: ',ntostr(gamma)+' +'+ntostr(gerrhigh[0])+' - '+ntostr(gerrlow[0])
  print,'Final r0:  ',ntostr(r0)+' +'+ntostr(r0errhigh[0])+' - '+ntostr(r0errlow[0])

  IF keyword_set(replace) THEN BEGIN 
      struct.gamma = gamma
      struct.r0 = r0
      struct.gammalow = gammalow
      struct.gammahigh = gammahigh
      struct.r0low = r0low
      struct.r0high = r0high

      IF tag_exist(struct,'xi_chisq') THEN struct.wgm_chisq=minchisq $
      ELSE struct = create_struct(struct,'xi_chisq',minchisq)

      IF tag_exist(struct,'xi_degfree') THEN struct.xi_degfree=degfree $
      ELSE struct = create_struct(struct,'xi_degfree',degfree)

      IF tag_exist(struct,'xi_gamvals') THEN struct.xi_gamvals=gamvals $
      ELSE struct = create_struct(struct,'xi_gamvals',gamvals)

      IF tag_exist(struct,'xi_r0vals') THEN struct.xi_r0vals=r0vals $
      ELSE struct = create_struct(struct,'xi_r0vals',r0vals)

      IF tag_exist(struct,'xi_chisq_surf') THEN struct.xi_chisq_surf=chisq_surf $
      ELSE struct = create_struct(struct,'xi_chisq_surf',chisq_surf)

  ENDIF 

END 
