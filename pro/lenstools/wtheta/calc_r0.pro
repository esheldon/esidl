PRO calc_r0, norm, normlow, normhigh, $
             power, powlow, powhigh, $
             rhobar, rhobarlow, rhobarhigh, $
             r0, r0errl, r0errh, r0low, r0high, level=level, doplot=doplot, silent=silent

  ;; norm in something/Mpc^2, rhobar in something/Mpc^3 -> returns r0 in Mpc

  IF n_params() LT 9 THEN BEGIN 
      print,'-Syntax: calc_r0, norm, normlow, normhigh, '
      print,'          power, powlow, powhigh, '
      print,'          rhobar, rhobarlow, rhobarhigh, '
      print,'          r0, r0errl, r0errh, r0low, r0high, level=level, doplot=doplot, silent=silent'
      return
  ENDIF 

  IF n_elements(level) EQ 0 THEN level=0

  ranges = [.683, .955, .997]
  halfranges = ranges/2.
  frac = halfranges[level]

  nrand = 100000

  delvarx, r0, r0errl, r0errh

  IF n_elements(level) EQ 0 THEN level=0

  gam = 1. + power

  f=gamma(0.5)*gamma(0.5*(gam-1))/gamma(0.5*gam)
  r0 = (norm/f/rhobar)^(1./gam)

  normerr = max( [normhigh-norm, $
                  norm-normlow] )
  powerr = max( [powhigh-power, $
                 power-powlow] )
  rhobarerr = max( [rhobarhigh-rhobar, $
                    rhobar-rhobarlow] )

  ;; generate gaussians in each
  rnorm = randomu(seed, nrand, /normal)*normerr + norm
  rpower = randomu(seed, nrand, /normal)*powerr + power
  rrhobar = randomu(seed, nrand, /normal)*rhobarerr + rhobar
  
  tgam = 1.+power
  tf = gamma(0.5)*gamma(0.5*(tgam-1))/gamma(0.5*tgam)

  rr0 = (rnorm/tf/rrhobar)^(1./tgam)
  mr0 = mean(r0)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find range around mean
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wlow  = where(rr0 LT r0, nlow)
  whigh = where(rr0 GE r0, nhigh)

  sl = sort(rr0[wlow])
  sl = wlow[sl]
  sl = reverse(sl)

  sh = sort(rr0[whigh])
  sh = whigh[sh]

  ;; fraction of objects below, as we add more and more
  lfrac = findgen(nlow)/nrand
  wflow = max(where(lfrac LE frac))
  r0low = rr0[sl[wflow]]

  ;; fraction of objects above, as we add more and more
  hfrac = findgen(nhigh)/nrand
  wfhigh = max(where(hfrac LE frac))
  r0high = rr0[sh[wfhigh]]

  r0errl = r0-r0low
  r0errh = r0high-r0

  IF NOT keyword_set(silent) THEN BEGIN 
      print
      print,ntostr(2.*frac*100., 4, /round)+'% region'
      print,'r0 = '+ntostr(r0)+' + '+ntostr(r0errh)+' - '+ntostr(r0errl)
  ENDIF 

  IF keyword_set(doplot) THEN BEGIN 
      maxnum = 1.e6

      plothist,rr0,xhist,yhist,bin=0.05,xtitle='r!D0!N'
      oplot,[r0,r0],[0,maxnum]
      oplot,[mr0,mr0],[0,maxnum],color=!red
      oplot, [r0low, r0low], [0, maxnum], color=!green
      oplot, [r0high, r0high], [0, maxnum], color=!green
  ENDIF 

return

END 
