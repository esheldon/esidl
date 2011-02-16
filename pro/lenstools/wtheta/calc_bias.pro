PRO calc_bias, gamma, $
               rgm, rgmerr, $
               rgg, rggerr, $
               bias, biaserr, $
               level=level

  IF n_params() LT 5 THEN BEGIN 
      print,'-Syntax: calc_bias, gamma, '
      print,'           rgm, rgmerr, '
      print,'           rgg, rggerr, '
      print,'           bias, biaserrlow, biaserrhigh, '
      print,'           level=level'
      print
      print,'gamma should be 3d power index'
      return
  ENDIF 

  ;; assume they are the same shape: scale independent bias
  ;; xi_{gm} = (r/r_{gm})^{-gamma}
  ;;
  ;; The result is b/r = (rgg/rgm)^gamma

  IF n_elements(level) EQ 0 THEN level=0

  ranges = [.683, .955, .997]
  halfranges = ranges/2.
  frac = halfranges[level]

  nrand = 1000000

  delvarx, bias, biaserrlow, biaserrhigh

  bias = (rgg/rgm)^gamma

  biaserr = bias*gamma*sqrt( (rgmerr/rgm)^2 + (rggerr/rgg)^2 )

  print,'b/r = '+ntostr(bias)+' '+!plusminus+' '+ntostr(biaserr)

  ;; now try a monte-carlo way

  nrand = 100000
  rrgm = randomu(seed, nrand, /normal)*rgmerr + rgm
  rrgg = randomu(seed, nrand, /normal)*rggerr + rgg

  rbias = (rrgg/rrgm)^gamma

  ;; recalculate the mean
  ;;mbias = mean(rbias)

  ;; use the mean
  mbias = bias

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find range around mean
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wlow=where(rbias LT mbias, nlow)
  whigh=where(rbias GE mbias, nhigh)

  sl = sort(rbias[wlow])
  sl = wlow[sl]
  sl = reverse(sl)

  sh = sort(rbias[whigh])
  sh = whigh[sh]

  ;; fraction of objects below, as we add more and more
  lfrac = findgen(nlow)/nrand
  wflow = max(where(lfrac LE frac))
  biaslow = rbias[sl[wflow]]

  ;; fraction of objects above, as we add more and more
  hfrac = findgen(nhigh)/nrand
  wfhigh = max(where(hfrac LE frac))
  biashigh = rbias[sh[wfhigh]]

  biaserrl = bias-biaslow
  biaserrh = biashigh-bias

  print
  print,ntostr(2.*frac*100., 4, /round)+'% region'
  print,'b/r = '+ntostr(mbias)+' + '+ntostr(biaserrh)+' - '+ntostr(biaserrl)

  maxnum = 1.e6
  plothist,rbias,xhist,yhist,bin=0.03,xtitle='b/r'
  oplot,[bias,bias],[0,maxnum]
  oplot,[mbias,mbias],[0,maxnum],color=!red
  oplot,[biaslow,biaslow],[0,maxnum],color=!green
  oplot,[biashigh,biashigh],[0,maxnum],color=!green


END 
