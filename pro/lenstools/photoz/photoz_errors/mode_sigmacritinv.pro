PRO mode_sigmacritinv, zlens, meanzs, sigzs, nrand, sigcritinv_mode, $
                       doplot=doplot

  randzs = randomu(seed, nrand, /normal)*sigzs + meanzs
  siginv = sigmacritinv((randzs > zlens),zlens,$
                        /use_lambda,/pcmsolar) > 0.

;  sigma_clip, siginv, meansc, sigsc, $
;              nsig=4, niter=10, /silent

  meansc = mean(siginv)
  sigsc = sdev(siginv)

  ninsig = 20
  sbin = 2.*sigsc/ninsig
  plothist, siginv, sigcxhist, sigcyhist, bin=sbin, xcrange=xcrange,/noplot

  w = where(sigcyhist EQ max(sigcyhist))
  sigcritinv_mode = sigcxhist[w[0]]

  IF keyword_set(doplot) THEN BEGIN 
      plot, sigcxhist, sigcyhist, psym=10
      oplot, [sigcritinv_mode, sigcritinv_mode], [0, max(sigcyhist)], color=!red
  ENDIF 

END 
