PRO fit_dn_sigma, lensum, dn, dnarc, dnpix, cut=cut

;  !p.multi=[0,0,2]

  IF n_elements(dnpix) EQ 0 THEN BEGIN 
      thresh=20.75
      dnclr=1
      fit_dn, lensum, dn, dnarc, dnpix, thresh=thresh, clr=dnclr
      lensum.dn = dn
  ENDIF 

  IF keyword_set(cut) THEN BEGIN 
      mindn = 1.0
      maxdn = 3.5
      minpix = 4.0
  ENDIF ELSE BEGIN 
      mindn = 0.0
      maxdn = 100.0
      minpix = 4.
  ENDELSE 
  w=where(dn GT mindn AND dn LT maxdn AND dnpix GT minpix AND lensum.sigma_cc GT 0.,nw)
;  w=where(dn GT 0. AND lensum.sigma_cc GT 0.)

  ytag = 'sigma_cc'
  yerrtag = 'sigma_cc_err_p'
  xtag = 'dn'

  normrange=[100.,140.]
  powrange=[0.5,0.8]

;  fit_tag_vs_tag_nobin, lensum[w], xtag, ytag, yerrtag, normrange,powrange, $
;    fitstruct

  xlog=1
  ylog=1

  !p.charsize = 1.0
  charsize = 1.2
  !p.thick = 1
  !x.thick = 1
  !y.thick = 1
  !p.charthick=1

;  ymax = 1.2*max(lensum[w].sigma_cc)
;  ymin = 0.8*min(lensum[w].sigma_cc)
;  yrange=[ymin,ymax]
;  xmax = 1.2*max(dn[w])
;  xmin = 0.8*min(dn[w])
;  xrange=[xmin,xmax]
  xrange=[0.129502, 11.0838]
  yrange=[43.2000, 392.400]
  print,xrange
  print,yrange
  ploterror,dn[w],lensum[w].sigma_cc,lensum[w].sigma_cc_err_p,psym=4, $
    xlog=xlog,ylog=ylog,ytitle='sigma',xtitle='Dn',$
    xrange=xrange,yrange=yrange,xstyle=1,ystyle=1
;  plot,dn[w],lensum[w].sigma_cc,psym=4,xlog=xlog,ylog=ylog,$
;    ytitle='sigma',xtitle='Dn',xrange=xrange,yrange=yrange,xstyle=1,ystyle=1

  nx=1000
  xx = arrscl(findgen(nx), 0, 20.)
;  yy = fitstruct.bestnorm*xx^(fitstruct.bestpow)
;  oplot,xx,yy

  logx = alog10(lensum[w].dn)
  logy = alog10(lensum[w].sigma_cc)
  logyerr = alog10(lensum[w].sigma_cc_err_p)

  print
  print,'Using FITLIN'
  fitlin, logx, logy, logyerr, a, siga, b, sigb,chisq,/silent
  print,10.^(a),10.^(siga),b,sigb
  print,'chisq = ',ntostr(chisq)

  yy = (10.^(a))*(xx)^(b)
  oplot,xx,yy

  yy = (10.^(a))*(xx)^(0.75)
  oplot,xx,yy,line=4

;  key=get_kbrd(1)


  !p.multi=0

  return
END 
