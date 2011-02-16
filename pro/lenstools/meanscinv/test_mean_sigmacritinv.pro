
PRO test_mean_sigmacritinv, zlens, meanzs, sigzs, npts, sigcritinv_mean, $
                            silent=silent, double=double

  ;; This demonstrates that integrating P(z)*sigmacrit_inv gives same
  ;; answer as the mean over monte carlo trials, but other methods such
  ;; as using the value at the mean zs or trying to integrate the 
  ;; histogram in sigmacrit_inv (which has a delta function) do not
  ;; work as well

  zL = zlens
  mzs = meanzs
  szs = sigzs

  minzs = meanzs - 4.0*sigzs
  maxzs = meanzs + 4.0*sigzs

  nz = 2000
  zvals = arrscl(findgen(nz), minzs, maxzs)
  pofzs = gaussprob(zvals, meanzs, sigzs, double=double)

  siginv = sigmacritinv(zlens, zvals)*1.e4

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; integrate over the redshift probability distribution
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  func = pofzs*siginv
  sigcritinv_mean = qgauss(func, zvals, npts)

  sigcritinv_meanzs = sigmacritinv(zlens, meanzs)*1.e4 

  !p.multi=[0,0,4]
  pold = !p.charsize
  !p.charsize=2
  maxplot = 1.e6
  plot, zvals, pofzs, xtitle='zs', ytitle='P(z)'
  oplot,[zlens,zlens],[0,1.e4]
  plot, zvals, siginv,xtitle='zs',ytitle='siginv'
  oplot,zvals,siginv,color=!green
  plot, zvals, func,xtitle='zs',ytitle='P(z)*siginv'
  
;  plothist, siginv, bin=0.05
;  oplot, [sigcritinv_meanzs, sigcritinv_meanzs], [0, maxplot], color=!red
;  oplot, [sigcritinv_mean, sigcritinv_mean], [0, maxplot], color=!blue
 
  ;;;;;;;;;;;;;;;;;;;;;
  ;; Monte carlo it
  ;;;;;;;;;;;;;;;;;;;;;

  nrand = 500000
  randzs = randomu(seed, nrand,/normal)*sigzs + meanzs
  siginv = sigmacritinv(zlens, randzs)*1.e4
  mean_monte = mean(siginv)

  plothist, siginv, xsiginv, ysiginv, bin=0.001,$
            xtitle='siginv', ytitle='N';,/ylog,yrange=[1,3.e5]
;print,max(ysiginv)
  wzero=where(xsiginv EQ 0, nzero, comp=wnotzero)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Try the siginv histogram
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  norm = nzero + qgauss(ysiginv[wnotzero], xsiginv[wnotzero], 100)
  pofsiginv = ysiginv[wnotzero]/norm
  xsiginv = xsiginv[wnotzero]

  mean_sig = qgauss(pofsiginv*xsiginv, xsiginv, 100)

  print,'At mean zs: ',sigcritinv_meanzs
  print,'Mean integrated over z: ',sigcritinv_mean
  print,'Over monte carlo: ',mean_monte
  print,'Mean integrated over sig: ',mean_sig


  oplot, [sigcritinv_meanzs, sigcritinv_meanzs], [1, maxplot], color=!red
  oplot, [sigcritinv_mean, sigcritinv_mean], [1, maxplot], thick=5
  oplot, [mean_monte, mean_monte], [1, maxplot], color=!blue
  oplot, [mean_sig, mean_sig], [1, maxplot], color=!green


  legend,['At mean zs',$
          'over p(z)',$
          'over monte carlo', $
          'over p(sig)'],$
         line=[0,0,0,0],$
         color=[!red, !p.color, !blue, !green],/left,charsize=0.7

  !p.charsize=pold
  !p.multi=0


END 
