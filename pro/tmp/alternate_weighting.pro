PRO alternate_weighting,lensum, noweight=noweight

  ;; remove the Nsource weighting from lensing measurements
  ;; or if /noweight, remove 1/sigma_crit weighting as well

  IF n_elements(lensum) EQ 0 THEN $
    lensum=mrdfits('/sdss5/data0/lensout/stripe10/sublum/r/lum_zgal_gal_stripe10_r_lensum_N1.fit',1)

  nbin=n_elements(lensum[0].rsum)

  combine_zlensum,lensum,80.,20.,1000.,1.,sumst,shst
  tlensum=lensum

  FOR i=0L, nbin-1 DO BEGIN 

      ww=where(tlensum.npair[i] GT 0,nww)
      IF nww NE 0 THEN BEGIN 

          IF keyword_set(noweight) THEN wdiv = 1./tlensum[ww].npair[i]/tlensum[ww].scritinv^2 $
          ELSE wdiv = 1./tlensum[ww].npair[i]

          tlensum[ww].wsum[i] = tlensum[ww].wsum[i]*wdiv

          tlensum[ww].tansigsum[i] = tlensum[ww].tansigsum[i]*wdiv
          tlensum[ww].tansigerrsum[i] = tlensum[ww].tansigerrsum[i]*wdiv^2
          tlensum[ww].radsigsum[i] = tlensum[ww].radsigsum[i]*wdiv
          tlensum[ww].radsigerrsum[i] = tlensum[ww].radsigerrsum[i]*wdiv^2
          
          tlensum[ww].etansum[i] = tlensum[ww].etansum[i]*wdiv
          tlensum[ww].etanerrsum[i] = tlensum[ww].etanerrsum[i]*wdiv^2
          tlensum[ww].eradsum[i] = tlensum[ww].eradsum[i]*wdiv
          tlensum[ww].eraderrsum[i] = tlensum[ww].eraderrsum[i]*wdiv^2
      ENDIF 

  ENDFOR 
  combine_zlensum,tlensum,80.,20.,1000.,1.,tsumst,tshst
  w=lindgen(nbin-1)

  yt=!tsym.delta_cap+!tsym.sigma_cap+'(R)'
  xt='Projected Radius (h!U'+!tsym.minus+'1!N kpc)'
  
  erase & multiplot,[1,2]
  ploterror,tshst.meanr[w],tshst.sigma[w]/tshst.ssh,$
    tshst.sigmaerr[w]/tshst.ssh,psym=1,ytit=yt
  oploterror,shst.meanr[w]+10,shst.sigma[w]/shst.ssh,shst.sigmaerr[w]/shst.ssh,psym=4,$
    color=!blue,errcolor=!blue
  oplot,[0,10000],[0,0]
  legend,['new','orig'],psym=[1,4],colors=[!black,!blue],/right,$
    thick=[!p.thick,!p.thick]

  ;key=get_kbrd(1)

  yt=!tsym.delta_cap+!tsym.sigma_cap+'('+!tsym.ltequal+'R)'

  multiplot
  ploterror,tshst.rmax_act[w],tshst.tsigma[w]/tshst.ssh,$
    tshst.tsigmaerr[w]/tshst.ssh,psym=1,xtit=xt,ytit=yt
  oploterror,shst.rmax_act[w]+10,shst.tsigma[w]/shst.ssh,shst.tsigmaerr[w]/shst.ssh,psym=4,$
    color=!blue,errcolor=!blue



  multiplot,/reset

  ;; find mean luminosity
  lensave, lensum, 'lum', meanlum, lumerr, element=2

  IF keyword_set(noweight) THEN lw = replicate(1.0, n_elements(tlensum)) $
  ELSE lw=tlensum.scritinv^2

  lwsum=total(lw)
  sum=total( lw*tlensum.lum[2] )
  tmeanlum = sum/lwsum
  errsum=total( lw^2*(tlensum.lum[2]-meanlum)^2 )
  tlumerr = sqrt(errsum)/lwsum
  
  print,'Luminosity means'
  print,meanlum,lumerr,tmeanlum,tlumerr
  print
  legend, ['new <L> = '+ntostr(tmeanlum/1.e10,4), 'old <L> = '+ntostr(meanlum/1.e10,4)],/right,textcolors=[!black,!blue]

  colprint,tshst.meanr,shst.sigma,shst.sigmaerr,tshst.sigma,tshst.sigmaerr
  print
  colprint,tshst.meanr,shst.tsigma/shst.tsigmaerr,tshst.tsigma/tshst.tsigmaerr,$
    (shst.tsigma/shst.tsigmaerr)/(tshst.tsigma/tshst.tsigmaerr)
END 
