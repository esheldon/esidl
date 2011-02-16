PRO test_systematic, clr, rstruct, rrstruct, diff, differr

  !p.multi=[0,0,2]

  chytitle='Norm'
  chxtitle='Power'

  xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)'
  sigytitle = '!S'+!tsym.sigma_cap+'!R!A'+!tsym.minus+'!N ('+!tsym.ltequal+'R) '+$
    !tsym.minus+' !S'+!tsym.sigma_cap+'!R!A'+!tsym.minus+$
    '!N (R) (h M'+sunsymbol()+' pc!U'+!tsym.minus+'2!N)'
  xrange=[10,11000]
  yrange=[0.05,200]

  stripes = ['10','82']
  nend = ['N8.fit','N4.fit']
  title='Stripes '+stripes[0]+' '+stripes[1]

  stripes=stripes[0]
  nend = nend[0]
  title='Stripe '+stripes

  indir = '/sdss5/data0/lensout/stripe'+stripes+'/'

  combfiles = indir+'zgal_gal_stripe'+stripes+'_'+!colors[clr]+'_sum_'+nend
  rcombfiles = indir+'zrand_stripe'+stripes+'_'+!colors[clr]+'_sum_'+nend

;  combfiles = [indir+'zgal_gal_stripe'+stripes+'_r_sum_'+nend, $
;               indir+'zgal_gal_stripe'+stripes+'_g_sum_'+nend]
;  rcombfiles = [indir+'zrand_stripe'+stripes+'_r_sum_'+nend, $
;                indir+'zrand_stripe'+stripes+'_g_sum_'+nend]

  combine_zshear, combfiles, rstruct
  combine_zshear, rcombfiles, rrstruct

  frac_overdense, rstruct.npair,  rstruct.nlenses, $
    rrstruct.npair, rrstruct.nlenses, rfrac, rfracerr,$
    rndiff, rndifferr

  rcorr = 1. + (rfrac > 0.)
  applycorr, rstruct, rstruct.Ssh, rcorr, rcorr ;last is dummy

  ;; subtract signal in random points
  diff = rstruct.sigma-rrstruct.sigma
  differr = sqrt(rstruct.sigmaerr^2 + rrstruct.sigmaerr^2)

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Fits
  ;;;;;;;;;;;;;;;;;;;;;;;;

  siguess = [40., -.7]
  fitpower, rstruct.meanr/1000., rstruct.sigma, rstruct.sigmaerr, $
    siguess, yfit, aout, aouterr

  powrange=[-.75,-0.25]
  npow = 100
  normrange=[3.0,10.0]
  nnorm = 100
  pow_chisq_conf_gen, rstruct.meanr/1000., rstruct.sigma, rstruct.sigmaerr,$
    powrange, normrange, npow, nnorm, chisq_surf,$
    pmin, nmin, plow, phigh, nlow,  nhigh, yfit=yfit,$
    ytitle=chytitle,xtitle=chxtitle,aspect=1,$
    title=title
  range2error, plow, pmin, phigh, perrhigh, perrlow
  range2error, nlow, nmin, nhigh, nerrhigh, nerrlow

  nkeep=[5,5]
  mean_error_legend,['Norm: ','Pow: '], $
    [nmin,pmin],[nerrlow[0],perrlow[0]],[nerrhigh[0],perrhigh[0]],$
    nkeep=nkeep,/right,box=0
  print,'Norm: ',ntostr(nmin),' + ',ntostr(nerrhigh[0]),' - ',ntostr(nerrlow[0])
  print,'Pow: ',ntostr(pmin),' + ',ntostr(perrhigh[0]),' - ',ntostr(perrlow[0])

  fitpower, rstruct.meanr/1000., diff, differr, $
    siguess, diffyfit, diffaout, diffaouterr

  powrange=[-1.0,-0.5]
  npow = 100
  normrange=[2.0,6.0]
  nnorm = 100
  pow_chisq_conf_gen, rstruct.meanr/1000., diff, differr,$
    powrange, normrange, npow, nnorm, chisq_surf, $
    fpmin, fnmin, fplow, fphigh, fnlow,  fnhigh, yfit=diffyfit,$
    ytitle=chytitle,xtitle=chxtitle,aspect=1
  range2error, fplow, fpmin, fphigh, fperrhigh, fperrlow
  range2error, fnlow, fnmin, fnhigh, fnerrhigh, fnerrlow

  mean_error_legend,['Norm: ','Pow: '], $
    [fnmin,fpmin],[fnerrlow[0],fperrlow[0]],[fnerrhigh[0],fperrhigh[0]],$
    nkeep=nkeep,/right,box=0
  print,'Norm: ',ntostr(fnmin),' + ',ntostr(fnerrhigh[0]),' - ',ntostr(fnerrlow[0])
  print,'Pow: ',ntostr(fpmin),' + ',ntostr(fperrhigh[0]),' - ',ntostr(fperrlow[0])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plot data and random
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  aploterror, 1, rstruct.meanr, rstruct.sigma, rstruct.sigmaerr, $
    psym=1, /xlog, /ylog, xrange=xrange, yrange=yrange, ystyle=1,xstyle=1,$
    xtitle=xtitle,ytitle=sigytitle,$
    title=title
  
  oplot,rstruct.meanr,yfit

  oploterror, rrstruct.meanr, rrstruct.sigma, rrstruct.sigmaerr,$
    psym=4, color=!red, errcolor=!red
;  oplot,rrstruct.meanr,rrstruct.sigma,psym=4,color=!red

  ;; plot subtracted

  aploterror, 1, rstruct.meanr, diff, differr, $
    psym=1, /xlog, /ylog, xrange=xrange, yrange=yrange, ystyle=1,xstyle=1,$
    xtitle=xtitle,ytitle=sigytitle
  oplot,rstruct.meanr,diffyfit
 
  forprint,rstruct.meanr,rstruct.sigma,diff,rstruct.sigma/diff,diff/differr


  !p.multi=0

END 
