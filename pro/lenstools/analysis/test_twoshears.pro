

PRO test_twoshears, clr, struct1, rstruct1, struct2, rstruct2,type=type

  !p.multi=[0,0,2]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; stuff for fitting and plotting
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  chytitle='Norm'
  chxtitle='Power'

  xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)'
  powrange=[-1.5,-0.1]
  xrange=[20,2400]
  npow = 200L
  nnorm = 200L

  IF n_elements(type) EQ 0 THEN type=''

  ;; shear stuff
;  sigytitle='Shear'
;  yrange=[0.00003,0.02]
;  normrange=[5.e-5,6.e-4]
;  tag='shear'
;  errtag='shearerr'

  ;; sigma stuff
  sigytitle = '!S'+!tsym.sigma_cap+'!R!A'+!tsym.minus+'!N ('+!tsym.ltequal+'R) '+$
    !tsym.minus+' !S'+!tsym.sigma_cap+'!R!A'+!tsym.minus+$
    '!N (R) (h M'+sunsymbol()+' pc!U'+!tsym.minus+'2!N)'
  yrange=[0.05,200]
  normrange=[0.0,6.0]
  tag='sigma'
  errtag='sigmaerr'

  ;;;;;;;;;;;;;;;;;;;;;;
  ;; first data set
  ;;;;;;;;;;;;;;;;;;;;;;

  stripes = ['10','82']
  nend = ['N4.fit','N4.fit']
  title='Stripes '+stripes[0]+' '+stripes[1]+' '+!colorsp[clr]

  nuse=0
  stripes=stripes[nuse]
  nend = nend[nuse]
  title='Stripe '+stripes+' '+!colorsp[clr]

  indir = '/sdss5/data0/lensout/stripe'+stripes+'/'

  combfiles = indir+type+'zgal_gal_stripe'+stripes+'_'+!colors[clr]+'_sum_'+nend
  rcombfiles = indir+type+'zrand_stripe'+stripes+'_'+!colors[clr]+'_sum_'+nend

;  combfiles = [indir+'zgal_gal_stripe'+stripes+'_r_sum_'+nend, $
;               indir+'zgal_gal_stripe'+stripes+'_g_sum_'+nend]
;  rcombfiles = [indir+'zrand_stripe'+stripes+'_r_sum_'+nend, $
;                indir+'zrand_stripe'+stripes+'_g_sum_'+nend]

  combine_zshear, combfiles, struct1
  combine_zshear, rcombfiles, rstruct1

  correct_shear, struct1, rstruct1

  nr=n_elements(struct1.meanr)
  indtmp=lindgen(nr)
;  wuse1=where( (indtmp GT 1) AND (struct1.sigma GT 0.) )
  wuse1=(lindgen(nr))[1:nr-1]

  IF NOT tag_exist(struct1, tag, index=tg) THEN message,'what!'
  IF NOT tag_exist(struct1, errtag, index=etg) THEN message,'what!'

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; second data set
  ;;;;;;;;;;;;;;;;;;;;;;;

  stripes = ['10','82']
  nend = ['N6.fit','N6.fit']

  stripes=stripes[nuse]
  nend = nend[nuse]

  indir = '/sdss5/data0/lensout/stripe'+stripes+'/'

  combfiles = indir+type+'zgal_gal_stripe'+stripes+'_'+!colors[clr]+'_sum_'+nend
  rcombfiles = indir+type+'zrand_stripe'+stripes+'_'+!colors[clr]+'_sum_'+nend

;  combfiles = [indir+'zgal_gal_stripe'+stripes+'_r_sum_'+nend, $
;               indir+'zgal_gal_stripe'+stripes+'_g_sum_'+nend]
;  rcombfiles = [indir+'zrand_stripe'+stripes+'_r_sum_'+nend, $
;                indir+'zrand_stripe'+stripes+'_g_sum_'+nend]

  combine_zshear, combfiles, struct2
  combine_zshear, rcombfiles, rstruct2

  correct_shear, struct2, rstruct2

  nr=n_elements(struct2.meanr)
  indtmp=lindgen(nr)
;  wuse2=where( (indtmp GT 1) AND (struct2.sigma GT 0.) )
  wuse2=(lindgen(nr))[1:nr-1]

  forprint,struct1.meanr,struct1.sigma,struct1.sigmaerr,$
    struct2.sigma,struct2.sigmaerr

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Fits to first
  ;;;;;;;;;;;;;;;;;;;;;;;;

  pow_chisq_conf_gen, $
    struct1.meanr[wuse1]/1000., struct1.(tg)[wuse1], struct1.(etg)[wuse1],$
    powrange, normrange, npow, nnorm, chisq_surf,$
    pmin1, nmin1, plow1, phigh1, nlow1,  nhigh1, yfit=yfit1,$
    ytitle=chytitle,xtitle=chxtitle,aspect=1,$
    title=title
  range2error, plow1, pmin1, phigh1, perrhigh1, perrlow1
  range2error, nlow1, nmin1, nhigh1, nerrhigh1, nerrlow1

  nkeep=[6,5]
  mean_error_legend,['Norm1: ','Pow1: '], $
    [nmin1,pmin1],[nerrlow1[0],perrlow1[0]],[nerrhigh1[0],perrhigh1[0]],$
    nkeep=nkeep,/right,box=0
  print,'Norm1: ',ntostr(nmin1),' + ',ntostr(nerrhigh1[0]),' - ',ntostr(nerrlow1[0])
  print,'Pow1: ',ntostr(pmin1),' + ',ntostr(perrhigh1[0]),' - ',ntostr(perrlow1[0])

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Fits to second
  ;;;;;;;;;;;;;;;;;;;;;;;;

  pow_chisq_conf_gen, $
    struct2.meanr[wuse2]/1000., struct2.(tg)[wuse2], struct2.(etg)[wuse2],$
    powrange, normrange, npow, nnorm, chisq_surf,$
    pmin2, nmin2, plow2, phigh2, nlow2,  nhigh2, yfit=yfit2,$
    ytitle=chytitle,xtitle=chxtitle,aspect=1,$
    title=title
  range2error, plow2, pmin2, phigh2, perrhigh2, perrlow2
  range2error, nlow2, nmin2, nhigh2, nerrhigh2, nerrlow2

  mean_error_legend,['Norm2: ','Pow2: '], $
    [nmin2,pmin2],[nerrlow2[0],perrlow2[0]],[nerrhigh2[0],perrhigh2[0]],$
    nkeep=nkeep,/right,box=0
  print,'Norm2: ',ntostr(nmin2),' + ',ntostr(nerrhigh2[0]),' - ',ntostr(nerrlow2[0])
  print,'Pow2: ',ntostr(pmin2),' + ',ntostr(perrhigh2[0]),' - ',ntostr(perrlow2[0])
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plot first data and random
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  aploterror, 1, $
    struct1.meanr[wuse1], struct1.(tg)[wuse1], struct1.(etg)[wuse1], $
    psym=1, /xlog, /ylog, xrange=xrange, yrange=yrange, ystyle=1,xstyle=1,$
    xtitle=xtitle,ytitle=sigytitle,$
    title=title
  
  oplot,struct2.meanr[wuse2],struct2.(tg)[wuse2],psym=4,color=!blue
  oplot,struct1.meanr[wuse1],yfit1

  oploterror, $
    rstruct1.meanr[wuse1], rstruct1.(tg)[wuse1], rstruct1.(etg)[wuse1],$
    psym=4, color=!red, errcolor=!red
;  oplot,rstruct1.meanr,rstruct1.shear,psym=4,color=!red

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plot second data and random
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;IF display_type() EQ 'X' THEN key=get_kbrd(1)

  aploterror, 1, $
    struct2.meanr[wuse2], struct2.(tg)[wuse2], struct2.(etg)[wuse2], $
    psym=1, /xlog, /ylog, xrange=xrange, yrange=yrange, ystyle=1,xstyle=1,$
    xtitle=xtitle,ytitle=sigytitle,$
    title=title
  
  oplot,struct2.meanr[wuse2],yfit2

  oploterror, $
    rstruct2.meanr[wuse2], rstruct2.(tg)[wuse2], rstruct2.(etg)[wuse2],$
    psym=4, color=!red, errcolor=!red
;  oplot,rstruct2.meanr,rstruct2.shear,psym=4,color=!red


  !p.multi=0

END 
