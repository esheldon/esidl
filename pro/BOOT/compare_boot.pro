PRO compare_boot, clr, ext=ext

  ;; compare regular errors and bootstrap
  ;; do some fits

  IF n_elements(ext) EQ 0 THEN ext='N1.fit'

  bootfile = $
    '/sdss5/data0/lensout/stripe10/wthetalumweq_stripe10_bootstrap_'+$
    !colors[clr]+'w_'+ext
  file = '/sdss5/data0/lensout/stripe10/wthetalumweq_stripe10_'+$
    !colors[clr]+'w_'+ext
  rfile = '/sdss5/data0/lensout/stripe10/wthetarandlumweq_stripe10_'+$
    !colors[clr]+'w_'+ext
  
  boot=mrdfits(bootfile,1)
  t=mrdfits(file,1)
  tr=mrdfits(rfile,1)
  diff = t.meanlum-tr.meanlum
  differr = sqrt(t.meanlumerr^2 + tr.meanlumerr^2)
  
  w=where(t.meanr GT 200)

  CASE clr OF
      1: BEGIN 
          normrange=[1.2,1.75]
          powrange=[-0.85,-0.7]
      END 
      2: BEGIN
          normrange=[1.6,2.0]
          powrange=[-0.85,-0.7]
      END 
      3: BEGIN
          normrange=[1.75,2.25]
          powrange=[-0.85,-0.7]
      END 
      4: BEGIN
          normrange=[2.0,2.5]
          powrange=[-0.85,-0.7]
      END 
  ENDCASE 
  npow=200L
  nnorm=200L

  !p.multi=[0,0,2]
  pow_chisq_conf_gen, t.meanr/1000., diff, differr, powrange, normrange,$
    npow, nnorm, ch, dpmin, dnmin, dpowlow, dpowhigh, dnormlow, dnormhigh,$
    yfit=dyfit,/center,aspect=1, title='Regular errors',$
    xtitle='power index', ytitle='Normalization',$
    perrlow=dperrlow,perrhigh=dperrhigh,$
    nerrlow=dnerrlow,nerrhigh=dnerrhigh,minchisq=minchisq,degfree=degfree,wuse=w
  mean_error_legend, ['index','norm'],$
    [dpmin,dnmin],[dperrlow,dnerrlow],[dperrhigh,dnerrhigh],/right,box=0,$
    digit=[3,3],nkeep=[6,5],/clear
  mess=!tsym.chi+'!U2!N/'+!tsym.nu+' = '+ntostr(minchisq,5)+'/'+ntostr(degfree)+$
    ' = '+ntostr(minchisq/degfree, 4)
  legend,mess,/right,/bottom,box=0
  legend,!colorsp[clr],/left,box=0

  pow_chisq_conf_gen, t.meanr/1000., boot.meanlum, boot.covariance, $
    powrange, normrange,$
    npow, nnorm, $
    ch, bpmin, bnmin, bpowlow, bpowhigh, bnormlow, bnormhigh,$
    yfit=byfit,/covariance,/center,aspect=1, title='Bootstrap and covariance',$
    xtitle='power index', ytitle='Normalization',$
    perrlow=bperrlow,perrhigh=bperrhigh,$
    nerrlow=bnerrlow,nerrhigh=bnerrhigh,minchisq=minchisq,degfree=degfree,wuse=w
  mean_error_legend, ['index','norm'],$
    [bpmin,bnmin],[bperrlow,bnerrlow],[bperrhigh,bnerrhigh],/right,box=0,$
    digit=[3,3],nkeep=[6,5],/clear
  mess=!tsym.chi+'!U2!N/'+!tsym.nu+' = '+ntostr(minchisq,5)+'/'+ntostr(degfree)+$
    ' = '+ntostr(minchisq/degfree, 4)
  legend,mess,/right,/bottom,box=0
  legend,!colorsp[clr],/left,box=0

  !p.multi=0
  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  xtickadd = [200,300,400,500,2000]
  ytickadd = [2,3,4,5]
  aploterror, 1.0, t.meanr[w],diff[w],differr[w],$
    psym=1,xran=[150,2500],yran=[0.8,11],/xstyle,/ystyle,/xlog,/ylog,$
    xtit=!kpcxtitle,ytit=!colorsp[clr]+' '+!lumytitle,hat=0
  oploterror, t.meanr[w], boot.meanlum[w], boot.meanlumerr[w], $
    psym=1,color=!red,errc=!red
  add_labels, xtickv=xtickadd, ytickv=ytickadd

  legend,['Regular Errors','Bootstrap'],line=[0,0], colors=[!p.color,!red],$
    /right,box=0,thick=[!p.thick,!p.thick]

  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  aploterror, 1.0, t.meanr[w],boot.meanlum[w], boot.meanlumerr[w],$
    psym=1,xran=[150,2500],yran=[0.8,11],/xstyle,/ystyle,/xlog,/ylog,$
    xtit=!kpcxtitle,ytit=!colorsp[clr]+' '+!lumytitle,hat=0,$
    title='Bootstrap and covariance'
  oplot, t.meanr[w], byfit[w]
  add_labels, xtickv=xtickadd, ytickv=ytickadd

END 
