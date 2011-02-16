PRO plot_test_error_formulae, t

  setupplot,dtype
  nc=n_elements(t)

  GOTO,jump1
  FOR i=0L, nc-1 DO BEGIN 

      w=where(t[i].e1 NE 0.0,nw)

      plothist,t[i].e1[w],bin=.01
      oplot,[t[i].inpute1,t[i].inpute1],[0,100000]
      IF dtype EQ 'X' THEN key=get_kbrd(1)

  ENDFOR 

jump1:

  erase & multiplot,[1,2]
  xtitle='S/N'
  yrange=prange(t.meane1,t.mease1err)
  yrange[0] = yrange[0] < 0
  ploterror,t.s2n,t.meane1,t.mease1err,psym=1,$
    ytitle='e!D1!N',yrange=yrange
  oplot,t.s2n,t.inpute1
  legend,['meas','input'],psym=[1,0],/right,thick=[!p.thick,!p.thick]
  multiplot
  yrange=prange(t.meane2,t.mease2err)
  yrange[0] = yrange[0] < 0
  ploterror,t.s2n,t.meane2,t.mease2err,psym=1,$
    ytitle='e!D2!N',xtitle=xtitle,yrange=yrange
  oplot,t.s2n,t.inpute2
  multiplot,/reset

  IF dtype EQ 'X' THEN key=get_kbrd(1)

  erase & multiplot,[1,2]
  ;; Gary's formulae
  ytitle=!tsym.sigma+'!S!De!L1!R!UMeas!N / '+!tsym.sigma+'!S!De!L1!R!UAdap'
  ratio = t.mease1err/t.meanmomerr
  err=ratio*(t.meanmomerrerr/t.meanmomerr)

  yrange=prange(ratio,err)
  yrange[0] = yrange[0] < 0
  yrange[1] = yrange[1] > 1
  ploterror,t.s2n,ratio,err,psym=1,$
    ytitle=ytitle,yrange=yrange
  oplot,t.s2n,replicate(1.0,nc)

  multiplot
  ytitle=!tsym.sigma+'!S!De!L2!R!UMeas!N / '+!tsym.sigma+'!S!De!L2!R!UAdap'
  ratio = t.mease2err/t.meanmomerr
  err=ratio*(t.meanmomerrerr/t.meanmomerr)
  yrange=prange(ratio,err)
  yrange[0] = yrange[0] < 0
  yrange[1] = yrange[1] > 1
  ploterror,t.s2n,ratio,err,psym=1,$
    xtitle=xtitle,ytitle=ytitle,yrange=yrange
  oplot,t.s2n,replicate(1.0,nc)
  multiplot,/reset

  IF dtype EQ 'X' THEN key=get_kbrd(1)

  erase & multiplot,[1,2]
  ;; Our formula without covariance
  ytitle=!tsym.sigma+'!S!De!L1!R!UMeas!N / '+!tsym.sigma+'!S!De!L1!R!U1'
  ratio = t.mease1err/t.meane1err
  err=ratio*(t.meane1errerr/t.meane1err)
  yrange=prange(ratio,err)
  yrange[0] = yrange[0] < 0
  yrange[1] = yrange[1] > 1
  ploterror,t.s2n,ratio,err,psym=1,$
    ytitle=ytitle,yrange=yrange
  oplot,t.s2n,replicate(1.0,nc)

  multiplot
  ytitle=!tsym.sigma+'!S!De!L2!R!UMeas!N / '+!tsym.sigma+'!S!De!L2!R!U1'
  ratio = t.mease2err/t.meane2err
  err=ratio*(t.meane2errerr/t.meane2err)
  yrange=prange(ratio,err)
  yrange[0] = yrange[0] < 0
  yrange[1] = yrange[1] > 1
  ploterror,t.s2n,ratio,err,psym=1,$
    xtitle=xtitle,ytitle=ytitle,yrange=yrange
  oplot,t.s2n,replicate(1.0,nc)
  multiplot,/reset


;  IF dtype EQ 'X' THEN key=get_kbrd(1)

;  erase & multiplot,[1,2]
;  ;; Our formula without covariance
;  ytitle=!tsym.sigma+'!S!De!L1!R!UMeas!N / '+!tsym.sigma+'!S!De!L1!R!UFix'
;  fixfac = 2./sqrt(1.+2.*t.inpute1^2)
;  ratio = t.mease1err/t.meane1err/fixfac
;  err=ratio*(t.meane1errerr/t.meane1err)
;  yrange=prange(ratio,err)
;  yrange[0] = yrange[0] < 0
;  yrange[1] = yrange[1] > 1
;  ploterror,t.s2n,ratio,err,psym=1,$
;    ytitle=ytitle,yrange=yrange
;  oplot,t.s2n,replicate(1.0,nc)

;  multiplot
;  ytitle=!tsym.sigma+'!S!De!L2!R!UMeas!N / '+!tsym.sigma+'!S!De!L2!R!UFix'
;  fixfac = 2./sqrt(1.+2.*t.inpute2^2)
;  ratio = t.mease2err/t.meane2err/fixfac
;  err=ratio*(t.meane2errerr/t.meane2err)
;  yrange=prange(ratio,err)
;  yrange[0] = yrange[0] < 0
;  yrange[1] = yrange[1] > 1
;  ploterror,t.s2n,ratio,err,psym=1,$
;    xtitle=xtitle,ytitle=ytitle,yrange=yrange
;  oplot,t.s2n,replicate(1.0,nc)
;  multiplot,/reset



  IF dtype EQ 'X' THEN key=get_kbrd(1)

  erase & multiplot,[1,2]
  ;; Our formula with covariance
  ytitle=!tsym.sigma+'!S!De!L1!R!UMeas!N / '+!tsym.sigma+'!S!De!L1!R!U2'
  ratio = t.mease1err/t.meane1err2
  err=ratio*(t.meane1err2err/t.meane1err2)
  yrange=prange(ratio,err)
  yrange[0] = yrange[0] < 0
  yrange[1] = yrange[1] > 1
  ploterror,t.s2n,ratio,err,psym=1,$
    ytitle=ytitle,yrange=yrange
  oplot,t.s2n,replicate(1.0,nc)

  multiplot
  ytitle=!tsym.sigma+'!S!De!L2!R!UMeas!N / '+!tsym.sigma+'!S!De!L2!R!U2'
  ratio = t.mease2err/t.meane2err2
  err=ratio*(t.meane2err2err/t.meane2err2)
  yrange=prange(ratio,err)
  yrange[0] = yrange[0] < 0
  yrange[1] = yrange[1] > 1
  ploterror,t.s2n,ratio,err,psym=1,$
    xtitle=xtitle,ytitle=ytitle,yrange=yrange
  oplot,t.s2n,replicate(1.0,nc)
  multiplot,/reset

  IF dtype EQ 'X' THEN key=get_kbrd(1)

  ;; Measured Err e1/ Err e2
  ytitle=!tsym.sigma+'!S!De!L1!R!UMeas!N / '+!tsym.sigma+'!S!De!L2!R!UMeas'
  ratio = t.mease1err/t.mease2err
  aplot,!gratio,t.s2n,ratio,psym=1,$
    xtitle=xtitle,ytitle=ytitle
  oplot,t.s2n,replicate(1.0,nc)


END 
