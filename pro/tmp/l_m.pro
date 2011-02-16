PRO l_m, lm, noplot=noplot

  IF n_params() EQ 0 THEN BEGIN
      print,'-Syntax: l_m, lm, noplot=noplot'
      return
  ENDIF 

  IF NOT keyword_set(noplot) THEN noplot=0

  file = '/sdss3/usrdevel/esheldon/idl.lib/tmp/zams.par'

  rdyanny, file, lm

  IF NOT noplot THEN BEGIN
      
      mess=['0.1 Msun', '1.0 Msun', '10. Msun']
      w=where(lm.msun EQ 1.)
      w10 = where(lm.msun EQ 10)
      wp1 = where(lm.msun EQ .1)
      psym=[1,2,4]

      yt='Log(L/Lsun)'
      xt = 'Log(M/Msun)'
      t='ZAMS'
      aplot,1.,alog10(lm.msun),lm.loglsun,xtitle=xt,ytitle=yt,title=t
      oplot,[alog10(lm[wp1].msun)],[lm[wp1].loglsun],psym=psym[0]
      oplot,[alog10(lm[w].msun)],[lm[w].loglsun],psym=psym[1]
      oplot,[alog10(lm[w10].msun)],[lm[w10].loglsun],psym=psym[2]
      legend,mess,psym=psym

      key=get_kbrd(1)
      yt = 'Log(L/Lsun)'
      xt = 'log(Te) (10^6 k)'
      xrange = [1.7,0.4]

      aplot, 1., alog10(lm.t6), lm.loglsun,xrange=xrange,xtitle=xt,ytitle=yt,title=t,xstyle=1
      oplot,[alog10(lm[wp1].t6)],[lm[wp1].loglsun],psym=psym[0]
      oplot,[alog10(lm[w].t6)],[lm[w].loglsun],psym=psym[1]
      oplot,[alog10(lm[w10].t6)],[lm[w10].loglsun],psym=psym[2]
      legend,mess,psym=psym,position=[.7,5]


  ENDIF 

return 
END 
