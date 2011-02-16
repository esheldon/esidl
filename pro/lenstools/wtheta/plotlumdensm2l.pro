PRO plotlumdensm2l, t;, lumclr

  ;indir = '/sdss5/data0/lensout/mass2light/'
  ;file = indir+'lumdens_matchrN6_wthetalumweq_stripe10_sum_'+!colors[lumclr]+'w_N5_m2l_omega_onlyall.fit'

  ;t=mrdfits(file,1)

  ;help,t,/str

  myusersym, 'fill_circle'
  psym=8

  xtitle = !kpcxtitle2
  ytitle='M(< R) [10!U12!N h!U'+!tsym.minus+'1!N M'+sunsymbol()+']'
  aploterror,1, t.rmax,t.massap/1.e12,t.massaperr/1.e12,psym=psym,$
    /xlog,/ylog,xrange=[20,2400],yrange=[0.01,50.0],/ystyle,/xstyle, $
    xtitle=xtitle,ytitle=ytitle, xticklen=0.04, yticklen=0.04

  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  ytitle = 'L(< R) [10!U10!N h!U'+!tsym.minus+$
    '2!NL'+sunsymbol()+']'
  aplot,!gratio, t.rmax,t.lumap/1.e10,xtitle=xtitle,ytitle=ytitle,yminor=5

  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  ytitle='M/L(< R) [h M'+sunsymbol()+'/L'+sunsymbol()+']'
  aploterror, !gratio, t.rmax,t.massap/t.lumap,t.massaperr/t.lumap,psym=psym,$
              /xlog,xrange=[20,3000],xtitle=xtitle,ytitle=ytitle,/xstyle,$
              yminor=5,xticklen=0.05

  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  ytitle='M/L(<R) [h M'+sunsymbol()+'/L'+sunsymbol()+']'
  aploterror, !gratio, t.rmax,t.masstot/t.lumtot,t.masstoterr/t.lumtot,psym=psym,$
    xtitle=xtitle,ytitle=ytitle,xticklen=0.03
  oplot, t.rmax, t.m2lfunc*100.

  print,'Asymptotic: ',t.M0/t.L0*100.

END 
