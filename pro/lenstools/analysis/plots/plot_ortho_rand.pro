PRO plot_ortho_rand, lenses, random, dops=dops

  IF keyword_set(dops) THEN BEGIN 
      psfile = '~/plots/rand_and_ortho.eps'
      begplot,name=psfile,/encap,ysize=6,xsize=7
  ENDIF 

  setup_mystuff

  xt=!mpcxtitle2
  topyt=!rdeltaytitle
  botyt=!randdeltaytitle

  xoplot = [1.e-6, 1.e6]
  yoplot = [0,0]

  ff = 1./1000.

  xrange = [15, 17000]*ff
  topyrange = [-2.5,2.5]
  botyrange = [-2.5,2.5]

  plottwoerr,lenses.meanr*ff,lenses.orthosig,lenses.orthosigerr,$
             random.meanr*ff,random.sigma,random.sigmaerr,$
             /xlog,psym=8,$
             xrange=xrange, $
             topyrange=topyrange,botyrange=botyrange,$
             xoplot=xoplot,yoplot=yoplot,$
             xtit=xt,topytit=topyt,botyt=botyt, $
             xstyle=1+2, botystyle=1+2,topystyle=1+2, $
             xtickf='loglabels', xticklen=0.06

  oplot,xoplot,yoplot

  IF keyword_set(dops) THEN endplot,/trim_bbox

END 
