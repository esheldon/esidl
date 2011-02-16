PRO plot_lens_absmag, struct,rxhist,ryhist,dops=dops, encap=encap, color=color, early=early, red=red, log=log

  IF keyword_set(log) THEN BEGIN 
      ymin = 7 
      yfac = 1.5
      ytickf = 'loglabels'
  ENDIF ELSE BEGIN 
      xticklen=0.02
      ymin=0
      yfac = 1
  ENDELSE 
  IF keyword_set(early) THEN BEGIN 
      earlystr = '_early'
      w=where(struct.eclass LT -0.02)
      yrange = [ymin,7000*yfac]
  ENDIF ELSE IF keyword_set(red) THEN BEGIN  
      earlystr = '_red'
      gmr = struct.abscounts[1] - struct.abscounts[2]
      w=where(gmr GT 0.7 AND $
              gmr GT 0.1 AND $
              gmr LT 1.1)
      yrange = [ymin,5000*yfac]
  ENDIF ELSE BEGIN 
      earlystr = ''
      w=lindgen(n_elements(struct))
      yrange = [ymin,9000*yfac]
  ENDELSE 

  IF keyword_set(color) THEN BEGIN 
      clrstr = '_color'
  ENDIF ELSE BEGIN 
      clrstr = ''
  ENDELSE 
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = '~/plots/absmag_hist'+clrstr+earlystr+'.eps'
          begplot,name=psfile,/color,/encap,xsize=7,ysize=7
      ENDIF ELSE BEGIN 
          psfile = '~/plots/absmag_hist'+clrstr+earlystr+'.ps'
          begplot,name=psfile,/color
      ENDELSE 
  ENDIF 

  IF keyword_set(color) THEN BEGIN 
      colors = [!blue, !green, !red, !magenta, !p.color]
      ;;nes = [0,0,0,0,0]
      lines = [1, 2, 3, 4, 0]
  ENDIF ELSE BEGIN 
      ;colors = [!grey0, !grey40, !grey60, !grey80, !grey90]
      colors = replicate(!p.color, 5)
      lines = [1, 2, 3, 4, 0]
  ENDELSE

  setup_mystuff
  xtitle = 'M - 5 log(h)'
  IF keyword_set(log) THEN ytitle = 'log (N)' ELSE ytitle = 'N(M)'

  bin = 0.15
  xmin = -26
  xmax = -14
  xrange = [xmax,xmin]
  aplot,!gratio,[0],xrange=xrange,yrange=yrange,$
        xtitle=xtitle, ytitle=ytitle,xstyle=3,ystyle=1,$
        xticklen=xticklen, ylog=log, ytickf=ytickf

  plothist,struct[w].absmag[0],uxhist,uyhist,bin=bin,min=xmin,max=xmax,/noplot
  plothist,struct[w].absmag[1],gxhist,gyhist,bin=bin,min=xmin,max=xmax,/noplot
  plothist,struct[w].absmag[2],rxhist,ryhist,bin=bin,min=xmin,max=xmax,/noplot
  plothist,struct[w].absmag[3],ixhist,iyhist,bin=bin,min=xmin,max=xmax,/noplot
  plothist,struct[w].absmag[4],zxhist,zyhist,bin=bin,min=xmin,max=xmax,/noplot

  oplot, uxhist, uyhist, color=colors[0],line=lines[0]
  oplot, gxhist, gyhist, color=colors[1],line=lines[1]
  oplot, rxhist, ryhist, color=colors[2],line=lines[2]
  oplot, ixhist, iyhist, color=colors[3],line=lines[3]
  oplot, zxhist, zyhist, color=colors[4],line=lines[4]

  idit_absmin = -22.2
  idit_absmax = -18.9

  ww = where(struct[w].absmag[2] LE idit_absmax AND $
             struct[w].absmag[2] GE idit_absmin, nww)
  print,float(nww)/n_elements(w)

  rdiff = abs(rxhist-idit_absmin)
  wmin  = where( rdiff EQ min(rdiff))
  rdiff = abs(rxhist-idit_absmax)
  wmax  = where(rdiff EQ min(rdiff))

  oplot, [rxhist[wmin]], [ryhist[wmin]], psym=8, color=colors[2]
  oplot, [rxhist[wmax]], [ryhist[wmax]], psym=8, color=colors[2]

  oplot,[-1000,1000],[0,0]

  legend,!colors,line=lines,colors=colors,box=0,charsize=1,thick=replicate(!p.thick,5)

  IF keyword_set(dops) THEN endplot
  IF keyword_set(encap) THEN set_bbox,psfile,'%%BoundingBox: 25 15 480 345'

  outdir = '~/lensout/absmag/'

  ufile = outdir+'umag_dist.dat'
  gfile = outdir+'gmag_dist.dat'
  rfile = outdir+'rmag_dist.dat'  
  ifile = outdir+'imag_dist.dat'
  zfile = outdir+'zmag_dist.dat'

  colprint, uxhist, uyhist, file=ufile
  colprint, gxhist, gyhist, file=gfile
  colprint, rxhist, ryhist, file=rfile
  colprint, ixhist, iyhist, file=ifile
  colprint, zxhist, zyhist, file=zfile

  print,total(uyhist),total(gyhist),total(ryhist),total(iyhist),total(zyhist)
  print,total(uyhist)/total(ryhist),total(gyhist)/total(ryhist),total(ryhist)/total(ryhist),total(iyhist)/total(ryhist),total(zyhist)/total(ryhist)

END 
