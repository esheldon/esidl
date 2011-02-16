PRO plot_lens_eclass_gmr, struct

  psfile = '~/plots/eclass_gmr_hist.eps'
  begplot,name=psfile,xsize=7,ysize=7,/color,/encap

;  ext=!csym.minus+'ECLASS'
;  eyt='N('+!csym.minus+'ECLASS)'
  ext='ECLASS'
  eyt='N(ECLASS)'

  gmrxt = 'g'+!csym.minus+'r'
  gmryt = 'N(g'+!csym.minus+'r)'

  gmr_cut = 0.7
  eclass_cut = -0.06

  gmr = struct.abscounts[1] - struct.abscounts[2]
  eclass = struct.eclass

  wgmr = where(struct.z GT 0.02 AND struct.z LT 0.3 AND $
               gmr GT 0.1 AND gmr LT 1.1)
  we = where(struct.z GT 0.02 AND struct.z LT 0.3)

  medgmr = median(gmr[wgmr])
  mede = median(eclass[we])

  wgmrlow = where(gmr[wgmr] LT gmr_cut, comp=wgmrhigh, ngmrlow, ncomp=ngmrhigh) 
  wgmrlow = wgmr[wgmrlow]
  wgmrhigh = wgmr[wgmrhigh]

  wehigh = where(eclass[we] LT eclass_cut, comp=welow, nehigh, ncomp=nelow) 
  welow = we[welow]
  wehigh = we[wehigh]

  print,medgmr,mede
  print,gmr_cut,eclass_cut
  print,gmr_cut,ngmrlow,ngmrhigh
  print,eclass_cut,nelow,nehigh

;  !p.multi=[0,0,2]

  pxsize = 0.75                 ; in normalized coordinates
  d_pxsize = pxsize*!d.x_vsize  ; device coordinates
  d_pysize = d_pxsize             ; device coordinates
  pysize = d_pysize/!d.y_vsize  ; normalized coords

  slack = 0.05

  minx = 0.16
  maxx = minx + pxsize

  bminy = 0.085 
  bmaxy = bminy + (pysize/2.0 - slack)
  plothist,eclass[we],exhist,eyhist,bin=0.01,min=-0.3,max=0.7,xrange=[0.6,-0.35],xtit=ext,ytit=eyt,/xsty, $
           position=[minx,bminy,maxx,bmaxy],yrange=[0,9500],/ystyl,xticklen=0.04
  wvdis = where(struct.vel_dis GT 50 AND struct.vel_dis LT 400)

  plothist,eclass[wvdis],bin=0.01,min=-0.3,max=0.7,/overplot, color=!grey35
  oplot,[eclass_cut,eclass_cut],[0,1.e6],color=!grey50

  miny = bmaxy + 2.0*slack
  maxy = bminy + pysize
  plothist,gmr[wgmr],gxhist,gyhist,min=-0.5,max=1.2,bin=0.01,xrange=[0,1],ytit=gmryt, $
           position=[minx,miny,maxx,maxy],/noerase,xticklen=0.04
  oplot,[gmr_cut,gmr_cut],[0,1.e6],color=!grey50

  xx = (maxx-minx)/2.0 + 0.132
  yy = bminy + (pysize/2.0-slack/4.)
  xyouts,xx,yy,gmrxt,/normal
  print,xx,yy
  endplot
  set_bbox,psfile,'%%BoundingBox: 10 0 470 430'

END 
