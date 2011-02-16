PRO plot_lensw_vs_lum, sigc=sigc

  nclr=n_elements(!colors)
  setupplot,dtype

  yt='Lens weight'
  xt='Luminosity (10!U10!N L!DSun!N)'
  FOR i=1L, nclr-1 DO BEGIN 

      wclr=i
      make_wthetaconv_func_lumw, 2, wclr, lumsolar, lensw, sigcritinv, /retlum

      IF keyword_set(sigc) THEN lensw=sigcritinv^2*1.e9

      histogram_2d, lumsolar, lensw, hist
      range=[0,max(hist.map)/4.]
      rdis, hist.map, xrange=hist.xrange,yrange=hist.yrange,$
        xtit=xt,ytit=yt, /invbw,range=range

;      aplot, !gratio, lumsolar, lensw, xtit=xt,ytit=yt,$
;        psym=4,symsize=0.3
      legend, !colorsp[wclr], box=0, /right, charsize=1.8

      IF (dtype EQ 'X') AND (i LT nclr-1) THEN key=get_kbrd(1)

  ENDFOR 

END 
