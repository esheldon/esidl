PRO testzbuff, pst, fr, fg, fi, imtot2

  IF n_elements(pst) EQ 0 THEN read_spec1d, 494,90,pst

  IF n_elements(fr) EQ 0 THEN BEGIN 

      rgbfchart,pst.raobj, pst.decobj, fr=fr, fg=fg, fi=fi

  ENDIF 

  rgbview, fi, fr, fg, imtot=imtot, /nodisplay

  setupplot,'Z'
  device, set_resolution=[750,750]
  
  implot_setup, fr, xsize, ysize, px, py, xrng, yrng, /center
  pos = [px[0], py[0], px[1], py[1]]

  tv, congrid( reform(imtot[0,*,*]), xsize, ysize), px[0], py[0]
  plot, [0,0], [0,0], xstyle=1, ystyle=1, $
                title=title,xtitle=xtitle,ytitle=ytitle, subtitle=subtitle, $
                xrange=xrng, yrange=yrng, position=pos,$
                /noerase, /device, /nodata

  tmp = tvrd()
  tsz = size(tmp,/dim)
  imtot2 = bytarr(3, tsz[0], tsz[1])
  imtot2[0,*,*] = tmp

  tv, congrid( reform(imtot[1,*,*]), xsize, ysize), px[0], py[0]
  plot, [0,0], [0,0], xstyle=1, ystyle=1, $
                title=title,xtitle=xtitle,ytitle=ytitle, subtitle=subtitle, $
                xrange=xrng, yrange=yrng, position=pos,$
                /noerase, /device, /nodata
  tmp = tvrd()
  imtot2[1,*,*] = tmp

  tv, congrid( reform(imtot[2,*,*]), xsize, ysize), px[0], py[0]
  plot, [0,0], [0,0], xstyle=1, ystyle=1, $
                title=title,xtitle=xtitle,ytitle=ytitle, subtitle=subtitle, $
                xrange=xrng, yrange=yrng, position=pos,$
                /noerase, /device, /nodata
  tmp = tvrd()
  imtot2[2,*,*] = tmp

  setupplot,'X'

  tvimage,imtot2,/keep

  write_jpeg, '~/tmp/test.jpg', imtot2, true=1

END 
