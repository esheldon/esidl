PRO readlog, file, M_Power, M_Powererr, lind, Map, Maperr, Rmax, rad, raderr, nbins, area, nave

  IF n_params() LT 1 THEN BEGIN
      print,'Syntax: readlog, file, M_Power, M_Powererr, lind, Map, Maperr, Rmax, rad, raderr, nbins, area, nave'
      return
  ENDIF 

  fmt = 'F, F, F, F, F, F, F, F'

  readcol,file,format=fmt, $
          nbins, $
          binsize, $
          area, $
          Map, $
          Maperr, $
          rad, $
          raderr, $
          nave

  Rmax = binsize/2.
  pold = !p.multi
  !p.multi=[0,0,2]

  title='Run 752/756   r'
  ytitle='<Map^2>'
  xtitle='Arcminutes'
  xrange=[.1, 100]

  ploterr, Rmax*60., Map, Maperr, $
           psym=1, $
           title=title, xtitle=xtitle, ytitle=ytitle, /xlog, $
           xrange=xrange, xstyle=1

  key = get_kbrd(1)

  
  yrange=[-.0002, .0008]
  ytitle = 'Radial'
  ploterr, Rmax*60., rad, raderr, $
           psym=1, $
           title=title, xtitle=xtitle, ytitle=ytitle, /xlog, $
           xrange=xrange, xstyle=1, yrange=yrange

  key = get_kbrd(1)

  w=where(Map GT 0., nw)
  print,Rmax[w]

  L = 1.                        
  G = (1+L)*(2+L)^2/(1+2*L)/(3+2*L)

  ;; Factor is for using arcminutes instead of radians
  lind = 2.16e4/(60.*Rmax[w])
print,lind

  M_Power = !pi/G*Map[w]/lind^2
  M_Powererr = !pi/G*Maperr[w]/lind^2

;  M_Powererr[nw-1] = M_Powererr[nw-1]*.2

  !p.multi=pold
  
  ytitle='P(kappa) l'
  xtitle = 'l'
  yrange=[1.e-13,1.d-10]
 
  plot, lind, M_Power, $
           psym=1, $
           title=title, xtitle=xtitle, ytitle=ytitle, /xlog,/ylog, $
           yrange=yrange

  oploterr, lind, M_Power, M_Powererr, $
           psym=1



  !p.multi=pold

return
END 
