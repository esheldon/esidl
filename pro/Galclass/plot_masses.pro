PRO plot_masses, a, e, s, bin=bin

  nr=n_elements(a.rmax_act)
  IF n_elements(bin) EQ 0 THEN bin=nr-1


  IF (!d.flags AND 1) EQ 0 THEN doX=1 ELSE doX=0
  IF doX THEN BEGIN
      !p.charsize = 1.5
      charsize = .7
      !p.thick = 1
      !x.thick = 1
      !y.thick = 1
      !p.charthick=1
  ENDIF ELSE BEGIN
      charsize=.7
      !p.thick = 5
      !x.thick = 5
      !y.thick = 5
      !p.charthick = 4
  ENDELSE 
  !x.margin=[10,10]
  simpctable
  !p.background=!white
  !p.color=!black

  IF tag_exist(a,'sismass') THEN BEGIN 
      yt='M(r)  (10!U12!N h!U-1!N M!M!Ln!N)!X'
      yt2='M(r)!Dellip!N / M(r)!Dspiral!N'
      xt='Radius (h!U-1!N kpc)'
      xrange=[0,400]
      yrange=[0,10]
      ;;  yrange=[0,15]

      aploterror,1, e.rmax_act[0:bin], e.sismass[0:bin]/1.e12, e.sismasserr[0:bin]/1.e12, psym=3, ytitle=yt,xtitle=xt, $
        xrange=xrange,yrange=yrange
      oploterror, e.rmax_act[0:bin], e.sismass[0:bin]/1.e12, e.sismasserr[0:bin]/1.e12, color=!red, errcolor=!red, line=2
      oploterror, s.rmax_act[0:bin], s.sismass[0:bin]/1.e12, s.sismasserr[0:bin]/1.e12, color=!blue, errcolor=!blue, line=3
      oploterror, a.rmax_act[0:bin], a.sismass[0:bin]/1.e12, a.sismasserr[0:bin]/1.e12, line=0
      ;;  oploterror, brg.rmax_act[0:bin], brg.sismass[0:bin]/1.e12, brg.sismasserr[0:bin]/1.e12, line=4,color=!magenta,errcolor=!magenta

      ;;  mess=['All','Spiral','Elliptical','BRG']
      mess=['All','Spiral','Elliptical']
      ;;  legend, mess, line=[0, 3, 2,4], thick=[!p.thick, !p.thick, !p.thick,!p.thick],colors=[!black,!blue, !red,!magenta]
      legend, mess, line=[0, 3, 2], thick=[!p.thick, !p.thick, !p.thick],colors=[!black,!blue, !red]
      
  ENDIF ELSE BEGIN 


      yt = 'Tangential Shear'
      yrange=prange(e.shear[0:bin],s.shear[0:bin],e.shearerr[0:bin],s.shearerr[0:bin])

      aploterror,1, e.rmax_act[0:bin], e.shear[0:bin], e.shearerr[0:bin], psym=3, ytitle=yt,xtitle=xt, $
        xrange=xrange,yrange=yrange
      oplot,[0,10000],[0,0]
      oploterror, e.rmax_act[0:bin], e.shear[0:bin], e.shearerr[0:bin], color=!red, errcolor=!red, line=2
      oploterror, s.rmax_act[0:bin], s.shear[0:bin], s.shearerr[0:bin], color=!blue, errcolor=!blue, line=3
      oploterror, a.rmax_act[0:bin], a.shear[0:bin], a.shearerr[0:bin], line=0

      legend, mess, line=[0, 3, 2], thick=[!p.thick, !p.thick, !p.thick],colors=[!black,!blue, !red],/right

      key=get_kbrd(1)

      xt = 'Projected Radius (h!U-1!N kpc)'
      yt = '!S!7R!3!R!A-!N(!Ml!3r) - !S!7R!3!R!A-!N(r) (h M!M!Ln!N!3 pc!U-2!N)!X'
      yrange=prange(e.sigma[0:bin],s.sigma[0:bin],e.sigmaerr[0:bin],s.sigmaerr[0:bin])

      aploterror,1, e.rmax_act[0:bin], e.sigma[0:bin], e.sigmaerr[0:bin], psym=3, ytitle=yt,xtitle=xt, $
        xrange=xrange,yrange=yrange
      oplot,[0,10000],[0,0]
      oploterror, e.rmax_act[0:bin], e.sigma[0:bin], e.sigmaerr[0:bin], color=!red, errcolor=!red, line=2
      oploterror, s.rmax_act[0:bin], s.sigma[0:bin], s.sigmaerr[0:bin], color=!blue, errcolor=!blue, line=3
      oploterror, a.rmax_act[0:bin], a.sigma[0:bin], a.sigmaerr[0:bin], line=0

      mess=['All','Spiral','Elliptical']
      legend, mess, line=[0, 3, 2], thick=[!p.thick, !p.thick, !p.thick],colors=[!black,!blue, !red],/right

  ENDELSE 


  !p.charsize = 1
  charsize = .7
  !p.thick = 1
  !x.thick = 1
  !y.thick = 1
  !p.charthick=1

END 
