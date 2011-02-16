PRO manasseplot

  IF display_type() EQ 'X' THEN BEGIN 
      !p.color=!black
      !p.background=!white
  ENDIF 
  pclr=!skyblue
  pclr=!lightsteelblue

  nn = 100
  zc = arrscl(findgen(nn), 0.0, 5.0)
  
  denom = 4.0+3.0*zc
  omega_m = (2.0+zc)/denom
  omega_lambda = 2.0*(1.0+zc)/denom
  
  omega_m_old = 2.0/denom
  delta_omega_m = zc/denom
  
  baryon = 0.0367
  
  omega_dm = omega_m - baryon
  omega_dm_old = omega_m_old - baryon

  ;; Calculate the allowed zc values given
  ;; constraints on omega_dm_old


  min_om = 0.20
  cmb_om = 0.227
  
  ;; old max_om
  max_om = 0.35
  ;; new max_om
  max_om = 0.3

  omdiff = abs(omega_m_old - min_om)
  zcmin_ind = where( min(omdiff) EQ omdiff )
  omdiff = abs(omega_m_old - max_om)
  zcmax_ind = where( min(omdiff) EQ omdiff )
  omdiff = abs(omega_m_old - cmb_om)

  zcallow = [ zc[zcmax_ind], zc[zcmin_ind] ]

  zc_cmb = 1.2
print,min_om,max_om
print,zcallow,zc_cmb

  erase & multiplot, [1,2]
  
  lines = [0,2,3,4]
  yrange=[0,1]
  ytitle=!tsym.omega_cap+'!DX!N'
  xtitle='z!Dc!N'
  
  colors = [!black, !red, !magenta, !blue]

  plot, zc, omega_lambda, yrange=yrange, ytitle=ytitle, line=lines[0],$
    color=colors[0]
  oplot, zc, omega_m, line=lines[1],$
    color=colors[1]
  oplot, zc, omega_m_old, line=lines[2],$
    color=colors[2]
  oplot, zc, delta_omega_m, line=lines[3],$
    color=colors[3]
  
  message=[!tsym.omega_cap+'!D'+!tsym.lambda_cap+',0!N',$
           !tsym.omega_cap+'!Dm,0!N',$
           !tsym.omega_cap+'!Dm,old!N',$
           !tsym.delta_cap+!tsym.omega_cap+'!Dm!N']
  legend,message,line=lines,/right,thick=replicate(!p.thick,4),spacing=1.6,$
    box=0, colors=colors
  legend,'(a)',/left,box=0
;  oplot, [zc_cmb, zc_cmb], [0,1],line=1
  multiplot
  
  smax=1
  
  plot, zc, omega_lambda, yrange=yrange, ytitle=ytitle, xtitle=xtitle,$
    line=lines[0],$
    color=colors[0]
;  polyfill, [0, zcallow[0], zcallow[0], 0], $
;    [0, 0, smax, smax], color=pclr
;  polyfill, [zcallow[1], 5, 5, zcallow[1]], $
;    [0, 0, smax, smax], color=pclr

;polyfill, [0, zcallow[0], zcallow[0], 0], $
;          [0, 0, smax, smax], /line_fill, orientation=45
;polyfill, [zcallow[1], 5, 5, zcallow[1]], $
;          [0, 0, smax, smax], color=pclr, /line_fill, orientation=45
;  plot, zc, omega_lambda, yrange=yrange, ytitle=ytitle, xtitle=xtitle,line=lines[0]
  
  
  oplot, zc, omega_dm, line=lines[1],$
    color=colors[1]
  oplot, zc, omega_dm_old, line=lines[2],$
    color=colors[2]
  oplot, zc, delta_omega_m, line=lines[3],$
    color=colors[3]

;  oplot, [zc_cmb, zc_cmb], [0,1],line=1

  message=[!tsym.omega_cap+'!D'+!tsym.lambda_cap+',0!N',$
           !tsym.omega_cap+'!Ddm,0!N',$
           !tsym.omega_cap+'!Ddm,old!N',$
           !tsym.delta_cap+!tsym.omega_cap+'!Dm!N']
  legend,message,line=lines,/right,thick=replicate(!p.thick,4),spacing=1.6,$
    box=0, colors=colors
  legend,'(b)',/left,box=0

  multiplot,/reset

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now repeat the last one with only upper bound
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  aplot, !gratio, zc, omega_lambda, yrange=yrange, ytitle=ytitle, xtitle=xtitle,line=lines[0]
  ;polyfill, [0, zcallow[0], zcallow[0], 0], $
  ;  [0, 0, smax, smax], color=pclr
  polyfill, [zcallow[1], 5, 5, zcallow[1]], $
    [0, 0, smax, smax], color=pclr
  axis, xaxis=0
  axis, xaxis=1,xtickname=replicate(' ',6)
  axis, yaxis=0
  axis, yaxis=1,ytickname=replicate(' ',6)
  oplot, zc, omega_lambda, line=lines[0]

  oplot, zc, omega_dm, line=lines[1]
  oplot, zc, omega_dm_old, line=lines[2]
  oplot, zc, delta_omega_m, line=lines[3]

  ;oplot, [zc_cmb, zc_cmb], [0,1],line=1

  message=[!tsym.omega_cap+'!D'+!tsym.lambda_cap+',0!N',$
           !tsym.omega_cap+'!Ddm,0!N',$
           !tsym.omega_cap+'!Ddm,old!N',$
           !tsym.delta_cap+!tsym.omega_cap+'!Dm!N']
  legend,message,line=lines,/right,thick=replicate(!p.thick,4),spacing=1.6,$
    box=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now repeat the last one with both bounds
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  aplot, !gratio, zc, omega_lambda, yrange=yrange, ytitle=ytitle, xtitle=xtitle,line=lines[0]
    polyfill, [0, zcallow[0], zcallow[0], 0], $
    [0, 0, smax, smax], color=pclr
  polyfill, [zcallow[1], 5, 5, zcallow[1]], $
    [0, 0, smax, smax], color=pclr
  axis, xaxis=0
  axis, xaxis=1,xtickname=replicate(' ',6)
  axis, yaxis=0
  axis, yaxis=1,ytickname=replicate(' ',6)
  oplot, zc, omega_lambda, line=lines[0]

  oplot, zc, omega_dm, line=lines[1]
  oplot, zc, omega_dm_old, line=lines[2]
  oplot, zc, delta_omega_m, line=lines[3]

  ;oplot, [zc_cmb, zc_cmb], [0,1],line=1

  message=[!tsym.omega_cap+'!D'+!tsym.lambda_cap+',0!N',$
           !tsym.omega_cap+'!Ddm,0!N',$
           !tsym.omega_cap+'!Ddm,old!N',$
           !tsym.delta_cap+!tsym.omega_cap+'!Dm!N']
  legend,message,line=lines,/right,thick=replicate(!p.thick,4),spacing=1.6,$
    box=0
  ;legend,'(b)',/left,box=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now repeat the last one with both bounds and cmb
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  aplot, !gratio, zc, omega_lambda, yrange=yrange, ytitle=ytitle, xtitle=xtitle,line=lines[0]
    polyfill, [0, zcallow[0], zcallow[0], 0], $
    [0, 0, smax, smax], color=pclr
  polyfill, [zcallow[1], 5, 5, zcallow[1]], $
    [0, 0, smax, smax], color=pclr
  axis, xaxis=0
  axis, xaxis=1,xtickname=replicate(' ',6)
  axis, yaxis=0
  axis, yaxis=1,ytickname=replicate(' ',6)
  oplot, zc, omega_lambda, line=lines[0]

  oplot, zc, omega_dm, line=lines[1]
  oplot, zc, omega_dm_old, line=lines[2]
  oplot, zc, delta_omega_m, line=lines[3]

  oplot, [zc_cmb, zc_cmb], [0,1],line=1

  message=[!tsym.omega_cap+'!D'+!tsym.lambda_cap+',0!N',$
           !tsym.omega_cap+'!Ddm,0!N',$
           !tsym.omega_cap+'!Ddm,old!N',$
           !tsym.delta_cap+!tsym.omega_cap+'!Dm!N']
  legend,message,line=lines,/right,thick=replicate(!p.thick,4),spacing=1.6,$
    box=0
  ;legend,'(b)',/left,box=0

  IF display_type() EQ 'X' THEN BEGIN 
      !p.color=!white
      !p.background=!black
  ENDIF 

return
END 
