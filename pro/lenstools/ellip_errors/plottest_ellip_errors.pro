PRO plottest_e1e2err_dist, str, title=title

  n=n_elements(str.s2n)

  pcharold = !p.charsize

  IF !d.name EQ 'X' THEN !p.charsize = 2.0 $
  ELSE !p.charsize=1.5

  IF n NE 1 THEN !p.multi=[0,0,n]

  xtitle = 'e!D1!Ne!D2!Nerr From Measured Parameters'
  FOR i=0L, n-1 DO BEGIN 

      s2n = str[i].s2n
      IF s2n LE 40.0 THEN binsize = 0.0005 ELSE binsize = 0.00001

      wuse = where(str[i].me1 EQ str[i].me1)
      wuse2 = where(str[i].me1[wuse] GT -1.0 AND $
                    str[i].me1[wuse] LT  1.0,nuse)
      wuse = wuse[wuse2]

      xmin = str[i].f_e1e2err_mean - 3.*str[i].f_e1e2err_err
      xmax = str[i].f_e1e2err_mean + 3.*str[i].f_e1e2err_err

      min = str[i].f_e1e2err_mean - 5.*str[i].f_e1e2err_err
      max = str[i].f_e1e2err_mean + 5.*str[i].f_e1e2err_err

      plothist, str[i].f_e1e2err[wuse], bin=binsize,$
                title=title+' S/N: '+ntostr(str[i].s2n,5),$
                xtitle=xtitle, min=min, max=max,xrange=[xmin,xmax]
      
      val = str[i].e1e2err
      oplot, [val, val], $
             [0, 1.e6], color=!red

;      val = str[i].f_e1e2err_mean
;      oplot, [val, val], $
;             [0, 1.e6], color=!yellow,thick=2

;      val = str[i].f_e1e2err_mean+str[i].f_e1e2err_err
;      oplot, [val, val], $
;             [0, 1.e6], color=!yellow,thick=2,line=2

;      val = str[i].f_e1e2err_mean-str[i].f_e1e2err_err
;      oplot, [val, val], $
;             [0, 1.e6], color=!yellow,thick=2,line=2


      val = str[i].m_e1e2err
      oplot, [val, val], $
             [0, 1.e6], color=!blue

      val = str[i].m_e1e2err+str[i].m_e1e2err_err
      oplot, [val, val], $
             [0, 1.e6], color=!blue,line=2

      val = str[i].m_e1e2err-str[i].m_e1e2err_err
      oplot, [val, val], $
             [0, 1.e6], color=!blue,line=2

;      legend,['No Noise','fmean','bootstrap'],$
;             line=[0,0,0],$
;             colors=[!red,!yellow,!blue],$
;             charsize=1.0

      legend,['Input Pars','bootstrap'],$
             line=[0,0],$
             colors=[!red,!blue],$
             charsize=1.0,thick=[!p.thick,!p.thick]

  ENDFOR 

  IF n NE 1 THEN !p.multi=0
  !p.charsize=pcharold

END 

PRO plottest_ellip_errors, str

  posangles = str[rem_dup(str.posangle)].posangle
  ellip = str[rem_dup(str.e)].e
  s2n = str[rem_dup(str.s2n)].s2n

  s2nxr = [0.9*min(s2n), 1.1*max(s2n)]

  npos = n_elements(posangles)
  nell = n_elements(ellip)
  ns2n = n_elements(s2n)

  FOR ipos = 0L, npos-1 DO BEGIN 
      tpos = posangles[ipos]
      FOR iell = 0L, nell-1 DO BEGIN 
          tell = ellip[iell]
          
          w=where(str.posangle EQ tpos AND $
                  str.e EQ tell, nw)
          
          IF nw NE 0 THEN BEGIN
              
              ;; only plot versus S/N if there is more than
              ;; one value for s2n!!
              xtitle='S/N'
              title='posang: '+ntostr(tpos,5)+' e: '+ntostr(tell,5)
              
              IF nw GT 1 THEN BEGIN 
                  
                  f_e1e1err_mean = str[w].f_e1e1err_mean
                  f_e1e1err_err  = str[w].f_e1e1err_err
                  f_e1e2err_mean = str[w].f_e1e2err_mean
                  f_e1e2err_err  = str[w].f_e1e2err_err
                  f_e2e2err_mean = str[w].f_e2e2err_mean
                  f_e2e2err_err  = str[w].f_e2e2err_err
                  
                  
                  xtitle='S/N'
                  title='posang: '+ntostr(tpos)+' e: '+ntostr(tell)
;                  !p.multi=[0,0,3]
                  !p.multi=[0,0,2]

                  pcharold=!p.charsize
                  IF !d.name EQ 'X' THEN !p.charsize=2.0 $
                  ELSE !p.charsize = 1.5
                  ytitle='meas_e1e1err/formal_e1e1err'
                  
                  ;; err measured from sample
                  m_e1e1err = str[w].m_e1e1err
                  ;;m_e1e1err_err = str[w].m_e1e1err/sqrt(2.*str[w].nuse)
                  m_e1e1err_err = str[w].m_e1e1err_err
                  
                  ratio = f_e1e1err_mean/m_e1e1err
                  raterr = ratio*sqrt( (f_e1e1err_err/f_e1e1err_mean)^2 + $
                                       (m_e1e1err_err/m_e1e1err)^2 )
                  
                  yrange = prange(ratio,raterr,/slack)
                  yrange[0] = 0
                  ploterror, str[w].s2n,$
                             ratio, raterr,$
                             psym=1,xtitle=xtitle,ytitle=ytitle,title=title,$
                             xrange=s2nxr, yrange=yrange
                  oplot,[-1000,1000],[1,1]
                  
                  ;; err measured from sample
                  ;; Don't plot: mean is often zero
;                  ytitle='meas_e1e2err/formal_e1e2err'
;                  m_e1e2err = str[w].m_e1e2err
;                  ;;m_e1e2err_err = str[w].m_e1e2err/sqrt(2.*str[w].nuse)
;                  m_e1e2err_err = str[w].m_e1e2err_err
                  
;                  ratio = f_e1e1err_mean/m_e1e1err
;                  raterr = ratio*sqrt( (f_e1e2err_err/f_e1e2err_mean)^2 + $
                                       (m_e1e2err_err/m_e1e2err)^2 )
;                  yrange = prange(ratio,raterr,/slack)
;                  yrange[0] = 0
;                  ploterror, str[w].s2n,$
;                             ratio, raterr,$
;                             psym=1,xtitle=xtitle,ytitle=ytitle,title=title,$
;                             xrange=s2nxr,yrange=yrange
;                  oplot,[-1000,1000],[1,1]
                  
                  ;; err measured from sample
                  ytitle='meas_e2e2err/formal_e2e2err'
                  m_e2e2err = str[w].m_e2e2err
                  ;;m_e2e2err_err = str[w].m_e2e2err/sqrt(2.*str[w].nuse)
                  m_e2e2err_err = str[w].m_e2e2err_err
                  
                  ratio = f_e1e1err_mean/m_e1e1err
                  raterr = ratio*sqrt( (f_e2e2err_err/f_e2e2err_mean)^2 + $
                                       (m_e2e2err_err/m_e2e2err)^2 )
                  yrange = prange(ratio,raterr,/slack)
                  yrange[0] = 0
                  ploterror, str[w].s2n,$
                             ratio, raterr,$
                             psym=1,xtitle=xtitle,ytitle=ytitle,title=title,$
                             xrange=s2nxr,yrange=yrange
                  oplot,[-1000,1000],[1,1]
                  
                  !p.multi=0
                  !p.charsize=pcharold
                  
                  IF !d.name EQ 'X' THEN BEGIN 
                      key=get_kbrd(1)
                      IF key EQ 'q' THEN return
                  ENDIF 
                  
              ENDIF ;; If nw > 1
              
              ;; now plot e1e2err
              plottest_e1e2err_dist, str[w],title=title
              
              IF !d.name EQ 'X' THEN BEGIN 
                  key=get_kbrd(1)
                  IF key EQ 'q' THEN return
              ENDIF 
                            
          ENDIF 

      ENDFOR 
  ENDFOR 

END 
