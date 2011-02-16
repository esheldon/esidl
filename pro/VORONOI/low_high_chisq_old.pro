PRO low_high_chisq, plot_both=plot_both, plot_sig=plot_sig, $
                    psfile=psfile, phil=phil,nops=nops

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; This plots the chisq surface for all bands for
  ;; High, Low, and All density regions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  outdir='/sdss4/data1/esheldon/TMP/conv/'

  IF n_elements(psfile) EQ 0 THEN BEGIN
      outfile = outdir + 'chisq_all_high_low.ps'
  ENDIF ELSE BEGIN
      outfile = outdir + psfile
  ENDELSE 

  IF NOT keyword_set(nops) THEN nops = 0
  IF NOT nops THEN begplot, name=outfile, /encapsulated

  colors=['u','g','r','i','z']
  pos = [75., 105.]

  clr = [1,2,3]
  nclr=n_elements(clr)
  type = [0,1,2]


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; All regions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  erase & multiplot, [1, nclr]

  xt='Cutoff Radius (arcsec)'
;  yt='Velocity Dispersion (km/s)'
  yt='!7r!3!DV!N (km/s)'
  tt='All Regions'
  FOR i=0L, nclr-1 DO BEGIN

      rdwtheta_conv, clr[i], type[2], datax, data, dataerr, modelx, model, $
                     sigma, cutoff, phil=phil

      xtitle=''
      ytitle=''
      title=''
      IF i EQ nclr-1 THEN BEGIN
          xtitle=xt
          ytitle=yt
          title=''
      ENDIF ELSE IF i EQ 0 THEN BEGIN
          title=tt
      ENDIF 
        
      yrange=[min(sigma),max(sigma)]
      xrange=[min(cutoff),max(cutoff)]

      chisq_conf, datax, data, dataerr, modelx, model, cutoff, sigma, $
        chisq_surf, min1, min2, low1, high1, low2, high2, $
        plot_both=plot_both, plot_sig=plot_sig, $
        xtitle=xtitle, title=title,$
        yrange=yrange,xrange=xrange, $
        ystyle=1+4,xstyle=1

      nyticks=7
      ytickn=[' ','120','140','160','180','200',' ']
      axis, yaxis=0,yticks=nyticks-1,ytickn=ytickn, ytitle=ytitle
      axis, yaxis=1, yticks=nyticks-1, ytickn=[replicate(' ',nyticks)]

      xyouts, pos[0], pos[1], 'Cutoff '+'['$
        +ntostr(low1[1])+', '+ntostr(high1[1])+']', charsize=0.8
      xyouts, pos[0], pos[1]+7., 'Vel. '+'['$
        +ntostr(low2[1])+', '+ntostr(high2[1])+']', charsize=0.8
      xyouts, pos[0], pos[1]+14., '95% conf. levels.',charsize=0.8

      print,'cutoff'
      forprint,low1,high1
      print,'sigma'
      forprint,low2,high2

      IF i NE nclr-1 THEN multiplot
  ENDFOR 
  multiplot, /reset

  IF nops THEN BEGIN
      key=get_kbrd(1)
      IF key EQ 'q' THEN return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; High density regions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  erase & multiplot, [1, nclr]

  tt='High Density Regions'
  FOR i=0L, nclr-1 DO BEGIN

      rdwtheta_conv, clr[i], type[0], datax, data, dataerr, modelx, model, $
                     sigma, cutoff, phil=phil

      xtitle=''
      ytitle=''
      title=''
      IF i EQ nclr-1 THEN BEGIN
          xtitle=xt
          ytitle=yt
          title=''
      ENDIF ELSE IF i EQ 0 THEN BEGIN
          title=tt
      ENDIF 

      yrange=[min(sigma),max(sigma)]
      xrange=[min(cutoff),max(cutoff)]

      chisq_conf, datax, data, dataerr, modelx, model, cutoff, sigma, $
        chisq_surf, min1, min2, low1, high1, low2, high2, $
        plot_both=plot_both, plot_sig=plot_sig, $
        xtitle=xtitle, title=title,xstyle=1,$
        yrange=yrange,xrange=xrange, $
        ystyle=1+4,xstyle=1

      nyticks=7
      ytickn=[' ','120','140','160','180','200',' ']
      axis, yaxis=0,yticks=nyticks-1,ytickn=ytickn, ytitle=ytitle
      axis, yaxis=1, yticks=nyticks-1, ytickn=[replicate(' ',nyticks)]

      xyouts, pos[0], pos[1], 'Cutoff '+'['$
        +ntostr(low1[1])+', '+ntostr(high1[1])+']', charsize=0.8
      xyouts, pos[0], pos[1]+7., 'Vel. '+'['$
        +ntostr(low2[1])+', '+ntostr(high2[1])+']', charsize=0.8
      xyouts, pos[0], pos[1]+14., '95% conf. levels.',charsize=0.8

      IF i NE nclr-1 THEN multiplot
  ENDFOR 
  multiplot, /reset

  IF nops THEN BEGIN
      key=get_kbrd(1)
      IF key EQ 'q' THEN return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Low density regions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  erase & multiplot, [1, nclr]
  
  tt='Low Density Regions'
  FOR i=0L, nclr-1 DO BEGIN

      rdwtheta_conv, clr[i], type[1], datax, data, dataerr, modelx, model, $
                     sigma, cutoff, phil=phil

      xtitle=''
      ytitle=''
      title=''
      IF i EQ nclr-1 THEN BEGIN
          xtitle=xt
          ytitle=yt
          title=''
      ENDIF ELSE IF i EQ 0 THEN BEGIN
          title=tt
      ENDIF 
          
      yrange=[min(sigma),max(sigma)]
      xrange=[min(cutoff),max(cutoff)]

      chisq_conf, datax, data, dataerr, modelx, model, cutoff, sigma, $
        chisq_surf, min1, min2, low1, high1, low2, high2, $
        plot_both=plot_both, plot_sig=plot_sig, $
        xtitle=xtitle, title=title,xstyle=1,$
        yrange=yrange,xrange=xrange, $
        ystyle=1+4,xstyle=1

      nyticks=7
      ytickn=[' ','120','140','160','180','200',' ']
      axis, yaxis=0,yticks=nyticks-1,ytickn=ytickn, ytitle=ytitle
      axis, yaxis=1, yticks=nyticks-1, ytickn=[replicate(' ',nyticks)]

      xyouts, pos[0], pos[1], 'Cutoff '+'['$
        +ntostr(low1[1])+', '+ntostr(high1[1])+']', charsize=0.8
      xyouts, pos[0], pos[1]+7., 'Vel. '+'['$
        +ntostr(low2[1])+', '+ntostr(high2[1])+']', charsize=0.8
      xyouts, pos[0], pos[1]+14., '95% conf. levels.',charsize=0.8

      IF i NE nclr-1 THEN multiplot
  ENDFOR 
  IF NOT keyword_set(nops) THEN endplot,/noprint

  multiplot, /reset

  !p.charsize=0

return
END 
