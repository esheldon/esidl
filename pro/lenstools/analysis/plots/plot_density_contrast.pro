PRO plot_density_contrast, struct, $
                           logbin=logbin, $
                           miny=miny, maxy=maxy, $
                           aspect=aspect, $
                           overplot=overplot, $
                           center=center, $
                           wuse=wuse, $
                           yrange=yrange, xrange=xrange, $
                           noxtitle=noxtitle, noytitle=noytitle, $
                           mpc=mpc, $
                           rebin=rebin, $
                           axis_shear=axis_shear, $
                           sigma_comoving=sigma_comoving, _extra=_extra

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: plot_density_contrast, struct, /logbin, miny=, maxy=, aspect=, /overplot, /center, wuse=, yrange=, xrange=, /noxtitle, /noytitle, /axis_shear, /sigma_comoving, _extra='
      return
  ENDIF 

  IF keyword_set(rebin) THEN BEGIN 
      meanr = struct.meanr_rebin
      sigma = struct.sigma_rebin
      sigmaerr = struct.sigmaerr_rebin
  ENDIF ELSE BEGIN 
      meanr = struct.meanr
      sigma = struct.sigma
      sigmaerr = struct.sigmaerr
  ENDELSE 

  esheldon_setup
  myusersym, 'fill_circle'
  psym=8
  minfac=0.5
  maxfac=1.5
  IF keyword_set(noxtitle) THEN xtitle='' ELSE BEGIN
      IF keyword_set(mpc) THEN BEGIN 
          ;xtitle=!mpcxtitle2
          xtitle=textoidl('R [h^{-1} Mpc]')
          radfac = 1./1000.
      ENDIF ELSE BEGIN 
          ;xtitle=!kpcxtitle2
          xtitle=textoidl('R [h^{-1} kpc]')
          radfac = 1.0
      ENDELSE 
  ENDELSE 
  IF keyword_set(noytitle) THEN BEGIN
      sigytitle=''
  ENDIF ELSE BEGIN 
        sigytitle = textoidl('\Delta\Sigma [h M_{\odot} pc^{-2}]')
  ENDELSE 
  IF n_elements(miny) EQ 0 THEN miny=0.1 ;Msolar/pc^2
  IF n_elements(maxy) EQ 0 THEN maxy=10000.

  IF n_elements(wuse) EQ 0 THEN BEGIN
      IF keyword_set(logbin) THEN w=where(meanr GT 0.0) $
      ELSE w=lindgen(n_elements(meanr))
  ENDIF ELSE w=wuse


  IF keyword_set(logbin) THEN BEGIN 
      IF n_elements(xrange) EQ 0 THEN $
        xrange=[minfac*min(meanr[w]*radfac), maxfac*max(meanr[w]*radfac)]

      xticklen=0.04 & yticklen=0.04
      xlog=1 & ylog=1

      w2=where(sigma[w] GT 0.0,nw)
      IF n_elements(yrange) EQ 0 THEN BEGIN 
          yrange=prange(sigma[w[w2]], sigmaerr[w[w2]],slack=2.0)
          yrange[0] = yrange[0] > miny
          yrange[1] = yrange[1] < maxy
      ENDIF 


      xtickf = 'loglabels'
      ytickf = 'loglabels'

  ENDIF ELSE BEGIN 

      IF n_elements(yrange) EQ 0 THEN yrange=prange(sigma[w], sigmaerr[w],/slack)
      yrange[0] = yrange[0] < 0
      IF n_elements(xrange) EQ 0 THEN xrange=[0, maxfac*max(meanr[w]*radfac)]

  ENDELSE 
  
  xstyle = 1
  ystyle = 1
  IF keyword_set(axis_shear) THEN ystyle = ystyle + 8

  pplot, meanr[w]*radfac, sigma[w], yerr=sigmaerr[w],$
      overplot=overplot, $
      aspect=aspect, $
      xtitle=xtitle, ytitle=sigytitle, $
      yrange=yrange, xrange=xrange,$
      xlog=xlog,ylog=ylog,psym=psym, $
      ystyle=ystyle, xstyle=xstyle, center=center,$
      xticklen=xticklen, yticklen=yticklen, $
      xtickf=xtickf, ytickf=ytickf, _extra=_extra

;  IF n_elements(aspect) EQ 0 THEN BEGIN 
;      ploterror, meanr[w]*radfac, sigma[w], sigmaerr[w],$
;                 xtitle=xtitle, ytitle=sigytitle, yrange=yrange, xrange=xrange,$
;                 xlog=xlog,ylog=ylog,psym=psym, $
;                 ystyle=ystyle, xstyle=xstyle, center=center,$
;                 xticklen=xticklen, yticklen=yticklen, xtickf=xtickf, ytickf=ytickf, _extra=_extra
;  ENDIF ELSE BEGIN 
;      aploterror, aspect, meanr[w]*radfac, sigma[w], sigmaerr[w],$
;                  xtitle=xtitle, ytitle=sigytitle, $
;                  yrange=yrange, xrange=xrange,$
;                  xlog=xlog,ylog=ylog,psym=psym, $
;                  ystyle=ystyle, xstyle=xstyle, center=center,$
;                  xticklen=xticklen, yticklen=yticklen, xtickf=xtickf, ytickf=ytickf, _extra=_extra
;  ENDELSE
 

  IF keyword_set(axis_shear) THEN BEGIN 

      shear = struct.sigma*struct.meanscritinv
      xr = xrange
      yr = yrange*struct.meanscritinv

      ;; if sigma was converted to comoving, then we must convert the
      ;; shear back
      IF keyword_set(sigma_comoving) THEN yr = yr*(1+struct.zmean)^2

      position = [!x.window[0],!y.window[0],!x.window[1],!y.window[1]]

;      plot, meanr[w]*radfac, shear, psym=psym, $
;            position=position, xstyle=xstyle+4, ystyle=1+2+4, /noerase, $
;            xlog=xlog, ylog=ylog, xrange=xr, yrange=yr

      plot, [0], [0], /nodata, $
            position=position, $
            xrange=xr, yrange=yr, xstyle=xstyle+4, ystyle=ystyle+4, /noerase, $
            xlog=xlog, ylog=ylog

      axis, yaxis=1, yticklen=yticklen, ytickf=ytickf, ytitle=!shytitle, $
            ystyle=ystyle

  ENDIF 

END 
