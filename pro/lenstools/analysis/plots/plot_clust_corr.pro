PRO plot_clust_corr, struct, logbin=logbin, wuse=wuse, aspect=aspect, $
                     miny=miny, maxy=maxy, $
                     xrange=xrange, yrange=yrange, $
                     noxtitle=noxtitle, noytitle=noytitle, _extra=_extra

  IF n_elements(struct) EQ 0 THEN BEGIN 
      print,'-Syntax: plot_clust_corr, struct'
      return
  ENDIF 

  nrad = n_elements(struct.meanr)
  radfac = 1.0/1000.0
  IF n_elements(wuse) EQ 0 THEN BEGIN 
      wuse = lindgen(nrad)
  ENDIF 

  esheldon_setup
  myusersym, 'fill_circle'
  psym=8
  minfac=0.5
  maxfac=1.5

  IF NOT keyword_set(noytitle) THEN BEGIN
      ytitle = 'C(R)-1'
  ENDIF  
  IF NOT keyword_set(noxtitle) THEN BEGIN 
      xtitle = !mpcxtitle2
  ENDIF 

  IF n_elements(miny) EQ 0 THEN miny = 1.e-3

  meanr = struct.meanr[wuse]*radfac
  corr = struct.corr[wuse]
  correrr = struct.corr_err[wuse]

  IF keyword_set(logbin) THEN BEGIN 

      IF n_elements(xrange) EQ 0 THEN $
        xrange=[minfac*min(meanr), maxfac*max(meanr)]

      xticklen=0.04 & yticklen=0.04
      xlog=1 & ylog=1

      w=where(corr GT 0.0,nw)
      IF n_elements(yrange) EQ 0 THEN BEGIN 
          yrange=prange(corr[w]-1, correrr[w],slack=0.5)
          yrange[0] = yrange[0] > miny
          IF n_elements(maxy) NE 0 THEN yrange[1] = yrange[1] > maxy
      ENDIF 

      xtickf = 'loglabels'
      ytickf = 'loglabels'

  ENDIF ELSE BEGIN 

      IF n_elements(yrange) EQ 0 THEN $
        yrange=prange(corr-1, correrr,/slack)
      yrange[0] = yrange[0] < 0
      IF n_elements(xrange) EQ 0 THEN $
        xrange=[0, maxfac*max(meanr)]

  ENDELSE 



  pplot, $
    meanr, corr-1.0, yerr=correrr, $
    xtitle=xtitle, ytitle=ytitle, $
    xrange=xrange, yrange=yrange, $
    xstyle=3, ystyle=3, $
    psym=8, xlog=xlog, ylog=ylog, $
    aspect=aspect, $
    xtickf=xtickf, ytickf=ytickf, $
    _extra=_extra



END 
