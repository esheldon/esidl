PRO plot_primary_bounds, stripes, overplot=overplot, _extra=_extra

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: plot_primary_bounds, stripes, overplot=overplot, _extra=_extra'
      return
  ENDIF 
  nst = n_elements(stripes)

  primary_bound_multi, stripes, bound_array, overall_bound
  nbound = n_elements(bound_array)

  IF NOT keyword_set(overplot) THEN BEGIN 
      
      etamin = overall_bound.etamin
      etamax = overall_bound.etamax
      lammin = overall_bound.lammin
      lammax = overall_bound.lammax

      xrange = fltarr(2)
      yrange = fltarr(2)

      fac1 = 0.9
      fac2 = 1.1

      IF lammin LT 0 THEN xrange[0] = fac2*lammin ELSE xrange[0] = fac1*lammin
      IF lammax LT 0 THEN xrange[1] = fac1*lammax ELSE xrange[1] = fac2*lammax
      IF etamin LT 0 THEN yrange[0] = fac2*etamin ELSE yrange[0] = fac1*etamin
      IF etamax LT 0 THEN yrange[1] = fac1*etamax ELSE yrange[1] = fac2*etamax

      print,xrange,yrange
      print,etamin,etamin*fac2
      xtitle = !csym.lambda+'!Dc!N'
      ytitle = !csym.eta+'!Dc!N'
      plot, $
        [0], /nodata, $
        xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, $
        xtitle=xtitle, ytitle=ytitle

  ENDIF 
      
  FOR i=0L, nst-1 DO BEGIN 
      plot_box, $
        bound_array[i].lammin, $
        bound_array[i].lammax, $
        bound_array[i].etamin, $
        bound_array[i].etamax,$
        _extra=_extra
  ENDFOR 

;  primary_bound_multi, stripes, bound, overall_bound
;  plot_box, $
;    overall_bound.lammin, overall_bound.lammax, $
;    overall_bound.etamin, overall_bound.etamax,$
;    _extra=_extra

END 
