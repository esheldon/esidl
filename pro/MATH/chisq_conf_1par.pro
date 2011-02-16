PRO chisq_conf_1par, datax, data, dataerr, modelx, model, param, $
                     chisq_surf, min, conf_low, conf_high, $
                     chisq_diff=chisq_diff, noplotmin=noplotmin, $
                     plot_both=plot_both,wparam=wparam,aspect=aspect, center=center, $
                     _extra=extra

  IF n_params() LT 6 THEN BEGIN
      print,'-Syntax: '
      return
  ENDIF 
  nlevels = 3                   ;68, 95, 99%
  line = [0,2,1]

  IF n_elements(aspect) EQ 0 THEN aspect = 1.0

  default = 1.e10
  conf_low = replicate(default,nlevels)
  conf_high = replicate(default,nlevels)

  w=where(dataerr EQ 0., nw)
  IF nw NE 0 THEN BEGIN
      print,'bad0'
      return
  ENDIF 

  szpar=n_elements(param)

  sz_model = size(model)
  sz_modelx = n_elements(modelx)

  IF sz_model[2] NE sz_modelx THEN BEGIN 
      print,'Bad1'
      return
  ENDIF 

  index = lindgen(szpar)

  chisq_surf = fltarr(szpar)

  errsquared=dataerr^2
  FOR ipar=0L, szpar-1 DO BEGIN
      
          tmpmod = model[ipar, *]
          
          modfunc=interpol(tmpmod, modelx, datax)
          
          ff= (data-modfunc)^2/errsquared
          chisq_surf[ipar] = total(ff)

  ENDFOR 

  minchisq = min(chisq_surf, /nan)
  

  chisq_diff = chisq_surf - minchisq
  w=where(chisq_surf EQ minchisq, nw)
  wparam=w

  min = (param[w])[0]

  ndata = n_elements(data)
  degfree = ndata-1

  wL = where(chisq_surf/degfree LE 8.,nwL)
  IF nwL EQ 0 THEN wL = lindgen(n_elements(param))
;  plot, param[wL], chisq_surf[wL]/degfree,_extra=extra, $
;    yrange=[0,min([8., max(chisq_surf/degfree)])]
;  IF NOT keyword_set(noplotmin) THEN $
;    oplot, [param[w]], [chisq_surf[w]/degfree], psym=7

  aplot, aspect, param[wL], chisq_diff[wL],_extra=extra, center=center
  IF NOT keyword_set(noplotmin) THEN $
    oplot, [param[w]], [chisq_diff[w]], psym=7

    print,'Min Chisq: ',minchisq,'/',ntostr(long(degfree)),' = ',ntostr(minchisq/degfree)

  FOR i=0L, nlevels-1 DO BEGIN

      w = where(chisq_diff LE !siglevels1[i], nw)

      IF nw NE 0 THEN BEGIN 
          conf_low[i] = min( param[w] )
          conf_high[i] = max( param[w] )

          oplot, [ conf_low[i], conf_low[i] ], [-10000., 10000], line=line[i]
          oplot, [ conf_high[i], conf_high[i] ], [-10000., 10000], line=line[i]

      ENDIF 
  ENDFOR 

  return
END 
