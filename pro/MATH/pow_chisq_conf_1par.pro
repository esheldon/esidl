PRO pow_chisq_conf_1par, datax, data, dataerr, normvals, power, $
                         chisq_surf, nmin, normlow, normhigh, $
                         normallow=normallow, $
                         chisq_per=chisq_per, $
                         chisq_diff=chisq_diff, noplotmin=noplotmin, $
                         aspect=aspect, center=center, $
                         _extra=extra, xtick_get=xtick_get, ytick_get=ytick_get,$
                         yfit=yfit

  ;; right now only changes norm

  IF n_params() LT 5 THEN BEGIN
      print,'pow_chisq_conf_1par, datax, data, dataerr, normvals, power, '
      print,'   chisq_surf, nmin, normlow, normhigh, '
      print,'   normallow=normallow, '
      print,'   chisq_per=chisq_per, '
      print,'   chisq_diff=chisq_diff, noplotmin=noplotmin, '
      print,'   aspect=aspect, center=center, '
      print,'   _extra=extra, xtick_get=xtick_get, ytick_get=ytick_get,'
      print,'   yfit=yfit'
      return
  ENDIF 

  nlevels = 3                   ;68, 95, 99.78%
  line = [0,2,1]

  default = 1.e10

  normlow = replicate(default,nlevels)
  normhigh = replicate(default,nlevels)

  w=where(dataerr EQ 0., nw)
  IF nw NE 0 THEN BEGIN
      print,'Some zero error vals'
      return
  ENDIF 

  nnorm = n_elements(normvals)

  chisq_surf = fltarr( nnorm )
  modfunc = data

  index = lindgen(nnorm)

  errsquared=dataerr^2

  FOR in=0L, nnorm-1 DO BEGIN
      
      norm = normvals[in]
      modfunc[*] = norm*(datax)^power
          
      ff= (data-modfunc)^2/errsquared
      chisq_surf[in] = total(ff)
      
  ENDFOR 

  minchisq = min(chisq_surf, /nan)
  
  chisq_diff = chisq_surf - minchisq

  w=where(chisq_surf EQ minchisq, nw)
  wparam=w
  nmin=normvals[w[0]]
  ndata = n_elements(data)
  degfree = ndata-1

  wL = where(chisq_surf/degfree LE 8.,nwL)
  IF nwL EQ 0 THEN wL = lindgen(n_elements(param))

  IF n_elements(aspect) EQ 0 THEN aspect = 1.0

  aplot, aspect, normvals[wL], chisq_diff[wL], _extra=extra, center=center
  IF NOT keyword_set(noplotmin) THEN $
    oplot, [normvals[w]], [chisq_diff[w]], psym=7

  chisq_per = minchisq/degfree
  print,'Min Chisq: ',minchisq,'/',ntostr(long(degfree)),' = ',ntostr(chisq_per)


  FOR i=0L, nlevels-1 DO BEGIN

      w = where(chisq_diff LE !siglevels1[i], nw)

      IF nw NE 0 THEN BEGIN 
          normlow[i] = min( normvals[w] )
          normhigh[i] = max( normvals[w] )

          oplot, [ normlow[i], normlow[i] ], [-10000., 10000], line=line[i]
          oplot, [ normhigh[i], normhigh[i] ], [-10000., 10000], line=line[i]

      ENDIF 
  ENDFOR 

  yfit = nmin*datax^power

  return
END 
