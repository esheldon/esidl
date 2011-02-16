PRO fit_const, data, covariance, crange, nc, $
               chisq, cmin, clow, chigh, wuse=wuse, $
               cerrhigh=cerrhigh, cerrlow=cerrlow, $
               verbose=verbose, $
               plot=plot, _extra=_extra

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: fit_const, data, errors/covariance, range, num, '
      print,'              chisq, cmin, clow, chigh, '
      print,'              wuse=, cerrhigh=, cerrlow=, /verbose, /plot'
      return
  ENDIF 

  ndata = n_elements(data)

  cvals = arrscl( findgen(nc), crange[0], crange[1] )
  chisq = fltarr(nc)

  ndim = size(covariance, /n_dim)

  IF ndim EQ 2 THEN BEGIN 
      cinv = invert(covariance)
      FOR i=0L, nc-1 DO BEGIN 
          diff = data - cvals[i]
          
          chisq[i] = transpose(diff)#(reform(cinv##diff))
          
      ENDFOR 
  ENDIF ELSE BEGIN 
      ;; Covariance is just errors.
      FOR i=0L, nc-1 DO BEGIN 
          diff = data - cvals[i]
          
          chisq[i] = total( diff^2/covariance^2 )
          
      ENDFOR
  ENDELSE 

  minchisq = min(chisq)
  wmin = where(chisq EQ minchisq)
  cmin = cvals[wmin[0]]

  degfree = ndata - 1

  IF keyword_set(verbose) THEN BEGIN 
      print,'Min Chisq: ',minchisq,'/',ntostr(long(degfree)),' = ',ntostr(minchisq/degfree)
  ENDIF 

  redchisq = chisq/degfree
  chisq_diff = chisq - minchisq

  w=where(chisq_diff LT 12)


  IF keyword_set(plot) THEN plot, cvals[w], chisq_diff[w], _extra=_extra

  chigh = fltarr(3)
  clow = fltarr(3)
  cerrhigh = fltarr(3)
  cerrlow = fltarr(3)

  lines = [0,2,1]
  FOR i=0,2 DO BEGIN 
      w=where(chisq_diff LT !siglevels1[i])
      tchigh = max(cvals[w], min=tclow)
      chigh[i] = tchigh
      clow[i] = tclow

      cerrhigh[i] = tchigh - cmin
      cerrlow[i] = cmin - tclow
      IF keyword_set(plot) THEN BEGIN 
          oplot, [tclow,tclow], [0, 1.e6],line=lines[i]
          oplot, [tchigh,tchigh], [0, 1.e6],line=lines[i]
      ENDIF 

  ENDFOR 

  IF keyword_set(plot) THEN BEGIN 
      mess2 = $
        !csym.chi+'!U2!N/'+!csym.nu+' = '+$
        ntostr(minchisq,4,/round)+'/'+ntostr(degfree)+$
        ' = '+ntostr(minchisq/degfree,4,/round)
      legend,mess2,/right,box=0
  ENDIF 

END 
