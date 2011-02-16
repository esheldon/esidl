PRO pow_chisq_conf_1par_gen, datax, data, covariance, $
                             normrange, nnorm, power, $
                             chisq_surf, nmin, normlow, normhigh, $
                             yfit=yfit,$
                             normallow=normallow, $
                             chisq_per=chisq_per, $
                             chisq_diff=chisq_diff, noplotmin=noplotmin, $
                             aspect=aspect, center=center, $
                             _extra=extra, xtick_get=xtick_get, $
                             minchisq=minchisq, degfree=degfree, $
                             ytick_get=ytick_get, wuse=wuse

  ;; right now only changes norm

  IF n_params() LT 6 THEN BEGIN
      print,'-Syntax: pow_chisq_conf_1par_gen, datax, data, covariance, '
      print,'           normrange, nnorm, power, '
      print,'           chisq_surf, nmin, normlow, normhigh, '
      print,'           yfit=yfit,'
      print,'           normallow=normallow, '
      print,'           chisq_per=chisq_per, '
      print,'           chisq_diff=chisq_diff, noplotmin=noplotmin, '
      print,'           aspect=aspect, center=center, '
      print,'           _extra=extra, xtick_get=xtick_get, '
      print,'           minchisq=minchisq, degfree=degfree, '
      print,'           ytick_get=ytick_get, wuse=wuse'
      return
  ENDIF 
  nlevels = 3                   ;68, 95, 99.78%
  line = [0,2,1]

  errcode=['Successful completion',$
           'Singular Covariance Matrix', $
           'Small pivot element in Covariance Matrix: Accuracy Compromised']



  normvals = arrscl( findgen(nnorm), normrange[0], normrange[1] )

  default = 1.e10

  normlow = replicate(default,nlevels)
  normhigh = replicate(default,nlevels)

  nnorm = n_elements(normvals)

  chisq_surf = fltarr( nnorm )
  modfunc = data

  index = lindgen(nnorm)

  ndata = n_elements(data)
  cs = size(covariance)
  IF (cs[0] EQ 2) THEN BEGIN ;; 2-D?
      IF cs[1] NE ndata THEN BEGIN
          message,'Covariance matrix must be NdataXNdata',/inf
          return
      ENDIF 
      docov = 1
  ENDIF ELSE BEGIN ;; its 1-D
      docov = 0
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; which data points do we use?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(wuse) EQ 0 THEN BEGIN
      ndata = n_elements(data)
      wuse = lindgen(ndata)
  ENDIF ELSE BEGIN 
      ndata = n_elements(wuse)
  ENDELSE 

  IF NOT docov THEN BEGIN 
      ;; Input covariance is really just dataerr
      errsquared=covariance^2

      FOR in=0L, nnorm-1 DO BEGIN ; loop over norms
          
          norm = normvals[in]
          modfunc[wuse] = norm*(datax[wuse])^power
          
          ff= (data[wuse]-modfunc[wuse])^2/errsquared[wuse]
          chisq_surf[in] = total(ff)
          
      ENDFOR 

  ENDIF ELSE BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; copy in block corresponding to wuse indices
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      tmpcov = replicate(covariance[0], ndata, ndata)
      FOR i=0L, ndata-1 DO BEGIN 
          FOR j=0L, ndata-1 DO BEGIN 
              tmpcov[i,j] = covariance[wuse[i], wuse[j]]
          ENDFOR 
      ENDFOR 

      cinv = invert(tmpcov, stat, /double)
      IF stat NE 0 THEN message,errcode[stat]

      FOR in=0L, nnorm-1 DO BEGIN ; loop over norms
          norm = normvals[in]

          modfunc[wuse] = norm*(datax[wuse])^power
          diff = data[wuse] - modfunc[wuse]

          chisq_surf[in] = transpose(diff)#(reform(cinv##diff))

      ENDFOR 

  ENDELSE 

  minchisq = min(chisq_surf, /nan)
  
  chisq_diff = chisq_surf - minchisq

  w=where(chisq_surf EQ minchisq, nw)
  wparam=w
  nmin=normvals[w[0]]
  degfree = ndata-1

  wL = where(chisq_surf/degfree LE 8.,nwL)
  IF nwL EQ 0 THEN wL = lindgen(n_elements(normvals))

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

  yfit = nmin*datax[wuse]^power

  return
END 
