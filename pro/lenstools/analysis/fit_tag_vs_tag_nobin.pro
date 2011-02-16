PRO pow_scale_fit_one_bandpass, datax, data, dataerr, power, norm, fitstruct, $
                                _extra=extra

  pow_chisq_conf, datax, data, dataerr, power, norm, $
    chisq_surf, bestpow, bestnorm, powlow, powhigh, normlow, normhigh, $
    xtitle='!7a!X', ytitle='A', _extra=extra
  clevel = 0
  
  IF n_elements(bestpow) EQ 0 THEN return

  highdiff = normhigh[clevel]-bestnorm & lowdiff = bestnorm - normlow[clevel]
  mess = 'A = '+ntostr(bestnorm,6)+$
    '!S!U+'+ntostr(highdiff,5)+'!R!D-'+ntostr(lowdiff,5)
  mess = [mess, '']
  highdiff = powhigh[clevel]-bestpow & lowdiff = bestpow - powlow[clevel]
  mess = [mess, $
          '!7a!X = '+ntostr(bestpow, 6) + $
          '!S!U+'+ntostr(highdiff,5)+'!R!D-'+ntostr(lowdiff,5)]
  legend, mess, /right, bcolor=bcolor, charsize=charsize
  
  print,'Norm: '+ntostr(bestnorm),' + ',$
    ntostr(normhigh[clevel]-bestnorm),' - ',ntostr(bestnorm-normlow[clevel])
  print,'Power: '+ntostr(bestpow),' + ',$
    ntostr(powhigh[clevel]-bestpow),' - ',ntostr(bestpow-powlow[clevel])

  fitstruct = create_struct('bestpow', bestpow, $
                            'bestnorm', bestnorm, $
                            'powlow', powlow, $
                            'powhigh', powhigh, $
                            'normlow', normlow, $
                            'normhigh', normhigh)

  return
END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Main
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO fit_tag_vs_tag_nobin, lensum, xtag, ytag, yerrtag, normrange, powrange, $
                          fitstruct, $
                          nnorm=nnorm,npower=npower, xrad=xrad, yrad=yrad

  IF n_params() LT 6 THEN BEGIN 
      print,'-Syntax: fit_tag_vs_tag_nobin, lensum, xtag, ytag, yerrtag, normrange, powrange, nnorm=nnorm,npower=npower'
      return
  ENDIF 

  IF n_elements(nnorm) EQ 0 THEN nnorm=50L
  IF n_elements(npower) EQ 0 THEN npower=50L
  IF n_elements(xrad) EQ 0 THEN xrad=0
  IF n_elements(yrad) EQ 0 THEN yrad=0

  ;; need tag, maxtag, or xrange for model, normrange, powrange,  nx, npower, nnorm

  ;; this is to test fitting to scatter mass plot using
  ;; my chisq fitter (extremely slow for this large data set)

  ;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;

  IF NOT tag_exist(lensum, xtag, index=wxt) THEN BEGIN 
      print,'Tag '+tag+' not found'
      return
  ENDIF 
  IF NOT tag_exist(lensum, ytag, index=wyt) THEN BEGIN 
      print,'Tag '+tag+' not found'
      return
  ENDIF 
  IF NOT tag_exist(lensum, yerrtag, index=wyerrt) THEN BEGIN 
      print,'Tag '+tag+' not found'
      return
  ENDIF 

  colors = ['u','g','r','i','z']

  IF (!d.flags AND 1) EQ 0 THEN doX=1 ELSE doX=0
  IF doX THEN BEGIN
      !p.charsize = 1.0
      charsize = 1.2
      !p.thick = 1
      !x.thick = 1
      !y.thick = 1
      !p.charthick=1
  ENDIF ELSE BEGIN
      charsize=1.0
      !p.charsize = 1.0
      !p.thick = 5
      !x.thick = 5
      !y.thick = 5
      !p.charthick = 4
  ENDELSE 

  simpctable
  !p.background=!white
  !p.color = !black

  time=systime(1)

  nl = n_elements(lensum)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; get tag/mass values for each lens
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  xtag_val = lensum.(wxt)[xrad]
  ytag_val = lensum.(wyt)[yrad]
  yerrtag_val = lensum.(wyerrt)[yrad]

  maxtag = max(xtag_val) & mintag = min(ytag_val)
  IF mintag LT 0. THEN xmin = 1.1*mintag ELSE xmin = 0.9*mintag
  IF maxtag LT 0. THEN xmax = 0.9*maxtag ELSE xmax = 1.1*maxtag

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; create the model
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  norm = arrscl( findgen(nnorm), normrange[0], normrange[1] )
  power = arrscl( findgen(npower), powrange[0], powrange[1] )
 

  pow_scale_fit_one_bandpass, xtag_val, ytag_val, yerrtag_val, power, norm, $
    fitstruct
      
  return

END 
