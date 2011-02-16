;+
; NAME:
;  GET_LIKE_CONF
;
;
; PURPOSE:
;  Use monte-carlo to convert the input likelihood to 1,2, and 3-sigma
;  values [68.3,95.5,99.7]
;
;
; CATEGORY:
;  Statistics
;
;
; CALLING SEQUENCE:
;  get_like_conf, vals, likelihood, vmax, vlow, vhigh, verrl, verrh, $
;                  silent=silent, domean=domean
;
;
; INPUTS:
;  vals: the x-values corresponding to the likelihood.
;  likelihood: the likelihood as a function of "vals".
;
;
; KEYWORD PARAMETERS:
;  /silent: shut up
;  /domean: Use the mean value rather than the maximum likelihood value.
;
;
; OUTPUTS:
;  vmax: maximum value (or mean if /domean)
;  vlow: 1,2,3 sigma lower limits.
;  vhigh: 1,2,3 sigma upper limits.
;  verrl: 1,2,3 sigma errors on low side.
;  verrh: 1,2,3 sigma errors on high side.
;
; EXAMPLE:
;  
;
;
; MODIFICATION HISTORY:
;  Created ??-Jan-2003: Erin Sheldon UofChicago
;
;-


PRO get_like_conf, vals, likelihood, vmax, vlow, vhigh, verrl, verrh, $
                   silent=silent, domean=domean, input_vmax=input_vmax

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: get_like_conf, vals, likelihood, vmax, vlow, vhigh, verrl, verrh, silent=silent, domean=domean'
      print
      print,'Default is vmax=maximum value. if /domean, then the expectation value is calculated'
      return
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Uses a monte carlo method to get the 1,2,3-sigma regions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ranges = [.683, .955, .997]
  halfranges = ranges/2.

  wgood=where( (likelihood EQ likelihood) AND (likelihood NE 0), ngood)
  IF ngood EQ 0 THEN message,'No good values in likelihood'

  IF n_elements(input_vmax) NE 0 THEN BEGIN
      vmax = input_vmax 
  ENDIF ELSE IF NOT keyword_set(domean) THEN BEGIN 
      ml = max(likelihood[wgood], /NAN)
      wmax = where(likelihood[wgood] EQ ml, nwmax)
      
      IF nwmax EQ 0 THEN message,'Bad values in likelihood'

      vmax = vals[ wgood[wmax[0]] ]
  ENDIF ELSE BEGIN 
      npts = ngood*2L
      x1 = min( vals[wgood], max=x2)
      
      gauleg, x1, x2, npts, XXi, WWi
      
      funcvals = interpol(likelihood[wgood], vals[wgood], XXi)
      norm = total(funcvals*WWi)

      vmax = total( funcvals*XXi*WWi )/norm

  ENDELSE 

  ss = sort(vals[wgood])
  ss = wgood[ss]

  nrand = 100000L
  genrand, likelihood[ss], vals[ss], nrand, randvals, /quadratic, /double

  vlow = dblarr(3)
  vhigh = vlow
  verrl = vlow
  verrh = vlow

  wlow  = where(randvals LE vmax, nlow)
  IF nlow NE 0 THEN BEGIN 
      lfrac = findgen(nlow)/nrand

      sl = sort(randvals[wlow])
      sl = wlow[sl]
      sl = reverse(sl)
  ENDIF ELSE BEGIN 
      vlow = replicate(vmax, 3)
      verrl = replicate(1.e10, 3)
  ENDELSE 

  whigh = where(randvals GE vmax, nhigh)
  IF nhigh NE 0 THEN BEGIN 
      hfrac = findgen(nhigh)/nrand

      sh = sort(randvals[whigh])
      sh = whigh[sh]
  ENDIF ELSE BEGIN 
      vhigh = replicate(vmax, 3)
      verrh = replicate(1.e10, 3)
  ENDELSE 

  FOR level=0,2 DO BEGIN 

      frac = halfranges[level]

      IF nlow NE 0 THEN BEGIN 
          ;; fraction of objects below, as we add more and more
          wflow = max(where(lfrac LE frac))
          vlow[level] = randvals[sl[wflow]]
          ;; convert ranges to errors
          verrl[level] = vmax-vlow[level]
      ENDIF 

      IF nhigh NE 0 THEN BEGIN 
          ;; fraction of objects above, as we add more and more
          wfhigh = max(where(hfrac LE frac))
          vhigh[level] = randvals[sh[wfhigh]]
          ;; convert ranges to errors
          verrh[level] = vhigh[level]-vmax
      ENDIF 

  ENDFOR 

  randvals = 0

  IF NOT keyword_set(silent) THEN BEGIN 
      print
      print,ntostr(2.*halfranges[0]*100., 4, /round)+'% region'
      print,'val = '+ntostr(vmax)+' + '+ntostr(verrh[0])+' - '+ntostr(verrl[0])
  ENDIF 

END 
