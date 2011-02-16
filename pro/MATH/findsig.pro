PRO findsig, array, bin, sig, nsig=nsig, med=med, zero=zero, _extra=extra

  IF n_params() LT 2 THEN BEGIN
      print,'-Syntax: findsig, array, bin, sig, nsig=nsig, med=med, zero=zero, _extra=extra'
      return
  ENDIF 

  IF n_elements(nsig) EQ 0 THEN nsig = 3

  s = dblarr(nsig)
  sig = dblarr(nsig) & sig[*] = -1.

  FOR i=0, nsig-1 DO BEGIN 
    
      s[i] = errorf( float(i+1)/sqrt(2.) )
  ENDFOR 

  title='Image Histogram'
  ytitle='Number'
  plothist, array, xhist, yhist, bin=bin, title=title, ytitle=ytitle, $
            _extra=extra

  n = total(yhist)

  IF keyword_set(med) THEN BEGIN 
      m = median(xhist)
      w=where(xhist EQ m)
      print,'Median is at: ',xhist[w]
      print,'Median is : ',yhist[w]
  ENDIF ELSE IF keyword_set(zero) THEN BEGIN 
      m = 0.
      w = where(xhist EQ m)
      print,'Value at 0: ',yhist[w]
  ENDIF ELSE BEGIN 
      m = max(yhist) 
      w = where(yhist EQ m, nw)
      IF nw GT 1 THEN BEGIN
          print,'More that one peak'
          return
      ENDIF 
      print,'Max is at: ',xhist[w]
      print,'Max is : ',yhist[w]
  ENDELSE 

  oplot,[xhist[w]],[yhist[w]],psym=4,symsize=2.
  w = w[0]

  step=0
  nx = n_elements(xhist)

  fix=0
  check = 1
  FOR i=0, nsig-1 DO BEGIN 
      istr = ntostr(i+1)
      continue = 1
      WHILE continue DO BEGIN
          
          step = step+1
          i1 = w-step
          i2 = w+step
          IF (i1 LT 0) OR (i2 GT nx-1) THEN BEGIN
              IF check AND (i NE 0) THEN BEGIN
                  fix = fix+1
                  print,'Cannot find the ',istr,' sigma point, extrapolating'
              ENDIF
              IF i NE 0 THEN BEGIN
                  useval = sig[i-1]/i*float(i+1)
                  IF check THEN BEGIN
                      print,'Might try : ',ntostr(i+1),'*1sig = ',$
                        ntostr((i+1)*sig[0])
                      print,'Using: n/(n-1)*sig[n-1] = ',ntostr(useval)
                  ENDIF 
                  sig[i] = useval
              ENDIF ELSE BEGIN
                  print,'Cannot characterize data with this binning'
                  return
              ENDELSE 
              continue = 0
              check=0
          ENDIF ELSE BEGIN 

              t = total( yhist[ i1:i2 ] )
    
              IF float(t)/n GE s[i] THEN BEGIN
                  x1 = xhist[i1]
                  x2 = xhist[i2]
                  y1 = yhist[i1]
                  y2 = yhist[i2]
                  sig[i] = abs( (x2 - x1)/2. )
                  oplot,[x1,x2],[y1,y2],psym=2,symsize=2.
                  continue = 0
              ENDIF 
          ENDELSE 
      ENDWHILE 
      print,istr,' sigma = ',ntostr(sig[i])
  ENDFOR 

  bad = nsig - fix
  leg = strarr(nsig)
  FOR jj = 0, nsig-1 DO BEGIN
      IF jj GE bad THEN BEGIN 
          leg[jj] = ntostr(jj+1)+'sig ~  '+ntostr(rnd(sig[jj],4), 6)
      ENDIF ELSE BEGIN 
          leg[jj] = ntostr(jj+1)+'sig = '+ntostr(rnd(sig[jj],4), 6)
      ENDELSE 
  ENDFOR 
  legend,leg
  xmin = min(xhist)
  xmax = max(xhist)
  x = arrscl(findgen(1000),xmin,xmax)
  f = yhist[w]*exp(-(x-xhist[w])^2/2./sig[0]^2)
  oplot,x,f
print,xhist[w]
print,max(x)
  return
END 


