PRO teste, me1, fwhm, err,name=name

e1=.2
n=11
start=0.5

fwhm = fltarr(n) + start
wt = fltarr(n)
me1 = fltarr(n)
err = fltarr(n)

step = 0.05
FOR i=0, n-1 DO BEGIN 

  fwhm[i] = start + i*step
  extractsim, 5, cat, fwhm=fwhm[i]

  ;; find stars
  w=where(cat.e1_ad NE -10.0 AND cat.stype EQ 2 AND cat.mag_best LT 20.0)
  c=cat[w]

  me1[i] = median(c.e1_ad)
  m = moment(c.e1_ad)
  err[i] = sqrt(m[1])
  ;; do weighted averages
  ;wt = 1.0/c.e1_aderr^2
  ;wtsum = total(wt)
  ;me1[i] = total( wt*c.e1_ad )/wtsum
  ;err[i] = sqrt(1.0/wtsum)
  
ENDFOR

IF keyword_set(name) THEN BEGIN
  dir = '/sdss4/data1/esheldon/'
  fname = dir+name
  makeps,fname
ENDIF
plot, fwhm, me1, xtitle='fwhm (arcsec)',ytitle='e1 psf',$
         title = 'Measurement Bias',psym=4

oploterr, fwhm, me1, err
oplot,[0,10], [e1,e1]

IF keyword_set(name) THEN ep

return
END

