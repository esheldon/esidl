PRO erange, cat

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Look at a range of e1 and e2 valuew
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

COMMON seed,seed

s=systime(1)

n = 40
nfields = 5

fwhm = fltarr(n)+1.5

;;; Use uniform sample of e1 and e2 in [-.5,.5]
tn = 1000

continue = 1
WHILE continue DO BEGIN 
  continue = 0
  e1 = randomu(seed,tn)-.5
  e2 = randomu(seed,tn)-.5
  w=where( e1^2 + e2^2 LT 1.0,nw)
  IF nw LT n  THEN BEGIN
    print,'Not enough good ones'
    return
  ENDIF

  e1 = e1[ w[0:n-1] ]
  e2 = e2[ w[0:n-1] ]

  plothist2, e1, bin=.05
  plothist2, e2, bin=.05,/overplot,linestyle=2
  
  print,'Are these ok (y/n)?  q to quit'
  key = get_kbrd(20)
  IF key EQ 'n' THEN continue = 1
  IF key EQ 'q' THEN return
ENDWHILE

findabtheta, e1, e2, aratio, posangle

posangle=posangle*180.0/!pi

print,e1
print
print,e2
print
name = '/sdss4/data1/esheldon/test_ell_5.fit'
print
print,'Output directed to : ',name

FOR i=0, n-1 DO BEGIN
  print,'Doing ',strtrim(string(i+1),2),' of ',$
         strtrim(string(n),2),' ',strtrim(string(nfields),2),' Fields each'
  cat=0
  extractsim, nfields, cat,aratio=aratio[i],posangle=posangle[i],fwhm=fwhm[i]
  cat.typeflag = i

  IF i EQ 0 THEN mwrfits, cat, name,/create ELSE mwrfits, cat, name

ENDFOR

print,strtrim(string( (systime(1)-s)/3600.0 ), 2),' Total Hours'

return
END





