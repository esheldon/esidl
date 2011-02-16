pro fitlin,xx,yy,sigy,a,siga,b,sigb,chisq,silent=silent

  On_error, 2                   ;Return to Caller

  if N_params() lt 3 then begin   
    print,'Syntax - fitlin, xx, yy, sigy [, a, siga, b, sigb,chisq,silent=silent]' 
    return
  endif

  IF NOT keyword_set(silent) THEN silent = 0

  x = double(xx)                ;Keep input X and Y vectors unmodified
  y = double(yy)
  rn = N_elements(x)

  if rn LT 2 then message,'Input X and Y vectors must contain at least 2 data points'

  if rn NE N_elements(y) then message,'Input X and Y vectors must contain equal number of data points'

; Compute averages and sums
; Taken from Numerical Recipes in c, page  662/663

  S = total(1./sigy^2)
  Sx = total(x/sigy^2)
  Sy = total(y/sigy^2)
  Sxx = total(x^2/sigy^2)
  Sxy = total(x*y/sigy^2)

  delta = S*Sxx - Sx^2
  
  a = (Sxx*Sy - Sx*Sxy)/delta
  b = (S*Sxy  - Sx*Sy)/delta
  siga = sqrt(Sxx/delta)
  sigb = sqrt(S/delta)

  chisq = total( (yy - (a + b*xx))^2/sigy^2 )

  IF NOT silent THEN BEGIN
;      print,'      a               siga              b                sigb           chisq'
      print,'  a         siga           b         sigb           chisq'
      print,ntostr(a),'  ',ntostr(siga),'   ',ntostr(b),'   ',ntostr(sigb),'   ',ntostr(chisq)
  ENDIF 
  return
end
