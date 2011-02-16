pro fitlin_old,xx,yy,sigy,a,siga,b,sigb,silent=silent

;	fitlin, xx, yy, a, siga, b, sigb   
  On_error, 2                   ;Return to Caller

  if N_params() lt 3 then begin   
    print,'Syntax - fitlin, xx, yy, sigy [, a, siga, b, sigb,silent=silent]' 
    return
  endif

  IF NOT keyword_set(silent) THEN silent = 0

  x = double(xx)                ;Keep input X and Y vectors unmodified
  y = double(yy)
  rn = N_elements(x)

  if rn LT 2 then $ message,'Input X and Y vectors must contain at least 2 data points'

  if rn NE N_elements(y) then message,'Input X and Y vectors must contain equal number of data points'

; Compute averages and sums

  xxx=total(x/(sigy)^2)
  yyy=total(y/(sigy)^2)
  ss=total(1./(sigy)^2)
  xavg =xxx/ss

  x=(x-xavg)/(sigy)
  totx=total(x*x)
  b=total(x*y/sigy)/totx
  a=(yyy-xxx*b)/ss

  siga=sqrt((1+xxx*xxx/(ss*totx))/ss)
  sigb=sqrt(1./totx)

  IF NOT silent THEN print,a,siga,b,sigb,n_elements(x)
  
  return
end
