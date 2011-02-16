PRO testslength, maxs2n
  

  IF n_elements(maxs2n) EQ 0 THEN maxs2n = 1000 $
  ELSE print,'Using max S/N: ',ntostr(maxs2n)
      

  f='/sdss4/data1/esheldon/GRIDSHEAR/KAPPA/s120names'
  readcol,f,format='A',names

  nn = n_elements(names)

  t=mrdfits(names[0],/silent)
  s=size(t)

  sx=s[1]
  sy=s[2]
  
  index = lindgen(sx*sy)

  x = index MOD sx
  y = index/sx

  cenx = (sx-1)/2.
  ceny = (sy-1)/2.

  noise = .046
  FOR i=0, nn-1 DO BEGIN

      t=mrdfits(names[i],/silent)
      max = max(t)

      IF max/noise GE maxs2n THEN BEGIN 
          w = where(t EQ max)

          print,'Max: ',ntostr(max/noise),'   N = ',ntostr(i)
          mx = x[w[0]]-cenx
          my = y[w[0]]-ceny
          print,'At    [',ntostr(mx,4),', ',ntostr(my,4),']'

          dist = sqrt( mx^2 + my^2)
          print,'Total Offset: ',ntostr(dist)
          print
      ENDIF 
  ENDFOR 

  return
END 

