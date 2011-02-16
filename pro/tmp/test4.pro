PRO test4, nn

  n=long(nn)

  nmid = n-d > 0
  IF n EQ 1 THEN nlast=0 ELSE nlast=1

  dim = 2L^n

  arr = intarr(dim,dim)

  elem = intarr(n, 4)
  FOR i=0L, n-1 DO elem[i,*] = [1,1,1,-1]
  elem = [1,1,1,-1]
  elem = replicate(elem,n)


  fi=0L
  midi = lonarr(nmid)
  lasti=0L
  FOR i=0L, dim-1 DO BEGIN 
      FOR j=0L, dim-1 DO BEGIN 

          fval = elem[fi]

          val = fval
          FOR mi=0L, nmid-1 DO BEGIN 
              val=val*elem[midi[mi]]
          ENDFOR 
          IF nlast NE 0 THEN val=val*elem[lasti]

      ENDFOR 
  ENDFOR 
      

          

  ENDFOR 

return
END 
