PRO rand_rtest, v1, v2, nrand, r, s, rr, ss

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: rand_rtest, v1, v2, nrand, r, s, r_rand, s_rand'
      return
  ENDIF 

  t=systime(1)

  COMMON seed,seed
  n=n_elements(v1)

  rtest, v1, v2, r, s

  rr=dblarr(nrand)
  ss=rr

  print,'Doing ',ntostr(nrand),' Random Tests'
  FOR i=0L, nrand-1 DO BEGIN

      a1 = randomu(seed, n)
      a2 = randomu(seed, n)

      srt1 = sort(a1)
      srt2 = sort(a2)

      sv1 = v1[srt1]
      sv2 = v2[srt2]

      rtest, sv1, sv2, tr, ts

      rr[i] = tr
      ss[i] = ts

  ENDFOR

  ptime,systime(1) - t

  
  return
END 
