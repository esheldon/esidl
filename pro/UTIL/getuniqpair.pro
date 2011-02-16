PRO getuniqpair, v1, v2, uv1, uv2, num

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: getuniqpair, v1, v2, uv1, uv2, num'
      return
  ENDIF 

  ;; Only good for arrays with many duplicates

  nv1 = n_elements(v1)
  nv2 = n_elements(v2)

  IF nv1 NE nv2 THEN return

  rm=rem_dup(v1)
  u1 = v1[rm]
  n1 = n_elements(u1)

  delvarx,uv1,uv2

  FOR i=0L, n1-1 DO BEGIN 

      w=where(v1 EQ u1[i], nw)

      u=v2[w[rem_dup(v2[w])]]
      nu=n_elements(u)
      IF n_elements(uv1) EQ 0 THEN BEGIN 
          uv1 = replicate(u1[i], nu)
          uv2 = u
      ENDIF ELSE BEGIN 
          uv1 = [uv1, replicate(u1[i], nu)]
          uv2 = [uv2, u]
      ENDELSE 

  ENDFOR 

  nuniq = n_elements(uv1)
  num = lonarr(nuniq)
  FOR i=0L, nuniq-1 DO BEGIN 
      w=where(v1 EQ uv1[i] AND v2 EQ uv2[i],nw)
      num[i] = nw
  ENDFOR 

  return
END 
