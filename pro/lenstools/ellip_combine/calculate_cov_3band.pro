PRO calculate_cov_3band, e1, e2, e1e1err, e1e2err, e2e2err, cove, cove1, cove2, $
                         fixdiag=fixdiag

  siz=size(e1)
  IF siz[1] NE 3 THEN message,'Must enter arrays of [3,n]'

  nobj = siz[2]
  print,'Number of objects: ',nobj
  
  defval = -9999.
  rcut = 0.6
  errcut = 0.4

  g=0
  r=1
  i=2

  rn = 1./(nobj-1.0)

  ;; shot noise for diag term screws up diagonal

  ;; 1st column

  v1g1g = total( e1[g,*]*e1[g,*] )*rn
  v1r1g = total( e1[r,*]*e1[g,*] )*rn
  v1i1g = total( e1[i,*]*e1[g,*] )*rn
  v2g1g = total( e2[g,*]*e1[g,*] )*rn
  v2r1g = total( e2[r,*]*e1[g,*] )*rn
  v2i1g = total( e2[i,*]*e1[g,*] )*rn
  ;; diag term
  IF keyword_set(fixdiag) THEN BEGIN
      v1g1g = v1g1g - mean(e1e1err[g,*]^2)
  ENDIF 

  ;; 2nd column
  v1g1r = v1r1g
  v1r1r = total( e1[r,*]*e1[r,*] )*rn
  v1i1r = total( e1[i,*]*e1[r,*] )*rn
  v2g1r = total( e2[g,*]*e1[r,*] )*rn
  v2r1r = total( e2[r,*]*e1[r,*] )*rn
  v2i1r = total( e2[i,*]*e1[r,*] )*rn
  ;; diag term
  IF keyword_set(fixdiag) THEN BEGIN
      v1r1r = v1r1r - mean(e1e1err[r,*]^2)
  ENDIF 


  ;; 3rd column
  v1g1i = v1i1g
  v1r1i = v1i1r
  v1i1i = total( e1[i,*]*e1[i,*] )*rn
  v2g1i = total( e2[g,*]*e1[i,*] )*rn
  v2r1i = total( e2[r,*]*e1[i,*] )*rn
  v2i1i = total( e2[i,*]*e1[i,*] )*rn
  ;; diag term
  IF keyword_set(fixdiag) THEN BEGIN
      v1i1i = v1i1i - mean(e1e1err[i,*]^2)
  ENDIF 

  ;; 4th column
  v1g2g = v2g1g
  v1r2g = v2g1r
  v1i2g = v2g1i
  v2g2g = total( e2[g,*]*e2[g,*] )*rn
  v2r2g = total( e2[r,*]*e2[g,*] )*rn
  v2i2g = total( e2[i,*]*e2[g,*] )*rn
  ;; diag term
  IF keyword_set(fixdiag) THEN BEGIN
      v2g2g = v2g2g - mean(e2e2err[g,*]^2)
  ENDIF 

  ;; 5th column
  v1g2r = v2r1g
  v1r2r = v2r1r
  v1i2r = v2r1i
  v2g2r = v2r2g
  v2r2r = total( e2[r,*]*e2[r,*] )*rn
  v2i2r = total( e2[i,*]*e2[r,*] )*rn
  ;; diag term
  IF keyword_set(fixdiag) THEN BEGIN
      v2r2r = v2r2r - mean(e2e2err[r,*]^2)
  ENDIF 

  ;; 6th column
  v1g2i = v2i1g
  v1r2i = v2i1r
  v1i2i = v2i1i
  v2g2i = v2i2g
  v2r2i = v2i2r
  v2i2i = total( e2[i,*]*e2[i,*] )*rn
  ;; diag term
  IF keyword_set(fixdiag) THEN BEGIN
      v2i2i = v2i2i - mean(e2e2err[i,*]^2)
  ENDIF 

  ;;v1g2g = 0.0 & v2g1g = 0.0
  ;;v1g2r = 0.0 & v2g1r = 0.0
  ;;v1g2i = 0.0 & v2g1i = 0.0

  ;;v1r2g = 0.0 & v2r1g = 0.0
  ;;v1r2r = 0.0 & v2r1r = 0.0
  ;;v1r2i = 0.0 & v2r1i = 0.0

  ;;v1i2g = 0.0 & v2i1g = 0.0
  ;;v1i2r = 0.0 & v2i1r = 0.0
  ;;v1i2i = 0.0 & v2i1i = 0.0

  ;; goes e1,e1,e1  e2,e2,e2
  cove = [ [v1g1g, v1g1r, v1g1i, v1g2g, v1g2r, v1g2i], $
           [v1r1g, v1r1r, v1r1i, v1r2g, v1r2r, v1r2i], $
           [v1i1g, v1i1r, v1i1i, v1i2g, v1i2r, v1i2i], $
           [v2g1g, v2g1r, v2g1i, v2g2g, v2g2r, v2g2i], $
           [v2r1g, v2r1r, v2r1i, v2r2g, v2r2r, v2r2i], $
           [v2i1g, v2i1r, v2i1i, v2i2g, v2i2r, v2i2i] ]
  
  
  cove1 = [ [v1g1g, v1g1r, v1g1i], $
            [v1r1g, v1r1r, v1r1i], $
            [v1i1g, v1i1r, v1i1i] ]

  cove2 = [ [v2g2g, v2g2r, v2g2i], $
            [v2r2g, v2r2r, v2r2i], $
            [v2i2g, v2i2r, v2i2i] ]



END 
