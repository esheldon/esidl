PRO calculate_cov_3band2, e1, e2, e1e1err, e1e2err, e2e2err, cove, cove1, cove2, $
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

  sigma_clip, e1[g,*]*e1[g,*], v1g1g, tsig, nsig=3, niter=4,/silent
  sigma_clip, e1[r,*]*e1[g,*], v1r1g, tsig, nsig=3, niter=4,/silent
  sigma_clip, e1[i,*]*e1[g,*], v1i1g, tsig, nsig=3, niter=4,/silent
  sigma_clip, e2[g,*]*e1[g,*], v2g1g, tsig, nsig=3, niter=4,/silent
  sigma_clip, e2[r,*]*e1[g,*], v2r1g, tsig, nsig=3, niter=4,/silent
  sigma_clip, e2[i,*]*e1[g,*], v2i1g, tsig, nsig=3, niter=4,/silent

  ;; diag term
  IF keyword_set(fixdiag) THEN BEGIN
      sigma_clip, e1e1err[g,*]^2, me1gerr, nsig=3, niter=4, /silent
      v1g1g = v1g1g - me1gerr
  ENDIF 

  ;; 2nd column
  v1g1r = v1r1g
  sigma_clip, e1[r,*]*e1[r,*], v1r1r, tsig, nsig=3, niter=4,/silent
  sigma_clip, e1[i,*]*e1[r,*], v1i1r, tsig, nsig=3, niter=4,/silent
  sigma_clip, e2[g,*]*e1[r,*], v2g1r, tsig, nsig=3, niter=4,/silent
  sigma_clip, e2[r,*]*e1[r,*], v2r1r, tsig, nsig=3, niter=4,/silent
  sigma_clip, e2[i,*]*e1[r,*], v2i1r, tsig, nsig=3, niter=4,/silent

  ;; diag term
  IF keyword_set(fixdiag) THEN BEGIN
      sigma_clip, e1e1err[r,*]^2, me1rerr, nsig=3, niter=4, /silent
      v1r1r = v1r1r - me1rerr
  ENDIF 


  ;; 3rd column
  v1g1i = v1i1g
  v1r1i = v1i1r
  sigma_clip, e1[i,*]*e1[i,*], v1i1i, tsig, nsig=3, niter=4,/silent
  sigma_clip, e2[g,*]*e1[i,*], v2g1i, tsig, nsig=3, niter=4,/silent
  sigma_clip, e2[r,*]*e1[i,*], v2r1i, tsig, nsig=3, niter=4,/silent
  sigma_clip, e2[i,*]*e1[i,*], v2i1i, tsig, nsig=3, niter=4,/silent

  ;; diag term
  IF keyword_set(fixdiag) THEN BEGIN
      sigma_clip, e1e1err[i,*]^2, me1ierr, nsig=3, niter=4, /silent
      v1i1i = v1i1i - me1ierr
  ENDIF 

  ;; 4th column
  v1g2g = v2g1g
  v1r2g = v2g1r
  v1i2g = v2g1i
  sigma_clip, e2[g,*]*e2[g,*], v2g2g, tsig, nsig=3, niter=4,/silent
  sigma_clip, e2[r,*]*e2[g,*], v2r2g, tsig, nsig=3, niter=4,/silent
  sigma_clip, e2[i,*]*e2[g,*], v2i2g, tsig, nsig=3, niter=4,/silent


  ;; diag term
  IF keyword_set(fixdiag) THEN BEGIN
      sigma_clip, e2e2err[g,*]^2, me2gerr, nsig=3, niter=4, /silent
      v2g2g = v2g2g - me2gerr
  ENDIF 

  ;; 5th column
  v1g2r = v2r1g
  v1r2r = v2r1r
  v1i2r = v2r1i
  v2g2r = v2r2g
  sigma_clip, e2[r,*]*e2[r,*], v2r2r, tsig, nsig=3, niter=4,/silent
  sigma_clip, e2[i,*]*e2[r,*], v2i2r, tsig, nsig=3, niter=4,/silent

  ;; diag term
  IF keyword_set(fixdiag) THEN BEGIN
      sigma_clip, e2e2err[r,*]^2, me2rerr, nsig=3, niter=4, /silent
      v2r2r = v2r2r - me2rerr
  ENDIF 

  ;; 6th column
  v1g2i = v2i1g
  v1r2i = v2i1r
  v1i2i = v2i1i
  v2g2i = v2i2g
  v2r2i = v2i2r
  sigma_clip, e2[i,*]*e2[i,*], v2i2i, tsig, nsig=3, niter=4,/silent

  ;; diag term
  IF keyword_set(fixdiag) THEN BEGIN
      sigma_clip, e2e2err[i,*]^2, me2ierr, nsig=3, niter=4, /silent
      v2i2i = v2i2i - me2ierr
  ENDIF 

  v1g2g = 0.0 & v2g1g = 0.0
  v1g2r = 0.0 & v2g1r = 0.0
  v1g2i = 0.0 & v2g1i = 0.0

  v1r2g = 0.0 & v2r1g = 0.0
  v1r2r = 0.0 & v2r1r = 0.0
  v1r2i = 0.0 & v2r1i = 0.0

  v1i2g = 0.0 & v2i1g = 0.0
  v1i2r = 0.0 & v2i1r = 0.0
  v1i2i = 0.0 & v2i1i = 0.0

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
