PRO calc_ratio_cov, x, covx, y, covy, ratio, ratiocov

  IF n_params() LT 4 THEN BEGIN  
      print,'-Syntax: calc_ratio_cov, x, covx, y, covy, ratio, ratiocov'
      print,'Compute ratio x/y and the covariance'
      print,'If input covariances are vectors, then ratioerr is returned'
      return
  ENDIF 

; calculate the covariance matrix for the ratio x/y, given the 
; covariances on x and y

  delvarx, ratio, ratiocov

  nx = n_elements(x)
  ny = n_elements(y)

  IF nx NE ny THEN BEGIN
      print,'x must be same size as y'
      return
  ENDIF 

  ;; vectors were entered for covx, covy: interpret as errors
  ncovx = n_elements(covx) & ncovy = n_elements(covy)
  IF ncovx EQ nx AND ncovy EQ nx THEN BEGIN 
      ratio = x/y
      ratiocov = sqrt( (covx/x)^2 + (covy/y)^2 )
      return
  ENDIF 

  IF ncovx NE nx^2 OR ncovy NE ny^2 THEN BEGIN 
      print,'input covariances must either be square matricies '
      print,'with size nxXnx, nyXny or vector error arrays'
      return
  ENDIF 

  ratiocov = covx
  ratiocov[*] = 0
  ratio = x
  ratio[*] = 0

  FOR i=0L, nx-1 DO BEGIN 
      xi = x[i]
      yi = y[i]

      ratioi = xi/yi
      FOR j=0L, nx-1 DO BEGIN 
          xj = x[j]
          yj = y[j]

          ratioj = xj/yj

          covxij = covx[i,j]
          covyij = covy[i,j]

          IF i EQ j THEN ratio[i] = xi/yi

          ratiocov[i,j] = ratioi*ratioj*( covxij/xi/xj + covyij/yi/yj )
          
      ENDFOR 

  ENDFOR 

return

END 
