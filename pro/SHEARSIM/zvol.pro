FUNCTION  zvol, zmax, zmin, nbins, h=h, omega=omega

  IF n_params() EQ 0 THEN BEGIN
    print,'-Syntax: result = zvol(zmax, zmin, nbins)
    return, 0.
  ENDIF 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Phil's formula (per square degree)
  ; dv=3.0461742e-4*c/H0*dlum(z)^2/(sqrt(1.+omega*z)*(1.+z)^6)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nmax = n_elements(zmax)
  nmin = n_elements(zmin)
  IF nmax NE nmin THEN BEGIN
    print,'zmax and zmin must be same size'
    return,0.
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(h) EQ 0 THEN h=.7
  IF n_elements(omega) EQ 0 THEN omega = 1.0
  c = 2.9979e5                  ;km/s
  H0 = h*100.                   ;km/s/Mpc

  v = dblarr(nmax)

  FOR i=0, nmax-1 DO BEGIN 
    z = findgen(nbins)/(nbins-1)*zmax[i] + zmin[i]
    dlum = lumdist(z)

    dv = 1.0/H0*dlum^2/( sqrt(1.+omega*z)*(1.+z)^6 )
    v[i] = int_tabulated(z, dv)
  ENDFOR 

  return, v

END 


