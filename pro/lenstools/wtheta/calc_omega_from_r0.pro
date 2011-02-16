PRO calc_omega_from_r0, norm, normerr, power, powerr, r0, r0err, omega, omega_errl, omega_errh

  ;; norm in Msolar/Mpc^2, r0 in Mpc -> returns Omega (unitless)

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax:  calc_omega_from_r0, norm, normerr, power, powerr, r0, r0err, omega, omega_errl, omega_errh'
      return
  ENDIF 

  rho_crit = double(2.77545e-7) ;Msolar/pc^3
  rho_crit = rho_crit*(10d)^18. ; Msolar/Mpc^3

  gam = 1. + power

  f=gamma(0.5)*gamma(0.5*(gam-1))/gamma(0.5*gam)

  omega = (norm/rho_crit)/(f * r0)
  
  ;; calculate extremes to get error (hard analytically)
  normext = [norm-normerr, norm+normerr]
  powext = [power-powerr, power+powerr]
  r0ext = [r0-r0err, r0+r0err]

  FOR nn=0,1 DO BEGIN 
      FOR np=0,1 DO BEGIN 
          FOR nr=0,1 DO BEGIN 
              
              tgam = 1.+powext[np]
              tf=gamma(0.5)*gamma(0.5*(tgam-1))/gamma(0.5*tgam)

              tomega = (normext[nn]/rho_crit)/(tf*r0ext[nr])

              add_arrval, tomega, tomega_ext

          ENDFOR 
      ENDFOR 
  ENDFOR 

  min_omega = min(tomega_ext, max=max_omega)
  
  omega_errh = max_omega-omega
  omega_errl = omega-min_omega

  print,omega,' + ',omega_errh,' - ',omega_errl

END 
