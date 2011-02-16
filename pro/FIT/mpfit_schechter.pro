
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    mpfit_schechter
;       
; PURPOSE:
;    returns Schechter function with parameters p at values X
;          p[0] = nstar
;          p[1] = xstar (lstar)
;          p[2] = alpha
;
;   func = nstar*(L/Lstar)^alpha * exp(-(L/Lstar)) / Lstar
;
; EXAMPLE:
;	
; Use with mpfit:
; IDL> MYFUNCT = 'mpfit_schechter'
; IDL> START_PARMS = [nstar_guess, lstar_guess, alpha_guess]
; IDL> parms = MPFITFUN(MYFUNCT, xvals, yvals, errvals, start_parms, $
; IDL>                  AUTODERIVATIVE=1)
;                          
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION mpfit_schechter, x, p

  ;; the model
  ymodel = p[0] * (x/p[1])^p[2] * exp(-(x/p[1])) / p[1]
  return, ymodel

END 
