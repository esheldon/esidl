FUNCTION mpfit_power_law, x, p

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    mpfit_power_law
;       
; PURPOSE:
;    returns model power law with parameters p at values X
;
; EXAMPLE:
;	
; Use with mpfit:
; IDL> MYFUNCT = 'mpfit_power_law'
; IDL> START_PARMS = [norm_guess, index_guess]
; IDL> parms = MPFITFUN(MYFUNCT, xvals, yvals, errvals, start_parms, $
; IDL>                  AUTODERIVATIVE=1)
;                          
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; the model
  ymodel = p[0]*x^p[1]

  return, ymodel

END 
