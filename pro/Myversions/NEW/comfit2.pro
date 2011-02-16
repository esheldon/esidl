;$Id: comfit2.pro,v 1.1.1.1 2004/11/11 20:17:22 esheldon Exp $
;
; Copyright (c) 1994-1998, Research Systems, Inc.  All rights reserved.
;       Unauthorized reproduction prohibited.
;+
; NAME:
;       COMFIT
;
; PURPOSE:
;       This function fits the paired data {X(i), Y(i)} to one of six common
;       types of appoximating models using a gradient-expansion least-squares
;       method. The result is a vector containing the model parameters.
;
; CATEGORY:
;       Statistics.
;
; CALLING SEQUENCE:
;       Result = COMFIT(X, Y, A)
;
; INPUTS:
;       X:    An n-element vector of type integer, float or double.
;
;       Y:    An n-element vector of type integer, float or double.
;
;       A:    A vector of initial estimates for each model parameter. The length
;             of this vector depends upon the type of model selected.
;
; KEYWORD PARAMETERS:
;       EXPONENTIAL:    If set to a non-zero value, the parameters of the 
;                       exponential model are computed. Y = a0  * a1^x + a2.
;
;         GEOMETRIC:    If set to a non-zero value, the parameters of the
;                       geometric model are computed.  Y = a0 * x^a1 + a2.
;    
;          GOMPERTZ:    If set to a non-zero value, the parameters of the
;                       Gompertz model are computed.  Y = a0 * a1^(a2*x) + a3.
;
;        HYPERBOLIC:    If set to a non-zero value, the parameters of the
;                       hyperbolic model are computed.  Y = 1./(a0 + a1*x)
;
;          LOGISTIC:    If set to a non-zero value, the parameters of the
;                       logistic model are computed.  Y = 1./(a0 * a1^x + a2)
;
;         LOGSQUARE:    If set to a non-zero value, the parameters of the
;                       logsquare model are computed.
;			Y = a0 + a1*alog10(x) + a2 * alog10(x)^2
;
;             SIGMA:    Use this keyword to specify a named variable which 
;                       returns a vector of standard deviations for the computed
;                       model parameters.
;                
;           WEIGHTS:    An n-element vector of weights. If the WEIGHTS vector 
;                       is not specified, an n-element vector of 1.0s is used.
;
;              YFIT:    Use this keyword to specify a named variable which 
;                       returns n-element vector of y-data corresponding to the 
;                       computed model parameters.
;
; EXAMPLE:
;       Define two n-element vectors of paired data.
;         x = [2.27, 15.01, 34.74, 36.01, 43.65, 50.02, 53.84, 58.30, 62.12, $
;             64.66, 71.66, 79.94, 85.67, 114.95]
;         y = [5.16, 22.63, 34.36, 34.92, 37.98, 40.22, 41.46, 42.81, 43.91, $
;             44.62, 46.44, 48.43, 49.70, 55.31]
;       Define a 3-element vector of initial estimates for the logsquare model.
;         a = [1.5, 1.5, 1.5]
;       Compute the model parameters of the logsquare model, a(0), a(1), & a(2).
;         result = comfit(x, y, a, sigma = sigma, yfit = yfit, /logsquare)
;       The result should be the 3-element vector:
;         [1.42494, 7.21900, 9.18794]
;
; REFERENCE:
;       APPLIED STATISTICS (third edition)
;       J. Neter, W. Wasserman, G.A. Whitmore
;       ISBN 0-205-10328-6
;
; MODIFICATION HISTORY:
;       Written by:  GGS, RSI, September 1994
;-

pro exp_func, x, a, f, pder
  f = a[0] * a[1]^x + a[2]
  if n_params() ge 4 then $
    pder = [[a[1]^x], [a[0] * x * a[1]^(x-1)], [replicate(1., n_elements(x))]]
end

pro geo_func, x, a, f, pder
  f = a[0] * x^a[1] + a[2]
  if n_params() ge 4 then $
    pder = [[x^a[1]], [a[0] * alog(x) * x^a[1]], [replicate(1., n_elements(x))]]
end

pro gom_func, x, a, f, pder
  f = a[0] * a[1]^(a[2]*x) + a[3]
  if n_params() ge 4 then $
    pder = [[a[1]^(a[2]*x)], [a[0] * a[2] * x * a[1]^(a[2]*x-1)], $
            [a[0] * x * alog(a[1]) * a[1]^(a[2]*x)], [replicate(1., n_elements(x))]]
end

pro hyp_func, x, a, f, pder
  f = 1.0 / (a[0] + a[1] * x)
  if n_params() ge 4 then $
    pder = [[-1.0 / (a[0] + a[1] * x)^2], [-x / (a[0] + a[1] * x)^2]]
end

pro log_func, x, a, f, pder
  f = 1.0 / (a[0] * a[1]^x + a[2])
  if n_params() ge 4 then begin
      denom =  -1.0/(a[0] * a[1]^x + a[2])^2
      pder = [[a[1]^x*denom], [a[0] * x * a[1]^(x-1)*denom], [denom]]
  endif
end

pro logsq_func, x, a, f, pder
  b = alog10(x)
  b2 = b^2
  f = a[0] + a[1] * b + a[2] * b2
  if n_params() ge 4 then $
    pder = [[replicate(1., n_elements(x))], [b], [b2]]
end

;; added 4/24/00 Erin Scott Sheldon
PRO powerlaw_func, x, a, f, pder
  f = a[0]*x^a[1]
  IF n_params() GE 4 THEN BEGIN
      pder = [ [x^a[1]], [a[0]*alog(x) * x^a[1]] ]
  ENDIF
END 

;; added 5/18/00 Erin Scott Sheldon
PRO sigdiffsis_trunc_func, x, a, f, pder


  ;; a[0] = sigma_v in units of 1000. km/s
  ;; a[1] = outer scale in units of Mpc.
  ;; x should be in Mpc
  ;; y should be in units of 1.15e14 Msolar/Mpc^2

  ;; core in Mpc
  core = 1./1000.

  term1 =  1./sqrt(core + x^2)
  term2 =  2.*a[1]*term1^2
  term3 = -1./sqrt( x^2 + a[1]^2 )
  term4 =  2.*a[1]^2*term1^2*term3
  termsum_f = term1 + term2 + term3 + term4 

  f = a[0]^2*termsum_f 

;  print,'Here is x: ',x
;  print,'Here is f: ',f

  IF n_params() GE 4 THEN BEGIN 

;      print,'.',format='(a,$)'

      term1 = 2./(core + x^2)
      term2 = -4.*a[1]/(core + x^2)/sqrt(x^2 + a[1]^2)
      term3 = 2.*a[1]^2*sqrt(x^2 + a[1]^2 )/(core + x^2)
      term4 = a[1]*sqrt(x^2 + a[1]^2 )
      termsum_pder = term1 + term2 + term3 + term4

      pder = [ [ 2.* a[0]^2 *termsum_f     ], $
               [ 2.* a[0] *termsum_pder   ] ]

  ENDIF 
return
END 

PRO lognormal_func, xin, a, f

  x = double(abs(xin))

  nx = n_elements(x)

  amp = abs(a[0])
  alpha = abs(a[1])
  sig = a[2]
  ;;alpha = a[2]
  ;;sig = a[3]

  w = where(x GT alpha,nw)
  F = fltarr(nx)

  IF nw EQ 0 THEN BEGIN
      message,'No good x values',/inf
      return
  ENDIF 

  F[w] = amp*1./sqrt(2.*!pi*sig^2)*exp(-(alog(x[w]-alpha))^2/2./sig^2 )/x[w]

return
END

function comfit2, x, y, a, weights = weights, sigma = sigma, yfit = yfit, $
                  exponential = exponential, geometric = geometric, $
                  gompertz = gompertz, hyperbolic = hyperbolic, $
                  logistic = logistic, logsquare = logsquare, $
                  powerlaw=powerlaw, sigdiff=sigdiff, $
                  lognormal=lognormal, $
                  chisq = chisq, iter=iter, itmax=itmax, tol=tol, $
                  noderivative=noderivative
  
  on_error, 2

  fcn_names = ['exp_func', 'geo_func', 'gom_func', 'hyp_func', 'log_func', $
               'logsq_func', 'powerlaw_func', 'sigdiffsis_trunc_func','lognormal_func']
  fcn_npar = [ 3, 3, 4, 2, 3, 3, 2, 2, 3]

  nx = n_elements(x)
  wx = n_elements(weights)
  if nx ne n_elements(y) then $
    message, 'x and y must be vectors of equal length.'
  
  if wx eq 0 then weights = replicate(1.0, nx) $
  else if wx ne nx then $
    message, 'x and weights must be vectors of equal length.'
  
  a1 = A			;Copy initial guess
  
  if keyword_set(exponential) then i = 0 $
  else if keyword_set(geometric) then i = 1 $
  else if keyword_set(gompertz) then i = 2 $
  else if keyword_set(hyperbolic) then i = 3 $
  else if keyword_set(logistic) then i = 4 $
  else if keyword_set(logsquare) then i = 5 $
  ELSE IF keyword_set(powerlaw) THEN i = 6 $
  ELSE IF keyword_set(sigdiff) THEN i = 7 $
  ELSE IF keyword_set(lognormal) THEN i=8 $
  else message, 'Type of model must be supplied as a keyword parameter.'
  
  if n_elements(a1) ne fcn_npar[i] then $
    message, 'A must be supplied as a '+strtrim(fcn_npar[i], 2) + $
    '-element initial guess vector.'
  yfit = curvefit(x, y, weights, a1, sigma, function_name = fcn_names[i], chisq=chisq, $
                  iter=iter, itmax=itmax, tol=tol, noderivative=noderivative)

  return, a1
end
