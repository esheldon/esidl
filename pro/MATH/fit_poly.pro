pro  fit_poly, x, f, order, index=index,fitfunc=fitfunc,coeff=coeff,_extra=extra

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; NAME:  fit_poly
;       
; PURPOSE: wrapper for svdfit.  Just does polynomials.
;	
;
; CALLING SEQUENCE: pro  fit_poly, x, f, order, index=index, ylog=ylog,
;	fitfunc=fitfunc,coeff=coeff
;      
;                 
;
; INPUTS: x: abcissa  f: the function
;         order: the order of polynomial fit
;       
; OUTPUTS: plots the fitted function and original function together
;
; OPTIONAL OUTPUT ARRAYS: fitfunc:  array containing values of fitted func
;       coeff:  array containing the coefficients
;
; INPUT KEYWORD PARAMETERS:
; 
; PROCEDURE: 
;	
;	
;
; REVISION HISTORY:
;	
;       
;                                      
;                                        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if N_params() eq 0 then begin
	print,'-Syntax: fit_poly, x, f, order, index=index,fitfunc=fitfunc,coeff=coeff, _extra=extra'
	return
  endif

  n=n_elements(x)
  if (keyword_set(index)) then begin
	print,'Fitting over domain: ',x[min(index)],x[max(index)]
  endif else begin
	index=indgen(n)
  endelse

  fitfunc=fltarr(n)

  coeff=svdfit(x(index),f(index),order+1,chisq=chisq)
  print,'Chi Squared: ',chisq
  print,'Polynomial coefficients.'
  print,coeff

  for i=0,order do begin

    fitfunc = fitfunc + coeff(i)*x^i

  endfor

  plot,x,f,psym=1, _extra=extra
  oplot,x,fitfunc
  legend,['Input','Fitted'],psym=[1,0]

  return
  end




















