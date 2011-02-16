pro  sis_deflect, zlens, zsource, vdisp, alpha, reduced_alpha, h=h, omegamat=omegamat, silent=silent

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; NAME: sis_deflect
;       
; PURPOSE: calculate the deflection angle for a singular isothermal sphere
;	
;
; CALLING SEQUENCE:
;      
;                 
;
; INPUTS: 
;       
; OUTPUTS: 
;
; OPTIONAL OUTPUT ARRAYS:
;
; INPUT KEYWORD PARAMETERS:
; 
; PROCEDURE: The procedure below calculates it exactly, but the angle 
;     is nearly equivalent to the following simple formula:
;
;     alpha = (1.4'')*(vdisp/220 km/s)^2 
;      
;     alpha_reduced = (Dls/Dos)*alpha
;	
;
; REVISION HISTORY:
;	
;       
;                                      
;                                        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if N_params() eq 0 then begin
	print,'Syntax: sis_deflect, zlens, zsource, vdisp, alpha, reduced_alpha, h=h, omegamat=omegamat, silent=silent'
        print,'vdisp: velocity dispersion in km/s'
	return
  endif

  ppi = 3.14159
  c=3.0e5  ;km/sec
  if (not keyword_set(omegamat) ) then omegamat = 0.3
  if (not keyword_set(h) ) then h=0.7

  zdist, zlens, dang=dol, h=h, omegamat=omegamat,silent=1
  zdist, zsource, dang=dos, h=h, omegamat=omegamat,silent=1

  ;;;  This is wrong!! But I don't know the right answer yet.  Good
  ;;;  within a factor of 2.
  dls = dos - dol
  alpha = 4*ppi*vdisp^2/c^2

  ;;;;;;  convert to arcseconds  Use 2.8e-4 degrees/arcsec
  alpha = alpha*180.0/ppi/2.8e-4  
  alpha_reduced = alpha*dls/dos

  if (not keyword_set(silent) ) then begin
        print,'-----------------------------------'
        print,'Dol = ',strtrim(string(dol),2),' Mpc'
        print,'Dos = ',strtrim(string(dos),2),' Mpc'
        print,'Dls = ',strtrim(string(dls),2),' Mpc'
	print,'Alpha: ', strtrim(string(alpha),2), ' arcsec'
        print,'Reduced Alpha: ', strtrim(string(alpha_reduced),2), ' arcsec'
        print,'-----------------------------------'
  endif

return
end 
  



