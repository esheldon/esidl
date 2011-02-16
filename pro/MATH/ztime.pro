pro  ztime, z, time, h=h, omega=omega, silent=silent

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; NAME: zdist
;       
; PURPOSE: calculate time as a function of redshift in the matter dominated
;          era.
;	
;
; CALLING SEQUENCE: zdist, z, time, h=h, omega=omega,silent=silent
;      
; INPUTS:  z: redshift
;
; OPTIONAL INPUTS: h: hubble parameter
;                  omega: density parameter
;       
; OUTPUTS: dist:
;
; REVISION HISTORY: Erin Scott Sheldon 2/24/99
;	
;       
;                                      
;                                        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if N_params() eq 0 then begin
	print,'-Syntax: ztime, z, time, h=h, omega=omega, silent=silent'
	print,'Returns time in seconds.  Matter dominated only, z < 1100'
	return
  endif
;;;;;;;;;;;;;;;;;;  check keywords  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (not keyword_set(h)) then h=0.7

if (not keyword_set(omega)) then omega=0.3

;;;;;;;; some parameters  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

c = 3.0e10           ;;  speed of light in cm/s
c = c/3.1e24        ;;  speed of light in Mpc/s

;;;;;;;; calculate time  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

time = 3000.0/sqrt(omega*h^2)/(1.0 + z)^(3.0/2.0)/3.0/c
years = time/3.1536e7

if (not keyword_set(silent) ) then begin
  print,'Using h = ',strtrim(string(h),2),'  omega  = ',strtrim(string(omega),2)
  print,''
  print,'Time in seconds:',strtrim(string(time),2)
  print,'Time in years:',strtrim(string(years),2)
endif

return
end
	
















