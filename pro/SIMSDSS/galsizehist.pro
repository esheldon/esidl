pro  galsizehist, yhist, xhist, plt=plt

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; NAME: galsizehist
;       
; PURPOSE: generate the distribution of galaxies as a function of size
;	
;
; CALLING SEQUENCE: galsizehist,yhist, xhist, plt=plt
;
; OUTPUTS: yhist: histogram, number of galaxies per frame vs. size
;          xhist:  the abscissae
;
;  size = sqrt( ixx + iyy )
;
; INPUT KEYWORD PARAMETERS:  plt: if plt is set a plot of rhist is shown.
;
; REVISION HISTORY:  Erin Scott Sheldon  2/20/99 
;	
;                                        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if N_params() eq 0 then begin
    print,'Syntax: galsizehist, yhist, xhist, plt=plt'
    return
  endif

  file = '~/idl.lib/SIMSDSS/FIT/galsize.fit'
  galsize = mrdfits(file,1,hdr,/silent)

  xhist = galsize.xhist
  yhist = galsize.yhist


  IF keyword_set(plt) THEN plot,xhist,yhist,psym=4,xtitle='Sqrt( ixx + iyy )',$
               ytitle='Number',title='Galaxy Size Histogram'

return
end










