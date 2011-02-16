PRO shearfit, files, shguess, yrange=yrange, names=names, str=str, shfit=shfit

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    SHEARPFIT
;       
; PURPOSE:
;    Fits up to three different bandpasses to a model.  
;
; CALLING SEQUENCE:
;    
;
; INPUTS: 
;    
;
; OPTIONAL INPUTS:
;    
;
; KEYWORD PARAMETERS:
;    
;       
; OUTPUTS: 
;    
;
; OPTIONAL OUTPUTS:
;    
;
; CALLED ROUTINES:
;    
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;  
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() LT 1 THEN BEGIN 
     print,'-Syntax: shearfit, files [, shguess, yrange=yrange, names=names, str=str, shfit=shfit, _extra=extra]'
     print,''
     print,'shguess~'
     print,'Use doc_library,""  for more help.'  
     return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nfile = n_elements(files)
  IF nfile gt 3 THEN BEGIN
      print,'At most 3 bandpasses'
      return
  ENDIF 

  IF n_elements(names) EQ 0 THEN names = strarr(nfile)
  IF n_elements(shguess) EQ 0 THEN shguess=[.05, -.7]

  ;; number of elements in strings
  nstring = 7
  ;; number from last x position to print xyouts
  xn = 4
  ;; Shear polarizeability correction
  Ssh = 1.13

  ytitle='Tangental Shear'
  xtitle = 'Projected Radius (arcsec)'

  pold=!p.multi
  !p.multi=0

  erase & multiplot, [1, nfile]

  FOR i=0, nfile-1 DO BEGIN 
 
      rdobjshear, files[i], str,/silent
      
      str.shear = str.shear*Ssh
      n=n_elements(str.meanr)
      str = str[0:n-2]          ;last bin is usually bad

      IF i EQ 0 THEN yrange = prange(str.shear, str.shearerr)
      
      fitpower, str.meanr, str.shear, str.shearerr, shguess,$
        yfit, a, asigma, P, convergence

      ;; a[1] (the power) will be negative.   Convert to positive
      ;; since we fit model a[0]*r^-beta
      a[1] = -a[1]

      CASE i OF
          nfile-1: ploterr,str.meanr,str.shear,str.shearerr,psym=1, $
                           xtitle=xtitle,ytitle=ytitle,yrange=yrange
          0:       ploterr, str.meanr, str.shearerr,str.shearerr,psym=1, $
                            yrange=yrange, _extra=extra
          ELSE:    ploterr, str.meanr,str.shearerr,str.shearerr,psym=1,$
                            yrange=yrange
      ENDCASE 

      ;; Overplot the model fits
      oplot, str.meanr, yfit
      oplot,[0,10000],[0,0]

      ;; annotate
      x=str[n-xn-1].meanr
      y1=str[0].shear
      y2 = .8*y1
      y3 = 1.2*y1
      xyouts,x,y1,'a = '+ntostr(a[0],nstring) + ' +/- '+$
             ntostr(asigma[0],nstring)
      xyouts,x,y2,'b = '+ntostr(a[1],nstring) + ' +/- '+$
             ntostr(asigma[1],nstring)
      xyouts,x,y3,names[i]
      IF i NE nfile-1 THEN multiplot
  ENDFOR 
  multiplot, /reset
  
  shfit = yfit
  meanr = str.meanr
  print
  print,'This used MULTIPLOT, you may need to use ERASE'
  !p.multi=pold

  return 
END 
