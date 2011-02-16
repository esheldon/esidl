PRO aplotderror, aspect, x, y, xlow, xhigh, ylow, yhigh, _extra=extra, center=center

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    APLOTDERROR
;       
; PURPOSE:
;    Wrapper for PLOTDERR that forces a user defined aspect ratio.
;
; CALLING SEQUENCE:
;    aplotderr, aspect, [x,] y, xlow, xhigh, ylow, yhigh, 
;        center=center, _extra=extra
;
; INPUTS: 
;    aspect: xsize/ysize
;    x: x values
;    xlow, xhigh: high and low x values
;    y: y values
;    yerr: y error values.
;
;
; KEYWORD PARAMETERS:
;    /center: if set then center up the display.
;    _extra:  plotting keywords.
;       
; CALLED ROUTINES:
;    PLOTERR
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    Author: Erin Scott Sheldon  UofMich 11/17/99  
;       
;
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  On_error,2

  np=n_params()
  IF (np LT 7) THEN BEGIN 
     print,'-Syntax: aplotderror, aspect, x, y, xlow, xhigh, ylow, yhigh, _extra=extra, center=center'
     print,''
     print,'Must be called with 7 parameters'
     print,' aspect = xsize/ysize'
     print,'Use doc_library,"aploterr"  for more help.'  
     return
  ENDIF 

  IF !p.multi[1] EQ 0 THEN !p.multi[1] = 1
  IF !p.multi[2] EQ 0 THEN !p.multi[2] = 1

  plot, [0,1], [0,1], /nodata, xstyle=4, ystyle=4
  px = !x.window*!d.x_vsize
  py = !y.window*!d.y_vsize

  xsize = px[1] - px[0]
  ysize = py[1] - py[0]

  a0 = xsize/ysize
  CASE 1 OF
      (aspect LE a0): xsize = xsize*(aspect/a0) ;shrink xsize
      (aspect GT a0): ysize = ysize*(a0/aspect) ;shrink ysize
  ENDCASE 

  px[1] = px[0] + xsize
  py[1] = py[0] + ysize

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; center up the display 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  IF keyword_set(center) THEN dcenter, xsize, ysize, px, py

  position = [ [px(0), py(0)], [px(1), py(1)] ]

  plotderror, x, y, xlow, xhigh, ylow, yhigh, position=position, $
    /device, _extra=extra

return
END 
