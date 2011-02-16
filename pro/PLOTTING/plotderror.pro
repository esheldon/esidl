PRO plotderror, x, y, xlow, xhigh, ylow, yhigh, $
                NOHAT=hat, HATLENGTH=hln, ERRTHICK=eth, $
                ERRSTYLE=est, TYPE=itype, XRANGE = xrange, $
                XLOG=xlog, YLOG=ylog, $
                NSKIP = nskip, NOCLIP = noclip, $
                ERRCOLOR = ecol, YRANGE = yrange, $
                NSUM = nsum, _EXTRA = pkey, ANONYMOUS_ = DUMMY_
  
 On_error, 2
 np = N_params()
 IF (np LT 6) THEN BEGIN
        print, "PLOTDERROR must be called with 6  parameters."
        print, "Syntax: plotderror, x, y, xlow, xhigh, ylow, yhigh"
        RETURN
 ENDIF

; Error bar keywords (except for HATLENGTH; this one will be taken care of 
; later, when it is time to deal with the error bar hats).

 IF (keyword_set(hat)) THEN hat = 0 ELSE hat = 1
 if not keyword_set( THICK ) then thick = !P.THICK
 if (n_elements(eth) EQ 0) THEN eth = thick
 IF (n_elements(est) EQ 0) THEN est = 0
 IF (n_elements(ecol) EQ 0) THEN ecol = !P.COLOR
 if N_elements( NOCLIP ) EQ 0 then noclip = 0
 if not keyword_set(NSKIP) then nskip = 1

;				Other keywords.

 IF (keyword_set(itype)) THEN BEGIN
	CASE (itype) OF
		   1 :  ylog = 1	; X linear, Y log
		   2 :  xlog = 1	; X log, Y linear
		   3 :  BEGIN		; X log, Y log
			xlog = 1
			ylog = 1
			END
		ELSE : 
	ENDCASE
 ENDIF
 if not keyword_set(XLOG) then xlog = 0
 if not keyword_set(YLOG) then ylog = 0
;			If no x array has been supplied, create one.  Make
;			sure the rest of the procedure can know which parameter
;			is which.


;			Determine the number of points being plotted.  This
;			is the size of the smallest of the three arrays
;			passed to the procedure.  Truncate any overlong arrays.

 xx=x
 yy=y

 n = N_elements(xx) < N_elements(yy)

 IF n LT 2 THEN BEGIN 
     ;message,'Not enough points to plot.'
     xx = [xx]
     yy = [yy]
 ENDIF 

 xx = xx[0:n-1]
 yy = yy[0:n-1]
 xlo = xlow[0:n-1]
 ylo = ylow[0:n-1]
 xhi = xhigh[0:n-1]
 yhi = yhigh[0:n-1]

; If no y-range was passed via keyword or system variable, force one large 
; enough to display all the data and the entire error bars.     If NSKIP is
; set, then need consider only the error bars actually plotted.
; If a reversed y-range was passed, switch ylo and yhi.

 if not keyword_set( YRANGE ) then yrange = !Y.RANGE
 IF yrange[0] EQ yrange[1] THEN BEGIN
	if keyword_set( XRANGE ) then  begin
		good = where( (xx GT min(xrange)) and (xx LT max(xrange)) )
		yrange = [min(ylo[good]),max(yhi[good])]
	endif else yrange = [min(ylo), max(yhi)]
 ENDIF ;ELSE IF yrange[0] GT yrange[1] THEN BEGIN
;	ylo = yy + yerror*ierr
	;yhi = yy - yerror*ierr
; ENDIF

;        Similarly for x-range

 if not keyword_set( XRANGE ) then xrange = !X.RANGE

 IF xrange[0] EQ xrange[1] THEN xrange = [min(xlo), max(xhi)]


;			Plot the positions.

 plot, xx, yy, XRANGE = xrange, YRANGE = yrange, XLOG = xlog, YLOG = ylog, $
   _EXTRA = pkey, NOCLIP = noclip

;	Plot the error bars.   Compute the hat length in device coordinates
;       so that it remains fixed even when doing logarithmic plots.

 data_low = convert_coord(xx,ylo,/TO_DEVICE)
 data_hi = convert_coord(xx,yhi,/TO_DEVICE)
 x_low = convert_coord(xlo,yy,/TO_DEVICE)
 x_hi = convert_coord(xhi,yy,/TO_DEVICE)
 
 ycrange = !Y.CRANGE   &  xcrange = !X.CRANGE
 if ylog EQ 1 then ylo = ylo > 10^ycrange[0]
 if (xlog EQ 1) then xlo = xlo > 10^xcrange[0]
 sv_psym = !P.PSYM & !P.PSYM = 0
 
 FOR i = 0L, (n-1), Nskip DO BEGIN     
     
     plots, [xx[i],xx[i]], [ylo[i],yhi[i]], LINESTYLE=est,THICK=eth,  $
       NOCLIP = noclip, COLOR = ecol
     
;                                                         Plot X-error bars 
     plots, [xlo[i],xhi[i]],[yy[i],yy[i]],LINESTYLE=est, $
       THICK=eth, COLOR = ecol, NOCLIP = noclip
     IF (hat NE 0) THEN BEGIN
         IF (N_elements(hln) EQ 0) THEN hln = !D.X_VSIZE/100. 
         exx1 = data_low[0,i] - hln/2.
         exx2 = exx1 + hln
         plots, [exx1,exx2], [data_low[1,i],data_low[1,i]],COLOR=ecol, $
           LINESTYLE=est,THICK=eth,/DEVICE, noclip = noclip
         plots, [exx1,exx2], [data_hi[1,i],data_hi[1,i]], COLOR = ecol,$
           LINESTYLE=est,THICK=eth,/DEVICE, noclip = noclip
         
;                                                        Plot Y-error bars
         
         IF (N_elements(hln) EQ 0) THEN hln = !D.Y_VSIZE/100.
         eyy1 = x_low[1,i] - hln/2.
         eyy2 = eyy1 + hln
         plots, [x_low[0,i],x_low[0,i]], [eyy1,eyy2],COLOR = ecol, $
           LINESTYLE=est,THICK=eth,/DEVICE, NOCLIP = noclip
         plots, [x_hi[0,i],x_hi[0,i]], [eyy1,eyy2],COLOR = ecol, $
           LINESTYLE=est,THICK=eth,/DEVICE, NOCLIP = noclip
     ENDIF
     NOPLOT:
 ENDFOR
 !P.PSYM = sv_psym
;
 RETURN
 END
