pro voronoi_density,x,y,d,plot=plot,rect=rect,_extra=ex,b=b,trueaspect=trueaspect, xstyle=xstyle, ystyle=ystyle
;+
;create a voronoi tesselation and then the density
;by inverting the area of each voronoi polygon
;if plot is set it will plot the tesselation
;rect is the bounding rectangle [xmin,ymin,xmax,ymax]
;if rect not present it will make a reasonable one containing 
;all the points
;extra keywords are handed to plot for ploting the points
;make sure there are no identicle points in the list
;
; Use trueaspect keyword to force plots to have correct aspect ratio.
;-

if n_params() eq 0 then begin
	print,'-syntax voronoi_density,x,y,d,plot=plot,rect=rect,_extra=ex,b=b,trueaspect=trueaspect'
	return
endif

triangulate, x, y, tr, b ,CONN=c                   ;Triangulate it
if n_elements(rect) eq 0 then rect = 0
IF NOT keyword_set(trueaspect) THEN trueaspect=0
IF n_elements(xstyle) EQ 0 THEN xstyle=0
IF n_elements(ystyle) EQ 0 THEN ystyle=0

IF NOT (xstyle MOD 2) THEN addx = 1 ELSE addx = 0
IF NOT (ystyle MOD 2) THEN addy = 1 ELSE addy = 0

n=n_elements(x)
area=fltarr(n)

if keyword_set(plot) then BEGIN

    IF trueaspect THEN BEGIN 
        sx=max(x)-min(x)
        sy=max(y)-min(y)

        aspect=sx/sy
        aplot,aspect,x,y,psym=3,_extra=ex, $
              yrange=[min(y),max(y)],xrange=[min(x),max(x)], $
              ystyle=ystyle+addy,xstyle=xstyle+addx
    ENDIF ELSE BEGIN 
        plot,x,y,psym=3,_extra=ex, $
              yrange=[min(y),max(y)],xrange=[min(x),max(x)], $
              ystyle=ystyle+addy,xstyle=xstyle+addx
    ENDELSE 
    for i=0l, n-1 do begin
        voronoi, x, y, i, c, xp, yp, rect ;Get the ith polygon
        oplot,[xp,xp(0)],[yp,yp(0)] ;Draw it
        area(i)=poly_area(xp,yp)
    endfor
endif else begin
    for i=0l, n-1 do begin
        voronoi, x, y, i, c, xp, yp, rect ;Get the ith polygon
        area(i)=poly_area(xp,yp)
    endfor
endelse

d=1.0/area

return
end








