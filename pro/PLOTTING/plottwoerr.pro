PRO plottwoerr, x1, y1, y1err, x2, y2, y2err, _extra=extra, center=center, $
                padding=padding, frac1=frac1, $
                aspect=aspect, topaspect=topaspect, $
                xtitle=xtitle, topytitle=topytitle, botytitle=botytitle,$
                xrange=xrange, topyrange=topyrange, botyrange=botyrange,$
                topystyle=topystyle, botystyle=botystyle, $
                topylog=topylog,$
                botylog=botylog,$
                botxminor=botxminor,botyminor=botyminor,$
                topxminor=topxminor,topyminor=topyminor,$
                toppsym=toppsym, botpsym=botpsym,$
                toplinestyle=toplinestyle, botlinestyle=botlinestyle, $
                xoplot=xoplot, $
                yoplot=yoplot, $
                oploterr=oploterr,$
                oplotcolor=oplotcolor,$
                oplotsym=oplotsym, $
                oplotlinestyle=oplotlinestyle, $
                xticklen=xticklen,$
                position1=position1, position2=position2,$
                xwindow1=xwindow1, ywindow1=ywindow1,$
                xwindow2=xwindow2, ywindow2=ywindow2,$
                xtickformat=xtickformat, $
                ytickformat1=ytickformat1, $
                ytickformat2=ytickformat2, $
                title=title

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    APLOT
;       
; PURPOSE:
;    Wrapper for PLOT that forces a user defined aspect ratio.
;
; CALLING SEQUENCE:
;    aplot, aspect, [x,] y, center=center, _extra=extra
;
; INPUTS: 
;    aspect: xsize/ysize
;    y: y values
;
; OPTIONAL INPUTS:
;    x: optional x values
;
; KEYWORD PARAMETERS:
;    /center: if set then center up the display.
;    _extra:  plotting keywords.
;       
; CALLED ROUTINES:
;    PLOT
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    Author:  Erin Scott Sheldon  UofMich 11/17/99  
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
;  On_error,2

  np = n_params()
;  IF (np LT 2) OR (np GT 3) THEN BEGIN 
  IF np LT 6 THEN BEGIN 
     print,'-Syntax:plottwoerr, x1, y1, y1err, x2, y2, y2err, $'
     print,'    _extra=extra, center=center, $'
     print,'    padding=padding, frac1=frac1, $'
     print,'    aspect=aspect, topaspect=topaspect, $'
     print,'    xtitle=xtitle, topytitle=topytitle, botytitle=botytitle,$'
     print,'    title=title'
     print,''
     print,' aspect = xsize/ysize'
     print,'Use doc_library,"aplot"  for more help.'  
     return
  ENDIF 

  IF !p.multi[1] EQ 0 THEN !p.multi[1] = 1
  IF !p.multi[2] EQ 0 THEN !p.multi[2] = 1

  ;IF n_elements(xrange) EQ 0 THEN xrange=[0,0]
  ;IF n_elements(botyrange) EQ 0 THEN botyrange=[0,0]
  ;IF n_elements(topyrange) EQ 0 THEN topyrange=[0,0]

  IF n_elements(xrange) NE 0 THEN BEGIN 
      botxrange = xrange
      topxrange = xrange
  ENDIF 

  IF (n_elements(aspect) NE 0) AND (n_elements(topaspect) NE 0) THEN $
    message,'Do not send both aspect and topaspect'
  IF n_elements(padding) EQ 0 THEN padding = 0.02
  IF n_elements(frac1) EQ 0 THEN frac1 = 0.5

  plot, [0,1], [0,1], /nodata, xstyle=4, ystyle=4
  px = !x.window*!d.x_vsize
  py = !y.window*!d.y_vsize
  
  xsize = px[1] - px[0]
  ysize = py[1] - py[0]

  ;; aratio of plot region
  a0 = xsize/ysize
  
  ;; user wants different aspect ratio for plot region
  IF n_elements(aspect) NE 0 THEN BEGIN
      CASE 1 OF
          (aspect LE a0): xsize = xsize*(aspect/a0) ;shrink xsize
          (aspect GT a0): ysize = ysize*(a0/aspect) ;shrink ysize
      ENDCASE
  ENDIF 

  px[1] = px[0] + xsize
  py[1] = py[0] + ysize

  f1 = frac1
  f2 = 1.-frac1

  ;; top plot
  xsize1 = xsize
  ysize1 = ysize*f1

  ;; bottom plot
  xsize2 = xsize1
  ysize2 = ysize*f2

  px1=px & py1=py
  ysize1 = ysize1-padding*!d.y_vsize/2.
  py1[0] = py[1] - ysize1

  px2=px & py2=py
  ysize2 = ysize2-padding*!d.y_vsize/2.
  py2[1] = py[0] + ysize2

  IF n_elements(topaspect) NE 0 THEN BEGIN 

      a1=xsize1/ysize1
      CASE 1 OF 
          (topaspect LE a1): BEGIN ;; shrink x
              txsize = xsize1*(topaspect/a1)
              diff = xsize1-txsize
              px1[1] = px1[1] - diff
              px2[1] = px2[1] - diff
              xsize1=txsize
              xsize2=txsize
          END 
          (topaspect GT a1): BEGIN ;; shrink y
              tysize = ysize1*(a1/topaspect)
              diff = ysize1 - tysize
              py1[1] = py1[1] - diff
              ysize1=tysize
          END 
      ENDCASE 
  ENDIF 

  IF keyword_set(center) THEN BEGIN 
      tpx = px1
      tpy = [py2[0], py1[1]]
      dcenter, xsize1, ysize1+ysize2+padding*!d.y_vsize, tpx, tpy

      px1 = tpx
      px2 = tpx
      
      ydiff = tpy[1] - py1[1]
      py1 = py1+ydiff
      py2 = py2+ydiff

  ENDIF 

  position1 = [ px1[0], py1[0], px1[1], py1[1] ]
  position2 = [ px2[0], py2[0], px2[1], py2[1] ]

  ploterror, x1, y1, y1err, position=position1, /device, /noerase, _extra=extra, $
             xrange=topxrange, yrange=topyrange, $
             ylog=topylog, psym=toppsym, line=toplinestyle, $
             xtickname=replicate(' ', 30), title=title, ytitle=topytitle,$
             xticklen=xticklen,xminor=topxminor,yminor=topyminor, $
             ytickf=ytickformat1, ystyle=topystyle

  ;; are we overplotting anything?
  IF (n_elements(xoplot) NE 0) AND (n_elements(yoplot) NE 0) THEN BEGIN 
      IF n_elements(oploterr) NE 0 THEN BEGIN 
          oploterror, xoplot, yoplot, oploterr, psym=oplotsym, linestyle=oplotlinestyle, color=oplotcolor, errc=oplotcolor
      ENDIF ELSE BEGIN
          oplot, xoplot, yoplot, psym=oplotsym, linestyle=oplotlinestyle, color=oplotcolor
      ENDELSE 
  ENDIF 

  ;; want x ticklenghts to be same absolute size
  IF n_elements(xticklen) NE 0 THEN BEGIN
      xticklen2=xticklen*(ysize1/ysize2)
  ENDIF ELSE xticklen2=!x.ticklen*(f1/f2)
  ploterror, x2, y2, y2err, position=position2, /device, /noerase, _extra=extra, $
             xtitle=xtitle, ytitle=botytitle, $
             ylog=botylog, psym=botpsym, linestyle=botlinestyle, $
             xrange=botxrange, yrange=botyrange, xticklen=xticklen2,$
             xminor=botxminor,yminor=botyminor, $
             xtickformat=xtickformat, ytickf=ytickformat2, ystyle=botystyle

  ;; return normalized coords
  position1 = [ px1[0]/!d.x_vsize, py1[0]/!d.y_vsize, $
                px1[1]/!d.x_vsize, py1[1]/!d.y_vsize ]
  position2 = [ px2[0]/!d.x_vsize, py2[0]/!d.y_vsize, $
                px2[1]/!d.x_vsize, py2[1]/!d.y_vsize ]

  xwindow1 = [ position1[0], position1[2] ]
  xwindow2 = [ position2[0], position2[2] ]

  ywindow1 = [ position1[1], position1[3] ]
  ywindow2 = [ position2[1], position2[3] ]

return

END 
