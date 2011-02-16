PRO circ_radec, photostr, index, objx, objy, ra, dec, radius, $
                color=color,$
                linestyle=linestyle,$
                box=box, $
                cross=cross, holeradius=holeradius, $
                clr=clr, silent=silent, $
                pos_ang=pos_ang, rmax=rmax, rmin=rmin, $
                order=order, ny=ny, xshift=xshift, yshift=yshift

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: 
;    circ_radec
;    THIS ROUTINE IS OBSOLETE.  USE THE MARK_RADEC PROCEDURE.
;       
; PURPOSE: 
;    Circles an ra and dec position on the current display.  We need 
;    a photo structure, an index and the position of the object with that
;    index inorder to complete the mapping.	
;
; CALLING SEQUENCE: 
;    circ_radec, photostr, index, objx, objy, ra, dec, 
;    nodisplay=nodisplay, radius=radius,color=color,box=box
;      
;                 
;
; INPUTS: photostr: a photo structure
;         index: the index of an object in that structure
;         objx, objy: the x and y positions of index object.
;         ra, dec:  the ra and dec to be circled.
;
; INPUT KEYWORD PARAMETERS:
;         radius=: radius of the circle. 
;         pos_ang=, rmax=, rmin=: instead of just circle, overplot an ellipse
;         clr=: color index to make the mapping from.
;
;         order=: order=0 (default) means images was plotted 0,0 in bottom
;                 left,  order=1 is top left.  Note objx,objy
;                 is still given in the original coord. system. User must input
;                 ny for the transformation.
;         ny=: number of pixels in the y-direction.
;         /box: use a box instead of a circle.  radius will be side of box
;         /silent: Shut off the messages except errors.
;
;       
; OUTPUTS: none
;
; CALLED ROUTINES: 
;                  CONVERT2XY
;                  POLYWARP
;                  KMAP
;                  (RDIS_SETUP)
;                  (RDIS)
;                  TVCIRCLE
;                  (TVBOX)
;
; 
; PROCEDURE: 
;	convert the ra and dec to x y coordinates
;	
;
; REVISION HISTORY:
;	Author: Erin Scott Sheldon Umich 5/25/99
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() LT 7 THEN BEGIN 
      print,'-Syntax: circ_radec, photostr, index, objx, objy, ra, dec, radius, color=color, box=box, clr=clr, ny=ny, order=, /silent]'
      print,''
      print,'Use doc_library,"circ_radec"  for more help.'  
      return
  ENDIF  

  IF NOT keyword_set(silent) THEN silent = 0

  IF n_elements(xshift) EQ 0 THEN xshift=0.0
  IF n_elements(yshift) EQ 0 THEN yshift=0.0

;;; calculate difference between image position and photo image position

  ; mapping will use position of this color
  IF n_elements(clr) EQ 0 THEN clr = 2

  offsetx = objx - photostr[index].colc[clr]
  offsety = objy - photostr[index].rowc[clr]

  field = photostr[index].field
  id = photostr[index].id
  w=where( photostr.field EQ field, count)

  IF (count LT  3) THEN BEGIN
    print,'Not enough found!'
    count=ntot
    w=indgen(ntot)
  ENDIF 

  map_order=1
                                ;so the number of elements is same
  photomap, photostr[w], map_order, ra, dec, row, col, clr=clr

  xc = col + offsetx
  yc = row + offsety

  IF n_elements(order) NE 0 THEN BEGIN 
      IF order[0] EQ 1 THEN BEGIN 
          IF n_elements(ny) EQ 0 THEN message,'Must input ny if order= is sent'
      
          yc = (ny-1) - yc
      ENDIF 
  ENDIF 

  IF (n_elements(color) EQ 0) THEN BEGIN
      IF !d.n_colors GT 10000000 THEN BEGIN
          simpctable
          color=!white 
      ENDIF ELSE color = !d.n_colors
  ENDIF

  IF n_elements(pos_ang) NE 0 THEN BEGIN 
      IF (n_elements(rmax) EQ 0) OR (n_elements(rmin) EQ 0) THEN $
        message,'Must give rmax/rmin if pos_ang is given'

      tvellipse, rmax, rmin, xc+xshift, yc+yshift, pos_ang, color, /data

  ENDIF ELSE BEGIN 

      linest=0
      IF n_elements(linestyle) NE 0 THEN linest=linestyle[0]

      IF keyword_set(box) THEN BEGIN 
          tvbox, radius, xc+xshift, yc+yshift, color[0], $
            /data,linestyle=linest, noclip=0
      ENDIF ELSE IF keyword_set(cross) THEN BEGIN 
          plot_cross, radius, xc+xshift, yc+yshift, $
            color=color[0], $
            holeradius=holeradius
      ENDIF ELSE BEGIN 
          tvcircle, radius, xc+xshift, yc+yshift, color[0], $
            /data,linestyle=linest, noclip=0
      ENDELSE 

  ENDELSE 

return
end
























