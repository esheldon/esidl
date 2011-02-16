pro  fchart, pstruct, index, radius, clr, fchart, dir=dir, $
             objx=objx, objy=objy, impos=impos, silent=silent, $
             maxsize=maxsize, nonoise=nonoise
             
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME: 
;    FCHART
;       
; PURPOSE: 
;    Create a finding chart for an SDSS object.
;	
; CALLING SEQUENCE: 
;
;    fchart, pstruct, index, radius, clr, fchart, dir=dir, 
;            silent=silent, objx=objx, objy=objy, nonoise=nonoise
;
; COMMENTS:  The wrapper obj_info is a very convenient way to use FCHART
;            since it makes a nice display.
;
; INPUTS: 
;    
;    pstruct:  A photo structure.
;              IMPORTANT NOTE: This MUST have objc_rowc as a continuously added
;              added number, not restarted at the beginning of each frame.
;              Such boundaries are artificial and make it difficult to make 
;              finding charts centered on the object.  
;              This is done automatically by READ_TSOBJ
;
;    index:    The index of the object in the photo structure which needs
;              a finding chart.
;    radius:   The box which the object is in will have sides 2*radius
;              unless radius is too big. (e.g. bigger that 2048 pixels)
;    clr:      The color used to make the finding chart. This follows the
;              photo convention: 0,1,2,3,4 -> u,g,r,i,z
;
; INPUT KEYWORD PARAMETERS:
;         dir:  The directory to find the atlas images.
;         clr:  The color used to make the finding chart. This follows the
;               photo convention: 0,1,2,3,4 -> u,g,r,i,z
;         silent: If silent is set, nothing is printed except some errors
;         nonoise: if set, no noise is added to image.
;         maxsize: maximum size for atlas images. Default is [500,500].  If 
;               your images are clipped off (as with large galaxies) increase
;               maxsize.
;         
; OUTPUTS: 
;
;    fchart:  The image of the finding chart. 
;
; OPTIONAL OUTPUTS: 
;    objx: The x-position of the input object in pixels in image
;               coordinates (as opposed to photo objc_rowc, etc)
;    objy: Same for y.
;    impos: Absolute position of left hand corner of image.
;
; CALLED ROUTINES: 
;
;    GET_FAMILY
;    GET_ATLAS
; 
; PROCEDURE: 
;  Create a finding chart around an input object from all the objects nearby.
;  The trick is to only use the complete images, not the ones with pieces
;  cut out of them.  This program uses the get_family procedure to figure
;  out which are the good ones:
;
;    orphans:  Those objects with no siblings or parents.
;    child of a bad parent: A bad parent will have one child which is 
;             a fainter version of itself.
;    good parents: The clean parent image. (from which good children 
;             are made)
;    grandparents: Sometimes there is a grandparent, from which comes
;             only one child.  This child has children of its own.
;	
;  From here, the atlas image  of each good neighbor is the proper color is 
;  found and placed within the appropriate box.  
;  Note this box may not center on the main object if it is less that 'radius'
;  from either the edge of the frame in the column direction. This is also, 
;  true if it is near either end of the series of frames read into pstruct.
;
;
; REVISION HISTORY: 
;     Author  Erin Scott Sheldon UM 03/21/99
;     Dave Johnston - was adding way too much noise
;        to background in some cases , now it just adds a 
;        trace amount of noise to background 
;        sky rms sqrt(30)  5/15/99	
;     Now allows objects with center outside image to 
;     contribute light to the image. Object centers must
;     be within maxsize/2.  14-NOV-2000
;       
;                                      
;-                                        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if N_params() LT 4 then begin
	print,'-Syntax: fchart, pstruct, index, radius, clr, fchart, dir=dir,'
        print,'   objx=objx, objy=objy, silent=silent, nonoise=nonoise'
        print,''
        print,'Use  doc_library,"fchart" for more help.'
	return
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; check some keywords, set parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  chartmin = 20  ;;; minimum value in final chart
  colors=['u','g','r','i','z']

  if (clr lt 0 or clr gt 4) then begin
     print,'Color index must be in [0,4]'
     return
  endif
  if n_elements(dir) EQ 0 then dir=''
  if (not keyword_set(silent)) then silent=0
  IF NOT keyword_set(nonoise) THEN nonoise=0
  IF NOT keyword_set(maxsize) THEN maxsize=[500,500]

  useind = lindgen(n_elements(pstruct))

  ;; average differences of rowc, colc between bandpasses and r-band

  cdiff = median( pstruct.colc[2] - pstruct.colc[clr] )
  rdiff = median( pstruct.rowc[2] - pstruct.rowc[clr] )
  difftol = 3.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; make sure picture area as defined by radius does not extend ;;;
  ;; beyond the bounds of sloan frame in the column direction    ;;;
  ;; Also, since all we have is pstruct, we must make sure that  ;;;
  ;; rows are within bounds set by pstruct                       ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  object = pstruct[index]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; This method of calculating the minrow/maxrow is not that good, it
  ;; will definitely be too large/small.  But it should only mattter
  ;; in the first or last frame where radius will go out of the bounds
  ;; Note the 64 is due to overlaps
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  tmp1=min(pstruct.objc_rowc)
  tmp2=max(pstruct.objc_rowc) 
  f=object.field
  
  IF clr NE 2 THEN BEGIN 
      offcheck = object.rowc[2] - object.rowc[clr]
      IF abs(offcheck - rdiff) GT difftol THEN BEGIN 
          object.rowc[clr] = object.rowc[2] - rdiff
      ENDIF 
      offcheck = object.colc[2] - object.colc[clr]
      IF abs(offcheck - cdiff) GT difftol THEN BEGIN 
          object.colc[clr] = object.colc[2] - cdiff
      ENDIF 
  ENDIF 

  CASE 1 OF
    (f eq min(pstruct.field)) : BEGIN
         offset=object.rowc[clr] - (object.objc_rowc - tmp1)
         if (not silent) then print,''
         if (not silent) then print,'OFFSET AT BEGINNING: ',$
                 strtrim(string(offset),2)
         minrow=tmp1 - offset
         maxrow=tmp2
       END
    (f eq max(pstruct.field)) : BEGIN
         offset = ((1489 - 1) - object.rowc[clr]) - $
                (tmp2 - object.objc_rowc)
         if (not silent) then print,''
         if (not silent) then print,'OFFSET AT END: ',$
                 strtrim(string(offset),2)
         maxrow = tmp2 + offset
         minrow=tmp1
       END
    else: BEGIN 
         minrow=tmp1
         maxrow=tmp2
       END
  ENDCASE

  maxcol = 2048-1  ;; note mincol is just 0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; may have to trim down radius if its too big
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  radius = long(radius)
  r2 = 2*radius
  diffrow=maxrow-minrow
  if (r2 ge diffrow) then begin
     print,'Radius too big.  Trimming'
     radius=diffrow/2
     r2=2*radius
  endif
  if (r2 ge maxcol) then begin
     print,'Radius*2 is larger than field width.  Trimming.'
     radius=maxcol/2 + 1
     r2=2*radius
  endif
  
  ;;; absolute position of object  ;;;
  objpos = [object.objc_colc, object.objc_rowc]

  ;;; objects initial relative position in image       ;;;
  ;;; using floating point for programs like tvellipse ;;;
  objx = radius-.5
  objy = radius-.5

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find center                                                
  ;; Check columns(x) first. Note center and relative position  
  ;; of object may change                                        
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  center=[ objpos[0],objpos[1] ]
  diff=objpos[0]-radius
  if (diff lt 0) then begin
     center[0]=center[0] - diff
     objx = objx + diff
  endif else begin
     diff=maxcol - (objpos[0] + radius)
     if (diff lt 0) then begin
        center[0]=center[0] + diff
        objx = objx - diff
     endif
  endelse

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Check rows (y) 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  diff=(objpos[1]-radius) - minrow
  if (diff lt 0) then begin
     center[1]=center[1] - diff
     objy=objy + diff
  endif else begin
     diff=maxrow - (objpos[1] + radius)
     if (diff lt 0) then begin
        center[1]=center[1] + diff
        objy = objy - diff
     endif
  endelse
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Bottom left corner in absolute frame. This  will be the origin 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  lcorner = [center[0] - radius, center[1] - radius]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; the region covered by the image
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  col = [center[0] - radius, center[0] + radius]
  row = [center[1] - radius, center[1] + radius]
  impos = [col[0], row[0]]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Print out a bunch of cool stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  id = strtrim(string(object.id),2)
  if (not silent) then begin
     print,''
     print,'-----------------------------'
     print,'Building Finding Chart From ',$
            colors[clr],' Images For Object id ',id
     strr2=strtrim(string(r2),2)
     print,'Chart size: [',strr2,',',strr2,']'
     print,'Object center is at: [',strtrim(string(long(objx)),2),',',$
            strtrim(string(long(objy)),2),']'
     print,'-----------------------------'

  endif

  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Define the finding chart
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fchart=intarr(r2,r2)
  fchart[*,*] = 1000

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;  select stuff within box of side 2*radius = r2. Centers can be
  ;;  maxsize/2 away from that edge.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  w=where( (pstruct[useind].objc_colc le (col[1]+maxsize[0]/2.) ) and $
           (pstruct[useind].objc_colc ge (col[0]-maxsize[0]/2.) ) and $
           (pstruct[useind].objc_rowc le (row[1]+maxsize[0]/2.) ) and $
           (pstruct[useind].objc_rowc ge (row[0]-maxsize[0]/2.) ), count )

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;  Find the atlas images without stuff cut out of them
  ;;  Speed is limited by all the if statements and size of pstruct
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if (count ne 0) then begin

    good = [-33]
    for i=0, count-1 do begin
       kk = w[i]
       parent = -33
       children = -33
       get_family, pstruct, useind[kk], parent, children, $
         dir=dir,/nodisplay,/silent
   
       ;;;;; is it good? ;;;;;;
       nchild = n_elements(children)
       goodid = -33
       if (parent[0] eq -33 and nchild eq 1) then begin
          if (children[0] eq -1) then begin ;;;;;; orphan  ;;;;;;;
             goodid = useind[kk]
          endif else begin ;;;;; bad parent or grandparent  ;;;;;
             goodid = children[0]
          endelse
       endif else if (parent[0] eq -33 and nchild gt 1) then begin
          goodid = useind[kk]          ;;;; good parent
       endif
       if (goodid ne -33) then begin
           if (good[0] eq -33) then good = goodid else good = [good,goodid]
       endif
   
    endfor

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;  Get the atlas images and place them in chart
    ;;  NOTE: Speed of this is limitied by get_atlas 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    max = n_elements(good)
    if (not silent) then begin
        printstr = 'Checking '+strtrim(string(max),2)+' neighbors'
      print,printstr
    endif
  
    strmax=strtrim(string(max),2)
    if (good[0] eq -33) then good = index
    for i=0, max-1 do begin
       IF ( ( (i EQ 0) OR (i EQ max-1) OR  ((i+1) mod 40) eq 0 ) $
            and not silent  ) then print,'.',format='(a,$)'
       gi=good[i]
       get_atlas, pstruct, gi, dir=dir, clr=clr, imtot=im, maxsize=maxsize,$
         /nodisplay,/noprompt,/silent, col0=col0, row0=row0

       IF n_elements(im) EQ 0 THEN GOTO,jump

       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       ;; [col0, row0] is position of bottom left corner of atlas 
       ;; in its field (col0 same for field and absolute coord) 
       ;; Convert row0 to absolute coordinates  (Uses fact that 
       ;; rowc are not added up by read_tsobj)   
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

       rowold = row0
       colold = col0

       row0 = row0 + ( pstruct[gi].objc_rowc - pstruct[gi].rowc[clr] )
       col0 = col0 + ( pstruct[gi].objc_colc - pstruct[gi].colc[clr] )

       IF clr NE 2 THEN BEGIN 
           rowdiff = pstruct[gi].rowc[2] - pstruct[gi].rowc[clr]
           coldiff = pstruct[gi].colc[2] - pstruct[gi].colc[clr]

           IF ( ( abs(rowdiff - rdiff) GT difftol )  OR $
                ( abs(coldiff - cdiff) GT difftol ) ) THEN BEGIN 

               row0 = rowold + ( pstruct[gi].objc_rowc - $
                                 (pstruct[gi].rowc[2]-rdiff )  )
               col0 = colold + ( pstruct[gi].objc_colc - $
                                 (pstruct[gi].colc[2]-cdiff) )

           ENDIF 
       ENDIF 

       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       ;; put row0 and col0 in our image coordinates 
       ;; lcorner is the absolute position of the 
       ;; bottom left hand corner of our image  
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

       im_col0 = round(col0 - lcorner[0])
       im_row0 = round(row0 - lcorner[1])

       sz = size(im)
       imcols = sz[1]
       imrows = sz[2]

       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       ;; The pixels of atlas in fchart coords 
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

       x = lindgen(imcols) + im_col0
       y = lindgen(imrows) + im_row0

       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       ;; Positions which are acutally in fchart 
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       wx = where(x gt 0 and x lt r2, nwx)
       wy = where(y gt 0 and y lt r2, nwy)
       IF (nwx NE 0) AND (nwy NE 0) THEN BEGIN 
         maxwx=max(wx)
         minwx=min(wx)
         maxwy=max(wy)
         minwy=min(wy)

         maxx=x[maxwx]
         minx=x[minwx]
         maxy=y[maxwy]
         miny=y[minwy]

         ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
         ;; Place the atlas images.  Sacrificing a bit of memory for speed
         ;; Should be fine for fcharts of this size ( < [2048,2048] )
         ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

         goodim = im[minwx:maxwx, minwy:maxwy]
         goodfchart = fchart[minx:maxx,miny:maxy]
         goodfchart = goodim + $
              (goodfchart ne 1000)*(goodfchart - 1000)
         fchart[minx:maxx,miny:maxy] = goodfchart
       ENDIF 
       delvarx,im
       jump:
    endfor
    if (not silent) then BEGIN
        print
        print,'-----------------------------'
    ENDIF 

  endif else begin ;;;; if count was zero, then no objects were placed! ;;;
     print, 'No objects were found.  Finding Chart is empty!'
  endelse

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Add noise (if requested) and correct saturation, etc
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  w=where(fchart EQ 0, nw)
  IF nw NE 0 THEN fchart[w] = 1000

  nw=0L
  IF NOT nonoise THEN w = where(fchart eq 1000, nw) ;Places to add noise

  ;; Add noise for sky of about 30.
  IF (NOT nonoise) AND (nw NE 0) THEN BEGIN
      fchart[w] = fchart[w] + sqrt( 30.0 )*randomn(seed,nw)
  ENDIF 

  l = where(fchart LT 0, nl)
  IF nl NE 0 THEN fchart[l] = 32767

return
end























