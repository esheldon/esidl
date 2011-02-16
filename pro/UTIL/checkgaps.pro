PRO checkgaps, first, last, isgap, names=names_in

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: checkgaps, first, last, isgap, names=names'
      print,'Input up to 8 first/last pairs'
      return
  ENDIF 

  !p.multi=[0,1,2]

  isgap = 0

  nfirst = n_elements(first)
  nlast = n_elements(last)
  IF nfirst NE nlast THEN message,'Nlast must equal Nfirst'

  nnames = n_elements(names_in)
  IF nnames EQ 0 THEN BEGIN
      names=ntostr(lindgen(nfirst))
      nnames = nfirst
  ENDIF ELSE names=names_in
  IF nnames NE nfirst THEN message,'Nnames must equal Nfirst'

  miny=-0.5
  maxy= 0.5
  IF !d.name EQ 'PS' THEN simpctable
  colors = [!red, !green, !blue, !magenta, !cyan, !seagreen, !orange, !purple,$
            !slategrey, !navyblue, !firebrick, !hotpink, !turquoise]
  ncolors = n_elements(colors)

  xrange = prange(first,last,/noerr,/slack)

  plot,[0], /nodata, xrange=xrange, yrange=[-1,1]

  FOR i=0L, nfirst-1 DO BEGIN 

      plot_box, first[i], last[i], miny, maxy, $
                color=colors[i MOD ncolors]

  ENDFOR 

  legend,names,line=replicate(0,nfirst),colors=colors[lindgen(nfirst) MOD ncolors],$
         charsize=0.75

  ;; throw out ranges completely contained within another range
  ;; obviously not a gap issue there

  keep = lindgen(nfirst)
  plot,[0], /nodata, xrange=xrange, yrange=[-1,1]

  FOR i=0L, nfirst-1 DO BEGIN

      w=where( first LT first[i] AND last GT last[i], nw)
      IF nw NE 0 THEN BEGIN
          
          message,names[i]+' is fully contained',/inf
          add_arrval, i, throwout

      ENDIF ELSE BEGIN 
          plot_box, first[i], last[i], miny, maxy, $
                color=colors[i MOD ncolors]
      ENDELSE 

  ENDFOR 

  IF n_elements(throwout) NE 0 THEN remove, throwout, keep
  nkeep = n_elements(keep)

  legend, names[keep], line=replicate(0,nkeep), colors=colors[keep MOD ncolors],$
          charsize=0.75

  ;; look for gaps

  IF nkeep GT 1 THEN BEGIN 
      FOR i=0L,nkeep-1 DO BEGIN 
         
          ind = keep[i]

          ;; find ranges than end past the 
          ;; current range
          w=where(keep NE ind AND $
                  last[keep] GT last[ind], nw)
          IF nw NE 0 THEN BEGIN 
              minfirst = min(first[keep[w]])
          
              ;; gap?
              IF minfirst GT last[ind] THEN BEGIN 
                  gap = abs(minfirst - last[ind])
                  wnext = where(first EQ minfirst)
                  isgap = 1
                  message,'Gap of '+ntostr(gap)+' between '+names[ind]+' and '+names[wnext],/inf
              ENDIF 

          ENDIF 
      ENDFOR 


  ENDIF 

  IF NOT isgap THEN message,'No gaps found',/inf ELSE message,'Gaps were found',/inf

  !p.multi=0

END 
