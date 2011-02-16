PRO make_mask_sdss2d_findstripe, y, bound_arr, w, quit

  delvarx,w
  nstripe = n_elements(bound_arr)

  FOR ist=0L,  nstripe-1 DO BEGIN 
      
      IF ( (y LE bound_arr[ist].etamax) AND $
           (y GE bound_arr[ist].etamin) ) THEN BEGIN 
          
          w = ist
      ENDIF 
              
  ENDFOR 

  IF n_elements(w) EQ 0 THEN BEGIN 
      print,'The point was not found within any etamin/max bounds'
      quit = 0
  ENDIF 

END 

PRO make_mask_sdss2d, datax, datay, bound_arr, mask, bad, good, ox=ox, oy=oy, noprimary_bound_cut=noprimary_bound_cut

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: make_mask_sdss2d, datax, datay, bound_arr, mask, bad, good, ox=ox, oy=oy, /noprimary_bound_cut'
      return
  ENDIF 
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; primary bound cut is optional: cut on the primary top and bottom
  ;; bound
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  spacing = 0.2
  orient = 45

  ;; how many stripes are here?
  nstripe = n_elements(bound_arr)

  ;; absolute max and min in eta
  etamax = max(bound_arr.etamax)
  etamin = min(bound_arr.etamin)

  etamaxstr = string(etamax,format='(F0,:,X)')
  etaminstr = string(etamin,format='(F0,:,X)')

  ndata = n_elements(datax)

  simpctable

  maxx = max(datax, min=minx)
  maxy = max(datay, min=miny)

  xdiff = maxx-minx
  ydiff = maxy-miny

  addx = xdiff/10.0
  addy = ydiff/10.0

  xrange = [0.,0.]
  myrange = [0.,0.]

  xrange[0] = minx - addx
  xrange[1] = maxx + addx

  myrange[0] = miny - addy
  myrange[1] = maxy + addy

  IF n_elements(ox) NE 0 THEN BEGIN 
      plot, ox, oy, psym=3, xrange=xrange, yrange=myrange,ystyle=1,xstyle=1, $
           title = 'Input Data'

      oplot, datax, datay, psym=8, color=!blue

  ENDIF ELSE BEGIN 

      plot,datax, datay, psym=3, xrange=xrange, yrange=myrange,ystyle=1,xstyle=1, $
           title = 'Input Data'

  ENDELSE 

  continue=1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; begin with the edge of the allowed region
  ;; edgecut is allowed: set keyword for
  ;; stripes 10,76,82,86 only, and never for
  ;; multiple stripes
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT keyword_set(noprimary_bound_cut) THEN BEGIN 
      mask = $
        "( (datay GE "+etamaxstr+") OR " + $
        "  (datay LE "+etaminstr+") )"
      addstr = " OR "

      oplot, [-1000.,1000.],[etamin,etamin], color=!red
      oplot, [-1000.,1000.],[etamax,etamax], color=!red

  ENDIF ELSE BEGIN 
      mask = ""
      addstr = ""
  ENDELSE 
  linestyle = 0
  nbadtot = 0L
  WHILE continue DO BEGIN 

      print,'------------------------------------------------------------------------'

      quit=0
      WHILE NOT quit DO BEGIN 
          reply = ' '
          print,'What kind of mask is this? (1 for REGULAR 2 for ENDCUT q to QUIT)'
          read,reply
          CASE ntostr(reply) OF 
              '1': BEGIN & type='regular' & quit=1 & END 
              '2': BEGIN & type='branch' & quit=1 & END
              'q': return
              ELSE: print,'invalid response'
          ENDCASE 
      ENDWHILE 

      quit=0
      WHILE NOT quit DO BEGIN 
          reply = ' '
          print,'Allow edgecut on this mask? (y/n or q to QUIT)'
          read,reply
          CASE strlowcase(reply) OF 
              'y': BEGIN & allowedgecut=1 & quit=1 & END 
              'n': BEGIN & allowedgecut=0 & quit=1 & END
              'q': return
              ELSE: print,'invalid response'
          ENDCASE 
      ENDWHILE 

      ;;;;;;;;;;;;;;;;;;;;;;;;
      ;; get the points
      ;;;;;;;;;;;;;;;;;;;;;;;;

      IF type EQ 'regular' THEN BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; regular 2d masks
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          quit=0
          WHILE NOT quit DO BEGIN 
              quit=1
              print,'click on a minimum x position'
              cursor, xmin, y_xmin, /data, /down
              make_mask_sdss2d_findstripe, y_xmin, bound_arr, w_min, quit 
              IF NOT quit THEN GOTO,jump
              
              yrange = [bound_arr[w_min].etamin, bound_arr[w_min].etamax]
              xminstr = string(xmin,format='(F0,:,X)')
              print,'   xmin = ',xminstr
              oplot, [xmin,xmin],yrange,color=!red, linestyle=linestyle
              
              print,'click on a maximum x position'
              cursor, xmax, y_xmax, /data, /down
              make_mask_sdss2d_findstripe, y_xmax, bound_arr, w_max, quit
              IF NOT quit THEN GOTO,jump
              
              yrange = [bound_arr[w_max].etamin, bound_arr[w_max].etamax]
              xmaxstr = string(xmax,format='(F0,:,X)')
              print,'   xmax = ',xmaxstr
              oplot, [xmax,xmax],yrange,color=!red, linestyle=linestyle
              
              IF w_min NE w_max THEN BEGIN 
                  print,'Both eta points should be within same stripe'
                  quit=0
              ENDIF 
              
jump:
          ENDWHILE 
          
          ymax = bound_arr[w_min].etamax
          ymin = bound_arr[w_min].etamin

          yminstr = string(ymin,format='(F0,:,X)')
          ymaxstr = string(ymax,format='(F0,:,X)')

          ;; Check if edgecut can be applied to this mask
          IF allowedgecut THEN BEGIN 
              xminstr = '('+xminstr+' - DATA_EDGECUT)'
              xmaxstr = '('+xmaxstr+' + DATA_EDGECUT)'

              yminstr = '('+yminstr+' - DATA_EDGECUT)'
              ymaxstr = '('+ymaxstr+' + DATA_EDGECUT)'
          ENDIF 

          xmin_condition = ' GE '
          xmax_condition = ' LE '
          ymin_condition = ' GE '
          ymax_condition = ' LE '
          combinestr = ' AND '

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; if /noprimary_bound_cut then no cuts in y for
          ;; this mask (can be applied later with maxy=maxy
          ;; type cuts in apply_mask2d
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF NOT keyword_set(noprimary_bound_cut) THEN BEGIN 
              xx=[xmin,xmax,xmax,xmin]
              yy=[yrange[0],yrange[0],yrange[1],yrange[1]]
              polyfill,xx,yy,/line_fill,color=!red,$
                spacing=spacing,orient=orient
              oplot,[xmin,xmax],[ymax,ymax],color=!red
              oplot,[xmin,xmax],[ymin,ymin],color=!red
              
              maskadd = $
                '( (datax'+xmin_condition+xminstr+')'+combinestr+$
                '  (datax'+xmax_condition+xmaxstr+')'+combinestr+$
                '  (datay'+ymin_condition+yminstr+')'+combinestr+$
                '  (datay'+ymax_condition+ymaxstr+') )'
          ENDIF ELSE BEGIN 
              xx=[xmin,xmax,xmax,xmin]
              yy=[myrange[0],myrange[0],myrange[1],myrange[1]]
              polyfill,xx,yy,/line_fill,color=!red,$
                spacing=spacing,orient=orient
              oplot, [xmin, xmin], [myrange[0], myrange[1]], color=!red
              oplot, [xmax, xmax], [myrange[0], myrange[1]], color=!red

              maskadd = $
                '( (datax'+xmin_condition+xminstr+')'+combinestr+$
                '  (datax'+xmax_condition+xmaxstr+') )'
          ENDELSE 

      ENDIF ELSE BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; this is for edgecut in lambda at ends
          ;; of data
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          quit = 0
          WHILE NOT quit DO BEGIN 
              quit = 1
              print,'click on a minimum x position'
              cursor, xmin, y_xmin, /data, /down
              make_mask_sdss2d_findstripe, y_xmin, bound_arr, w_min, quit 
              IF NOT quit THEN GOTO,jump2

              yrange = [bound_arr[w_min].etamin, bound_arr[w_min].etamax]
              xminstr = string(xmin,format='(F0,:,X)')
              print,'   xmin = ',xminstr
              oplot, [xmin,xmin],yrange,color=!red, linestyle=linestyle
          
              print,'click on a maximum x position'
              cursor, xmax, y_xmax, /data, /down
              make_mask_sdss2d_findstripe, y_xmax, bound_arr, w_max, quit
              IF NOT quit THEN GOTO,jump2

              yrange = [bound_arr[w_max].etamin, bound_arr[w_max].etamax]
              xmaxstr = string(xmax,format='(F0,:,X)')
              print,'   xmax = ',xmaxstr
              oplot, [xmax,xmax],yrange,color=!red, linestyle=linestyle

              IF w_min NE w_max THEN BEGIN 
                  print,'Both eta points should be within same stripe'
                  quit=0
              ENDIF
jump2:
          ENDWHILE 

          ymax = bound_arr[w_min].etamax
          ymin = bound_arr[w_min].etamin

          yminstr = string(ymin,format='(F0,:,X)')
          ymaxstr = string(ymax,format='(F0,:,X)')

          ;; Check if edgecut can be applied to this mask
          IF allowedgecut THEN BEGIN 
              xminstr = '('+xminstr+' - DATA_EDGECUT)'
              xmaxstr = '('+xmaxstr+' + DATA_EDGECUT)'

              yminstr = '('+yminstr+' - DATA_EDGECUT)'
              ymaxstr = '('+ymaxstr+' + DATA_EDGECUT)'
          ENDIF 

          xmin_condition = ' LE '
          xmax_condition = ' GE '
          ymin_condition = ' GE '
          ymax_condition = ' LE '
          xcombinestr = ' OR '
          ycombinestr = ' AND '

          xx=[xrange[0],xmin,xmin,xrange[0]]
          yy=[yrange[0],yrange[0],yrange[1],yrange[1]]
          polyfill,xx,yy,/line_fill,color=!red,spacing=spacing,orient=orient

          xx=[xmax,xrange[1],xrange[1],xmax]
          yy=[yrange[0],yrange[0],yrange[1],yrange[1]]
          polyfill,xx,yy,/line_fill,color=!red,spacing=spacing,orient=orient

;          maskadd = $
;            '( (datax'+xmin_condition+xminstr+')'+xcombinestr+$
;            '  (datax'+xmax_condition+xmaxstr+') )'

;          maskadd = $
;            ( ((y ge ymin) and (y le ymax)) and $
;              ( (x LE xmin) OR (x GE xmax) ) )

          maskadd = $
            '( ( (datay GE '+yminstr+') AND (datay LE '+ymaxstr+') ) AND '+$
            '  ( (datax LE '+xminstr+') OR  (datax GE '+xmaxstr+') ) )'

      ENDELSE 

      print,'mask:  '
      print,'   ',maskadd
      mask = maskadd + addstr +  mask

      quit=0
      WHILE NOT quit DO BEGIN 
          reply = ' '
          read,'Want to add another mask? (y/n) ',reply
          CASE strlowcase(ntostr(reply)) OF
          'y': BEGIN & continue = 1 & quit=1 & END 
          'n': BEGIN & continue = 0 & quit=1 & END
          'q': return
              ELSE: print,'invalid response'
          ENDCASE 
      ENDWHILE 

      addstr = " OR "

  ENDWHILE 

END 
