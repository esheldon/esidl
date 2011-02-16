PRO make_mask, datax, datay, mask, bad, good, apply=apply, ox=ox, oy=oy

  ;;
  ;; create a mask on the x-coordinate of 2-d data that can be used 
  ;; in the following way:
  ;;       if mask then "don't use this point"
  ;; or
  ;;       bad = where(mask)
  ;;
  ;; "across branch" masks are used when crossing branch cuts in the
  ;; coordinate system. e.g. the [-180,180] point 
  ;; currently doesn't work for [0,360] however...
  ;; 

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: make_mask, datax, datay, mask, bad, good, apply=apply, ox=ox, oy=oy'
      return
  ENDIF 

  IF NOT keyword_set(apply) THEN apply=0
  IF apply THEN BEGIN
      pold = !p.multi
      !p.multi = [0,0,2]
  ENDIF 

  ndata = n_elements(datax)

  simpctable

  maxx = max(datax, min=minx)
  maxy = max(datay, min=miny)

  xdiff = maxx-minx
  ydiff = maxy-miny

  addx = xdiff/10.0
  addy = ydiff/10.0

  xrange = [0.,0.]
  yrange = [0.,0.]

  xrange[0] = minx - addx
  xrange[1] = maxx + addx

  yrange[0] = miny - addy
  yrange[1] = maxy + addy

  IF n_elements(ox) NE 0 THEN BEGIN 
      plot, ox, oy, psym=3, xrange=xrange, yrange=yrange,ystyle=1,xstyle=1, $
           title = 'Input Data'

      oplot, datax, datay, psym=8, color=!blue

  ENDIF ELSE BEGIN 

      plot,datax, datay, psym=8, xrange=xrange, yrange=yrange,ystyle=1,xstyle=1, $
           title = 'Input Data'

  ENDELSE 
  continue=1

  mask = ""
  addstr = ""

  linestyle = 0
  nbadtot = 0L
  WHILE continue DO BEGIN 

      print,'------------------------------------------------------------------------'

      quit=0
      WHILE NOT quit DO BEGIN 
          reply = ' '
          print,'What kind of mask is this? (1 for regular 2 for branch q to quit)'
          read,reply
          CASE ntostr(reply) OF 
              '1': BEGIN & type='regular' & quit=1 & END 
              '2': BEGIN & type='without' & quit=1 & END
              'q': return
              ELSE: print,'invalid response'
          ENDCASE 
      ENDWHILE 

      quit=0
      WHILE NOT quit DO BEGIN 
          reply = ' '
          print,'Allow edgecut on this mask? (y/n or q to quit)'
          read,reply
          CASE strlowcase(reply) OF 
              'y': BEGIN & allowedgecut=1 & quit=1 & END 
              'n': BEGIN & allowedgecut=0 & quit=1 & END
              'q': return
              ELSE: print,'invalid response'
          ENDCASE 
      ENDWHILE 



      print,'click on a minimum x position'
      cursor, xmin, y, /data, /down
      xminstr = string(xmin,format='(F0,:,X)')
      print,'   xmin = ',xminstr
      oplot, [xmin,xmin],yrange,color=!red, linestyle=linestyle

      print,'click on a maximum x position'
      cursor, xmax, y, /data, /down
      xmaxstr = string(xmax,format='(F0,:,X)')
      print,'   xmax = ',xmaxstr
      oplot, [xmax,xmax],yrange,color=!red, linestyle=linestyle

      IF type EQ 'regular' THEN BEGIN 

          ;; Check if edgecut can be applied to this mask
          IF allowedgecut THEN BEGIN 
              xminstr = '('+xminstr+' - DATA_EDGECUT)'
              xmaxstr = '('+xmaxstr+' + DATA_EDGECUT)'
          ENDIF 

          xmin_condition = ' GE '
          xmax_condition = ' LE '
          combinestr = ' AND '

          xx=[xmin,xmax,xmax,xmin]
          yy=[yrange[0],yrange[0],yrange[1],yrange[1]]
          polyfill,xx,yy,/line_fill,color=!red

      ENDIF ELSE BEGIN 

          ;; Check if edgecut can be applied to this mask
          IF allowedgecut THEN BEGIN 
              xminstr = '('+xminstr+' + DATA_EDGECUT)'
              xmaxstr = '('+xmaxstr+' - DATA_EDGECUT)'
          ENDIF 

          xmin_condition = ' LE '
          xmax_condition = ' GE '
          combinestr = ' OR '

          xx=[xrange[0],xmin,xmin,xrange[0]]
          yy=[yrange[0],yrange[0],yrange[1],yrange[1]]
          polyfill,xx,yy,/line_fill,color=!red

          xx=[xmax,xrange[1],xrange[1],xmax]
          yy=[yrange[0],yrange[0],yrange[1],yrange[1]]
          polyfill,xx,yy,/line_fill,color=!red
      ENDELSE 

      maskadd = '( (data'+xmin_condition+xminstr+')'+combinestr+$
                              '(data'+xmax_condition+xmaxstr+') )'

      print,'mask:  '
      print,'   ',maskadd
      mask = maskadd + addstr +  mask

      IF apply THEN BEGIN 
          apply_mask, maskadd, datax, tbad, tgood
          nbad = ndata-n_elements(tgood)
          nbadtot = nbadtot+nbad
          print
          print,ntostr(nbad),' elements thrown out by this mask'
      ENDIF 

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

      addstr = ' OR '

  ENDWHILE 

  IF apply THEN BEGIN 
      apply_mask, mask, datax, bad, good
      IF good[0] NE -1 THEN BEGIN 
          plot,datax[good], datay[good], psym=3, $
            xrange=xrange, yrange=yrange,ystyle=1,xstyle=1, $
            title = 'Good Data'
      ENDIF 
      !p.multi=pold
      print
      print,ntostr(nbadtot),' total elements thrown out'
  ENDIF 

END 
