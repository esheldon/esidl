;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: 
;    DISPLAY_FCHART
;       
; PURPOSE: display a photo finding chart and circle the object that was used to
;          make the chart.  Label axis with arcminutes from the center.  Output
;          inmages, .jpg, .png, or .ps if requested.  Gif no longer supported
;          by idl.
;	
;
; CALLING SEQUENCE:  
;      
;                 
;
; INPUTS: fchart:  The finding chart
;         object:  The photo object struct used to make the finding chart.
;         objx:    x-position of object in the finding chart.
;         objy:    y-position
;         clr:     color index the finding chart was made from 
;                                          used in the title.
;
; INPUT KEYWORD PARAMETERS:
;         order=: Order to draw images. default is order=0 is 0,0 at bottom
;                 left. for order=1, 0,0 is top left.
;         /box:  Draw a box instead of a circle.
;         /jpeg: Write a jpeg from display.
;         /png: Write a png from display.
;         /ps:   If set a postscript files are created for the requested items.
;         fnamepng: name of png file. Default fchart.png
;         fnameps: name of ps file.  Default is fchart.ps
;         fnamegif: name if gif file.  Default is fchart.gif
;         silent:  If set, program will be silent except for error messages.
;         nodisplay: If set, there will be No x-window display of the chart.  
;               -May still make ps files.-
;         hideradec: Don't show the ra and dec on the plots.
;         circ_rad: radius of circle.
;         nocirc:  if set, no circle is drawn.
;         _extra: any extra plotting command keywords.
;       
; OUTPUTS: 
;
; OPTIONAL OUTPUTS: May output .ps or .gif files if requested.
;
; CALLED ROUTINES:
;
;                  RADECSTR
;                     RADEC
;                  SIGMA_CLIP
;                  BEGPLOT
;                  ENDPLOT
;                  RDIS_SETUP
;                  RDIS:
;                      TVIM2_SCL:
;                               TVIM2
;                  TVCIRCLE
;                  TVBOX
;
;
; REVISION HISTORY: 
;      Author: Erin Scott Sheldon  Umich 5/25/99
;      Changed order in title: run-rerun-camcol-field-id 
;      Fixed circle color problem.
;                                               02-Jul-2002
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO fchart_display_defradius, image, markStruct, markRadius, $
                              pixels=pixels, arcseconds=arcseconds, $
                              arcminutes=arcminutes

  COMMON markstruct_block, $
    defMarkRA, defMarkDEC, defMarkObjx, defMarkObjy, $
    defMarkType, defMarkHoleFraction, defMarkLineStyle, $
    defMarkRadius, defMarkColor, $
    defaultMarkStruct

  IF keyword_set(pixels) THEN BEGIN 
      cfac = 1.0
  ENDIF ELSE IF keyword_set(arcseconds) THEN BEGIN 
      cfac = 1.0/0.4
  ENDIF ELSE BEGIN ;; Default is arcminutes
      cfac = 60.0/0.4
  ENDELSE 

  ;; Some statistics
  nMark = n_elements(markStruct.ra)
  imsize = size(image)
  ny = imsize[2]
  ss = imsize[1]

  defrad = ss/10.0

  ;; radius
  w=where(markStruct.radius EQ defMarkRadius,nw)
  IF nw NE 0 THEN BEGIN
      IF nw EQ n_elements(markStruct.radius) THEN BEGIN 
          IF (nMark-1) GT 0 THEN BEGIN 
              markRadius = [defrad, replicate(default/2.0,nMark-1)]
          ENDIF ELSE BEGIN 
              markRadius = defrad
          ENDELSE 
      ENDIF
  ENDIF ELSE BEGIN 
      ;; Convert radius to pixels if needed
      markRadius = markStruct.radius*cfac
  ENDELSE 

END 

PRO fchart_display_getcolor, circ_color
  IF !d.name EQ 'X' THEN BEGIN 
      defsysv, '!white', exist=wexist
      IF NOT wexist THEN BEGIN 
          max = float(!d.n_colors)-1.
          circ_color = 1.*max
      ENDIF ELSE circ_color=!white
  ENDIF ELSE BEGIN 
      circ_color=!p.color
  ENDELSE 
END 

PRO fchart_display_markstruct_checktags, markStruct

  tags = tag_names(markStruct)

  match, tags, ['OBJY', 'OBJX'], mtags, mradec

  IF n_elements(mtags) NE 2 THEN BEGIN 
      message,'markStruct must contain objx and objy'
  ENDIF 

END 

PRO fchart_display, fchart, $
                    object=object, $
                    maguse=maguse,$
                    $
                    clr=clr, $
                    $
                    markStruct=markStruct, $
                    extra_markStruct=extra_markStruct, $
                    $           ; Units for all lenghts
                    pixels=pixels, $
                    arcminutes=arcminutes, $
                    arcseconds=arcseconds, $
                    $
                    radec=radec, $
                    directions=directions,$
                    $
                    order=order, $
                    $
                    png=png, jpeg=jpeg, gif=gif, ps=ps, $
                    fnameps=fnameps, $
                    fnamegif=fnamegif, $
                    $
                    silent=silent,$
                    nodisplay=nodisplay,$
                    $
                    _extra=extra

  IF N_params() LT 1 THEN BEGIN 
     print,'-Syntax: '
     print,''
     print,'Use doc_library,"display_fchart"  for more help.'  
     return
  ENDIF 

  on_error, 2

  IF n_elements(fchart) EQ 0 THEN BEGIN 
      print,'No fchart!'
      return
  ENDIF 

  IF n_elements(markStruct) NE 0 THEN BEGIN 
      fchart_display_markstruct_checktags, markStruct
      
      markstruct_copy_defaults, markStruct
  ENDIF 

  IF n_elements(extra_markStruct) NE 0 THEN BEGIN 
      fchart_display_markstruct_checktags, extra_markStruct
      
      markstruct_copy_defaults, extra_markStruct
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Some initializations
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT keyword_set(gif)       THEN gif = 0
  IF NOT keyword_set(ps)        THEN ps=0
  IF NOT keyword_set(fnameps)   THEN fnameps = 'fchart.ps'
  IF NOT keyword_set(fnamegif)  THEN fnamegif = 'fchart.gif'
  IF NOT keyword_set(silent)    THEN silent = 0
  IF NOT keyword_set(nodisplay) THEN nodisplay=0

  IF NOT keyword_set(nocirc) THEN nocirc = 0

  IF n_elements(clr) NE 0 THEN BEGIN 
      IF clr[0] LT 0 OR clr[0] GT 4 THEN BEGIN 
          print,'clr must be in [0,4]'
          clrstr = "?"
      ENDIF ELSE BEGIN 
          clrstr = colors[clr[0]]
      ENDELSE 
  ENDIF 

  ;; This is needed to fit all the info on.

  colors=['u','g','r','i','z'] 
 
  IF n_elements(object) NE 0 THEN BEGIN 

      tags = tag_names(object)

      ;; Is all the id info there?
      match,tags,['RUN','RERUN','CAMCOL','FIELD','ID'], mtags, mid

      IF n_elements(mtags) EQ 5 THEN BEGIN 

          run = run2string(object[0].run)
  
          field = field2string(object[0].field)
          camcol = strtrim(string(object[0].camcol),2)
          id = strn(object[0].id,length=5,padchar='0')
          rerun='-'+strtrim(string(object[0].rerun), 2)

          object_name = run+rerun+'-'+camcol+'-'+field+'-'+id

          title = object_name

          IF n_elements(clr) NE 0 THEN dclr = clr[0] ELSE dclr=2

          angle = $
            angle_rowcol2radec(object.run, $
                               object.camcol, $
                               object.field, dclr)

      ENDIF ELSE title = ''

      IF keyword_set(radec) THEN BEGIN 
          radecstr, object.ra, object.dec, rastr, decstr
          IF title EQ '' THEN title=rastr+' : ' + decstr $
          ELSE title = title + '  ' +rastr+' : ' + decstr
      ENDIF

      wmag = sdss_maguse(object, maguse=maguse, silent=silent)

      IF wmag NE -1 THEN BEGIN 
          
          c=object[0].(wmag)
          xtitle=''
          FOR kk=0, 4 DO BEGIN
              IF (xtitle NE '') THEN xtitle = xtitle+'  '
              mag = strmid( strtrim(string(c[kk]),2), 0, 5)
              xtitle = xtitle + colors[kk]+'='+mag
          ENDFOR
          
          addtitle = ''
          FOR kk=0,3 DO BEGIN
              diff = strmid(strtrim(string(c[kk] - c[kk+1]),2), 0, 5)
              addtitle = addtitle + '  '+colors[kk]+'-'+colors[kk+1]+'='+diff
          ENDFOR
          xtitle = xtitle + addtitle  
      ENDIF ELSE BEGIN 
          xtitle = 'Offset (arcminutes)'
      ENDELSE 


  ENDIF ELSE BEGIN 
      title = ''
      xtitle = 'Offset (arcminutes)'
  ENDELSE 

  IF n_elements(clr) NE 0 THEN BEGIN 
      IF title EQ '' THEN title = 'Filter: '+clrstr $
      ELSE title = title  + '   Filter: '+clrstr
  ENDIF 

  ytitle='Offset (arcminutes)'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Now, draw the circle around our object and then display and save.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
  ;;;; set appropriate scaling of the circle ;;;;

  ss = size(fchart)
  sx = ss[1]*.4/60.0
  sy = ss[2]*.4/60.0
  IF (NOT keyword_set(circ_rad)) THEN circ_rad = min([ss[1],ss[2]])/20.0

  ;;;; make x and y tick labels
  xlabels = [ntostr(-2*sx/4.0, 5), $
             ntostr(-sx/4.0, 5)  , $
             '0.0', $
             ntostr(sx/4.0,4), $
             ntostr(2*sx/4.0,4)]
  ylabels = xlabels

  pold=!p.multi
  !p.multi=[0,0,1]

  ;; Postscript output
  IF (ps) THEN BEGIN 

      begplot, name=fnameps

      fchart_display_getcolor, circ_color

      tvasinh, fchart, /noframe, title=title, _extra=extra, order=order

      IF (NOT nocirc) THEN BEGIN 
          IF (NOT keyword_set(box) ) THEN $
            tvcircle, circ_rad, objx, objy, circ_color, /data $
          ELSE tvbox, 2*circ_rad, objx, objy, circ_color, /data

      ENDIF
      axis,xaxis=0,xticks=4,xtickn=xlabels, xtitle=xtitle
      axis,xaxis=1,xticks=4,xtickn=[' ',' ',' ',' ',' ']
      axis,yaxis=0,yticks=4,ytickn=ylabels, ytitle=ytitle
      axis,yaxis=1,yticks=4,ytickn=[' ',' ',' ',' ',' ']

      IF keyword_set(directions) AND n_elements(angle) NE 0 THEN BEGIN 
          plot_ne_arrows, angle, fracsize=0.05, offset_frac=offset_frac,$
            order=order
      ENDIF 

      endplot

  ENDIF 

  ;; X window and gif output
  IF (NOT nodisplay) THEN BEGIN 

      tvasinh, fchart, /noframe, title=title, _extra=extra, order=order
      fchart_display_getcolor, circ_color
      
      IF n_elements(markStruct) NE 0 THEN BEGIN 
          markxy, $
            markRadius[i], markx, marky, $
            mark_type = markStruct.type[i], $
            holeFraction=markStruct.hole_fraction[i], $
            color=color, linestyle=markStruct.linestyle[i], $
            $
            order=order, ny=ny
      ENDIF 

;      IF (NOT nocirc) THEN BEGIN
;          IF (NOT keyword_set(box) ) THEN $
;            tvcircle, circ_rad, objx, objy, circ_color, /data $
;          ELSE tvbox, 2*circ_rad, objx, objy, circ_color, /data
;      ENDIF

      IF n_elements(markStruct) NE 0 THEN BEGIN 
          
      ENDIF 

      axis,xaxis=0,xticks=4,xtickn=xlabels, xtitle=xtitle
      axis,xaxis=1,xticks=4,xtickn=[' ',' ',' ',' ',' ']
      axis,yaxis=0,yticks=4,ytickn=ylabels, ytitle=ytitle
      axis,yaxis=1,yticks=4,ytickn=[' ',' ',' ',' ',' ']

      IF keyword_set(directions) AND n_elements(angle) NE 0 THEN BEGIN 
          plot_ne_arrows, angle, fracsize=0.05, offset_frac=offset_frac,$
            order=order
      ENDIF 

      IF (gif) THEN BEGIN
          wshow, 0, iconic=0
          write_gif, fnamegif, tvrd()
      ENDIF
  ENDIF 
  !p.multi=pold


return
end











