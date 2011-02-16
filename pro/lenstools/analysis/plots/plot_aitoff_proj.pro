PRO plot_aitoff_proj, struct, dops=dops, pscolor=pscolor, $
                      slam=slam, seta=seta, sx=sx, sy=sy

  IF n_elements(slam) NE 0 AND n_elements(sx) EQ 0 THEN BEGIN 
      csurvey2eq, slam, seta, sra, sdec
      myaitoff, sra, sdec, sx, sy
  ENDIF 

  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(pscolor) THEN fend = '_color.eps' ELSE fend = '.eps'
      IF n_elements(sx) NE 0 THEN psfile = '~/plots/spec_source_aitoff'+fend $
      ELSE psfile = '~/plots/spec_aitoff'+fend
      begplot, name=psfile, /color, /encapsulated, xsize = 7, ysize=6
  ENDIF 

  IF !d.name EQ 'PS' THEN BEGIN  
      IF keyword_set(pscolor) THEN ocolor = !red ELSE ocolor = !grey75
  ENDIF ELSE BEGIN 
      ocolor = !red
  ENDELSE 
ocolor=!red
  myaitoff,struct.ra,struct.dec,x,y
  
  aplot, !gratio, [0], xrange = [-100, 100], yrange = [-10, 80], $
         /nodata, xtickname=[' ',' ',' ', ' ', ' '], $
         ytickname=[' ',' ',' ', ' ', ' ',' '],$
         xticklen=1.e-5,yticklen=1.e-5
  IF n_elements(sx) NE 0 THEN BEGIN 
      plotrand, sx, sy, psym=3, fracuse=0.05, /overplot
  ENDIF 
  oplot, x, y, psym=3, color=ocolor

  myaitoff_grid,/label,charsize=2,charthick=2,longrange=[90,240], latrange=[0,60]

;  !p.color = !black
;  myaitoff_grid,/label,charsize=1,longrange=[90,240], latrange=[45,60]
;  !p.color = !white

  IF !d.name EQ 'PS' THEN BEGIN 
      decpos = [-20, 33.9]
      rapos = [-80, 1]
  ENDIF ELSE BEGIN 
      decpos = [-9, 33.9]
      rapos = [-70, 1]
  ENDELSE 

  IF keyword_set(dops) THEN BEGIN 
      endplot
      set_bbox,psfile,"%%BoundingBox: 85 55 480 347"
  ENDIF 

;  xyouts, decpos[0], decpos[1], 'DEC = ',charsize=1
;  xyouts, rapos[0], rapos[1], 'RA = ',charsize=1

END 
