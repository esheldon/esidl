PRO getgreys, data, colors, reverse=reverse

  IF n_params() LT 1 THEN BEGIN  
      print,'-Syntax: getgreys, data, colors, reverse=reverse'
      return
  ENDIF 

  IF !d.n_colors LE 256 THEN BEGIN 

      ;; assume simpctable has been run
      max = !grey75
      min = !grey0

      colors = long( arrscl(data, min, max) )

      IF keyword_set(reverse) THEN colors=reverse(temporary(colors))

  ENDIF ELSE BEGIN 
      
      max = 255
      min = 75

      R = long(arrscl( data, min, max ))
      IF keyword_set(reverse) THEN R=reverse(temporary(R))

      G = R
      B = R

      colors = R + 256L*(G+256L*B)
  ENDELSE 

END 
