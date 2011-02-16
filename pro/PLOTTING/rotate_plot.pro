PRO rotate_plot, image, gridx, gridy, $
                 step=step, shade=shade, psfile=psfile, _extra=e

  IF n_params() EQ 0 THEN BEGIN 
      print,'-Syntax: rotate_plot, z [, x, y, step=step, shade=shade, psfile=psfile, _extra=e ]'
      return
  ENDIF 

  IF n_elements(step) EQ 0 THEN step=10 ;degrees

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Build the command string
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(shade) THEN command='shade_surf' ELSE command='surface'
  string = command+', image'
  IF n_elements(gridx) NE 0 THEN string = string+', gridx'
  IF n_elements(gridy) NE 0 THEN string = string+', gridy'
  string = string + ', ax=ax, az=az, _extra=e'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up the keys
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  up = '[A'
  down = '[B'
  left = '[D'
  right = '[C'

  ang0 = 30
  z=0
  x=0
  continue=1
  WHILE continue DO BEGIN 
      ax = ang0 + x*step
      az = ang0 + z*step

;      call_procedure, command, z, ax=ax, az=az, _extra=e
      r=execute(string)
      go=1
      WHILE go DO BEGIN
          getseq,key
          CASE key OF
              up:    BEGIN & x=x-1 & go=0 & END 
              down:  BEGIN & x=x+1 & go=0 & END 
              left:  BEGIN & z=z+1 & go=0 & END 
              right: BEGIN & z=z-1 & go=0 & END 
              '':    BEGIN & continue=0 & go=0 & END 
              ELSE: go=1
          ENDCASE 
      ENDWHILE 
  ENDWHILE 
      
  IF n_elements(psfile) NE 0 THEN BEGIN
      makeps, psfile, /noland
;      call_procedure, command, z, ax=ax, az=az, _extra=e
      r=execute(string)
      ep
  ENDIF 

  return
END 
