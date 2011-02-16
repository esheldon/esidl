PRO shapelet_open_png, struct, clr

  ;; open the png files output by the shapelet
  ;; wrapper code
  nstruct = n_elements(struct)
  message = "Hit a key: space: next  p: previous  q: quit"
  i=0
  WHILE i LE nstruct-1 DO BEGIN 

      w=where(!run_status.run EQ struct[i].run, nm)
      IF nm EQ 0 THEN message,'Cannot find this run!'
      
      stripe = !run_status[w[0]].stripe
      sstr = stripe2string(stripe)
      dir = '~/shapelet_outputs/stripe'+sstr+'_images/'
      
      name = shapelet_decomp_png_name(struct[i], clr, struct[i].dtype, $
                                          struct[i].nmax)
      
      file = dir+name
      
      im=read_image(file)
      
      erase
      tv,im

      IF nstruct GT 1 THEN BEGIN 

          key=prompt_kbrd(message)
          
          key = strlowcase(key)
          CASE key OF
              'q': return
              'p': i= (i-1) > 0
              ELSE: i=i+1
          ENDCASE 
      ENDIF ELSE i=i+1
  ENDWHILE 


END 
