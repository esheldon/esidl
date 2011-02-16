PRO display_pixeljack, jackknife_ids, pixnums, pngfile=pngfile

  IF n_params() NE 2 THEN BEGIN 
      print,'-Syntax: display_pixeljack, jackknife_ids, pixnums, pngfile='
      return
  ENDIF 

  IF n_elements(pngfile) NE 0 THEN BEGIN 
      setupplot,'z'
      device, resolution = [1024,768]
  ENDIF 

  nj = n_elements(jackknife_ids) & np = n_elements(pixnums)
  IF nj NE np THEN message,'jackknife_ids and pixnums must be same size'

  simpctable, rct, gct, bct, colorlist=colorlist
  nc = n_elements(colorlist)

  rmd = rem_dup(pixnums)

  display_pixel, pixnums[rmd], /iso, resolution=256

  h = histogram(jackknife_ids[rmd], rev=rev, min=0)
  nh = n_elements(h)
  FOR i=0L, nh-1 DO BEGIN 

      IF rev[i] NE rev[i+1] THEN BEGIN 

          w = rev[ rev[i]:rev[i+1]-1 ]

          display_pixel, pixnums[rmd[w]], $
            color=colorlist[ i MOD nc ], /over_plot, resolution=256

      ENDIF 

  ENDFOR 

  IF n_elements(pngfile) NE 0 THEN BEGIN 
      write_png, pngfile, tvrd(), rct, gct, bct
      setupplot,'X'
  ENDIF 

END 
