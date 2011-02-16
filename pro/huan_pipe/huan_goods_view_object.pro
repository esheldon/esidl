PRO huan_goods_view_object, image, cat, w

  OBJECT2_AMOMENT_FAINT =  '200000'X ; too faint for adaptive moments 
  OBJECT2_AMOMENT_UNWEIGHTED = '200000'X ; failed so tried unweighted mom
  OBJECT2_AMOMENT_SHIFT =  '400000'X ; centre moved too far while
                                ;determining adaptive moments 
  OBJECT2_AMOMENT_MAXITER = '800000'X ; Too many iterations while

  FAINT_STR = 'FAINT'
  UNW_STR = 'UNWHT'
  SHIFT_STR = 'SHIFT'
  MAXITER_STR = 'MAXITER'

  sz = size(image, /dim)
  nx = sz[0]
  ny = sz[1]

  nw = n_elements(w)

  subim_sz = 20

  FOR i=0L, nw-1 DO BEGIN 

      ii = w[i]

      wh = cat[ii].whyflag
      print
      IF (wh AND OBJECT2_AMOMENT_FAINT) NE 0 THEN print,FAINT_STR
      ;;IF (wh AND OBJECT2_AMOMENT_UNWEIGHTED) NE 0 THEN print,UNW_STR
      IF (wh AND OBJECT2_AMOMENT_SHIFT) NE 0 THEN print,SHIFT_STR
      IF (wh AND OBJECT2_AMOMENT_MAXITER) NE 0 THEN print,MAXITER_STR


      x = cat[ii].x_sect - 1.0
      y = cat[ii].y_sect - 1.0

      xmin = (x-subim_sz) > 0
      xmax = (x+subim_sz) < (nx-1)

      ymin = (y-subim_sz) > 0
      ymax = (y+subim_sz) < (ny-1)
      

      tim = image[ xmin:xmax, ymin:ymax ]

      title = 'i_AB = '+ntostr(cat[ii].mag_auto)
      title = title + ' S/N = '+ntostr(cat[ii].flux_auto/cat[ii].fluxerr_auto)
      tvasinh, tim, title=title

      key=prompt_kbrd("Hit a key")
      IF strlowcase(key) EQ 'q' THEN return

  ENDFOR 


END 
