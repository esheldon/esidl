PRO class_fchart, cat, radius=radius, contrast=contrast, saturation=saturation

  make_clflag_struct, cl
  cl.failed = 'N'
  cl.ellip = 'Y'


  clflag_select, cat, cl, good
  ngood = n_elements(good)
  print
  print,'Found '+ntostr(ngood)+' with good classifications'
  print
  ;ngood=n_elements(cat)
  ;good = lindgen(ngood)

  ELLIP = 2L^0
  ELLIP_LIKELY = 2L^1
  SPIRAL_LIKELY = 2L^2
  SPIRAL = 2L^3

  IF display_type() EQ 'PS' THEN noframe=1 ELSE noframe=0
  FOR i=0L, ngood-1 DO BEGIN 
      
      ii=good[i]
      
      class = cat[ii].classification
      CASE 1 OF
          ((class AND ELLIP) NE 0): typestr = 'ELLIP'
          ((class AND ELLIP_LIKELY) NE 0): typestr = 'ELLIP_LIKELY'
          ((class AND SPIRAL) NE 0): typestr = 'SPIRAL'
          ((class AND SPIRAL_LIKELY) NE 0): typestr = 'SPIRAL_LIKELY'
          ELSE: typestr= '????'
      ENDCASE 

      ;rgbfchart, cat[ii].ra, cat[ii].dec, radius=radius, contrast=contrast

      fetch_dir, cat[ii].run, cat[ii].camcol, cat[ii].rerun, dir, atldir
      get_atlas, cat, ii, dir=atldir, img=img, imr=imr, imi=imi, /nodisplay

      isig   =    8.07816
      rsig   =    6.71133
      gsig   =    6.03404
      isky = 1000.
      rsky = 1000.
      gsky = 1000.

      rgbview2, imi, imr, img, $
                contrast=contrast, /noprompt, saturation=saturation,$
                rsky=isky, gsky=rsky, bsky=gsky,$
                rsig=isig, gsig=rsig, bsig=gsig,noframe=noframe

      ;legend, typestr, /right, box=0,color=!white
      ;xyouts, 0.9,0.9, typestr, /normal
      print
      print,'Class: '+typestr
      print

      IF display_type() EQ 'X' THEN BEGIN
          key=get_kbrd(1)
          IF key EQ 'q' THEN return
      ENDIF 
  ENDFOR 

END 
