PRO huan_goods_addnoise, type, seeing=seeing

  IF n_elements(seeing) NE 0 THEN BEGIN 
      seestr = '_'+ntostr(seeing, 4, /round)
  ENDIF ELSE BEGIN 
      seestr = ''
  ENDELSE 

  dir = "~/Huan/goods/"
  imdir = dir + 'images/'
  catdir = dir + 'cat/'

  image_sections = [22, 23, 33, 34]
  nsec = n_elements(image_sections)

  ga = obj_new('goods_addnoise',type, seeing=seeing)
  FOR i=0L, nsec-1 DO BEGIN 

      imsec  = image_sections[i]

      catfile = catdir + 'h_goods_si_sect'+ntostr(imsec)+'_r1.0z_cat_mod.fit'
      imfile = imdir + 'h_si_sect'+ntostr(imsec)+'_v1.0_drz_img.fits'

      out_catfile = $
        catdir + type + seestr + $
        '_h_goods_si_sect'+ntostr(imsec)+'_r1.0z_cat_mod.fit'
      out_imfile = $
        imdir +  type + seestr + $
        '_h_si_sect'+ntostr(imsec)+'_v1.0_drz_img.fits'

;      print
;      print,'Input Files: '
;      print,'  '+catfile
;      print,'  '+imfile

      print
      print,'Output Files: '
      print,'  '+out_catfile
      print,'  '+out_imfile
      print

      ga->read_cat, catfile
      ga->read_image, imfile

      IF n_elements(seeing) NE 0 THEN BEGIN 
          ga->convolve
      END 

      ga->addnoise

      print,'Writing outputs'
      ga->write_cat, out_catfile
      ga->write_image, out_imfile

      print,'Freeing memory'
      ga->free

  ENDFOR 

END 
