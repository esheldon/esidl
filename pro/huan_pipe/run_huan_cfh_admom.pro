PRO run_huan_cfh_admom

  info_file = '~/Huan/cfh/image_info.dat'

  readcol, info_file, imlist, catlist, sky, sky_rms, zeropoint, $
           format='A,A,F,F,F', /silent

  nfile = n_elements(imlist)
  ;;nfile = 1

  FOR i=0L, nfile-1 DO BEGIN 

      huan_cfh_admom, imlist[i], catlist[i], sky[i], sky_rms[i]

  ENDFOR 
  
 

END 
