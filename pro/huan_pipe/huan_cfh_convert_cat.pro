PRO huan_cfh_convert_cat_read, file, struct

  ss=create_struct('number', 0L, $
                   'mag_auto', 0.0, $
                   'magerr_auto', 0.0, $
                   'kron_radius', 0.0, $
                   'x_image', 0.0, $
                   'y_image', 0.0, $
                   'a_image', 0.0, $
                   'b_image', 0.0, $
                   'theta_image', 0.0, $
                   'fwhm_image', 0.0, $
                   'flags', 0L, $
                   'class_star', 0.0, $
                   'x_world', 0.0, $
                   'y_world', 0.0, $
                   'mag_aper', fltarr(5), $
                   'magerr_aper', fltarr(5))

  
  nhead = 16

  nl = numlines(file)
  nobj = nl - nhead
  struct = replicate(ss, nobj)

  openr, lun, file, /get_lun

  ;; skip lines
  tmp = ' '
  FOR i=0L, nhead-1 DO BEGIN 
      readf, lun, tmp
      print,tmp
      tmp = ' '
  ENDFOR 
  readf, lun, struct

  free_lun, lun

END 

PRO huan_cfh_convert_cat

  ;; read the file
  indir = '~/Huan/cfh/cat/'
  outdir = indir

  file = indir + '0920dpIAB1.cat'

  huan_cfh_convert_cat_read, file, struct

  ;; First output a fits version of this catalog
  outfile = outdir + '0920dpIAB1.fit'
  print,'Outputting file: ',outfile
  mwrfits, struct, outfile, /create

  ;; The id is actually chip*10000 + id
  ;; Use this to split catalog into individual catalogs
  ;; And output as fits

  ;; The number offset is 1 in sextractor

  FOR chip=0,11 DO BEGIN 

      outfile = outdir + '0920dpIAB1_'+strn(chip,length=2,padchar='0')+'.fit'

      addnum = chip*10000L
      lownum = 1L+addnum
      highnum = 10000L+addnum
      
      print,lownum,highnum

      w=where(struct.number GE lownum AND struct.number LE highnum, ncat)

      print
      print,'number range = ['+ntostr(lownum)+', '+ntostr(highnum)+']'
      print,'Found '+ntostr(ncat)+' in this range'
      print,'range of numbers: ['+ntostr(min(struct[w].number))+', '+ntostr(max(struct[w].number))+']'
      print,'Outputting file: ',outfile
      mwrfits, struct[w], outfile, /create

  ENDFOR 

END 
