PRO addsdsshdr2, run, camcol, rerun, front=front

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: addsdsshdr, run, camcol, rerun, front=front'
      return
  ENDIF 

  ;; This is different from addsdsshdr in that it deals with
  ;; case where headers are too big to use modfits

  ;; add in sdss header stuff to adatc files
  IF n_elements(front) EQ 0 THEN front='adat'
  IF front EQ 'adat' THEN nchar = 24
  fixfront = front+'new'

  print
  print,'Outputs files will start with ',fixfront
  print

  fetch_dir, run, camcol, rerun, dir, corrdir=corrdir

  fetch_file_list, dir, tsfiles, tsfnums
  fetch_file_list, corrdir, corrfiles, fnums, front=front, nchar=nchar
  rename_tsobj, tsfiles, corrdir, fixfront, fixfiles

  nf = n_elements(fnums)

  time=systime(1)

  FOR i=0, nf-1 DO BEGIN 
      
      field = fnums[i]

      corrf = corrfiles[i]
      infile = tsfiles[i]

      hdr0=headfits(infile, exten=0)
      IF i EQ 0 THEN BEGIN
          print
          print, 'RUN ',      SXPAR(hdr0, 'RUN')
          print, 'CAMCOL ',   SXPAR(hdr0, 'CAMCOL')
          print, 'RERUN ',    SXPAR(hdr0, 'RERUN')
          print, 'STRIPE ',   SXPAR(hdr0, 'STRIPE')
          print, 'STRIP ',    SXPAR(hdr0, 'STRIP')
          print, 'PHOTO_VER', SXPAR(hdr0, 'PHOTO_VER')
          print, 'GAIN ',     SXPAR(hdr0, 'GAIN')
          print
      ENDIF 
 
      hdr1 = headfits(infile, exten=1)
      fxhclean, hdr1
      IF (i MOD 20) EQ 0 THEN BEGIN
          print,'------------------------------'
          print,'Field: ',ntostr(field)
          print,'FIELD ', SXPAR(hdr1, 'FIELD')
          print, 'SEEING_R ', SXPAR(hdr1, 'SEEING_R')
          print,'------------------------------'
      ENDIF 

      openr, lun, corrf, /get_lun
      IF i EQ 0 THEN BEGIN 
          tmp = mrdfits3(lun, 1, 0, /silent)
      ENDIF ELSE BEGIN
          tmp = mrdfits3(lun, 1, 0, /silent, /deja_vu)
      ENDELSE 
      free_lun, lun
      
      outfile = fixfiles[i]
      mwrfits2, tmp, outfile, hdr1, /create, hdr0=hdr0

      delvarx,tmp
      
  ENDFOR 
  ptime,systime(1)-time

  return
END 
