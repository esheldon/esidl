PRO addsdsshdr, run, camcol, rerun, front=front

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: addsdsshdr, run, camcol, rerun, front=front'
      return
  ENDIF 

  ;; add in sdss header stuff to adatc files
  IF n_elements(front) EQ 0 THEN front='adatc'

  fetch_dir, run, camcol, rerun, dir, corrdir=corrdir

  fetch_file_list, dir, tsfiles, tsfnums
  fetch_file_list, corrdir, corrfiles, fnums, front=front

  nf = n_elements(fnums)

  time=systime(1)

  FOR i=0, nf-1 DO BEGIN 
      
      field = fnums[i]

      outfile = corrfiles[i]
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
      IF (i MOD 20) EQ 0 THEN BEGIN
          print,'------------------------------'
          print,'Field: ',ntostr(field)
          print,'FIELD ', SXPAR(hdr1, 'FIELD')
          print, 'SEEING_R ', SXPAR(hdr1, 'SEEING_R')
          print,'------------------------------'
      ENDIF 

      MODFITS,   outfile, data, hdr0, exten=0 

      FXHMODIFY, outfile, 'FIELD',    SXPAR(hdr1, 'FIELD'), exten=1
      FXHMODIFY, outfile, 'SEEING_U', SXPAR(hdr1, 'SEEING_U'), exten=1
      FXHMODIFY, outfile, 'SEEING_G', SXPAR(hdr1, 'SEEING_G'), exten=1
      FXHMODIFY, outfile, 'SEEING_R', SXPAR(hdr1, 'SEEING_R'), exten=1
      FXHMODIFY, outfile, 'SEEING_I', SXPAR(hdr1, 'SEEING_I'), exten=1
      FXHMODIFY, outfile, 'SEEING_Z', SXPAR(hdr1, 'SEEING_Z'), exten=1
      
  ENDFOR 
  ptime,systime(1)-time

  return
END 
