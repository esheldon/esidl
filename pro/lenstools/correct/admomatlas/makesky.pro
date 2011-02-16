PRO makesky, run, rerun, camcol

  IF n_params() NE 3 THEN BEGIN
      print,'-Syntax: makesky, run, rerun, camcol'
      return
  ENDIF 

  rstr = ntostr(run)
  cstr = ntostr(camcol)
  colors = ['u','g','r','i','z']
  nclr = 5

  fetch_dir, run, camcol, rerun, dir, atldir,corratldir=skydir
  fetch_file_list,dir,ff,fnums, fieldmax=fieldmax, fieldmin=fieldmin

  sdssidl_setup

  nfield = n_elements(fnums)

  s=create_struct('frame', 0, $
                  'sky_int', fltarr(5), $
                  'sky_raw', fltarr(5), $
                  'sky_slope', fltarr(5), $
                  'sky_sig', fltarr(5), $
                  'skypluserr', fltarr(5), $
                  'skyminerr', fltarr(5), $
                  'gain', fltarr(5) )
  struct = replicate(s, nfield)
  print,'Run: ',ntostr(run), ' Rerun: ',ntostr(rerun),' Camcol: ',$
        ntostr(camcol)
  FOR i=0, nfield-1 DO BEGIN

      field = fnums[i]
      IF ( (i MOD 20) EQ 0 ) OR (i EQ nfield-1)  THEN BEGIN
          print,'field = '+ntostr(field)
      ENDIF 
      psfield_name, run, camcol, field, fname

      t=mrdfits(atldir+fname, 6, /silent)

      ; copy stuff we want for each bandpass
      struct[i].frame = field
      struct[i].sky_int = t.sky
      struct[i].sky_raw = t.sky
      struct[i].sky_slope = t.skyslope
      struct[i].sky_sig = t.skysig
      struct[i].skypluserr = t.sky+t.skyerr
      struct[i].skyminerr = t.sky - t.skyerr
      struct[i].gain = t.gain

  ENDFOR 
  message = ' frame      sky_int  sky_raw    sky_slope   sky_sig    sky+error   sky-error  gain'
  fmt='(I7, 7F15.5)'
  FOR clr=0, nclr-1 DO BEGIN 
      skyfile = skydir+'Sky-'+rstr+'-'+cstr+'-'+colors[clr]+'.dat'
      openw, lun, skyfile, /get_lun
      !textunit = lun

      forprint,struct.frame, struct.sky_int[clr], struct.sky_raw[clr], $
               struct.sky_slope[clr], struct.sky_sig[clr], $
               struct.skypluserr[clr], struct.skyminerr[clr],$
               struct.gain[clr], TEXT=5, F=fmt, /silent

      free_lun, lun
  ENDFOR 

  return
END           
