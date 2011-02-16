PRO makeskyfit, run, rerun, camcol

  IF n_params() NE 3 THEN BEGIN
      print,'-Syntax: makeskyfit, run, rerun, camcol'
      return
  ENDIF 

  rstr = ntostr(run)
  cstr = ntostr(camcol)
  colors = ['u','g','r','i','z']
  nclr = 5

  sdssidl_setup
  fetch_dir, run, camcol, rerun, dir, atldir,corratldir=skydir
  fetch_file_list,dir,ff,fnums, fieldmax=fieldmax, fieldmin=fieldmin

  skyname, run, camcol, skyfile
  skyfile = skydir + skyfile
  print
  print,'Output file name: ',skyfile
  print

  nfield = n_elements(fnums)

  s=create_struct('field', 0, $
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

      IF datatype(t) NE 'INT' THEN BEGIN  
          struct[i].field = field
          struct[i].sky_int = t.sky
          struct[i].sky_raw = t.sky
          struct[i].sky_slope = t.skyslope
          struct[i].sky_sig = t.skysig
          struct[i].skypluserr = t.sky+t.skyerr
          struct[i].skyminerr = t.sky - t.skyerr
          struct[i].gain = t.gain
      ENDIF ELSE print,'File '+fname+' is missing!'

  ENDFOR 

  mwrfits, struct, skyfile, /create

  return
END           
