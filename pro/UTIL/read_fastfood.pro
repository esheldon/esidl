;+
;  read_fastfood, filename [, timestep, infostruct=infostruct, x=x, y=y, z=z,
;                 vx=vx, vy=vy, vz=vz]'
;
; Files are f77_unformatted.  Must be opened with that keyword to openr
; Layout of this file (ignoring the extra record info which is skipped
; for use by using openr with the /f77_unformatted keyword):
;
;   header: 
;     integer record array length 5
;        ? 
;        number of particles 
;        ? 
;        ? 
;        ?
;     float record array length 9
;        size of box
;        8 numbers all zero
;
;     float record length 1
;        zoutput
;
;   data:
;     6 float records, length np
;        x
;        y
;        z
;        vx
;        vy
;        vz
;
; This is repeated for each timestep
;
;-

PRO rff_skip2timestep, lun, timestep

  ;; Note, this will not skip the first timestep
  FOR i=0, timestep-1 DO BEGIN 
      ;; Read the first header chunk
      idat = rff_readlongrecord(lun, 5)
      nparticles2skip = idat[1]
      ;; skip next header chunks
      rff_skip4byterecord, lun, 9
      ;; skip zoutput
      rff_skip4byterecord, lun, 1
      
      ;; skip records
      FOR i=1,6 DO rff_skip4byterecord, lun, nparticles2skip
  ENDFOR 

END 
PRO rff_skip4byterecord, lun, np
  ;; Get current position in the file
  nbytes = 4ULL
  point_lun, -lun, pos

  ;; 2 is for the leading int and trailing int (what are they?)
  record_length = (2ULL+ulong64(np))*nbytes
  point_lun, lun, pos + record_length
END 
FUNCTION rff_readfloatrecord, lun, np
  tmp = fltarr(np)
  readu, lun, tmp
  return,tmp
END 
FUNCTION rff_readlongrecord, lun, np
  tmp = longarr(np)
  readu, lun, tmp
  return,tmp
END 

PRO read_fastfood, file, timestep, infostruct=infostruct, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: read_fastfood, filename [, timestep, infostruct=infostruct, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz]'
      print
      print,'Just send the keywords you wish to be read from the file'
      print,'Send timestep to retrieve a particular output timestep'
      print
      message,'Halting'
  ENDIF 

  xe = arg_present(x)
  ye = arg_present(y)
  ze = arg_present(z)
  vxe = arg_present(vx)
  vye = arg_present(vy)
  vze = arg_present(vz)

  num_present = xe+ye+ze+vxe+vye+vze

  print,'Opening file: ',file
  openr, lun, file, /get_lun, /f77_unformatted


  ;; Header data
  idat = lonarr(5)
  fdat = fltarr(9)


  ;; Are we skipping timesteps?
  nt = n_elements(timestep)
  IF nt NE 0 THEN BEGIN 
      IF nt GT 1 THEN message,'Can only read one time step'
      IF timestep LT 0 THEN message,'Timestep must be positive'
      IF timestep GT 0 THEN BEGIN 
          print,'Skipping to timestep '+ntostr(timestep)
          rff_skip2timestep, lun, timestep      
      ENDIF 
  ENDIF

  ;; Reading this time step
  readu, lun, idat
  readu, lun, fdat

  colprint,idat
  print
  colprint,fdat
  
  zoutput = 0.
  readu, lun, zoutput

  ;; The max val for x/y/z. Note, this
  ;; is not the physical max val, which
  ;; is box_size
  maxval = idat[0]
  np  = idat[1]
  box_size = fdat[0]
  print,'  nparticles = ',np
  print,'  Box size = ',box_size
  print,'  zoutput = ',zoutput

  infostruct = {nparticles: np, $
                box_size: box_size, $
                zoutput: zoutput}



  IF num_present NE 0 THEN BEGIN 

      scale = (box_size/maxval)
      IF NOT xe THEN BEGIN 
          rff_skip4byterecord, lun, np
      ENDIF ELSE BEGIN 
          print,'Reading x'
          x = rff_readfloatrecord(lun, np)
          x = x*scale
      ENDELSE 
      IF NOT ye THEN BEGIN 
          rff_skip4byterecord, lun, np
      ENDIF ELSE BEGIN 
          print,'Reading y'
          y = rff_readfloatrecord(lun, np)
          y = y*scale
      ENDELSE 
      IF NOT ze THEN BEGIN 
          rff_skip4byterecord, lun, np
      ENDIF ELSE BEGIN 
          print,'Reading z'
          z = rff_readfloatrecord(lun, np)
          z = z*scale
      ENDELSE 


      IF NOT vxe THEN BEGIN 
          rff_skip4byterecord, lun, np
      ENDIF ELSE BEGIN 
          print,'Reading vx'
          vx = rff_readfloatrecord(lun, np)
      ENDELSE 
      IF NOT vye THEN BEGIN 
          rff_skip4byterecord, lun, np
      ENDIF ELSE BEGIN 
          print,'Reading vy'
          vy = rff_readfloatrecord(lun, np)
      ENDELSE 
      IF NOT vze THEN BEGIN 
          rff_skip4byterecord, lun, np
      ENDIF ELSE BEGIN 
          print,'Reading vz'
          vz = rff_readfloatrecord(lun, np)
      ENDELSE 

  ENDIF 

  free_lun, lun
  
END 
