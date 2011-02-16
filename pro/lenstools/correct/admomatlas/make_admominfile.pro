
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; Create adatin input file for admomatlas.c
;
; Inputs:  color_index: bandpass to use for galaxy/star selection.
;          run, rerun, camcol:
;          nchar: number of characters in output file name. (optional)
;          nframes:  Optional parameter which tells how many field to read
;          start:    Beginning field
;          maxmag:   maximum magnitude for galaxies in color_index
;                    DEFAULT = 25.
;          stmaxmag: maximum mag for stars.  Should be size=3 for g,r,i
;                    DEFAULT = [21,20,20] A bit too faint for flexibility
; Outputs: fits files containing adaptive moments, with tsObj-* relpaced
;          with adat-*
;
; Author:  Phil Fischer
; Date: 1/20/99
; Erin Scott Sheldon   2/18/2000
;       Adapted to new file formats
;       Now does all bandpasses at once
;       Now extracts stars and sets flag: starflag
;           add 2b^cindex for star in that color
;           This is just a first time through, stars can be selected later
;           with different cuts. Then starflag can be reset. 
;           Sped up, cleaned up.
; 26-OCT-2000 E.S.S. Changed the way sky is found.
;                    Now uses read_tsobj since it is fast enough now.
;                    Simplified structure creation.
;                    Added rotation tag
; 04-JAN-2001 E.S.S. Now outputs adat files without ixx,iyy, etc measured
;                    and creates input file for admomatlas.c
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO admomatlas_select, select_clr, pstruct, maxmag, _ss

  make_flag_struct,fs
  fs.satur='N'
  fs.satur_center='N'
  fs.bright = 'N'
  fs.nopetro_big = 'N'
  flag_select,pstruct,fs,select_clr,_ss1
  IF (_ss1[0] EQ -1) THEN BEGIN
      print,'No objects passed flag cuts!'
      GOTO,jump
  ENDIF 
  
  make_flag_struct, fs
  fs.moved='N'
  fs.deblended_as_psf = 'N'
  flag_select,pstruct[_ss1],fs,select_clr,_ss2,/objc
  IF (_ss2[0] EQ -1) THEN BEGIN
      print,'No objects passed flag cuts!'
      _ss1 = _ss2
      GOTO,jump
  ENDIF 
  _ss1 = _ss1[_ss2]
  
  _ss2 = where( ( pstruct[_ss1].petrocounts[select_clr] $
                  - pstruct[_ss1].reddening[select_clr] ) LE  maxmag, nss2)
  IF nss2 EQ 0 THEN BEGIN
      print,'No objects passed flag cuts!'
      _ss1 = _ss2
      GOTO,jump       
  ENDIF 
  _ss1 = _ss1[_ss2]
  
  jump:
  ;; also include any object which has primtarget set
  _ss3 = where(pstruct.primtarget NE 0,nss3)

  ;; any with primtarget set?
  IF nss3 EQ 0 THEN _ss = _ss1 ELSE BEGIN 
      ;; any pass flag cuts?
      IF _ss1[0] EQ -1 THEN _ss=_ss3 ELSE BEGIN 
          ;; combine and remove duplicates
          _ss = [_ss1, _ss3]
          _ss = _ss[rem_dup(_ss)]
      ENDELSE 
  ENDELSE 

  return
END 

pro make_admominfile, select_clr, run, rerun, camcol, nchar, taglist=taglist, $
                      psFieldrerun=psFieldrerun, $
                      nframes=nframes,start=start,maxmag=maxmag,$
                      outdir=outdir, status=status

  ;; status is 1 unless we reach end
  status=1

  IF n_params() LT 4 THEN BEGIN ;Help message
      print,'-syntax admomatlas, select_clr, run, rerun, camcol, nchar, psFieldrerun=psFieldrerun, taglist=taglist, nframes=nframes,start=start,maxmag=maxmag,stmaxmag=stmaxmag, outdir=outdir, status=status'
      return 
  ENDIF 
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(psFieldrerun) EQ 0 THEN psFieldrerun=rerun

  setup_mystuff

  colors=['u','g','r','i','z']
  rstr = ntostr(run)
  cstr = ntostr(camcol)

  ;; color indices for bandpasses
  gg = 1
  rr = 2
  ii = 3

  IF n_elements(maxmag) EQ 0 THEN maxmag=22.5 ;galaxy maxmag in select_clr

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Setup the column for read_tsobj
  ;; and get corrdir
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  make_admomtags, taglist, /default
  fetch_dir, run, camcol, rerun, dir, atldir, corrdir=corrdir, $
    corratldir=skydir           ;We keep sky files in corratldir
  fetch_file_list, dir, files, fnums, start=start, nframes=nframes, $
    fieldmin=fieldmin, fieldmax=fieldmax

  nfields = n_elements(files)

  IF n_elements(outdir) NE 0 THEN corrdir=outdir

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; set up admomatlas.c input file name and open file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  admomatlas_infile, run, rerun, camcol, adfile
  adfile = corrdir + adfile

  openw, lun, adfile, /get_lun

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; get psField directory. If psFieldrerun = rerun, then use atldir
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF psFieldrerun EQ rerun THEN psField_dir = atldir $
  ELSE fetch_dir, run, camcol, psFieldrerun, ddd, psField_dir

  ;; print some stuff

  print,'-----------------------------------------------------'
  print,' column '+ntostr(camcol)
  print,'-----------------------------------------------------'
  print,'max '+colors[select_clr]+'-band mag: ',ntostr(maxmag)
  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; define fields to skip
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  badrun = [745]
  badcol  = [6]
  badfield = [247]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over fields
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time=systime(1)
  FOR ic = 0L, nfields-1 DO BEGIN

      infile = files[ic]
      field = fnums[ic]

      fstr = ntostr(field)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Read the tsObj file for this field
      ;; IMPORTANT: use /noadd: we don't want to sum objc_rowc!!
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      read_tsobj, dir, pstruct, start=field, nframes=1, $
        taglist=taglist, tsobjstr=tsobjstr, $
        /noadd,verbose=0
          
      ;; read_tsobj returns nothing for empty file
      nps = n_elements(pstruct)
      IF (nps NE 0) THEN BEGIN  
          _ss=replicate(-1,1)

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Make sure sky/psf can be measured
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          psField_name, run, camcol, field, fname
          skystr = mrdfits(psField_dir+fname, 6, /silent)
          IF (size(skystr))[0] NE 0 THEN BEGIN

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Need to address flags.
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              admomatlas_select, select_clr, pstruct, maxmag, _ss
              IF _ss[0] NE -1 THEN BEGIN 

                  nps = n_elements(_ss)
                  pstruct = temporary(pstruct[_ss])
      
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; print info
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  print,'Run: ',rstr,' Camcol: ',cstr,' Field: ',fstr,$
                        ' After Cuts: ',ntostr(nps)

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; print to admomatlas file
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


                  FOR oi=0L, nps-1 DO BEGIN ;loop over objects
                      
                      printf, lun, $
                              pstruct[oi].run, pstruct[oi].rerun, $
                              pstruct[oi].camcol, $
                              pstruct[oi].field, pstruct[oi].id,$
                              pstruct[oi].colc[gg], pstruct[oi].rowc[gg], $
                              pstruct[oi].colc[rr], pstruct[oi].rowc[rr], $
                              pstruct[oi].colc[ii], pstruct[oi].rowc[ii], $
                              pstruct[oi].petrorad[gg], $
                              pstruct[oi].petrorad[rr], $
                              pstruct[oi].petrorad[ii],$
                              skystr.skysig[gg], skystr.skysig[rr], $
                              skystr.skysig[ii], $
                              format='(5(I0,:,1X),12(F0,:,1X))'

                  ENDFOR ;; end loop over objects

              ENDIF ELSE print,'No objects passed cuts'
          ENDIF ELSE print,'psField file '+psField_dir+fname+' not found or corrupt'
      ENDIF ;; not empty

      delvarx, pstruct
      skystr=0 & _ss=0 & _ss2=0 & stars=0

  ENDFOR                        ;End loop over fields

  ;; close output file
  free_lun, lun

  ptime,systime(1)-time
  status = 0

  return
END 
