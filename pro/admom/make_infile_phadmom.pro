PRO phadmom_select, select_clr, pstruct, maxmag, _ss

  make_flag_struct,fs
  fs.satur='N'
  fs.satur_center='N'
  fs.bright = 'N'
  fs.nopetro_big = 'N'
  flag_select,pstruct,fs,select_clr,_ss1
  IF (_ss1[0] EQ -1) THEN BEGIN
      print,'No objects passed flag cuts!'
      GOTO,jump:
  ENDIF 
  
  make_flag_struct, fs
  fs.blended='N'
  fs.moved='N'
  fs.deblended_as_psf = 'N'
  flag_select,pstruct[_ss1],fs,select_clr,_ss2,/objc
  IF (_ss2[0] EQ -1) THEN BEGIN
      print,'No objects passed flag cuts!'
      _ss1 = _ss2
      GOTO,jump:
  ENDIF 
  _ss1 = _ss1[_ss2]
  
  _ss2 = where( ( pstruct[_ss1].petrocounts[select_clr] $
                  - pstruct[_ss1].reddening[select_clr] ) LE  maxmag, nss2)
  IF nss2 EQ 0 THEN BEGIN
      print,'No objects passed flag cuts!'
      _ss1 = _ss2
      GOTO,jump:       
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


PRO make_infile_phadmom, run, rerun, camcol, start=start,nframes=nframes

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: make_infile_phadmom, run, rerun, camcol, start=start,nframes=nframes'
      return
  ENDIF 
  g=1 & r=2 & i=3
  rstr = run2string(run)
  rrstr = ntostr(rerun)
  cstr = ntostr(camcol)
  

  fetch_dir, run, camcol, rerun, dir,atldir,corrdir=corrdir
  fetch_file_list, dir, files, fnums, start=start, nframes=nframes

  nf = n_elements(files)
  f1str = field2string(fnums[0])
  f2str = field2string(fnums[nf-1])

  outfile = corrdir + 'admom_'+rstr+'-'+cstr+'-'+rrstr+'-'+f1str+'-'+f2str+'.dat'
  print
  print,'Output file: ',outfile

  openw, lun, outfile, /get_lun

  taglist = ['run','rerun','camcol','field','id', 'colc', 'rowc', 'petrorad']

  FOR fi=0L, nf-1 DO BEGIN 

      field=fnums[fi]
      ;; sky structure
      psField_dir = atldir
      psField_name, run, camcol, field, fname
      skystr = mrdfits(psField_dir+fname, 6, /silent)

      ;; photo structure
      str=0
      read_tsobj, dir, str, start=field, nframes=1, tsobjstr=tsobjstr

      phadmom_select, 2, str, 25., _ss

;      get_atlas, str(_ss), 0, clr=1, dir=atldir, /nodisplay,/noprompt,$
;        row0=row0,col0=col0
;      print,'xcen = '+ntostr(str[_ss[0]].colc[1]-col0)+' ycen = '+ntostr(str[_ss[0]].rowc[1]-row0)

      nobj = n_elements(_ss)
      FOR oi=0L, nobj-1 DO BEGIN 

          ind = _ss[oi]
          printf, lun, $
                  str[ind].run, str[ind].rerun, str[ind].camcol, $
                  str[ind].field, str[ind].id,$
                  str[ind].colc[g], str[ind].rowc[g], $
                  str[ind].colc[r], str[ind].rowc[r], $
                  str[ind].colc[i], str[ind].rowc[i], $
                  str[ind].petrorad[g], str[ind].petrorad[r], str[ind].petrorad[i],$
                  skystr.skysig[g], skystr.skysig[r], skystr.skysig[i], $
                  format='(5(I0,:,1X),12(F0,:,1X))'
      ENDFOR 

  ENDFOR 

  jump:
  close,lun
  free_lun,lun

END 
