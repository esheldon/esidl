pro fitmom, run, rerun, camcol, maxmag=maxmag,$
            psFieldrerun=psFieldrerun, start=start, nframes=nframes,$
            status=status

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;    FITMOM
;
; PURPOSE:
;    Finds moments of the psf at the position of each object by
;    using the eigentemplates in the psField files
; 
; 
; INPUTS:  
;    run, rerun, camcol:
;    maxmag=maxmag: the maximum magnitude of objects for which to compute
;                   the psf info.  default is 25.
; Outputs: Fits file containing the relevant info.
;
; Author:  Phil Fischer
; Date: 1/29/99
; Converted to new file formats. Added comments and cleaned up
;                                     2/20/2000 E.S.S.
; Converted to using eigentemplates and outputting fits file
;                                     26-OCT-2000 E.S.S.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; status is 1 unless we reach end
  status=1
  
  IF n_params() LT 3 THEN BEGIN ;Help message
      print,'-syntax fitmom2, run, rerun, camcol, maxmag=maxmag, psFieldrerun=psFieldrerun,start=start, nframes=nframes, status=status'
      return 
  ENDIF 
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(psFieldrerun) EQ 0 THEN psFieldrerun=rerun

  IF n_elements(maxmag) EQ 0 THEN maxmag=25.

  newfront = 'adat'             ;The front of files with adaptive moments
;  newfront = 'adatc'

  setup_mystuff
  bindir = !MYBINDIR
  colors=['u','g','r','i','z']

  bands = [1,2,3]               ; we will process g,r,i
  nband = n_elements(bands)
  ofiles = strarr(nband)
  
  runstr = ntostr(run)
  rrstr = ntostr(rerun)
  cstr = ntostr(camcol)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Setup the column
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  fetch_dir, run, camcol, rerun, dir, atldir, corrdir=corrdir,$
    corratldir = fitdir
  fetch_file_list, dir, files, fnums, start=start, nframes=nframes, $
                   fieldmin=fieldmin, fieldmax=fieldmax

  nfields = n_elements(files)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; set up output structures
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  arrval=fltarr(5)
  newstr = create_struct('run', 0,$
                         'rerun',0,$
                         'camcol',0,$
                         'field',0,$
                         'id',0,$
                         'psfixx',arrval,$
                         'psfiyy',arrval,$
                         'psfixy',arrval,$
                         'psfrho4',arrval)
                         
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up new filenames. 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  rename_tsobj, files, corrdir, newfront, adatfiles, nchar

  print,'-----------------------------------------------------'
  print,' Run: ',runstr,' Rerun: ',rrstr,' Camcol: ',cstr
  print,'-----------------------------------------------------'
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; get psField directory. If psFieldrerun = rerun, then use atldir
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF psFieldrerun EQ rerun THEN psField_dir = atldir $
  ELSE fetch_dir, run, camcol, psFieldrerun, ddd, psField_dir

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over fields
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;; read in admom info
  admomatlas_infile, run, rerun, camcol, fname_in, fname_out
  fname_out = repstr(fname_out, '.dat', '.fit')
  admomstruct = mrdfits(corrdir + fname_out, 1)

  time = systime(1)
  FOR ic=0, nfields-1 DO BEGIN 
      infile = adatfiles[ic]
      field = fnums[ic]
      fstr=ntostr(field)

      ;; Get the psfield info (g,r,i)
      psfield_name, run, camcol, field, fname
      psfield = psField_dir+fname
      psp=ptrarr(5)
      tmp1 = mrdfits(psfield, 2, /silent)
      tmp2 = mrdfits(psfield, 3, /silent)
      tmp3 = mrdfits(psfield, 4, /silent)
      IF (  ( (size(tmp1))[0] EQ 0 ) OR $
            ( (size(tmp2))[0] EQ 0 ) OR $
            ( (size(tmp3))[0] EQ 0 ) )     THEN BEGIN
          print,'psfield file '+psfield+' missing or corrupted'
          tmp1=0 & tmp2=0 & tmp3=0
          GOTO,jump
      ENDIF 
      psp[1] = ptr_new( tmp1, /no_copy ) ;deletes tmp1
      psp[2] = ptr_new( tmp2, /no_copy )
      psp[3] = ptr_new( tmp3, /no_copy )
      
      pstruct=0
      ;; Fields have already been trimmed of overlap.  Don't need to
      ;; use read_tsobjmin
      openr, lun, infile, /get_lun
      IF ic EQ 0 THEN BEGIN 
          pstruct = mrdfits3(lun, 1, 0, /silent)
      ENDIF ELSE BEGIN
          pstruct = mrdfits3(lun, 1, 0, /silent, /deja_vu)
      ENDELSE 
      free_lun, lun
      IF ( size(pstruct) )[0] EQ 0 THEN GOTO,jump

      wad = where( admomstruct.field EQ field, ng)

      IF ng NE 0 THEN BEGIN 
          wad2=where( (admomstruct[wad].ixx[1] NE 0.) OR $
                      (admomstruct[wad].ixx[2] NE 0.) OR $
                      (admomstruct[wad].ixx[3] NE 0.) ,ng)
          IF ng NE 0 THEN wad = wad[wad2]
      ENDIF 
      print
      print,'Run: ',runstr,' Camcol: ',cstr,' Field: ',fstr,' Nobj: ',ntostr(ng)

      IF ng NE 0. THEN BEGIN 

          ;; match up the structs
          sphoto_match, admomstruct[wad], pstruct, m_ad, wg

          psfstr = replicate(newstr, ng)
          psfstr.run = pstruct[wg].run
          psfstr.rerun = pstruct[wg].rerun
          psfstr.camcol = pstruct[wg].camcol
          psfstr.field = pstruct[wg].field
          psfstr.id = pstruct[wg].id

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Loop over objects
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          
          FOR i=0L, ng-1 DO BEGIN 

              ind = wg[i]
              ind2 = wad[ m_ad[i] ]
              psf_rec,psp,pstruct[ind].rowc,pstruct[ind].colc,psf,bands

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Loop over bandpasses
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              FOR iband = 0, nband-1 DO BEGIN 

                  cindex = bands[iband]

                  ;; only measure psf if this object is
                  ;; measured
                  IF admomstruct[ind2].ixx[cindex] NE 0.0 THEN BEGIN 
                      ixx = 0.
                      iyy = 0.
                      ixy = 0.
                      momerr = 0.
                      rho4 = 0.
                      whyflag = 0L
                      nrow = 51L
                      ncol = 51L
                      tmp = call_external(value = [0B,0B,0B,0B,0B,0B,0B,0B,0B],$
                                          bindir+'psfadmom.so','psfadmom', $
                                          *psf[iband], nrow, ncol,$
                                          ixx, iyy, ixy, momerr, rho4, whyflag)
                      psfstr[i].psfixx[cindex] = ixx
                      psfstr[i].psfiyy[cindex] = iyy
                      psfstr[i].psfixy[cindex] = ixy
                      psfstr[i].psfrho4[cindex] = rho4
                  ENDIF 
              ENDFOR ;Loop over bandpasses     
              ptrarr_free, psf
              IF (i MOD 20) EQ 0 THEN print,'.',format='($,a)'
          ENDFOR                ;Loop over objects

          psfname2, run, rerun, camcol, field, outfile
          mwrfits, psfstr, fitdir + outfile, /create

          psfstr = 0
      ENDIF 
jump:
      ptrarr_free, psp
      pstruct = 0
  ENDFOR                        ;Loop over fields


  ptime,systime(1)-time
  status=0

END 
