pro fitmom_test, run, rerun, camcol, maxmag=maxmag,$
                    psFieldrerun=psFieldrerun, start=start, nframes=nframes
  
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
  
  IF n_params() LT 3 THEN BEGIN ;Help message
      print,'-syntax fitmom2, run, rerun, camcol, maxmag=maxmag, psFieldrerun=psFieldrerun,start=start, nframes=nframes'
      return 
  ENDIF 
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  COMMON fitmomseed,seed

  IF n_elements(psFieldrerun) EQ 0 THEN psFieldrerun=rerun

  IF n_elements(maxmag) EQ 0 THEN maxmag=25.

  ;; we deleted adatc
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
  print,corrdir+fname_out
  admomstruct = mrdfits(corrdir + fname_out, 1)

  ;; make round gaussian to send to
  ;; admom
  makegauss,psf,[51,51],sqrt(3.),counts=200000.

  time = systime(1)
  FOR ic=0, nframes-1 DO BEGIN 
      infile = adatfiles[ic]
      field = fnums[ic]
      fstr=ntostr(field)

      wad = where( admomstruct.field EQ field, ng)

      IF ng NE 0 THEN BEGIN 
          wad2=where( (admomstruct[wad].ixx[1] NE 0.) OR $
                      (admomstruct[wad].ixx[2] NE 0.) OR $
                      (admomstruct[wad].ixx[3] NE 0.) ,ng)
          IF ng NE 0 THEN wad = wad[wad2]
      ENDIF 

      IF ng NE 0. THEN BEGIN 
          print
          print,'Run: ',runstr,' Camcol: ',cstr,' Field: ',fstr,' Nobj: ',ntostr(ng)
          
          psfstr = replicate(newstr, ng)
          psfstr.run = admomstruct[wad].run
          psfstr.rerun = admomstruct[wad].rerun
          psfstr.camcol = admomstruct[wad].camcol
          psfstr.field = admomstruct[wad].field
          psfstr.id = admomstruct[wad].id

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Loop over objects
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          
          FOR i=0L, ng-1 DO BEGIN 

              ind = wad[i]

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Loop over bandpasses
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              FOR iband = 0, nband-1 DO BEGIN 

                  cindex = bands[iband]

                  ;; only measure psf if this object is
                  ;; measured
                  IF admomstruct[ind].ixx[cindex] NE 0.0 THEN BEGIN 
                      momerr = 0.
                      rho4 = 0.
                      whyflag = 0L
                      nrow = 51L
                      ncol = 51L

                      ;; using weight of object
                      ixx = admomstruct[ind].ixx[cindex]
                      iyy = admomstruct[ind].iyy[cindex]
                      ixy = admomstruct[ind].ixy[cindex]

                      tmp = call_external(value = [0B,0B,0B,0B,0B,0B,0B,0B,0B],$
                                          bindir+'psfadmom_retmom.so',$
                                          'psfadmom_retmom', $
                                          psf, nrow, ncol,$
                                          ixx, iyy, ixy, momerr, rho4, whyflag)

                      psfstr[i].psfixx[cindex] = ixx
                      psfstr[i].psfiyy[cindex] = iyy
                      psfstr[i].psfixy[cindex] = ixy
                      psfstr[i].psfrho4[cindex] = rho4
                  ENDIF 
              ENDFOR ;Loop over bandpasses     

              IF (i MOD 20) EQ 0 THEN print,'.',format='($,a)'
          ENDFOR                ;Loop over objects

          psfname2, run, rerun, camcol, field, outfile
outfile=outfile+'_test'
          mwrfits, psfstr, fitdir + outfile, /create

          psfstr = 0
      ENDIF 
jump:

  ENDFOR                        ;Loop over fields


  ptime,systime(1)-time

END 
