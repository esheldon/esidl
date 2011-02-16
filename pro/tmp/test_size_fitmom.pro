pro test_size_fitmom, run, rerun, camcol, thresh, iorder, nframes=nframes, start=start, stmaxmag=stmaxmag
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;    FITMOM
;
; PURPOSE:
;    Fits polynmial to star moments
; 
; 
; INPUTS:  
;    run, rerun, camcol:
;    thresh: Clipping threshold in sigmas
;    iorder: Order of the fit
;    nframes:  Optional parameter which tells how many frames to read
;    start:    Beginning frame
;    stmaxmag: maximum mag for stars.  Should be size=3 for g,r,i
;                    DEFAULT = [20, 19, 19]
; Outputs: Ascii files for this camcol and each color,  containing fit info 
;          for each field and color.
;
; Author:  Phil Fischer
; Date: 1/29/99
; Converted to new file formats. Added comments and cleaned up
;                                     2/20/2000 E.S.S.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  IF n_params() LT 5 THEN BEGIN ;Help message
      print,'-syntax test_size_fitmom, run, rerun, camcol, thresh, iorder, nframes=nframes, start=start, stmaxmag=stmaxmag'
      print,'Typically, iorder=2, thresh = 3'
      return 
  ENDIF 
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  newfront = 'adat'             ;The front of files with adaptive moments
  newfront = 'adatc'            ;this is because running test

  setup_mystuff
  bindir = !MYBINDIR
  colors=['u','g','r','i','z']

  bands = [1,2,3]               ; we will process g,r,i
  nband = n_elements(bands)
  ofiles = strarr(nband)
  
  runstr = ntostr(run)
  rrstr = ntostr(rerun)
  cstr = ntostr(camcol)

  nmag=n_elements(stmaxmag)
  IF nmag NE nband THEN BEGIN 
      IF nmag EQ 0  THEN BEGIN
          stmaxmag = [20., 19., 19.]
      ENDIF ELSE BEGIN
          print,'stmaxmag must be size ',nband
      ENDELSE 
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Setup the column
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  fetch_dir, run, camcol, rerun, dir, atldir, corrdir=corrdir,$
    corratldir = fitdir
  fetch_file_list, dir, files, fnums, start=start, nframes=nframes, $
                   fieldmin=fieldmin, fieldmax=fieldmax
  ntot = fieldmax-fieldmin+1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; set up output structures
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;put this in
  arrval = dblarr(15)
  newstr = create_struct('field', 0L, $
                         'nstars', 0L, $
                         'szcoeff', arrval, $
                         'szrms', arrval, $
                         'e1coeff', arrval, $
                         'e1rms', arrval, $
                         'e2coeff', arrval, $
                         'e2rms', arrval, $
                         'rho4coeff', arrval, $
                         'rho4rms', arrval)
  gstr = replicate(newstr, nframes)
  rstr = gstr
  istr = gstr
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up new filenames. 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  rename_tsobj, files, corrdir, newfront, adatfiles, nchar

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Open the output files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;!!!!!!!!!!!CHANGE TO .FIT FILE!!!!!!!!!

;  luns = lonarr(nband)
  FOR i=0, nband-1 DO BEGIN
      
      cindex = bands[i]
      momname, run, camcol, cindex, iorder, tmpname
      tmpname=repstr(tmpname,'.txt','.fit')
      ofiles[i] = fitdir+'test'+tmpname

      ;; take this out
;      openw, lun, ofiles[i], /get_lun
;      luns[i] = lun
  ENDFOR 
  print
  colprint,ofiles
  print
  print,'-----------------------------------------------------'
  print,' Run: ',runstr,' Rerun: ',rrstr,' Camcol: ',cstr
  print,'-----------------------------------------------------'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over fields
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  time = systime(1)
  FOR ic=0, nframes-1 DO BEGIN 
      infile = adatfiles[ic]
      field = fnums[ic]
      fstr=ntostr(field)
      print,'Run: ',runstr,' Camcol: ',cstr,' Field: ',fstr
      
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

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Loop over bandpasses
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      

;;;;;;!!!!!!!! Altered to just find size from petrorad. Should
;;;;;;;;;;;;;;;ignore the output shape info.
      FOR iband = 0, nband-1 DO BEGIN 
          
          cindex = bands[iband]
          ;; altered to make sure petroraderr > 0
          stars=where(pstruct.starflag[cindex] EQ 1 AND $
                      pstruct.petroraderr[cindex] GT 0., nstars)
          print,'Band: ',colors[cindex],' Nstars: ',ntostr(nstars)
          a1=dblarr(4,15)
          rms=fltarr(4)
          thresh=double(thresh)

          IF nstars NE 0 THEN BEGIN 
              x=double(pstruct[stars].colc[cindex])
              y=double(pstruct[stars].rowc[cindex])

              ii=indgen(nstars)

              sigg=dblarr(4,nstars)
;              sigg(0,ii)=pstruct(stars[ii]).ixx[cindex]
;              sigg(1,ii)=pstruct(stars[ii]).iyy[cindex]
;              sigg(2,ii)=pstruct(stars[ii]).ixy[cindex]
;              sigg(3,ii)=pstruct(stars[ii]).rho4[cindex]
              ;; going to make it fit size from petrorad
              sigg(0,ii)=(pstruct(stars[ii]).petrorad[cindex])^2/2.
              sigg(1,ii)=(pstruct(stars[ii]).petrorad[cindex])^2/2.
              sigg(2,ii)=0.     ;pstruct(stars[ii]).ixy[cindex]
              sigg(3,ii)=2.

              ntmp = nstars
              CASE 1 OF 
                  (ntmp LT 10): iord=0
                  (ntmp LT 20): iord = iorder < 1
                  (ntmp LT 30): iord = iorder < 2
                  (ntmp LT 50): iord = iorder < 3
                  ELSE: iord = iorder
              ENDCASE 

              iord=long(iord)
              IF (ntmp GT 1)THEN BEGIN 
                  ff=call_external(bindir+'fitmomsdss.so',$
                                   'fitmom_',x,y,sigg,a1,iord,thresh,ntmp,rms)
              ENDIF ELSE BEGIN 
                  ntmp=0
              ENDELSE 
          ENDIF ELSE BEGIN
              print,'No stars for field: ',fstr
              ntmp = 0
          ENDELSE 

          ;; take this out
;          FOR k=0,3,1 DO BEGIN 
;              printf, luns[iband], ofiles[iband], field, $
;                a1[k,0],a1[k,1],a1[k,2],a1[k,3],a1[k,4], $ 
;                a1[k,5],a1[k,6],a1[k,7],a1[k,8],a1[k,9],$
;                a1[k,10],a1[k,11],a1[k,12], $
;                a1[k,13],a1[k,14],rms[k],ntmp,$
;                format='(a20,i5,15e16.8,f10.3,i5)'
;          ENDFOR 
             
          ;; put this in 
          CASE cindex OF 
              1: BEGIN          ;g-band
                  gstr[ic].field = field
                  gstr[ic].nstars = ntmp
                  gstr[ic].szcoeff = a1[0,*]
                  gstr[ic].szrms   = rms[0]
                  gstr[ic].e1coeff = a1[1,*]
                  gstr[ic].e1rms   = rms[1]
                  gstr[ic].e2coeff = a1[2,*]
                  gstr[ic].e2rms   = rms[2]
                  gstr[ic].rho4coeff = a1[3,*]
                  gstr[ic].rho4rms = rms[3]
              END
              2: BEGIN        ;r-band 
                  rstr[ic].field = field
                  rstr[ic].nstars = ntmp
                  rstr[ic].szcoeff = a1[0,*]
                  rstr[ic].szrms   = rms[0]
                  rstr[ic].e1coeff = a1[1,*]
                  rstr[ic].e1rms   = rms[1]
                  rstr[ic].e2coeff = a1[2,*]
                  rstr[ic].e2rms   = rms[2]
                  rstr[ic].rho4coeff = a1[3,*]
                  rstr[ic].rho4rms = rms[3]
              END 
              3: BEGIN          ;i-band
                  istr[ic].field = field
                  istr[ic].nstars = ntmp
                  istr[ic].szcoeff = a1[0,*]
                  istr[ic].szrms   = rms[0]
                  istr[ic].e1coeff = a1[1,*]
                  istr[ic].e1rms   = rms[1]
                  istr[ic].e2coeff = a1[2,*]
                  istr[ic].e2rms   = rms[2]
                  istr[ic].rho4coeff = a1[3,*]
                  istr[ic].rho4rms = rms[3]
             END 
              ELSE: BEGIN 
                  print,'What!'
                  return
              END 
          ENDCASE 

      ENDFOR                    ;Loop over bandpasses
  ENDFOR                        ;Loop over fields
  
  ptime,systime(1)-time

  ;; take this out
;  FOR iband = 0, nband-1 DO BEGIN
;      free_lun, luns[iband]
;  ENDFOR 
  ;; put this in
  mwrfits, gstr, ofiles[0], /create
  mwrfits, rstr, ofiles[1], /create
  mwrfits, istr, ofiles[2], /create

END 
