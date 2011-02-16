PRO calc_corrshapes, pstruct, status=status

  ;; status is 1 unless we reach end
  status=1

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: calc_corrshapes, pstruct, status=status'
      return
  ENDIF 

  ;; bandpasses g,r,i
  bands = [1,2,3]
  stmaxmag = [21.,20.,20.]
  nband = n_elements(bands)

  ;; loop over bandpass and do corrected shape (not including dilution)

  FOR iband = 0L, nband-1 DO BEGIN 
      
      cindex = bands[iband]

      psfsize = pstruct.psfixx[cindex] + pstruct.psfiyy[cindex]
      objsize = pstruct.ixx[cindex] + pstruct.iyy[cindex]

      ;; seeing in arcseconds
      pstruct.seeing[cindex] = sqrt(psfsize/2.)*2.35*0.4

      ;; now correct shapes for well-measured objects
      ww=where(objsize GT 0.0 AND psfsize GT 0., nww)
      IF nww NE 0 THEN BEGIN 

          psfsize = psfsize[ww]
          objsize = objsize[ww]

          psfixx  = pstruct[ww].psfixx[cindex]
          psfiyy  = pstruct[ww].psfiyy[cindex]
          psfixy  = pstruct[ww].psfixy[cindex]
          psfrho4 = pstruct[ww].psfrho4[cindex]

          cor1 = psfsize/objsize*(4./psfrho4-1.)/(4./pstruct[ww].rho4[cindex]-1.)

          psfe1 = (psfixx - psfiyy)/psfsize
          psfe2 = 2.*psfixy/psfsize

          obje1 = ( pstruct[ww].ixx[cindex] - pstruct[ww].iyy[cindex] )/objsize
          obje2 = 2.*pstruct[ww].ixy[cindex]/objsize

          pstruct[ww].e1[cindex] = obje1 - psfe1*cor1 
          pstruct[ww].e2[cindex] = obje2 - psfe2*cor1
          pstruct[ww].R[cindex]  = cor1

          print,'Band: '+!colors[cindex]+' Ngood: '+ntostr(nww)+'  ',format='($,a)'

          ;; extract stars
          phextract_stars, pstruct, cindex, stars, max_mag=stmaxmag[iband], $
                           /silent
          ;; 2^0 u-gand, 2^1 g-band, 2^2 r-band, etc.
          IF stars[0] NE -1 THEN BEGIN
              pstruct[stars].starflag = pstruct[stars].starflag + 2b^cindex
              print,'Nstars:  ',ntostr(n_elements(stars))
          ENDIF ELSE print

      ENDIF ELSE print,'No objects measured!'

      setzero, psfsize, objsize, psfixx, psfiyy, psfixy, psfrho4, cor1,$
               psfe1, psfe2, obje1, obje2
  ENDFOR 

  status=0
END 

PRO admom_psfield, psfname, pstruct, status=status

  ;; status is 1 unless we reach end
  status=1
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: admom_psfield, psfname, pstruct'
      return
  ENDIF 

  bands = [0,1,2,3,4]               ; we will process u,g,r,i,z
  nband = n_elements(bands)

  bindir = !mybindir
  ;; psf info for each bandpass
  tmp0 = mrdfits(psfname, 1, /silent)
  tmp1 = mrdfits(psfname, 2, /silent)
  tmp2 = mrdfits(psfname, 3, /silent)
  tmp3 = mrdfits(psfname, 4, /silent)
  tmp4 = mrdfits(psfname, 5, /silent)

  IF (  ( (size(tmp0))[0] EQ 0 ) OR $
        ( (size(tmp1))[0] EQ 0 ) OR $
        ( (size(tmp2))[0] EQ 0 ) OR $
        ( (size(tmp3))[0] EQ 0 ) OR $
        ( (size(tmp4))[0] EQ 0 ) )     THEN BEGIN
      print,'psfield file '+psfield+' missing or corrupted'
      tmp0=0 & tmp1=0 & tmp2=0 & tmp3=0 & tmp4=0
      return
  ENDIF 

  ;; make pointers
  psp=ptrarr(5)

  psp[0] = ptr_new( tmp0, /no_copy )
  psp[1] = ptr_new( tmp1, /no_copy ) ;deletes tmp1
  psp[2] = ptr_new( tmp2, /no_copy )
  psp[3] = ptr_new( tmp3, /no_copy )
  psp[4] = ptr_new( tmp4, /no_copy )

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; loop over the objects
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nobj = n_elements(pstruct)
  FOR ind=0L, nobj-1 DO BEGIN 

      psf_rec, psp, pstruct[ind].rowc, pstruct[ind].colc, psf, bands
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Loop over bandpasses
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      FOR iband = 0, nband-1 DO BEGIN 

          cindex = bands[iband]
          
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
          pstruct[ind].psfixx[cindex] = ixx
          pstruct[ind].psfiyy[cindex] = iyy
          pstruct[ind].psfixy[cindex] = ixy
          pstruct[ind].psfrho4[cindex] = rho4

      ENDFOR                    ;Loop over bandpasses     
      ptrarr_free, psf
      IF (ind MOD 20) EQ 0 THEN print,'.',format='($,a)'

  ENDFOR ;; Loop over objects

  print
  ptrarr_free, psp

  status=0
END 

PRO make_corrected_files, run, rerun, camcol, $
                          psFieldrerun=psFieldrerun, start=start, nframes=nframes,$
                          status=status, outdir=outdir, addphotomom=addphotomom

  ;; status is 1 unless we reach end
  status =1

  IF n_params() LT 3 THEN BEGIN ;Help message
      print,'-syntax make_corrected_files, run, rerun, camcol, psFieldrerun=psFieldrerun,start=start, nframes=nframes, outdir=outdir, status=status, addphotomom=addphotomom'
      return 
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  setup_mystuff

  newfront = 'adatc'             ;The front of files with adaptive moments

  colors=['u','g','r','i','z']
  runstr = ntostr(run)
  rerunstr=ntostr(rerun)
  cstr = ntostr(camcol)

  make_admomtags, taglist, /default, addphotomom=addphotomom ;don't include extra stuff here, add later
  fetch_dir, run, camcol, rerun, dir, atldir, corrdir=corrdir, $
    corratldir=skydir           ;We keep sky files in corratldir
  fetch_file_list, dir, files, fnums, start=start, nframes=nframes, $
    fieldmin=fieldmin, fieldmax=fieldmax
  make_lensstruct, lensstr

  nfields = n_elements(files)

  IF n_elements(psFieldrerun) EQ 0 THEN psFieldrerun=rerun
  fetch_dir, run, camcol, psFieldrerun, ddd, psField_dir

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; read the admomatlas.c output file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  admomatlas_infile, run, rerun, camcol, fname_in, fname_out
  fname_out = repstr(fname_out, '.dat', '.fit')
  admomstruct = mrdfits(corrdir + fname_out, 1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; allow user to put output files anywhere 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(outdir) NE 0 THEN corrdir=outdir

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  rename_tsobj, files, corrdir, newfront, outfiles, outnchar

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Get the files that contain the rotation angle for each field
  ;; rotation between (row,col) and (lambda,eta)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  rotname = skydir + 'surveyrot_'+run2string(run)+'_'+ntostr(camcol)+'_'+$
    ntostr(rerun)+'.fit'
  rotstruct = mrdfits(rotname,1,/silent)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over fields and find psf moments and correct shapes
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Outputs files are '+newfront+'-*'
  print

  FOR ic = 0L, nfields-1 DO BEGIN 

      infile = files[ic]
      outfile = outfiles[ic]
      field = fnums[ic]
      fstr=ntostr(field)

      ;; read headers
      hdr0=headfits(infile, exten=0)
      hdr1=headfits(infile, exten=1)
      fxhclean, hdr1

      ;; Read tsObj file
      read_tsobj, dir, pstruct, start=field, nframes=1, $
                  taglist=taglist, tsobjstr=tsobjstr, ex_struct=lensstr, $
                  /noadd,verbose=0

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Rotation (one per field for now)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      nps = n_elements(pstruct)
      print
      print,'Run: ',runstr,' Camcol: ',cstr,' Field: ',fstr,' Nobj: ',ntostr(nps)
      IF nps NE 0 THEN BEGIN 

          wrot = where(rotstruct.field EQ field,nrot)
          IF nrot NE 0 THEN BEGIN 
              pstruct.rotation = rotstruct[wrot].angle
              
              ;; Match to admomstruct
              aa = where(admomstruct.field EQ field,naa)
              IF naa NE 0 THEN BEGIN 
                  
                  photo_match, pstruct.run, pstruct.rerun, pstruct.camcol, $
                               pstruct.field, pstruct.id,$
                               admomstruct[aa].run, admomstruct[aa].rerun,  $
                               admomstruct[aa].camcol,  $
                               admomstruct[aa].field,  admomstruct[aa].id, $
                               mobj, tmpadmom
                  
                  IF mobj[0] NE -1 THEN BEGIN 
                      madmom = aa[tmpadmom]
                      
                      ;; only keep matches
                      pstruct = temporary(pstruct[mobj])
                      
                      ;; get psf info
                      psfield_name, run, camcol, field, psfname
                      psfname = psField_dir + psfname
                      print,'Calculating PSF moments'
                      admom_psfield, psfname, pstruct
                      
                      ;; copy in moments
                      pstruct.ixx = admomstruct[madmom].ixx
                      pstruct.ixy = admomstruct[madmom].ixy
                      pstruct.iyy = admomstruct[madmom].iyy
                      pstruct.momerr = admomstruct[madmom].momerr
                      pstruct.rho4 = admomstruct[madmom].rho4
                      pstruct.whyflag = admomstruct[madmom].whyflag
                      
                      ;; make corrected shapes
                      print,'Correcting shapes'
                      calc_corrshapes, pstruct
                      
                  ENDIF ELSE BEGIN 
                      print,'No matches'
                      delvarx, pstruct
                  ENDELSE ;; no matches
              ENDIF ;; field not in admomastruct 
          ENDIF ;; no rotation info       
          
      ENDIF ;; empty pstruct

      ;; write file (may be empty)
      mwrfits2, pstruct, outfile, hdr1, /create, /destroy, hdr0=hdr0

  ENDFOR 

  status=0

END 
