PRO makecor, run, rerun, camcol, start=start, nframes=nframes, stmaxmag=stmaxmag,$
             status=status

  ;; status is 1 unless we reach end
  status=1

  IF n_params() LT 3 then begin
      print,'-syntax makecor, run, rerun, camcol, start=start, nframes=nframes, stmaxmag=stmaxmag, status=status'
      return
  ENDIF

  time=systime(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  newfront = 'adat'             ;The front of files with adaptive moments
  outfront = 'adatc'

  colors=['u','g','r','i','z']
  bands = [1,2,3]               ; we will process g,r,i
  nband = n_elements(bands)
  
  nmag=n_elements(stmaxmag)
  IF nmag NE nband THEN BEGIN 
      IF nmag EQ 0  THEN BEGIN
          stmaxmag = [21.,20.,20.]
      ENDIF ELSE BEGIN
          print,'stmaxmag must be size ',nband
      ENDELSE 
  ENDIF 

  rstr = ntostr(run)
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

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up new filenames. 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  rename_tsobj, files, corrdir, newfront, adatfiles, nchar
  rename_tsobj, files, corrdir, outfront, outfiles, outnchar

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read in file containing the shape measurements
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  admomatlas_infile, run, rerun, camcol, fname_in, fname_out
  fname_out = corrdir + repstr(fname_out, '.dat', '.fit')
  print,'Reading admom file: ',fname_out
  admomstruct = mrdfits(fname_out, 1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over fields and find corrected e1, e2 and R
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'-----------------------------------------------------'
  print,' Run: ',rstr,' Rerun: ',rrstr,' Camcol: ',cstr
  print,'-----------------------------------------------------'

  print,'Outputs files are '+outfront+'-*'
  print

  FOR ic = 0L, nfields-1 DO BEGIN 
      infile = adatfiles[ic]
      outf   = outfiles[ic]

      IF outf EQ infile THEN BEGIN
          print,'Error: infile = outfile'
          print,infile,outf
          return
      ENDIF 

      field = fnums[ic]
      fstr = ntostr(field)

      ;; psf fit name
      psfname2, run, rerun, camcol, field, psffitname
      psffitname = fitdir + psffitname

      ;; Fields have already been trimmed of overlap.  Don't need to
      ;; use read_tsobjmin
      openr, lun, infile, /get_lun, ERROR=error1
      IF error1 NE 0 THEN BEGIN
          print,!ERR_STRING
          free_lun, lun
          GOTO,jump2
      ENDIF 
      openr, lun2, psffitname, /get_lun, ERROR=error2
      IF ic EQ 0 THEN BEGIN 
          pstruct = mrdfits3(lun, 1, 0, hdr1, /silent)
          IF error2 EQ 0 THEN psfstr = mrdfits4(lun2, 1, 0, /silent)
          IF tag_exist(pstruct,'seeing') THEN doseeing = 1 ELSE doseeing=0
      ENDIF ELSE BEGIN
          pstruct = mrdfits3(lun, 1, 0, hdr1, /silent, /deja_vu)
          IF error2 EQ 0 THEN psfstr = mrdfits4(lun2, 1, 0, /silent, /deja_vu )
      ENDELSE 
      free_lun, lun
      free_lun, lun2

      fxhclean, hdr1
      hdr0 = headfits(infile, exten=0)

      IF error2 NE 0 THEN BEGIN 
          print,!ERR_STRING
          GOTO,jump
      ENDIF 

      IF ( size(pstruct) )[0] EQ 0 THEN GOTO,jump
      
      ;; match up with admomstruct
      aa = where(admomstruct.field EQ field,naa)
      IF naa EQ 0 THEN BEGIN 
          print,'Field '+fstr+' not in admom struct!'
          GOTO,jump
      ENDIF 

      photo_match, pstruct.run, pstruct.rerun, pstruct.camcol, $
        pstruct.field, pstruct.id,$
        admomstruct[aa].run,  admomstruct[aa].rerun,  admomstruct[aa].camcol,  $
        admomstruct[aa].field,  admomstruct[aa].id, $
        tmobj, tmpadmom

      ;; copy in the shape measurements
      pstruct[tmobj].ixx = admomstruct[aa[tmpadmom]].ixx
      pstruct[tmobj].ixy = admomstruct[aa[tmpadmom]].ixy
      pstruct[tmobj].iyy = admomstruct[aa[tmpadmom]].iyy
      pstruct[tmobj].momerr = admomstruct[aa[tmpadmom]].momerr
      pstruct[tmobj].rho4 = admomstruct[aa[tmpadmom]].rho4
      pstruct[tmobj].whyflag = admomstruct[aa[tmpadmom]].whyflag
      
      ;; we will use our own seeing measurements
      IF doseeing THEN pstruct.seeing = 0.

      ;; match up with psf struct
      photo_match, pstruct.run, pstruct.rerun, pstruct.camcol, $
        pstruct.field, pstruct.id,$
        psfstr.run,  psfstr.rerun,  psfstr.camcol,  $
        psfstr.field,  psfstr.id, $
        tmobj, tmpsf

      nmatch = n_elements(tmobj)

      print,'Run: ',rstr,' Camcol: ',cstr,' Field: ',fstr,' Use: ',ntostr(nmatch)

      FOR iband = 0, nband-1 DO BEGIN 
         
          mpsf=tmpsf
          mobj=tmobj
          cindex = bands[iband]

          psfsize = psfstr[mpsf].psfixx[cindex] + psfstr[mpsf].psfiyy[cindex]
          objsize = pstruct[mobj].ixx[cindex] + pstruct[mobj].iyy[cindex]

          IF doseeing THEN BEGIN 
              ;; seeing in arcseconds
              seeing = sqrt(psfsize/2.)*2.35*0.4
              pstruct[mobj].seeing[cindex] = seeing
              seeing=0
          ENDIF 

          ;; Bad measurements
          ww = where(objsize GT 0. AND psfsize GT 0.,nww)
          IF nww EQ 0 THEN BEGIN
              print
              print,'No objects measured!'
              print
              GOTO,jump
          ENDIF 
          mpsf = mpsf[ww]
          mobj = mobj[ww]
          psfsize = psfsize[ww]
          objsize = objsize[ww]

          psfixx = psfstr[mpsf].psfixx[cindex]
          psfiyy = psfstr[mpsf].psfiyy[cindex]
          psfixy = psfstr[mpsf].psfixy[cindex]
          psfrho4 = psfstr[mpsf].psfrho4[cindex]

          cor1 = psfsize/objsize*(4./psfrho4-1.)/(4./pstruct[mobj].rho4[cindex]-1.)

          psfe1 = (psfixx - psfiyy)/psfsize

          psfe2 = 2.*psfixy/psfsize
      
          obje1 = ( pstruct[mobj].ixx[cindex] - pstruct[mobj].iyy[cindex] )/objsize
          obje2 = 2.*pstruct[mobj].ixy[cindex]/objsize
      
          ;; can remove when not using previously made adatc files!
;          pstruct.e1[cindex] = 1.e10
;          pstruct.e2[cindex] = 1.e10
;          pstruct.R[cindex] = -1.

          pstruct[mobj].e1[cindex] = obje1 - psfe1*cor1 
          pstruct[mobj].e2[cindex] = obje2 - psfe2*cor1
          pstruct[mobj].R[cindex]  = cor1
          
          print,'Band: '+colors[cindex]+' Ngood: '+ntostr(nww)+'  ',format='($,a)'

          ;; extract stars
          phextract_stars, pstruct, cindex, stars, max_mag=stmaxmag[iband], $
                           /silent
          ;; 2^0 u-gand, 2^1 g-band, 2^2 r-band, etc.
          IF stars[0] NE -1 THEN BEGIN
              pstruct[stars].starflag = pstruct[stars].starflag + 2b^cindex
              print,'Nstars:  ',ntostr(n_elements(stars))
          ENDIF ELSE print
          
      ENDFOR                    ;End loop over bandpasses

jump:
      mwrfits2, pstruct, outf, hdr1, /create, hdr0=hdr0
jump2:
                                ;Free up the memory
      setzero, pstruct, psfstr, cor1, psfsize, psfixx, psfiyy, psfixy, $
        psfrho4, psfe1, psfe2, obje1, obje2

      
  ENDFOR                        ;End loop over fields
   
  ptime,systime(1)-time

  status=0

  return
END 

