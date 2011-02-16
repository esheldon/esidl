PRO make_scat53, run, rerun, clr, $
                 nsig=nsig, $
                 typecut=typecut,$
                 primary=primary, $
                 outdir=outdir

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: 
;   MAKE_SCAT
;       
; PURPOSE: 
;   For a given color and run, combine all "SOURCE" galaxies into one big file
;   with a limited number of tags.  
;	
; CALLING SEQUENCE:
;      make_scat, run, rerun, clr, $
;                     nsig=nsig, $
;                     typecut=typecut,$
;                     primary=primary,$
;                     onefile=onefile,$
;                     outdir=outdir
; INPUTS:
;  run: run in integer form
;  rerun: rerun in integer form
;  clr: bandpass in integer form (g=1, r=2, i=3)
;
; OPTIONAL INPUTS:
;  nsig=nsig: the number of sigma obove stellar locus to cut for galaxies.
;        default is 4 sigma (I no longer do this cut)
;  /typecut: cut on photo type
;  /primary: cut on status=primary (not as good as cut by ra/dec/seeing)
;  outdir=outdir: output directory
;
; OUTPUTS:
;  One large combined file of all six columns.
;
; CALLED ROUTINES:
;  SDSSIDL_SETUP
;  FETCH_DIR
;  FETCH_FILE_LIST
;  PSFNAME
;  SXPAR
;  SXADDPAR
;  FXHCLEAN
;  EQ2SURVEY
;  (MAKE_STATUS_STRUCT)
;  (STATUS_SELECT)
;  MWRFITS
;
;
; PROCEDURE: 
;	
;	
;
; REVISION HISTORY:
;   1-June-2000  Erin Scott Sheldon UofMichigan
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() LT 3 THEN BEGIN 
      print,'-Syntax: make_survey_scat, run, rerun, clr, '
      print,'          nsig=nsig, '
      print,'         typecut=typecut,'
      print,'         primary=primary,'
      print,'         outdir=outdir'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time = systime(1)

  IF n_elements(nsig) EQ 0 THEN nsig = 4.
  IF NOT keyword_set(typecut) THEN typecut=0
  IF NOT keyword_set(primary) THEN primary=0
  
  ;; set up status flag select stuff 
  make_status_struct,statusstr
  statusstr.primary='Y'

  newfront = 'adatc'
  selectclr = 2.                ; select based on r-band magnitude

  edef = -9999.                 ; default value for e1 and e2 if not meas.
  rotdef = 9999.                ; default rotation.
 
  rcut = 0.8                    ; smear polarizeability cut.

  galtype = 3

  rr = ntostr(run)
  rrr = ntostr(rerun)
                       
  colmin = 1                    ; What columns to use
  colmax = 6
                                
  colors = ['u','g','r','i','z']

  maxmag = 22.
  minmag = 16.

  ;; Arrays to hold first/last lambda's
  flam_arr = dblarr(6)
  llam_arr = dblarr(6)
  
  crflag = 0

  typ = 'blah1'+ntostr(long(systime(1)))

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Catalog of corrected shapes and positions.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  val = -9999d
  struct = create_struct(name=typ,$
                         'run', 0L, $
                         'rerun', 0, $
                         'camcol', 0b, $
                         'field', 0, $
                         'id', 0L, $
                         'e1',0.0, $
                         'e2',0.0, $
                         'e1_recorr',0.0,$
                         'e2_recorr',0.0,$
                         'e1e1err', 0.0, $
                         'e1e2err', 0.0, $
                         'e2e2err', 0.0, $
                         'ra', val, $
                         'dec', val, $
                         'r',0.0, $
                         'rpetro',0.0,$
                         'ipetro',0.0, $
                         'grmodel',0.0, $
                         'rimodel',0.0, $
                         'seeing', 0.0,$
                         'rotation', 0.0, $
                         'photoz_z', 0.0, $
                         'photoz_zerr', 0.0, $
                         'photoz_quality', 0) 

  sdssidl_setup, /silent
  setup_mystuff
  
  IF n_elements(outdir) EQ 0 THEN BEGIN 
      outdir=sdssidl_config('shapecorr_dir')+$
        'corr'+rr+'/'+ntostr(rerun)+'/combined/'
  ENDIF 
  outname = $
    outdir + 'run'+rr+'_'+rrr+'_srcgal_'+colors[clr]+'.fit'

  print
  print,'Output file = ',outname

  FOR camcol = colmin, colmax DO BEGIN

      ci = camcol-1
      cstr = ntostr(camcol)
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Setup the column
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      tmp = 0
      fetch_dir, run, camcol, rerun, dir, atldir, $
        corrdir=corrdir, corratldir=fitdir

      fetch_file_list,corrdir, files, fnums,front=newfront, $
        fieldmin=fieldmin, fieldmax=fieldmax
      nfields = n_elements(fnums)

      IF camcol EQ colmin THEN BEGIN
          ;; allow for some fields missing in
          ;; some column.  Screwed if first/last
          ;; missing from column 1
          tnfields = fieldmax-fieldmin+1
          numlist = lonarr(6, tnfields)
          ptrlist = ptrarr(6, tnfields)
          ntotal = lonarr(6)
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; correct for slope in egal vs epsf?
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      corshape_names, run, rerun, camcol, clr,$
                      psfile, corshapefile
      print
      print,'Reading corshape file for slope correction: ',corshapefile

      slopestruct=mrdfits(corshapefile,1)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; read in corrected fields
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      print
      print,'reading from ',corrdir
      print,'outputting to ',outdir

      first=1
      FOR fi=0, nfields-1 DO BEGIN 

          field = fnums[fi]
          fstr = ntostr(field)
          ind = fi*4L

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; get corrected file
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          read_tsobj, corrdir, tmp, start=field, tsobjstr=tsobjstr, $
                      front=newfront, verbose=0, status=status

          ;; is file non-empty?
          IF status EQ 0 THEN BEGIN 

              ;; get info from header
              IF first THEN BEGIN 
                  
                  ;; get header info
                  print
                  print,'Getting HDR info'
                  outhdr = headfits(files[fi])
                  IF (sxpar(outhdr, 'RUN') EQ 0) OR $
                    (sxpar(outhdr,'STRIPE') EQ 0 ) THEN BEGIN 
                      print
                      print,'No stripe info found in hdr'
                  ENDIF  
                  
                  first=0
              ENDIF 

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; make cut on status (primary)
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              IF keyword_set(primary) THEN BEGIN 
                  status_select, tmp, statusstr, gind
                  IF gind[0] EQ -1 THEN ng=0 ELSE ngal=n_elements(gind)
              ENDIF ELSE BEGIN
                  ngal = n_elements(tmp)
                  gind = lindgen(ngal)
              ENDELSE 
              
              IF ngal NE 0 THEN BEGIN 
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; Select objects by r-band magnitude.  This is probably not 
                  ;; a good idea for clusters because you will get a lot of 
                  ;; contamination.
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  ;; Reddening correct
                  tmp.petrocounts = tmp.petrocounts-tmp.reddening
                  tmp.counts_model = tmp.counts_model-tmp.reddening

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; Choose based on selectclr magnitude
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  w = where( tmp[gind].petrocounts[selectclr] LT maxmag AND $
                             tmp[gind].petrocounts[selectclr] GE minmag, ngal)
                  IF ngal NE 0 THEN BEGIN 

                      gind = gind[w]

                      ;; cut galaxies based on size relative to PSF
                      ;; the r > 0 is good enough to ensure e1/e2
                      ;; have been measured
                      ;; Also, rotation must have been measured for
                      ;; this field

                      w = where( (tmp[gind].m_r[clr] LT rcut) AND $
                                 (tmp[gind].m_r[clr] GT   0.) AND $
                                 (tmp[gind].rotation[clr] LT rotdef), ngal )
                      
                      IF ngal NE 0 THEN BEGIN 
                          gind=gind[w]
                          
                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                          ;; make cut on PHOTO type if requested
                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                          IF typecut THEN BEGIN 
                              w=where( (tmp[gind].type[1] EQ galtype AND $
                                        tmp[gind].type[2] EQ galtype) OR $
                                       (tmp[gind].type[2] EQ galtype AND $
                                        tmp[gind].type[3] EQ galtype),ngal)
                          ENDIF ELSE w=lindgen(ngal)
                              
                          IF ngal NE 0 THEN BEGIN 

                              gind=gind[w]
                              
                              
                              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                              ;; seeing cut if requested
                              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                              meanseeing=mean_check(tmp[gind].seeing[clr])
                              IF n_elements(maxseeing) NE 0 THEN BEGIN 
                                  w=where(tmp[gind].seeing[clr] LT $
                                          maxseeing, ngal)
                              ENDIF ELSE w=lindgen(ngal)
                              IF ngal NE 0 THEN BEGIN 

                                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                  ;; Final cuts have been made. Copy 
                                  ;; into structure
                                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                  gind=gind[w]
          
                                  srctmp = replicate(struct, ngal)
                                  srctmp.seeing = tmp[gind].seeing[clr]

                                  ;; dilution correction
                                  corr = 1./(1.- tmp[gind].m_r[clr])
                                  srctmp.e1e1err = $
                                    tmp[gind].m_e1e1err[clr]*corr
                                  srctmp.e1e2err = $
                                    tmp[gind].m_e1e2err[clr]*corr
                                  srctmp.e2e2err = $
                                    tmp[gind].m_e2e2err[clr]*corr

                                  e1 = tmp[gind].m_e1_corr[clr]
                                  e2 = tmp[gind].m_e2_corr[clr]

                                  ;; Correct for slope
                                  ;; Done before 
                                  ;; blurring correction below

                                  psfe1 = tmp[gind].m_e1_psf[clr]
                                  psfe2 = tmp[gind].m_e2_psf[clr]
                                  
                                  correct_eslope, slopestruct, $
                                                  e1, $
                                                  e2, $
                                                  psfe1, psfe2,$
                                                  e1_recorr, e2_recorr

                                  ;; correct for dilution 
                                  e1 = e1*corr
                                  e2 = e2*corr
                                  e1_recorr = e1_recorr*corr
                                  e2_recorr = e2_recorr*corr

                                  srctmp.e1 = e1
                                  srctmp.e2 = e2
                                  srctmp.e1_recorr = e1_recorr
                                  srctmp.e2_recorr = e2_recorr
                                  srctmp.rotation = tmp[gind].rotation[clr]

                                  ;; CHANGED: just convert ra/dec to
                                  ;;          lam/eta later if needed

                                  srctmp.ra = tmp[gind].ra
                                  srctmp.dec = tmp[gind].dec
                                  srctmp.r = tmp[gind].m_r[clr]
                                  srctmp.rpetro = tmp[gind].petrocounts[2]
                                  srctmp.ipetro = tmp[gind].petrocounts[3]
                                  srctmp.grmodel = $
                                    tmp[gind].counts_model[1]-$
                                    tmp[gind].counts_model[2]
                                  srctmp.rimodel = $
                                    tmp[gind].counts_model[2]-$
                                    tmp[gind].counts_model[3]

                                  ;; photoz info
                                  srctmp.photoz_z = tmp[gind].photoz_z
                                  srctmp.photoz_zerr = tmp[gind].photoz_zerr
                                  srctmp.photoz_quality = $
                                    tmp[gind].photoz_quality
                                  
                                  srctmp.run = tmp[gind].run
                                  srctmp.rerun = tmp[gind].rerun
                                  srctmp.camcol = tmp[gind].camcol
                                  srctmp.field = tmp[gind].field
                                  srctmp.id = tmp[gind].id

                                  ;; copy into pointers
                                  numlist[ci, fi] = ngal
                                  ntotal[ci] = ntotal[ci] + ngal
                                  ptrlist[ci, fi] = ptr_new(srctmp)
                                  
                                  IF fi MOD 20 EQ 0 THEN BEGIN 
                                      print,'Run: ',rr,' Camcol: ',cstr,$
                                            ' Field: ',fstr,'  ',$
                                            colors[clr]+$
                                            '-band Src Galaxies: '+$
                                            ntostr(ngal)
                                      print,'Seeing['+colors[clr]+'] = ',$
                                            ntostr(meanseeing)+'"',$
                                            ' Mean e1: ',ntostr(mean(e1))
                                  ENDIF 

                                  ;; free some memory
                                  tmp=0 & srctmp=0 & corr=0
                                  e1rot=0 & e2rot=0
                              ENDIF ELSE BEGIN
                                  message,'No galaxies passed seeing cut of '+$
                                          ntostr(maxseeing),/inf
                              ENDELSE 
                          ENDIF ELSE BEGIN
                              message,'No galaxies passed typecut '+$
                                      'in field '+fstr,/inf
                          ENDELSE 
                      ENDIF ELSE BEGIN
                          message,'No galaxies passed size '+$
                                  'and rotation cuts in field '+fstr,/inf
                      ENDELSE 
                  ENDIF ELSE BEGIN
                      message,'No galaxies in magnitude range '+$
                              colors[clr]+' = ['+ntostr(minmag)+', '+$
                              ntostr(maxmag)+']'+' in field '+fstr,/inf
                  ENDELSE 
              ENDIF ELSE BEGIN
                  message,'No objects passed status cut in field '+$
                          fstr,/inf
              ENDELSE 
          ENDIF ELSE BEGIN
              message,'File '+files[fi]+' is empty',/inf
          ENDELSE 
      ENDFOR ;; loop over fields

      rotstruct=0
      print,'Run: ',rr,' Camcol: ',cstr,' total ',colors[clr],$
            '-band galaxies: ',$
            ntostr(ntotal[ci])
      print

  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Combine all columns together
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Combining columns into larger structures'
  print

  ntotal = long( total(ntotal) )

  sources = replicate(struct, ntotal )
  beg=0L
  FOR camcol = colmin, colmax DO BEGIN 
      ci = camcol-1
      FOR fi=0, nfields-1 DO BEGIN 

          ;; Don't use fields flagged as bad
          IF (numlist[ci, fi] NE 0) THEN BEGIN 
              sources(beg:beg+numlist[ci, fi]-1)=*ptrlist[ci, fi]
          ENDIF 
          ptr_free,ptrlist[ci, fi]   ;kill the heap variables associated 
                                ;with the pointer
          beg=beg+numlist[ci, fi]

      ENDFOR 
  ENDFOR 

  ;; get first clambda, last clambda
  eq2csurvey, sources.ra, sources.dec, clam, ceta
  firstlam = min(clam,max=lastlam)

  sxaddpar, outhdr, 'RCUT', rcut, ' Maximum Rsmear cut for star-galaxy separation'
  sxaddpar, outhdr, 'TYPECUT', typecut, ' Was a typecut made?'
  sxaddpar, outhdr, 'PRIMARY', primary, ' Was a primary cut made?'

  sxaddpar, outhdr, 'FIRSTLAM', firstlam
  sxaddpar, outhdr, 'LASTLAM', lastlam

  print,'FIRSTLAM: ',firstlam
  print,'LASTLAM: ',lastlam
  print

  print,'Outfile: ',outname
  mwrfits2, temporary(sources), outname, /destroy, /create, hdr0=outhdr

  print
  print,'Total ',colors[clr],'-band galaxies: ',ntostr(ntotal)

  ptime,systime(1)-time

  return
END 
