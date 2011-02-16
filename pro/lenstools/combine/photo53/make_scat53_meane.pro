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
;    make_scat53_meane, run, rerun, $
;                       hardrcut=hardrcut, $
;                       err_cut=err_cut,minerr=miner, $
;                       besterror=besterror,$
;                       wmomerror=wmomerror,$
;                       outdir=outdir
;
; INPUTS:
;  run: run in integer form
;  rerun: rerun in integer form
;  clr: bandpass in integer form (g=1, r=2, i=3)
;
; OPTIONAL INPUTS:
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

PRO make_scat53_meane_applyrot, rotation, nclrs, e1, e2, e1_recorr, e2_recorr, $
                                e1e1err, e1e2err, e2e2err


  FOR clr=0,nclrs-1 DO BEGIN 

      w=where(e1e1err[clr,*] GT 0.0,nw)

      IF nw NE 0 THEN BEGIN 
          rotate_e1e2error, rotation[clr,w], e1[clr,w], e2[clr,w], $
                            e1e1err[clr,w],e1e2err[clr,w],e2e2err[clr,w],$
                            e1out,e2out,e1e1errout,e1e2errout,e2e2errout

          e1[clr,w] = e1out
          e2[clr,w] = e2out
          e1e1err[clr,w] = e1e1errout
          e1e2err[clr,w] = e1e2errout
          e2e2err[clr,w] = e2e2errout
          
          rotate_e1e2, rotation[clr,w], e1_recorr[clr,w], e2_recorr[clr,w],$
                       e1out, e2out
          e1_recorr[clr,w] = e1out
          e2_recorr[clr,w] = e2out
      ENDIF 
  ENDFOR 

END 

PRO make_scat53_meane_flagselect, pstruct, flag_struct, statusstr, clrs, keep, nkeep

  nclr = n_elements(clrs)

  status_select, pstruct, statusstr, wtmp
  IF wtmp[0] EQ -1 THEN BEGIN 
      nkeep = 0
      keep = -1
      return
  ENDIF 

  ;; flag select on each bandpass
  FOR i=0L, nclr-1 DO BEGIN 
      flag_select, pstruct, flag_struct, clrs[i], keep, input_index=wtmp
      IF keep[0] EQ -1 THEN BEGIN 
          nkeep = 0
          return
      ENDIF ELSE IF i NE nclr-1 THEN wtmp = temporary(keep)
  ENDFOR 
  nkeep = n_elements(keep)

  return
END 

PRO make_scat53_meane_rselect_varcut, tmp, gind, clrs, $
                                      hardrcut, hardprobcut, $
                                      grcutstruct,rrcutstruct,ircutstruct, $
                                      gslopestruct,rslopestruct,islopestruct,$
                                      e1, e2, e1e1err, e1e2err, e2e2err, $
                                      e1_recorr, e2_recorr, $
                                      smear, corr, $
                                      good, ngood

  ;; keep and nkeep are the union of good measurements in
  ;; all bandpasses "clrs"

  COMMON make_scat53_block, edef, e_errdef, smdef, rotdef

  ntotal = n_elements(tmp)

  nclrs = n_elements(clrs)

  e1 = reform( replicate(edef, nclrs, ntotal) )
  e2 = e1
  e1_recorr = e1
  e2_recorr = e2

  psfe1 = tmp.m_e1_psf[clrs]
  psfe2 = tmp.m_e2_psf[clrs]

  e1e1err = reform( replicate(e_errdef, nclrs,  ntotal) )
  e1e2err = e1e1err
  e2e2err = e1e1err

  smear = reform( replicate(smdef, nclrs, ntotal) )
  corr = smear

  rmag = tmp.petrocounts[2]

  ;; rotation?
  w = where( tmp[gind].rotation[2] LT rotdef, tngal )
  IF tngal NE n_elements(gind) THEN BEGIN 
      good = -1
      ngood = 0
      print,'No rotation found for this field'
      return
  ENDIF

  ngal = 0L
  FOR ii=0L, nclrs-1 DO BEGIN 
      tclr = clrs[ii]
      CASE tclr OF
          1: trcutstruct = grcutstruct
          2: trcutstruct = rrcutstruct
          3: trcutstruct = ircutstruct
          ELSE: message,'Bad clr' 
      ENDCASE 
      tsmear = tmp[gind].m_r[tclr]

      cut = interpol(trcutstruct.rsmear_cuts, $
                     trcutstruct.meanmag, rmag[gind]) < hardrcut

      tgind = where(tsmear LT cut AND $
                    tsmear GT 0.0 AND $
                    tmp[gind].objc_prob_psf GE 0.0 AND $
                    tmp[gind].objc_prob_psf LT (1.-hardprobcut), tngal)

      IF tngal NE 0 THEN BEGIN 
          
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; keep track of how many found in any bandpass
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          ind = gind[tgind]

          IF good[0] EQ -1 THEN BEGIN
              good = ind
              ngood = tngal
          ENDIF ELSE BEGIN
              good = [good, ind]
              good = good[uniq(good, sort(good))]
              ngood = n_elements(good)
          ENDELSE 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; copy in from the structure: those that failed
          ;; will get defval and be ignored later
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          
          e1[ii, ind] = tmp[ind].m_e1_corr[tclr]
          e2[ii, ind] = tmp[ind].m_e2_corr[tclr]
          e1e1err[ii, ind] = tmp[ind].m_e1e1err[tclr]
          e1e2err[ii, ind] = tmp[ind].m_e1e2err[tclr]
          e2e2err[ii, ind] = tmp[ind].m_e2e2err[tclr]
          smear[ii, ind] = tmp[ind].m_r[tclr]
          corr[ii, ind] = 1./(1.-smear[ii, ind])

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; now do slope correction
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          CASE clrs[ii] OF
              1: tslopestr = gslopestruct
              2: tslopestr = rslopestruct
              3: tslopestr = islopestruct
              ELSE: message,'Bad clr' 
          ENDCASE 

          correct_eslope, tslopestr, $
                          e1[ii,ind], $
                          e2[ii,ind], $
                          psfe1[ii,ind], psfe2[ii,ind],$
                          te1_recorr, te2_recorr
          e1_recorr[ii,ind] = te1_recorr
          e2_recorr[ii,ind] = te2_recorr

          e1[ii, ind] = e1[ii, ind]*corr[ii, ind]
          e2[ii, ind] = e2[ii, ind]*corr[ii, ind]
          e1_recorr[ii, ind] = e1_recorr[ii, ind]*corr[ii, ind]
          e2_recorr[ii, ind] = e2_recorr[ii, ind]*corr[ii, ind]
          e1e1err[ii, ind] = e1e1err[ii, ind]*corr[ii, ind]
          e1e2err[ii, ind] = e1e2err[ii, ind]*corr[ii, ind]
          e2e2err[ii, ind] = e2e2err[ii, ind]*corr[ii, ind]

      ENDIF 
  ENDFOR 

  IF n_elements(good) EQ 0 THEN BEGIN
      good = -1
      ngood = 0
  ENDIF ELSE BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; apply rotation (its ok to rotate on corrected shapes)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      rotation = tmp.rotation[clrs]
      make_scat53_meane_applyrot,rotation,nclrs,e1,e2, $
                                 e1_recorr,e2_recorr, $
                                 e1e1err, e1e2err, e2e2err
  ENDELSE 

END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; this one we apply a cut in probability
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO make_scat53_meane_rselect, tmp, gind, clrs, $
                               hardrcut, hardprobcut, $
                               gslopestruct,rslopestruct,islopestruct,$
                               e1, e2, e1e1err, e1e2err, e2e2err, $
                               e1_recorr, e2_recorr, $
                               smear, corr, $
                               good, ngood, hirata=hirata

  ;; keep and nkeep are the union of good measurements in
  ;; all bandpasses "clrs"

  COMMON make_scat53_block, edef, e_errdef, smdef, rotdef

  ntotal = n_elements(tmp)

  nclrs = n_elements(clrs)

  e1 = reform( replicate(edef, nclrs, ntotal) )
  e2 = e1
  e1_recorr = e1
  e2_recorr = e2

  psfe1 = tmp.m_e1_psf[clrs]
  psfe2 = tmp.m_e2_psf[clrs]

  e1e1err = reform( replicate(e_errdef, nclrs,  ntotal) )
  e1e2err = e1e1err
  e2e2err = e1e1err

  smear = reform( replicate(smdef, nclrs, ntotal) )
  corr = smear

  rmag = tmp.petrocounts[2]

  ;; rotation?
  w = where( tmp[gind].rotation[2] LT rotdef, tngal )
  IF tngal NE n_elements(gind) THEN BEGIN 
      good = -1
      ngood = 0
      print,'No rotation found for this field'
      return
  ENDIF

  ngal = 0L
  FOR iclr=0L, nclrs-1 DO BEGIN 

      IF keyword_set(hirata) THEN BEGIN 
          tsmear = tmp.m_r_h[clrs[iclr]]
      ENDIF ELSE BEGIN 
          tsmear = tmp.m_r[clrs[iclr]]
      ENDELSE 
      
      tgind = where(tsmear[gind] LT hardrcut AND $
                    tsmear[gind] GT 0.0 AND $
                    tmp[gind].objc_prob_psf GE 0.0 AND $
                    tmp[gind].objc_prob_psf LT (1.-hardprobcut), tngal)

      IF tngal NE 0 THEN BEGIN 
          
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; keep track of how many found in any bandpass
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          ind = gind[tgind]

          IF good[0] EQ -1 THEN BEGIN
              good = ind
              ngood = tngal
          ENDIF ELSE BEGIN
              good = [good, ind]
              good = good[uniq(good, sort(good))]
              ngood = n_elements(good)
          ENDELSE 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; copy in from the structure: those that failed
          ;; will get defval and be ignored later
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          
          IF keyword_set(hirata) THEN BEGIN 
              e1[iclr, ind] = tmp[ind].m_e1_corr_h[clrs[iclr]]
              e2[iclr, ind] = tmp[ind].m_e2_corr_h[clrs[iclr]]
          ENDIF ELSE BEGIN 
              e1[iclr, ind] = tmp[ind].m_e1_corr[clrs[iclr]]
              e2[iclr, ind] = tmp[ind].m_e2_corr[clrs[iclr]]
          ENDELSE 

          e1e1err[iclr, ind] = tmp[ind].m_e1e1err[clrs[iclr]]
          e1e2err[iclr, ind] = tmp[ind].m_e1e2err[clrs[iclr]]
          e2e2err[iclr, ind] = tmp[ind].m_e2e2err[clrs[iclr]]
          smear[iclr, ind] = tsmear[ind]
          corr[iclr, ind] = 1./(1.-smear[iclr, ind])

          ;; hiratas stuff already dilution corrected
          IF NOT keyword_set(hirata) THEN BEGIN 
              e1[iclr, ind] = e1[iclr, ind]*corr[iclr, ind]
              e2[iclr, ind] = e2[iclr, ind]*corr[iclr, ind]
          ENDIF 

          e1e1err[iclr, ind] = e1e1err[iclr, ind]*corr[iclr, ind]
          e1e2err[iclr, ind] = e1e2err[iclr, ind]*corr[iclr, ind]
          e2e2err[iclr, ind] = e2e2err[iclr, ind]*corr[iclr, ind]

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; now do slope correction (after dilution corr)
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          CASE clrs[iclr] OF
              1: tslopestr = gslopestruct
              2: tslopestr = rslopestruct
              3: tslopestr = islopestruct
              ELSE: message,'Bad clr' 
          ENDCASE 

          correct_eslope, tslopestr, $
                          e1[iclr,ind], $
                          e2[iclr,ind], $
                          psfe1[iclr,ind], psfe2[iclr,ind],$
                          te1_recorr, te2_recorr

          e1_recorr[iclr,ind] = te1_recorr
          e2_recorr[iclr,ind] = te2_recorr

      ENDIF 
  ENDFOR 

  IF n_elements(good) EQ 0 THEN BEGIN
      good = -1
      ngood = 0
  ENDIF ELSE BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; apply rotation (its ok to rotate on corrected shapes)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      rotation = tmp.rotation[clrs]
      make_scat53_meane_applyrot,rotation,nclrs,e1,e2, $
                                 e1_recorr,e2_recorr, $
                                 e1e1err, e1e2err, e2e2err
  ENDELSE 

END 

PRO make_scat53_meane, run, rerun, clrs, purity, $ ; purity not being used
                       hardrcut=hardrcut, $
                       hardprobcut=hardprobcut, $
                       err_cut=err_cut,minerr=minerr, $
                       besterror=besterror,$
                       wmomerror=wmomerror,$
                       outdir=outdir, hirata=hirata



  IF N_params() LT 3 THEN BEGIN 
      print,'-Syntax: make_scat53_meane, run, rerun, avgclrs, purity, '
      print,'          hardrcut=hardrcut, hardprobcut=hardprobcut, '
      print,'          err_cut=err_cut,'
      print,'          besterror=besterror, '
      print,'          wmomerror=wmomerror, '
      print,'          outdir=outdir, hirata=hirata'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time = systime(1)
  
  make_status_struct,statusstr
  statusstr.secondary='Y'

  make_flag_struct, flag_struct
  flag_struct.DEBLENDED_AS_PSF = 'N'

  newfront = 'adatc'
  rind = 2                ; select based on r-band magnitude
  nclrs = n_elements(clrs)      ; bandpasses to average

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Strings for identifying which bandpasses were used
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  clrstring = ''
  fclrstring = ''
  FOR i=0L, nclrs-1 DO BEGIN
      fclrstring = fclrstring+!colors[clrs[i]]
      clrstring = clrstring+!colors[clrs[i]]+' '
  ENDFOR 
  clrstring = ntostr(clrstring)

  COMMON make_scat53_block, edef, e_errdef, smdef, rotdef

  edef = -9999.                 ; default value for e1 and e2 if not meas.
  e_errdef = 9999.
  rotdef = 9999.                ; default rotation.
  smdef = 9999.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Cuts and floors
  ;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; smear polarizeability cut.
  IF n_elements(hardrcut) EQ 0 THEN hardrcut = !hardrcut 

  ;; cut in galaxy probability
  IF n_elements(hardprobcut) EQ 0 THEN hardprobcut = !hardprobcut

  ;; maximum error on _each_ ellipticity. 
  ;; A harsher cut can be made later on combined error
  IF n_elements(err_cut) EQ 0 THEN err_cut = 0.4*sqrt(nclrs)
  IF n_elements(minerr) EQ 0 THEN minerr = 0.01 ; A floor placed on error for each ellipticity

  IF keyword_set(hirata) THEN BEGIN
      hirstr = '_h' 
      hdr_hir = 'yes'
  ENDIF ELSE BEGIN 
      hirstr = ''
      hdr_hir = 'no'
  ENDELSE 

  ;; Strings for run/rerun
  rr = ntostr(run)
  rrr = ntostr(rerun)

  ;; Camcols to use
  colmin = 1                    ; What columns to use
  colmax = 6
                                
  colors = ['u','g','r','i','z']

  ;; r-band mag limits.  14.5 is lower limit for spectro sample
  maxmag = 22.
  minmag = 14.5

  ;; Arrays to hold first/last lambda's
  flam_arr = dblarr(6)
  llam_arr = dblarr(6)
  
  crflag = 0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; structures containing the rsmear_cut vs. rmag relation
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  read_rsmear_cuts, run, rerun, purity, 1, grcutstruct, status=gstatus
;  read_rsmear_cuts, run, rerun, purity, 2, rrcutstruct, status=rstatus
;  read_rsmear_cuts, run, rerun, purity, 3, ircutstruct, status=istatus
  
;  read_rsmear_cuts, run, rerun, purity, 1, grcutstruct_h, status=gstatus, $
;                    /hirata
;  read_rsmear_cuts, run, rerun, purity, 2, rrcutstruct_h, status=rstatus, $
;                    /hirata
;  read_rsmear_cuts, run, rerun, purity, 3, ircutstruct_h, status=istatus, $
;                    /hirata

;  IF (gstatus NE 0) OR (rstatus NE 0) OR (istatus NE 0) THEN return

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Catalog of corrected shapes and positions.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  typ = 'blah1'+ntostr(long(systime(1)))
  val = -9999d
  arrval = replicate(val, 5)
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
                         'combine_flags',0,$
                         'runcombine_flags',0,$
                         'nave',0b,$
                         'nrunave', 0b, $
                         'clambda', val, $
                         'ceta', val, $
                         'rsmear',0.0, $
                         'objc_prob_psf', 0.0, $
                         'upetro', 0.0, $
                         'gpetro', 0.0, $
                         'rpetro', 0.0, $
                         'ipetro', 0.0, $
                         'zpetro', 0.0, $
                         'rpetrorad', 0.0, $
                         'grmodel',0.0, $
                         'rimodel',0.0, $
                         'seeing', fltarr(nclrs),$
                         'rotation', 0.0, $
                         'photoz_z', 0.0, $
                         'photoz_type', 0.0, $
                         'photoz_zerr', 0.0, $
                         'photoz_covar_tz', 0.0, $
                         'photoz_typeerr', 0.0, $
                         'photoz_abscounts', arrval, $
                         'photoz_quality', 0)  

  sdssidl_setup, /silent
  setup_mystuff

  IF n_elements(outdir) EQ 0 THEN $
    outdir=sdssidl_config('shapecorr_dir')+'corr'+rr+'/'+rrr+'/combined/'
  outname = $
    outdir + 'run'+rr+'_'+rrr+'_srcgal_'+fclrstring+hirstr+'.fit'

  print
  print,'------------------------------------'
  print,'hard Rsmear cut: '+ntostr(hardrcut)
  print,'Error cut: '+ntostr(err_cut)
  print,'Minimum error floor: '+ntostr(minerr)
  print
  print,'Output file = ',outname
  print,'------------------------------------'
  print

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

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; structures containing the rsmear_cut vs. rmag relation
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;      read_rsmear_cuts, run, rerun, purity, 1, grcutstruct, $
;                        status=gstatus, camcol=camcol
;      read_rsmear_cuts, run, rerun, purity, 2, rrcutstruct, $
;                        status=rstatus, camcol=camcol
;      read_rsmear_cuts, run, rerun, purity, 3, ircutstruct, $
;                        status=istatus, camcol=camcol

;      IF (gstatus NE 0) OR (rstatus NE 0) OR (istatus NE 0) THEN BEGIN 
;          message,'Cannot find rsmear cut struct'
;      ENDIF 

;      read_rsmear_cuts, run, rerun, purity, 1, grcutstruct_h, $
;                        status=gstatus, camcol=camcol, /hirata
;      read_rsmear_cuts, run, rerun, purity, 2, rrcutstruct_h, $
;                        status=rstatus, camcol=camcol, /hirata
;      read_rsmear_cuts, run, rerun, purity, 3, ircutstruct_h, $
;                        status=istatus, camcol=camcol, /hirata

;      IF (gstatus NE 0) OR (rstatus NE 0) OR (istatus NE 0) THEN BEGIN 
;          message,'Cannot find rsmear cut struct'
;      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Use in correcting for slope in egal vs epsf
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      read_corshape, run, rerun, camcol, 1, gslopestruct, hirata=hirata
      read_corshape, run, rerun, camcol, 2, rslopestruct, hirata=hirata
      read_corshape, run, rerun, camcol, 3, islopestruct, hirata=hirata

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
              ;; make cut on status and flag cuts
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              make_scat53_meane_flagselect, tmp, flag_struct, statusstr, clrs, gind, ngal
 
              IF ngal NE 0 THEN BEGIN
 
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; Select objects by r-band magnitude.  
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  ;; Reddening correct
                  tmp.petrocounts = tmp.petrocounts-tmp.reddening
                  tmp.counts_model = tmp.counts_model-tmp.reddening
                  rmag = tmp.petrocounts[rind]

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; Choose based on r magnitude
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  w = where( rmag[gind] LT maxmag AND rmag[gind] GE minmag, ngal)
                  IF ngal NE 0 THEN BEGIN 

                      gind = gind[w]

                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                      ;; Choose based on the rsmear_cut vs. rmag  relationship
                      ;; do each of the requested bandpasses separately. Also
                      ;; check for rotation.  Apply recorrection for slope of
                      ;; egal vs. epsf and apply rotation.  Also, correct for
                      ;; dilution
                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;                      make_scat53_meane_rselect_varcut, tmp, gind, clrs, $
;                                                        hardrcut, hardprobcut, $
;                                                        grcutstruct, rrcutstruct, ircutstruct, $
;                                                        gslopestruct, rslopestruct, islopestruct, $
;                                                        e1, e2, e1e1err, e1e2err, e2e2err, $
;                                                        e1_recorr, e2_recorr, $
;                                                        smear, corr, $
;                                                        w, ngal
                      
                      make_scat53_meane_rselect, tmp, gind, clrs, $
                                                 hardrcut, hardprobcut, $
                                                 gslopestruct, rslopestruct, islopestruct, $
                                                 e1, e2, e1e1err, e1e2err, e2e2err, $
                                                 e1_recorr, e2_recorr, $
                                                 smear, corr, $
                                                 w, ngal, hirata=hirata

                      IF ngal NE 0 THEN BEGIN 
                          gind=gind[w]
                              
                          meanseeing=mean_check(tmp[gind].seeing[rind])

                                                            
                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                          ;; find best ellipticity (should use dilution 
                          ;; corrected shapes and errors)
                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                          
                          IF nclrs GT 1 THEN BEGIN 
                              combine_ellip_cove1e2, e1, e2, e1e1err, e1e2err, e2e2err, smear, $
                                                     new_e1, new_e2, $
                                                     new_e1e1err, new_e1e2err, new_e2e2err,$
                                                     new_smear, $
                                                     combine_flags, nave, gind, verbose=0, $
                                                     err_cut=err_cut, rcut=hardrcut,$
                                                     besterror=besterror, wmomerror=wmomerror
                              
                              combine_ellip_cove1e2, e1_recorr, e2_recorr, e1e1err, e1e2err, e2e2err, smear, $
                                                     new_e1_recorr, new_e2_recorr, $
                                                     tnew_e1e1err, tnew_e1e2err, tnew_e2e2err,$
                                                     tnew_smear, $
                                                     combine_flags_recorr, tnave, gind_recorr, verbose=0, $
                                                     err_cut=err_cut, rcut=hardrcut,$
                                                     besterror=besterror, wmomerror=wmomerror
                          ENDIF ELSE BEGIN 
                              new_e1 = e1 & new_e1_recorr = e1_recorr
                              new_e2 = e2 & new_e2_recorr = e2_recorr
                              new_e1e1err = e1e1err
                              new_e1e2err = e1e2err
                              new_e2e2err = e2e2err

                              new_smear = smear
                          ENDELSE 
                              
                          IF gind[0] NE -1 THEN BEGIN 
                              
                              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                              ;; Final cuts have been made. Copy 
                              ;; into structure
                              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                              ;; "gind" should be same for both
                              ngal = n_elements(gind)
                              ngal_recorr = n_elements(gind_recorr)
                              IF ngal NE ngal_recorr THEN message,'WHAT!!! ngal ne ngal_recorr'
                              
                              srctmp = replicate(struct, ngal)
                              srctmp.seeing = tmp[gind].seeing[clrs]
                              
                              srctmp.e1 = new_e1[gind]
                              srctmp.e2 = new_e2[gind]
                              srctmp.e1_recorr = new_e1_recorr[gind]
                              srctmp.e2_recorr = new_e2_recorr[gind]
                              srctmp.e1e1err = new_e1e1err[gind]
                              srctmp.e1e2err = new_e1e2err[gind]
                              srctmp.e2e2err = new_e2e2err[gind]
                              
                              srctmp.combine_flags = combine_flags[gind]
                              srctmp.nave = nave[gind]
                              
                              ;; close enough, most likely
                              srctmp.rotation = tmp[gind].rotation[rind]
                              
                              eq2csurvey, tmp[gind].ra,tmp[gind].dec,$
                                          clambda, ceta

                              srctmp.clambda = clambda
                              srctmp.ceta = ceta
                              
                              srctmp.rsmear = new_smear[gind]
                              srctmp.objc_prob_psf = tmp[gind].objc_prob_psf

                              srctmp.upetro = tmp[gind].petrocounts[0]
                              srctmp.gpetro = tmp[gind].petrocounts[1]
                              srctmp.rpetro = tmp[gind].petrocounts[2]
                              srctmp.ipetro = tmp[gind].petrocounts[3]
                              srctmp.zpetro = tmp[gind].petrocounts[4]
                              
                              srctmp.rpetrorad = tmp[gind].petrorad[2]

                              srctmp.grmodel = $
                                tmp[gind].counts_model[1]-$
                                tmp[gind].counts_model[2]
                              srctmp.rimodel = $
                                tmp[gind].counts_model[2]-$
                                tmp[gind].counts_model[3]
                                  
                              ;; photoz info
                              srctmp.photoz_z = tmp[gind].photoz_z
                              srctmp.photoz_type = tmp[gind].photoz_type
                              srctmp.photoz_zerr = tmp[gind].photoz_zerr
                              srctmp.photoz_covar_tz = $
                                tmp[gind].photoz_covar_tz
                              srctmp.photoz_typeerr = tmp[gind].photoz_typeerr
                              srctmp.photoz_quality = tmp[gind].photoz_quality
                              srctmp.photoz_abscounts = $
                                tmp[gind].photoz_abscounts         
           
                              srctmp.run = tmp[gind].run
                              srctmp.rerun = tmp[gind].rerun
                              srctmp.camcol = tmp[gind].camcol
                              srctmp.field = tmp[gind].field
                              srctmp.id = tmp[gind].id
                              
                              ;; copy into pointers
                              numlist[ci, fi] = ngal
                              ntotal[ci] = ntotal[ci] + ngal
                              ptrlist[ci, fi] = ptr_new(srctmp, /no_copy)
                              
                              IF fi MOD 20 EQ 0 THEN BEGIN 
                                  print,'Run: ',rr,' Camcol: ',cstr,$
                                        ' Field: ',fstr,'  ',$
                                        'Src Galaxies: '+$
                                        ntostr(ngal)
                                  print,'Seeing['+colors[rind]+'] = ',$
                                        ntostr(meanseeing)+'"',$
                                        ' Mean e1: ',ntostr(mean_check(new_e1[gind]))
                              ENDIF 
                              
                              ;; free some memory
                              tmp=0 & corr=0
                              e1=0 & e2=0 & e1_recorr=0 & e2_recorr=0
                              new_e1=0 & new_e2=0 & new_e1_recorr=0 & new_e2_recorr=0
                              e1e1err=0 & e1e2err=0 & e2e2err=0
                              new_e1e1err=0 & new_e1e2err=0 & new_e2e2err=0
                              combine_flags=0 & nave=0
                              rotation=0 & clambda=0 & ceta=0
                          ENDIF ELSE BEGIN 
                              message,'No galaxies passed final rcut and err_cut',/inf
                          ENDELSE 

                          
                      ENDIF ELSE BEGIN
                          message,'No galaxies passed size '+$
                                  'and rotation cuts in field '+fstr,/inf
                      ENDELSE 
                  ENDIF ELSE BEGIN
                      message,'No galaxies in magnitude range '+$
                              colors[rind]+' = ['+ntostr(minmag)+', '+$
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
      print,'Run: ',rr,' Camcol: ',cstr,' total src galaxies: ',$
            ntostr(ntotal[ci])
      print

  ENDFOR ;; loop over camcols

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
      wt = where(numlist[ci, *] NE 0, nfields)
      FOR ff=0, nfields-1 DO BEGIN 

          fi = wt[ff]

          ;; Don't use fields flagged as bad
          IF (numlist[ci, fi] NE 0) THEN BEGIN 
              sources[beg:beg+numlist[ci, fi]-1]=*ptrlist[ci, fi]
          ENDIF 
          ptr_free,ptrlist[ci, fi]   ;kill the heap variables associated 
                                ;with the pointer
          beg=beg+numlist[ci, fi]

      ENDFOR 
  ENDFOR 

  ;; get first clambda, last clambda
  firstlam = min(sources.clambda,max=lastlam)

  sxaddpar, outhdr, 'RCUT', hardrcut, ' Hard Rsmear cut in addition to magnitude dependent cut'
  sxaddpar, outhdr, 'ERR_CUT', err_cut, ' Maximum error cut for individual ellipticity measurements'
  sxaddpar, outhdr, 'MINERR', minerr, ' Floor on error for ellipticity measurements (not a cut)'
;  sxaddpar, outhdr, 'PURITY', purity, ' Purity cut for sg separation'
  sxaddpar, outhdr, 'HIRATA', hdr_hir, ' Did we use the Hirata method?'

  sxaddpar, outhdr, 'FIRSTLAM', firstlam, ' First clambda in this data'
  sxaddpar, outhdr, 'LASTLAM', lastlam, ' Last clambda in this data'

  ;; add CMETHOD to header: how did we calculate
  ;; the covariance matrix?

  IF keyword_set(besterror) THEN BEGIN 
      cmethod = 'best'
  ENDIF ELSE IF keyword_set(wmomerror) THEN BEGIN 
      cmethod = 'waverage'
  ENDIF ELSE BEGIN
      cmethod = 'default'
  ENDELSE 

  sxaddpar, outhdr, 'AVGCLRS', clrstring, ' Bandpasses used to get mean ellipticity'
  sxaddpar, outhdr, 'CMETHOD', cmethod, $
            ' Method used to calculate covariance matrix'

  print,'FIRSTLAM: ',firstlam
  print,'LASTLAM: ',lastlam
  print

  print,'Outfile: ',outname
  mwrfits2, sources, outname, /destroy, /create, hdr0=outhdr

  print
  print,'Total source galaxies: ',ntostr(ntotal)

  ptime,systime(1)-time

  return
END 
