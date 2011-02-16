PRO adatc_dbbuild_match2mass, pstruct, psc_cat, xsc_cat

  TWOMASS_EXT = 2b^3
  TWOMASS_POINT = 2b^4

  ;; code from Ryan's add2mass2adatc_byrun.pro
  tol = 2.0d/3600.0

  close_match_radec, pstruct.ra, pstruct.dec, $
                     psc_cat.ra, psc_cat.dec, $
                     pmatch, psc_match,       $
                     tol, 1, miss1, /silent

  IF pmatch[0] NE -1 THEN BEGIN 
      pstruct[pmatch].j_2mass    = psc_cat[psc_match].j
      pstruct[pmatch].j_2masserr = psc_cat[psc_match].dj
      pstruct[pmatch].h_2mass    = psc_cat[psc_match].h
      pstruct[pmatch].h_2masserr = psc_cat[psc_match].dh
      pstruct[pmatch].k_2mass    = psc_cat[psc_match].k
      pstruct[pmatch].k_2masserr = psc_cat[psc_match].dk

      pstruct[pmatch].value_flags = $
        pstruct[pmatch].value_flags + TWOMASS_POINT

      npsc = n_elements(pmatch)
  ENDIF ELSE npsc = 0

  close_match_radec, pstruct.ra, pstruct.dec, $
                     xsc_cat.ra, xsc_cat.dec, $
                     pmatch, xsc_match,       $
                     tol, 1, miss1, /silent

  IF pmatch[0] NE -1 THEN BEGIN 
      pstruct[pmatch].j_2mass    = xsc_cat[xsc_match].j
      pstruct[pmatch].j_2masserr = xsc_cat[xsc_match].dj
      pstruct[pmatch].h_2mass    = xsc_cat[xsc_match].h
      pstruct[pmatch].h_2masserr = xsc_cat[xsc_match].dh
      pstruct[pmatch].k_2mass    = xsc_cat[xsc_match].k
      pstruct[pmatch].k_2masserr = xsc_cat[xsc_match].dk

      pstruct[pmatch].value_flags = $
        pstruct[pmatch].value_flags + TWOMASS_EXT

      nxsc = n_elements(pmatch)
  ENDIF ELSE nxsc = 0

  print,'2MASS matches: '+ntostr(npsc+nxsc)

END 

PRO adatc_dbbuild_read2mass, run, psc_cat, xsc_cat

  ;; code from Ryan's add2mass2adatc_byrun.pro

  ;; get the stripe for this run

  wrst = where(!run_status.run EQ run AND $
               !run_status.stripe NE -1, nwrst)
  IF nwrst EQ 0 THEN message,'No stripe info for run '+ntostr(run)
  stripe = !run_status[wrst[0]].stripe
  
  infile = '/home/caruso1/scranton/2MASS/psc_sdss.fit'
  print,'Reading in 2MASS data and selecting stripe ',ntostr(stripe)
  print,infile
  psc_cat = mrdfits(infile, 1, hdr, /silent)
  
  eq2csurvey,psc_cat.ra,psc_cat.dec,lambda,eta
  
  psc_stripe = eta2stripenum(eta)
  
  w = where(psc_stripe GE stripe-1 AND psc_stripe LE stripe+1,n_psc)
  
  IF w[0] EQ -1 THEN BEGIN
      print,'Cannot find data in psc_sdss.fit for stripe ',ntostr(stripe)
      RETURN
  ENDIF ELSE BEGIN
      psc_cat = psc_cat(w)
      psc_stripe = 0 & lambda = 0 & eta = 0
  ENDELSE
  
  xsc_cat = mrdfits('/home/caruso1/scranton/2MASS/xsc_sdss.fit',1,$
                    hdr,/silent)
  
  eq2csurvey,xsc_cat.ra,xsc_cat.dec,lambda,eta
  
  xsc_stripe = eta2stripenum(eta)
  
  w = where(xsc_stripe GE stripe-1 AND xsc_stripe LE stripe+1,n_xsc)
  
  IF w[0] EQ -1 THEN BEGIN
      print,'Cannot find data in xsc_sdss.fit for stripe ',ntostr(stripe)
      RETURN
  ENDIF ELSE BEGIN
      xsc_cat = xsc_cat(w)
      xsc_stripe = 0 & lambda = 0 & eta = 0
  ENDELSE
  
END 

PRO adatc_dbbuild_matchphotoz, pstruct, photoz, wff

  PHOTOZ_FLAG = 2b^1

  ;; using default sky and first field
  unique_id,pstruct.run,pstruct.rerun,pstruct.camcol,$
            pstruct.field,pstruct.id,uid,/silent
  
  match, photoz[wff].unique_id, uid, wphotoz, wpstruct, /sort
  
  IF wphotoz[0] NE -1 THEN BEGIN 
      wphotoz = wff[wphotoz]
      
      ;;help,wff,wphotoz
      
      pstruct[wpstruct].photoz_z = photoz[wphotoz].z
      pstruct[wpstruct].photoz_zerr = photoz[wphotoz].zerr
      
      pstruct[wpstruct].photoz_type = photoz[wphotoz].t
      pstruct[wpstruct].photoz_typeerr = photoz[wphotoz].terr
      
      pstruct[wpstruct].photoz_covar_tt = photoz[wphotoz].covar_tt
      pstruct[wpstruct].photoz_covar_tz = photoz[wphotoz].covar_tz
      pstruct[wpstruct].photoz_covar_zz = photoz[wphotoz].covar_zz
      
      pstruct[wpstruct].photoz_chisq = photoz[wphotoz].chisq
      pstruct[wpstruct].photoz_quality = photoz[wphotoz].quality
      pstruct[wpstruct].photoz_dist_mod = photoz[wphotoz].distmod
      
      pstruct[wpstruct].photoz_kcorr[0] = photoz[wphotoz].kcorr_u
      pstruct[wpstruct].photoz_kcorr[1] = photoz[wphotoz].kcorr_g
      pstruct[wpstruct].photoz_kcorr[2] = photoz[wphotoz].kcorr_r
      pstruct[wpstruct].photoz_kcorr[3] = photoz[wphotoz].kcorr_i
      pstruct[wpstruct].photoz_kcorr[4] = photoz[wphotoz].kcorr_z
      
      pstruct[wpstruct].photoz_abscounts[0] = photoz[wphotoz].abscounts_u
      pstruct[wpstruct].photoz_abscounts[1] = photoz[wphotoz].abscounts_g
      pstruct[wpstruct].photoz_abscounts[2] = photoz[wphotoz].abscounts_r
      pstruct[wpstruct].photoz_abscounts[3] = photoz[wphotoz].abscounts_i
      pstruct[wpstruct].photoz_abscounts[4] = photoz[wphotoz].abscounts_z
      
      ;; set photoz flag
      wp2 = where(pstruct[wpstruct].photoz_quality GT 0 AND $
                  pstruct[wpstruct].photoz_quality NE 12, nwp2)
      IF nwp2 NE 0 THEN BEGIN 
          pstruct[wpstruct[wp2]].value_flags = $
            pstruct[wpstruct[wp2]].value_flags + PHOTOZ_FLAG
      ENDIF 

  ENDIF ELSE print,'No matching PHOTOZs for field: '+ntostr(pstruct[0].field)
    
END 

PRO adatc_dbbuild_getsmear, pstruct, indices, clr, corr

  corr = pstruct[indices].m_rr_cc_psf[clr]/pstruct[indices].m_rr_cc[clr]*(4./pstruct[indices].m_cr4_psf[clr] -1.)/(4./pstruct[indices].m_cr4[clr]-1.)

END 

PRO adatc_dbbuild_select, select_clr, pstruct, maxmag, ss, nss

  make_flag_struct,fs
  fs.satur='N'
  fs.satur_center='N'
  fs.bright = 'N'
  fs.nopetro_big = 'N'


  flag_select,pstruct,fs,select_clr,ss1
  IF (ss1[0] EQ -1) THEN BEGIN
      print,'No objects passed flag cuts!'
      GOTO,jump
  ENDIF 
  
  make_flag_struct, fs
  fs.deblended_as_moving='N'

  ;; Want the catalogs to be more general
  ;; fs.deblended_as_psf = 'N'
  flag_select,pstruct[ss1],fs,select_clr,ss2,/objc
  IF (ss2[0] EQ -1) THEN BEGIN
      print,'No objects passed flag cuts!'
      ss1 = ss2
      GOTO,jump
  ENDIF 
  ss1 = ss1[ss2]

;  ss2 = where( ( pstruct[ss1].petrocounts[select_clr] $
;                  - pstruct[ss1].reddening[select_clr] ) LE  maxmag, nss2)

  ss2 = where( ( pstruct[ss1].cmodel_counts[select_clr] $
                  - pstruct[ss1].reddening[select_clr] ) LE  maxmag, nss2)

  
  ;; Undo for an OR statement
;  mag = pstruct[ss1].petrocounts - pstruct[ss1].reddening

;  ss2 = where( ( mag[0,*] LE maxmag ) OR $
;               ( mag[1,*] LE maxmag ) OR $
;               ( mag[2,*] LE maxmag ) OR $
;               ( mag[3,*] LE maxmag ) OR $
;               ( mag[4,*] LE maxmag ), nss2)


;  mag = 0

  IF nss2 EQ 0 THEN BEGIN
      print,'No objects passed flag cuts!'
      ss1 = ss2
      GOTO,jump       
  ENDIF 
  ss1 = ss1[ss2]
  
  jump:
  ;; also include any object which has primtarget set
  ss3 = where(pstruct.primtarget NE 0,nss3)

  ;; any with primtarget set?
  IF nss3 EQ 0 THEN ss = ss1 ELSE BEGIN 
      ;; any pass flag cuts?
      IF ss1[0] EQ -1 THEN ss=ss3 ELSE BEGIN 
          ;; combine and remove duplicates
          ss = [ss1, ss3]
          ss = ss[rem_dup(ss)]
      ENDELSE 
  ENDELSE 

  IF ss[0] EQ -1 THEN nss=0L ELSE nss=n_elements(ss)

  return

END 

PRO adatc_dbbuild, run, camcol, rerun=rerun, $
                   status=status, outdir=outdir, $
                   addphotoz=addphotoz, addbayes=addbayes, $
                   add2mass=add2mass, $
                   photoz=photoz, psc_cat=psc_cat, xsc_cat=xsc_cat
  

  ;; status is 1 unless we reach end
  status =1

  IF n_params() LT 2 THEN BEGIN ;Help message
      print,'-syntax '
      return 
  ENDIF 

  time=systime(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  maxmag = 22.5                 ; 

  setup_mystuff                 ; set up my variables

  newfront = 'adatc'            ; The front name for output files

  bands = [0,1,2,3,4]           ; bandpasses to correct
  nband = n_elements(bands)

  colors=['u','g','r','i','z']  ; some strings for later use
  runstr = ntostr(run)
  rerunstr=ntostr(rerun)
  cstr = ntostr(camcol)

  ;; flags for value-added stuff
  BAYES_FLAG = 2b^0
  PHOTOZ_FLAG = 2b^1
  PHOTOQSO_FLAG = 2b^2
  TWOMASS_EXT = 2b^3
  TWOMASS_POINT = 2b^4
  FIRST_FLAG = 2b^5
  ROSAT_FLAG = 2b^6

  ;; One of the bayes flags set in compute_bayes_prob
  PROBFLAG_NOPROB = 2b^2

  make_flag_struct, fs          ;for selecting good moments
  fs.AMOMENT_FAINT = 'N'
  fs.AMOMENT_UNWEIGHTED = 'N'
  fs.AMOMENT_SHIFT = 'N'
  fs.AMOMENT_MAXITER = 'N'

  print
  print,'Using maxmag = '+ntostr(maxmag)
  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; versions of photoz and probgal to place in headers
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  bayes_version  = 'v1_5_0'
  photoz_version = 'v1_5_0'
  IF NOT keyword_set(addbayes) THEN bayes_version = 'NOBAYES'
  IF NOT keyword_set(addphotoz) THEN photoz_version = 'NOPHOTOZ'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; tags to get from tsObj, and the extra structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  make_admomtags, taglist, /default
  make_lensstruct, lensstr

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; read the file names and get rotation directory
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fetch_dir, run, camcol, rerun, dir, atldir, corrdir=corrdir, $
             corratldir=rotdir  ;We keep rotation files in corratldir          
  fetch_file_list, dir, files, fnums, start=start, nframes=nframes

  nfields = n_elements(files)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; should we attempt to add PHOTOZ?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(addphotoz) THEN BEGIN 
      photoz_dir = sdssidl_config('photoz_dir')

      IF n_elements(photoz) EQ 0 THEN BEGIN 
          phfile = photoz_dir + $
            'tsObj_ascii_'+runstr+'_'+rerunstr+'_'+cstr+'.dat'
          
          read_photoz_ascii, phfile, photoz, status=phstatus

          IF phstatus NE 0 THEN BEGIN
              print,'NOT DOING PHOTOZ!!!!'
              photoz_version = 'NOPHOTOZ'
              addphotoz = 0
          ENDIF 
      ENDIF 

  ENDIF 
  ;; did we find photoz?
  IF keyword_set(addphotoz) THEN BEGIN 

      ;; get the field from photoz struct
      extract_from_id, photoz.unique_id, 'field', phzfields
      
      ;; lots of objects! Use histogram on fields, from zero (so
      ;; we can subscript rev_ind with field) to the maximum 
      ;; possible field (no more than 2000)
      hist=histogram(phzfields, min=0L, max=2000L, reverse=rev_ind)
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; should we match to 2mass?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(add2mass) THEN BEGIN 
      IF n_elements(psc_cat) EQ 0 AND n_elements(xsc_cat) EQ 0 THEN BEGIN 
          adatc_dbbuild_read2mass, run, psc_cat, xsc_cat
      ENDIF ELSE BEGIN 
          print,'Using input 2mass catalogs'
      ENDELSE 
      IF n_elements(psc_cat) EQ 0 OR n_elements(xsc_cat) EQ 0 THEN BEGIN 
          print,'Returning'
          return
      ENDIF 

  ENDIF 

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

  rotname = rotdir + 'surveyrot_'+run2string(run)+'_'+ntostr(camcol)+'_'+$
    ntostr(rerun)+'.fit'
  rotstruct = mrdfits(rotname,1,/silent)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over fields and find psf moments and correct shapes
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Output files are '+newfront+'-*'

  FOR ic = 0L, nfields-1 DO BEGIN 

      infile = files[ic]
      outfile = outfiles[ic]
      field = fnums[ic]
      fstr=ntostr(field)

      ;; read headers
      hdr0=headfits(infile, exten=0)
      hdr1=headfits(infile, exten=1)
      fxhclean, hdr1

      ;; add version numbers
      fxaddpar, hdr0, 'BAYE_VER', bayes_version, $
                ' Version of bayesian s/g separation code used.'
      fxaddpar, hdr0, 'PHTZ_VER', photoz_version, $
                ' Version of template photoz code used.'

      ;; add date this pipeline run
      fxaddpar, hdr0, 'MCDATE', systime(), $
                ' Date adatc_dbbuild.pro was run'

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Read tsObj file
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      read_tsobj, dir, pstruct, start=field, nframes=1, $
                  taglist=taglist, tsobjstr=tsobjstr, ex_struct=lensstr, $
                  /noadd,verbose=0

      nps = n_elements(pstruct)
      print
      print,'Run: ',runstr,' Camcol: ',cstr,' Field: ',fstr,' Nobj: ',$
            ntostr(nps)
      IF nps NE 0 THEN BEGIN 

          make_cmodel_counts, pstruct, cmodel, cmodelerr
          
          pstruct.cmodel_counts = cmodel
          pstruct.cmodel_countserr = cmodelerr

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; matches to other catalogs
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          wfirst = where(pstruct.firstmatch GT 0, nfirst)
          IF nfirst NE 0 THEN BEGIN 
              pstruct[wfirst].value_flags = $
                pstruct[wfirst].value_flags + FIRST_FLAG
          ENDIF 
          wrosat = where(pstruct.rosatmatch GT 0, nrosat)
          IF nrosat NE 0 THEN BEGIN 
              pstruct[wrosat].value_flags = $
                pstruct[wrosat].value_flags + ROSAT_FLAG
          ENDIF 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; don't need these tags any more.  Create a struct without
          ;; them for output
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF n_elements(repstruct) EQ 0 THEN BEGIN 
              tt = pstruct[0]
              repstruct = remove_tags(tt, ['FIRSTMATCH','ROSATMATCH'])
              delvarx, tt
          ENDIF 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Rotation (one per field for now)
          ;; do we have rotation? Use later
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          wrot = where(rotstruct.field EQ field,nrot)
          IF nrot EQ 0 THEN print,'No rotation for this field!!'

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; select objects
          ;; first of all, get rid of parents
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          wp = where(pstruct.nchild EQ 0, nwp)
          IF nwp NE 0 THEN BEGIN 
              pstruct = temporary(pstruct[wp])

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; now flag and magnitude selection
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              flags = pstruct.corrselect_flags

              adatc_dbbuild_select, 0, pstruct, maxmag, uss, unss
              IF unss NE 0 THEN flags[uss] = flags[uss] + 2b^0
                  
              adatc_dbbuild_select, 1, pstruct, maxmag, gss, gnss
              IF gnss NE 0 THEN flags[gss] = flags[gss] + 2b^1
                  
              adatc_dbbuild_select, 2, pstruct, maxmag, rss, rnss
              IF rnss NE 0 THEN flags[rss] = flags[rss] + 2b^2
                  
              adatc_dbbuild_select, 3, pstruct, maxmag, iss, inss
              IF inss NE 0 THEN flags[iss] = flags[iss] + 2b^3
                  
              adatc_dbbuild_select, 4, pstruct, maxmag, zss, znss
              IF znss NE 0 THEN flags[zss] = flags[zss] + 2b^4
                  
              ss = where(flags GT 0, nss)
              
              IF nss NE 0 THEN BEGIN

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; we made it this far: objects will be written
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  pstruct = temporary(pstruct[ss])
                  pstruct.corrselect_flags = flags[ss]

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; Match PHOTOZ struct
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          
                  IF keyword_set(addphotoz) THEN BEGIN
                      IF rev_ind[field] NE rev_ind[field+1] THEN BEGIN 
                          wff = rev_ind[ rev_ind[field]:rev_ind[field+1] -1 ]
                          ;;print,'Using field: ',phzfields[wff[0]]
                          adatc_dbbuild_matchphotoz, pstruct, photoz, wff
                      ENDIF ELSE print,'No PHOTOZs for field '+fstr
                  ENDIF 

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; make "corrected" survey coordinates
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  eq2csurvey, pstruct.ra, pstruct.dec, clambda, ceta
                  pstruct.clambda = temporary(clambda)
                  pstruct.ceta = temporary(ceta)

                  FOR iband=0L, nband-1 DO BEGIN 

                      clr = bands[iband]

                      ;; seeing in arcseconds
                      wsee = where(pstruct.m_rr_cc_psf[clr] GT 0.0,nsee)
                      IF nsee NE 0 THEN BEGIN 
                          pstruct[wsee].seeing[clr] = $
                            sqrt(pstruct[wsee].m_rr_cc_psf[clr]/2.)*2.35*0.4
                      ENDIF 

                      ;; objects with no good psf measurement: get seeing
                      ;; of nearest neighborz

                      IF nsee NE nss THEN BEGIN 
                          wseeb = lindgen(nss)
                          remove, wsee, wseeb

                          close_match_radec, $
                            pstruct[wseeb].ra, pstruct[wseeb].dec, $
                            pstruct[wsee].ra, pstruct[wsee].dec, $
                            mseeb, msee, 50.0/3600., 1, /silent
                          IF msee[0] NE -1 THEN BEGIN 
                              ;; copy in closest object's seeing
                              pstruct[wseeb[mseeb]].seeing[clr] = $
                                pstruct[wsee[msee]].seeing[clr]
                          ENDIF 
                      ENDIF 

                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                      ;; Only correct those with rotation
                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                      wrot = where(rotstruct.field EQ field,nrot)
                      IF nrot NE 0 THEN BEGIN 
                          pstruct.rotation = rotstruct[wrot].angle

                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                          ;; find good moment measurements
                          ;; in this bandpass
                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                          flag_select, pstruct, fs, clr, good
                          IF good[0] NE -1 THEN BEGIN

                              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                              ;; still need to check not -1000.
                              ;; PHOTO flags missed these!
                              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                              good2 = $
                                where(pstruct[good].m_e1[clr] NE -1000.,ngood2)
                              IF ngood2 NE 0 THEN BEGIN 
                                  good = good[good2]
                                  ngood = n_elements(good)

                                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                  ;; check for NAN
                                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                  good3 = where( (pstruct[good].m_e1e1err[clr] EQ $
                                                  pstruct[good].m_e1e1err[clr]) AND $
                                                 (pstruct[good].m_e1e2err[clr] EQ $
                                                  pstruct[good].m_e1e2err[clr]) AND $
                                                 (pstruct[good].m_e2e2err[clr] EQ $
                                                  pstruct[good].m_e2e2err[clr]),ngood3)

                                  IF ngood3 NE 0 THEN BEGIN 

                                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                      ;; Get smear polarizability, make sure
                                      ;; it is reasonable
                                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                      adatc_dbbuild_getsmear,pstruct,good,clr, $
                                                              corr

                                      good4 = where(corr GT 0.0, ngood4)
                                      IF ngood4 NE 0 THEN BEGIN 

                                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                          ;; correct shapes, copy in new fields
                                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                          good = good[good4]
                                          ngood = ngood4
                                          corr = corr[good4]
                                          
                                          pstruct[good].m_e1_corr[clr] = $
                                            pstruct[good].m_e1[clr] - $
                                            corr*pstruct[good].m_e1_psf[clr]
                                          
                                          pstruct[good].m_e2_corr[clr] = $
                                            pstruct[good].m_e2[clr] - $
                                            corr*pstruct[good].m_e2_psf[clr]
                                          
                                          pstruct[good].m_R[clr] = corr
                                          
                                      ENDIF ELSE ngood=0L ;corr > 0?

                                      ;; Chris Hirata's correction scheme
                                      compea4_struct, pstruct[good], clr, $
                                                      e1_out, e2_out, R_out, compflags
                                      pstruct[good].CompEA4flags[clr] = compflags
                                      good5 = where(compflags EQ 0L, ngood5)
                                      IF ngood5 NE 0 THEN BEGIN 
                                          tgood = good[good5]
                                          pstruct[tgood].m_e1_corr_h[clr] = e1_out[good5]
                                          pstruct[tgood].m_e2_corr_h[clr] = e2_out[good5]
                                          pstruct[tgood].m_R_h[clr] = R_out[good5]
                                      ENDIF 
                                  ENDIF ELSE ngood=0L ;NAN cut

                              ENDIF ELSE ngood=0L ;; e1 cut

                          ENDIF ELSE ngood=0L ;; AMOMENT flags

                          print,'Band: '+colors[clr]+' Ngood: '+ntostr(ngood)

                      ENDIF ;; rotation?

                  ENDFOR ;; loop over shape bandpasses

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; Bayesian Galaxy Probabilities, which
                  ;; depend on the seeing
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  IF keyword_set(addbayes) THEN BEGIN 
                      compute_bayes_prob, pstruct, probgal, probflags

                      pstruct.objc_prob_psf = 1. - probgal
                      pstruct.objc_prob_flags = probflags
                      wgood_prob = $
                        where((probflags AND PROBFLAG_NOPROB) EQ 0, ngood_prob)
                      IF ngood_prob NE 0 THEN BEGIN 
                          pstruct[wgood_prob].value_flags = $
                            pstruct[wgood_prob].value_flags + BAYES_FLAG
                      ENDIF 
                  ENDIF 

                  IF keyword_set(add2mass) THEN BEGIN 
                      adatc_dbbuild_match2mass, pstruct, psc_cat, xsc_cat
                  ENDIF 

              ENDIF ELSE BEGIN ;; flag objects passed cuts
                  print,'No objects written: Flag cuts'
                  delvarx, pstruct
              ENDELSE 

          ENDIF ELSE BEGIN  ;; no objects with nchild == 0
              print,'No objects written: All parents!' 
              delvarx, pstruct
          ENDELSE 
      ENDIF ELSE BEGIN ;; no objects in file
          print,'No objects written: empty file'
          delvarx,pstruct
      ENDELSE 

      nout = n_elements(pstruct)
      IF nout NE 0 THEN BEGIN 
          ;; copy into struct with extraneous tags removed
          outstruct = replicate(repstruct, nout)
          struct_assign, pstruct, outstruct, /nozero;, /verbose
      ENDIF 

      ;; write file (may be empty)
      mwrfits2, outstruct, outfile, hdr1, /create, /destroy, hdr0=hdr0

  ENDFOR 

  ptime,systime(1)-time

  status=0

END 
