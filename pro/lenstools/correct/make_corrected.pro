PRO make_corrected_getsmear, pstruct, indices, clr, corr

  corr = pstruct[indices].m_rr_cc_psf[clr]/pstruct[indices].m_rr_cc[clr]*(4./pstruct[indices].m_cr4_psf[clr] -1.)/(4./pstruct[indices].m_cr4[clr]-1.)

END 

FUNCTION make_corrected_select, select_clr, pstruct, adatc, maxmag, nused, $
                                input_index=input_index


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some basic flag cuts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  make_flag_struct,fs
  fs.satur='N'
  fs.satur_center='N'
  fs.bright = 'N'
  fs.nopetro_big = 'N'

  flag_select,pstruct,fs,select_clr,used1, input_index=input_index

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; If these are passed, do some objc cuts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF used1[0] NE -1 THEN BEGIN 

      ;; Objc flag cuts
      make_flag_struct, fs
      fs.deblended_as_moving='N'

      ;; Want the catalogs to be more general
      ;; fs.deblended_as_psf = 'N'
      flag_select,pstruct,fs,select_clr,used2,/objc, input_index=used1
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; If these cuts passed, do some magnitude cuts
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF used2[0] NE -1 THEN BEGIN 

          used1=used2

          ;; Magnitude cuts
          used2 = $
            where(adatc[used1].cmodel_counts_ext[select_clr] LE maxmag, nused2)

          IF nused2 NE -1 THEN BEGIN 
              used1 = used1[used2]
          ENDIF ELSE BEGIN 
              print,'No objects passed mag cuts!'
              used1 = used2
          ENDELSE 

      ENDIF ELSE BEGIN 
          print,'No objects passed objc flag cuts!'
          used1 = used2
      ENDELSE 

  ENDIF ELSE BEGIN 
      print,'No objects passed flag cuts in bandpass!'
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; also include any object which has primtarget set
  ;; This object would have to pass the flag cuts, but may
  ;; not pass the magnitude cut.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  used3 = where(pstruct.primtarget NE 0,nused3)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Combined the results
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF nused3 EQ 0 THEN used = used1 ELSE BEGIN 
      ;; any pass main cuts?
      IF used1[0] EQ -1 THEN used=used3 ELSE BEGIN 
          ;; combine and remove duplicates
          used = [used1, used3]
          used = used[rem_dup(used)]
      ENDELSE 
  ENDELSE 

  IF used[0] EQ -1 THEN nused=0L ELSE nused=n_elements(used)

  return,used

END 

FUNCTION make_corrected, pstruct, used, rotstruct=rotstruct, status=status


  ;; status is 1 unless we reach end
  status =1

  IF n_params() LT 1 THEN BEGIN ;Help message
      print,'-syntax adatc = make_corrected(pstruct, used, rotstruct=, status=)'
      return,-1
  ENDIF 

  setup_mystuff                 ; set up my variables

  maxmag = 22.5

  bands = [0,1,2,3,4]           ; bandpasses to correct
  nband = n_elements(bands)

  colors=['u','g','r','i','z']  ; some strings for later use

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

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; versions of probgal to place in headers
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  bayes_version  = 'v1_5_0'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; tags to get from tsObj, and the extra structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  newstruct = make_corrected_struct()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Get the files that contain the rotation angle for each field
  ;; rotation between (row,col) and (lambda,eta)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  run = pstruct[0].run
  camcol = pstruct[0].camcol
  rerun = pstruct[0].rerun
  field = pstruct[0].field

  wfield=where(pstruct.field NE field, nwfield)
  IF nwfield NE 0 THEN BEGIN 
      message,'You may only enter one field at a time',/inf
      return,-1
  ENDIF 

  newstruct.run = run
  newstruct.rerun = rerun
  newstruct.camcol = camcol
  newstruct.field = field

  nps = n_elements(pstruct)

  ;; Start with all, then trim down
  adatc = replicate(newstruct, nps)
  adatc.id = pstruct.id

  ;; The photoid
  adatc.photoid = photoid(pstruct)

  ;; Cmodel counts
  make_cmodel_counts, pstruct, cmodel, cmodelerr;, /defaults

  adatc.cmodel_counts = cmodel
  adatc.cmodel_countserr = cmodelerr

  ;; keep extinction corrected cmodel counts
  adatc.cmodel_counts_ext = cmodel - pstruct.reddening

  ;; make "corrected" survey coordinates
  eq2csurvey, pstruct.ra, pstruct.dec, clambda, ceta
  adatc.clambda = temporary(clambda)
  adatc.ceta = temporary(ceta)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; matches to other catalogs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wfirst = where(pstruct.firstmatch GT 0, nfirst)
  IF nfirst NE 0 THEN BEGIN 
      adatc[wfirst].value_flags = $
        adatc[wfirst].value_flags + FIRST_FLAG
  ENDIF 
  wrosat = where(pstruct.rosatmatch GT 0, nrosat)
  IF nrosat NE 0 THEN BEGIN 
      adatc[wrosat].value_flags = $
        adatc[wrosat].value_flags + ROSAT_FLAG
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Rotation (one per field for now)
  ;; do we have rotation? Use later
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(rotstruct) EQ 0 THEN BEGIN 
      ss = obj_new('sdss_survey_rot')
      rotstruct = ss->read(run,rerun=rerun)
      rotstruct = rotstruct[ where(rotstruct.camcol EQ camcol) ]
  ENDIF 

  wrot = where(rotstruct.field EQ field,nrot)
  IF nrot EQ 0 THEN print,'No rotation for this field!!'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; select objects
  ;; first of all, get rid of parents
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wp = where(pstruct.nchild EQ 0, nwp)

  IF nwp NE 0 THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; now flag and magnitude selection
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      flags = adatc.corrselect_flags

      uused=make_corrected_select(0, pstruct, adatc, maxmag, unused, $
                                input_index=wp)
      IF unused NE 0 THEN flags[uused] = $
        flags[uused] + corrselect_flag('GOODU')
      
      gused=make_corrected_select(1, pstruct, adatc, maxmag, gnused, $
                                input_index=wp)
      IF gnused NE 0 THEN flags[gused] = $
        flags[gused] + corrselect_flag('GOODG')
      rused=make_corrected_select(2, pstruct, adatc, maxmag, rnused, $
                                input_index=wp)
      IF rnused NE 0 THEN flags[rused] = $
        flags[rused] + corrselect_flag('GOODR')
                  
      iused=make_corrected_select(3, pstruct, adatc, maxmag, inused, $
                                input_index=wp)
      IF inused NE 0 THEN flags[iused] = $
        flags[iused] + corrselect_flag('GOODI')
                  
      zused=make_corrected_select(4, pstruct, adatc, maxmag, znused, $
                                input_index=wp)
      IF znused NE 0 THEN flags[zused] = $
        flags[zused] + corrselect_flag('GOODZ')
                  
      used = where(flags GT 0, nused)
              
      IF nused NE 0 THEN BEGIN

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; we made it this far so some objects will be corrected
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          adatc.corrselect_flags = flags

          FOR iband=0L, nband-1 DO BEGIN 

              clr = bands[iband]

              ;; seeing in arcseconds
              wsee = where(pstruct[used].m_rr_cc_psf[clr] GT 0.0,nsee, $
                           comp=wcompsee)

              IF nsee NE 0 THEN BEGIN 
                  wsee = used[wsee]
                  adatc[wsee].seeing[clr] = $
                    sqrt(pstruct[wsee].m_rr_cc_psf[clr]/2.)*2.35*0.4

                  ;; objects with no good psf measurement use median of
                  ;; the good ones

                  IF nsee NE nused THEN BEGIN 
                      wcompsee = used[wcompsee]
                      adatc[wcompsee].seeing[clr] = $
                        median(adatc[wsee].seeing[clr])
                  ENDIF 
              ENDIF 

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Only correct those with rotation
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              wrot = where(rotstruct.field EQ field,nrot)
              IF nrot NE 0 THEN BEGIN 
                  adatc.rotation = rotstruct[wrot].angle

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; find good moment measurements
                  ;; in this bandpass
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  flag_select, pstruct, fs, clr, good, input_index=used
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
                              
                              make_corrected_getsmear,pstruct,good,clr, corr

                              good4 = where(corr GT 0.0, ngood4)
                              IF ngood4 NE 0 THEN BEGIN 

                                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                  ;; correct shapes, copy in new fields
                                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                  good = good[good4]
                                  ngood = ngood4
                                  corr = corr[good4]
                                          
                                  adatc[good].m_e1_corr[clr] = $
                                    pstruct[good].m_e1[clr] - $
                                    corr*pstruct[good].m_e1_psf[clr]
                                          
                                  adatc[good].m_e2_corr[clr] = $
                                    pstruct[good].m_e2[clr] - $
                                    corr*pstruct[good].m_e2_psf[clr]
                                          
                                  adatc[good].m_R[clr] = corr
                                          
                              ENDIF ELSE ngood=0L ;corr > 0?

                              ;; Chris Hirata's correction scheme
                              compea4_struct, pstruct[good], clr, $
                                e1_out, e2_out, R_out, compflags
                              adatc[good].CompEA4flags[clr] = compflags

                              good5 = where(compflags EQ 0L, ngood5)
                              IF ngood5 NE 0 THEN BEGIN 
                                  tgood = good[good5]
                                  adatc[tgood].m_e1_corr_h[clr] = e1_out[good5]
                                  adatc[tgood].m_e2_corr_h[clr] = e2_out[good5]
                                  adatc[tgood].m_R_h[clr] = R_out[good5]
                              ENDIF 
                          ENDIF ELSE ngood=0L ;NAN cut

                      ENDIF ELSE ngood=0L ;; e1 cut
                      
                  ENDIF ELSE ngood=0L ;; AMOMENT flags

                  print,'Band: '+colors[clr]+' Ngood: '+ntostr(ngood)

              ENDIF ;; rotation?

          ENDFOR ;; loop over shape bandpasses

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Bayesian Galaxy Probabilities, which
          ;; depends on the seeing
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          compute_bayes_prob, pstruct, probgal, probflags
              
          adatc.objc_prob_gal = probgal
          adatc.objc_prob_flags = probflags
          wgood_prob = $
            where((probflags AND PROBFLAG_NOPROB) EQ 0, ngood_prob)
          IF ngood_prob NE 0 THEN BEGIN 
              adatc[wgood_prob].value_flags = $
                adatc[wgood_prob].value_flags + BAYES_FLAG
          ENDIF 

          ;; Only keep the good ones
          adatc = adatc[used]
          status = 0
          
      ENDIF ELSE BEGIN ;; flag objects passed cuts
          print,'No objects corrected: Flag and magnitude cuts'
          adatc = -1
      ENDELSE 

  ENDIF ELSE BEGIN  ;; no objects with nchild == 0
      print,'No objects corrected: All parents!' 
      adatc = -1
  ENDELSE 

  status = 0
  return,adatc

END 
