PRO vagc_getsmear, pstruct, indices, clr, corr

  corr = pstruct[indices].m_rr_cc_psf[clr]/pstruct[indices].m_rr_cc[clr]*(4./pstruct[indices].m_cr4_psf[clr] -1.)/(4./pstruct[indices].m_cr4[clr]-1.)

END 


PRO vagc_calculate_stuff, colstruct

  ;; Cmodel stuff

  fracpsf = colstruct.fracpsf
  devflux = colstruct.devflux
  devflux_ivar = colstruct.devflux_ivar
  expflux = colstruct.expflux
  expflux_ivar = colstruct.expflux_ivar

  cmodelflux = fracpsf*devflux + (1.0-fracpsf)*expflux
  cmodelflux_ivar = fracpsf*devflux_ivar + (1-fracpsf)*expflux_ivar

  colstruct.cmodelflux = cmodelflux
  colstruct.cmodelflux_ivar = cmodelflux_ivar

  calibflux2mag,cmodelflux,cmodel_counts,filter=5

  colstruct.cmodel_counts = cmodel_counts

  vagcdir = getenv('VAGC_REDUX')
  sdss_spec_dir = sdssidl_config('spec_dir')
  rotdir = sdss_spec_dir + 'blanton/survey_rot/'

  run = colstruct[0].run
  rerun = long(colstruct[0].rerun)
  camcol = colstruct[0].camcol

;  IF keyword_set(addphotoz) THEN BEGIN
;      IF rev_ind[field] NE rev_ind[field+1] THEN BEGIN 
;          wff = rev_ind[ rev_ind[field]:rev_ind[field+1] -1 ]
;          ;;print,'Using field: ',phzfields[wff[0]]
;          make_corrected_matchphotoz, pstruct, photoz, wff
;      ENDIF ELSE print,'No PHOTOZs for field '+fstr
;  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; make "corrected" survey coordinates
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  eq2csurvey, colstruct.ra, colstruct.dec, clambda, ceta
  colstruct.clambda = temporary(clambda)
  colstruct.ceta = temporary(ceta)
    
  ;; Add rotation

  rotname = rotdir + 'surveyrot_'+run2string(run)+'_'+ntostr(camcol)+'_'+$
    ntostr(rerun)+'.fit'
  rotstruct = mrdfits(rotname,1,/silent)


  ;; min=0 so can subscript with field
  IF n_elements(colstruct) EQ 1 THEN BEGIN 
      wf=where(rotstruct.field EQ colstruct.field, nwf)
      IF nwf NE 0 THEN BEGIN 
          colstruct.rotation = rotstruct[wf].angle
      ENDIF 
  ENDIF ELSE BEGIN 
      fh = histogram(colstruct.field, min=0, max=max(rotstruct.field),rev=frev)
      nf = n_elements(nfh)
      
      nf = n_elements(rotstruct.field)
      FOR fi=0L, nf-1 DO BEGIN 
          
          field = rotstruct[fi].field
          
          ;; Objects in this field?
          IF frev[field] NE frev[field+1] THEN BEGIN 
              
              wf = frev[ frev[field]:frev[field+1] -1 ]
              colstruct[wf].rotation = rotstruct[fi].angle
              
          ENDIF 
          
      ENDFOR 
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Correct adaptive moments (from make_corrected_files53)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  make_flag_struct, fs          ;for selecting good moments
  fs.AMOMENT_FAINT = 'N'
  fs.AMOMENT_UNWEIGHTED = 'N'
  fs.AMOMENT_SHIFT = 'N'
  fs.AMOMENT_MAXITER = 'N'

  nss = n_elements(colstruct)
  FOR clr=0,4 DO BEGIN 
      
      ;; seeing in arcseconds
      wsee = where(colstruct.m_rr_cc_psf[clr] GT 0.0,nsee)
      IF nsee NE 0 THEN BEGIN 
          colstruct[wsee].seeing[clr] = $
            sqrt(colstruct[wsee].m_rr_cc_psf[clr]/2.)*2.35*0.4
      ENDIF 

      ;; objects with no good psf measurement: get seeing
      ;; of nearest neighbor

      IF nsee NE nss THEN BEGIN 
          wseeb = lindgen(nss)
          remove, wsee, wseeb
          
          close_match_radec, colstruct[wseeb].ra, colstruct[wseeb].dec, $
            colstruct[wsee].ra, colstruct[wsee].dec, $
            mseeb, msee, 50.0/3600., 1, /silent
          IF msee[0] NE -1 THEN BEGIN 
              ;; copy in closest object's seeing
              colstruct[wseeb[mseeb]].seeing[clr] = $
                colstruct[wsee[msee]].seeing[clr]
          ENDIF 
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Only correct those with rotation
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      wrot = where(colstruct.rotation[clr] NE 9999.,nrot)
      IF nrot NE 0 THEN BEGIN 
          
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; find good moment measurements
          ;; in this bandpass
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          flag_select, colstruct, fs, clr, good
          IF good[0] NE -1 THEN BEGIN

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; still need to check not -1000.
              ;; PHOTO flags missed these!
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              good2 = $
                where(colstruct[good].m_e1[clr] NE -1000.,ngood2)
              IF ngood2 NE 0 THEN BEGIN 
                  good = good[good2]
                  ngood = n_elements(good)
                  
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; check for NAN
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  
                  good3 = where( (colstruct[good].m_e1e1err[clr] EQ $
                                  colstruct[good].m_e1e1err[clr]) AND $
                                 (colstruct[good].m_e1e2err[clr] EQ $
                                  colstruct[good].m_e1e2err[clr]) AND $
                                 (colstruct[good].m_e2e2err[clr] EQ $
                                  colstruct[good].m_e2e2err[clr]),ngood3)
                  
                  IF ngood3 NE 0 THEN BEGIN 
                      
                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                      ;; Get smear polarizability, make sure
                      ;; it is reasonable
                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                      vagc_getsmear,colstruct,good,clr, corr
                      
                      good4 = where(corr GT 0.0, ngood4)
                      IF ngood4 NE 0 THEN BEGIN 
                          
                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                          ;; correct shapes, copy in new fields
                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                          
                          good = good[good4]
                          ngood = ngood4
                          corr = corr[good4]
                          
                          colstruct[good].m_e1_corr[clr] = $
                            colstruct[good].m_e1[clr] - $
                            corr*colstruct[good].m_e1_psf[clr]
                          
                          colstruct[good].m_e2_corr[clr] = $
                            colstruct[good].m_e2[clr] - $
                            corr*colstruct[good].m_e2_psf[clr]
                          
                          colstruct[good].m_R[clr] = corr
                          
                      ENDIF ELSE ngood=0L ;corr > 0?
                      
                      ;; Chris Hirata's correction scheme
                      compea4_struct, colstruct[good], clr, $
                        e1_out, e2_out, R_out, compflags
                      colstruct[good].CompEA4flags[clr] = compflags
                      good5 = where(compflags EQ 0L, ngood5)
                      IF ngood5 NE 0 THEN BEGIN 
                          tgood = good[good5]
                          colstruct[tgood].m_e1_corr_h[clr] = e1_out[good5]
                          colstruct[tgood].m_e2_corr_h[clr] = e2_out[good5]
                          colstruct[tgood].m_R_h[clr] = R_out[good5]
                      ENDIF 
                  ENDIF ELSE ngood=0L ;NAN cut
                  
              ENDIF ELSE ngood=0L ;; e1 cut
              
          ENDIF ELSE ngood=0L ;; AMOMENT flags
          
          print,'   Band: '+!colors[clr]+' Ngood: '+ntostr(ngood)
          
      ENDIF ;; rotation?
      
  ENDFOR ;; loop over shape bandpasses
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Bayesian Galaxy Probabilities, which
  ;; depend on the seeing
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  IF keyword_set(addbayes) THEN BEGIN 
;      compute_bayes_prob, colstruct, probgal, probflags, $
;        /lognormal, /cmodel
;      colstruct.objc_prob_psf = 1. - probgal
;      colstruct.objc_prob_flags = probflags
;      wgood_prob = $
;        where((probflags AND PROBFLAG_NOPROB) EQ 0, ngood_prob)
;      IF ngood_prob NE 0 THEN BEGIN 
;          colstruct[wgood_prob].value_flags = $
;            colstruct[wgood_prob].value_flags + BAYES_FLAG
;      ENDIF 
;  ENDIF 


END 

PRO vagc_collate_spec, sp, im, kpetro, kmodel, slss, $
                       lss=lss, letter=letter, post=post, sample=sample

  ;; get the tags we want and add them to the galaxy catalog
  ;; Then output, but only by run since we have memory issues

  ;; Which sample to use: 
  ;; http://wassup.physics.nyu.edu/vagc/data/sdss/lss_letter.par

  IF n_elements(lss_vers) EQ 0 THEN BEGIN 
      lss_vers = sdssidl_config('lss_vers')
  ENDIF ELSE BEGIN 
      lss_vers = sample
  ENDELSE 
  catname = vagc_catname(lss=lss, letter=letter, post=post, sample=sample)

  sdss_spec_dir = sdssidl_config('spec_dir')
  sdss_lss_dir = sdssidl_config('lss_dir')
  sdss_vagc_dir = sdssidl_config('vagc_dir')
  sdss_vagc_vers = sdssidl_config('vagc_vers')
  sdss_shapecorr_dir = sdssidl_config('shapecorr_dir')

  outDir = sdss_spec_dir + 'blanton/gal_collated/byrun/'

  lssDir = sdss_lss_dir + lss_vers + "/" + letter + "/" + post + "/"
  vagcDir = sdss_vagc_dir + sdss_vagc_vers+'/'
  calibDir = vagcDir + 'sdss/parameters/'
  maskDir = sdss_shapecorr_dir + 'masks/'

  lssFile          = lssDir  + $
    'post_catalog.'+lss_vers +letter+post+'.fits'
  lssWindowFile    = maskDir + $
    'post_catalog.'+lss_vers +letter+post+'_inwindow.fit'

  vagcFile         = vagcDir + 'object_sdss_spectro.fits'
  imFile           = vagcDir + 'object_sdss_imaging.fits'

  kcorrPetroFile   = vagcDir + 'kcorrect/kcorrect.photoz.petro.z0.10.fits'
  kcorrModelFile   = vagcDir + 'kcorrect/kcorrect.photoz.model.z0.10.fits'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The LSS sample file if requested
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(lss) THEN BEGIN 
      IF n_elements(slss) EQ 0 THEN BEGIN 
          print
          print,'Reading LSS file: ',lssFile
          slss = mrdfits(lssFile,1)
          
          nlss = n_elements(slss)
          print,"Only Using "+ntostr(nlss)+" rows from the spec file"

          print
          print,'Reading lss Window File: ',lssWindowFile
          lssWindow = mrdfits(lssWindowFile,1)
          
      ENDIF 

      rows = slss.object_position
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The following files all match up 1-1. Indices of one are indices
  ;; of the other. Note for lss, since we set rows above, the lssWindow file
  ;; will also match up
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read the spectro information file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(sp) EQ 0 THEN BEGIN 
      spcolumns = ['mjd', 'plate', 'fiberid', 'objid', 'primtarget', $
                   'sectarget', 'class', 'subclass', 'z', 'z_err', 'rchi2', $
                   'dof', 'zwarning', 'vdisp', 'vdisp_err', $
                   'vdispchi2', 'vdispdof']
      print,'Reading vagc file: ',vagcfile
      sp = mrdfits(vagcfile, 1, columns=spcolumns, rows=rows)
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The imaging info file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(im) EQ 0 THEN BEGIN 
      imcolumns = ['RUN','RERUN','CAMCOL','FIELD','ID',$
                   'RESOLVE_STATUS']
      print,'Reading imaging file: ',imfile
      im = mrdfits(imfile, 1, columns=imcolumns, rows=rows)
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; K-corrections
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(kpetro) EQ 0 THEN BEGIN 
      print,'Reading kcorr petro file: ',kcorrPetroFile
      kpetro = mrdfits(kcorrPetroFile,1,rows=rows)
  ENDIF 

  IF n_elements(kmodel) EQ 0 THEN BEGIN 
      print,'Reading kcorr petro file: ',kcorrModelFile
      kmodel = mrdfits(kcorrModelFile, 1, rows=rows)
  ENDIF 

  ;; The output structure.
  ;; Stuff from the spectra
  defint = -9999
  deflon = -9999L
  defflt = -9999.
  spec_struct = create_struct('mjd', deflon, $
                              'plate', defint, $
                              'fiberid', defint, $
                              'objid', lonarr(5), $
                              'primtarget', deflon, $
                              'sectarget', deflon, $
                              'class', '', $
                              'subclass', '', $
                              'z', defflt, $
                              'z_err', defflt, $
                              'rchi2',defflt,$
                              'dof', deflon, $
                              'zwarning', deflon, $
                              'vdisp', defflt, $
                              'vdisp_err', defflt, $
                              'vdispchi2', defflt, $
                              'vdispdof', defflt)

  ;; Stuff from the imaging
  
  val  = -9999.
  val2 =  9999.
  val3 = -9999
  val4 = -9999d
  val5 =  0b
  val6 = -1b
  arrval  = replicate(-9999., 5)
  arrval2 = replicate( 9999., 5)
  arrval3 = replicate( 2b^6,  5)
  arrval8 = replicate( -9999., 8) ; for k-corrections
  lensstr = create_struct( 'VALUE_FLAGS',        val5, $
                           'CORRSELECT_FLAGS',   val5, $
                           'OBJC_PROB_FLAGS',    val5, $
                           'CMODELFLUX',       arrval, $ ; corrected model cnt 
                           'CMODELFLUX_IVAR',  arrval, $
                           'CMODEL_COUNTS',    arrval, $
                           'M_E1_CORR',        arrval, $ ; lensing stuff
                           'M_E2_CORR',        arrval, $
                           'M_R',              arrval, $
                           'M_E1_CORR_H',      arrval, $ ; Hirata's stuff
                           'M_E2_CORR_H',      arrval, $
                           'M_R_H',            arrval, $
                           'COMPEA4FLAGS',    arrval3, $
                           'SEEING',          arrval2, $
                           'ROTATION',        arrval2, $ ; shape rotation
                           'CLAMBDA',            val4, $ ; corrected coords.
                           'CETA',               val4, $
                           'kcorrpetro',      arrval8, $ ; kcorrections
                           'kcorrmodel',      arrval8, $ ; model kcorr
                           'kpetroflux',      arrval8, $ ; kcorrected flux
                           'kpetroflux_ivar', arrval8, $
                           'kmodelflux',      arrval8, $ ; model
                           'kmodelflux_ivar', arrval8, $
                           'abspetromag',     arrval8, $ ; petrosian
                           'absmodelmag',     arrval8  $ ; model
                         )

  ex_struct = create_struct( lensstr, spec_struct )

  ;; Mask info
  ex_struct = create_struct(ex_struct, $
                            'completeness', 0.0, $
                            'poly_id', 0L, $
                            'poly_area', 0d)

  taglist=['run', $
           'rerun', $
           'camcol', $
           'field', $
           'ID',              $
           'PARENT',          $
           'NCHILD',          $
           'OBJC_TYPE',       $
           'TYPE',            $
           'FLAGS',           $
           'FLAGS2',          $
           'OBJC_FLAGS',      $
           'OBJC_FLAGS2',     $
           'OBJC_ROWC',       $
           'OBJC_COLC',       $
           'ROWC',            $
           'COLC',            $
           'modelflux',       $
           'modelflux_ivar',  $
           'expflux',         $
           'expflux_ivar',   $
           'devflux',         $
           'devflux_ivar',   $
           'FRACPSF',         $
           'petroflux',       $
           'petroflux_ivar',  $
           'skyflux',         $
           'skyflux_ivar',    $
           'PETRORAD',        $
           'PETRORADERR',     $
           'PETROR50',        $
           'PETROR50ERR',     $
           'PETROR90',        $
           'PETROR90ERR',     $
           'psfflux',         $
           'psfflux_ivar',    $
           'resolve_status',  $
           'RA',              $
           'DEC',             $
           'PSF_FWHM',        $
           'EXTINCTION',      $
           'M_E1',            $
           'M_E2',            $
           'M_E1E1ERR',       $
           'M_E1E2ERR',       $
           'M_E2E2ERR',       $
           'M_RR_CC',         $
           'M_RR_CCERR',      $
           'M_CR4',           $
           'M_E1_PSF',        $
           'M_E2_PSF',        $
           'M_RR_CC_PSF',     $
           'M_CR4_PSF',       $
           'OBJC_PROB_PSF',   $
           'PROB_PSF',        $
           'tmass_j',         $
           'tmass_j_ivar',    $
           'tmass_h',         $
           'tmass_h_ivar',    $
           'tmass_k',         $
           'tmass_k_ivar'     $
          ]


  ;; grab the galaxies.  Will be unnessecary for the LSS samples
  nsp = n_elements(sp)
  wgal = where(strmatch(sp.class, 'GALAXY*') AND $
               sp.z GT 0.0, ngal, comp=comp, ncomp=ncomp)

  print
  print,'Found '+ntostr(ngal)+$
    ' galaxies with redshifts out of '+ntostr(nsp)+' objects'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; get the basic struct: create from one of the calibObj files plus
  ;; our extra stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  trun = im[wgal[0]].run
  trunstr = run2string(trun)
  tcalib_file = calibDir + 'calibObj-'+trunstr+'-1.fits'
  print
  print,'Creating output struct'
  tmpst = mrdfits(tcalib_file, 1, columns=taglist)

  tags_tmpst = tag_names(tmpst)
  IF n_elements(tags_tmpst) NE n_elements(taglist) THEN $
    message,'Some of the columns were not retrieved'

  tmpst = tmpst[0]
  zero_struct, tmpst
  ostruct = create_struct(tmpst, ex_struct)

  ;; output for reference later
;  mwrfits, ostruct, outDir+'ostruct.fits', /create

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; get the runs, camcols, fields, ids from the im struct
  ;; Note: runrev will subscript wgal!!
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  runs = im.run
  reruns = long(im.rerun)
  camcols = im.camcol
  fields = im.field
  ids = im.id

  ;; histogram runs
  print
  print,'Histogramming runs'
  runh = histogram(runs[wgal], rev=runrev)
  nrh = n_elements(runh)

  ;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over runs
  ;;;;;;;;;;;;;;;;;;;;;

  FOR ir=0L, nrh-1 DO BEGIN 

      ;; get the galaxies in this run
      IF runrev[ir] NE runrev[ir+1] THEN BEGIN 

          wrun = runrev[ runrev[ir]:runrev[ir+1]-1 ]
          wrun = wgal[wrun]
          nwrun = n_elements(wrun)

          rmd = rem_dup(runs[wrun])
          IF n_elements(rmd) GT 1 THEN message,'More than one run!'
          runstr = run2string(runs[wrun[0]])

          outstruct = replicate(ostruct, nwrun)

          ;; the output file
          outfile = outDir + 'run'+runstr+'-'+catname+'_collate.fits'
          print,'Run: '+runstr+' Outfile: ',outfile

          ;; histogram the camcols
          colh = histogram( camcols[wrun], rev=colrev )
          nch = n_elements(colh)
          
          FOR ic=0L, nch-1 DO BEGIN 
              ;; Get the gals in this camcol.  colrev will subscript wrun
              IF colrev[ic] NE colrev[ic+1] THEN BEGIN 

                  ;; wruncol will subscript outstruct, since outstruct
                  ;; is already the wrun subset
                  wruncol = colrev[ colrev[ic]:colrev[ic+1] -1 ]
                  wcol = wrun[wruncol]
                  nwcol = n_elements(wcol)

                  rmd = rem_dup(camcols[wcol])
                  IF n_elements(rmd) GT 1 THEN message,'More than one camcol!'
                  colstr = ntostr(camcols[wcol[0]])

                  ;; read the camcol
                  calib_file = calibDir + $
                    'calibObj-'+runstr+'-'+colstr+'.fits'


                  print,'  Camcol: '+colstr;+'  Infile: ',calib_file

                  IF n_elements(tmp) EQ 0 THEN BEGIN 
                      tmp = mrdfits_deja_vu(calib_file, /silent)
                  ENDIF ELSE BEGIN 
                      tmp = mrdfits_deja_vu(calib_file, /silent, /deja_vu)
                  ENDELSE 

                  ;; match them up
                  photo_match, runs[wcol], reruns[wcol], camcols[wcol], $
                    fields[wcol], ids[wcol], $
                    $
                    tmp.run, long(tmp.rerun), tmp.camcol, tmp.field, tmp.id, $
                    im_match, tmp_match, count=nmatch

                  IF nmatch NE nwcol THEN message,'Some did not match!'

                  colstruct = replicate(ostruct, nwcol)
                  
                  ;; copy in the stuff we need
                  im_match = wcol[im_match]
                  tmp = tmp[tmp_match]
                  sptmp = sp[im_match] ; sp and im match one-to-one

;                  print,'Copying from imaging'
                  struct_assign, tmp, colstruct, /nozero;, /verbose
;                  print,'Copying from spectra'
                  struct_assign, sptmp, colstruct, /nozero;, /verbose
;stop
                  ;; add cmodel, adaptive correction, rotation
                  vagc_calculate_stuff, colstruct

                  ;; K-corrections
                  colstruct.kcorrpetro = kpetro[im_match].kcorrect
                  colstruct.kcorrmodel = kmodel[im_match].kcorrect

                  colstruct.kpetroflux = kpetro[im_match].abmaggies
                  colstruct.kpetroflux_ivar = kpetro[im_match].abmaggies_ivar

                  colstruct.kmodelflux = kmodel[im_match].abmaggies
                  colstruct.kmodelflux_ivar = kmodel[im_match].abmaggies_ivar

                  colstruct.abspetromag = kpetro[im_match].absmag
                  colstruct.absmodelmag = kmodel[im_match].absmag

                  ;; Mask info
                  colstruct.completeness = lssWindow[im_match].completeness
                  colstruct.poly_id      = lssWindow[im_match].poly_id
                  colstruct.poly_area    = lssWindow[im_match].poly_area

;stop
                  ;; copy in these columns
;                  print,'Copying columns into main struct'
                  outstruct[wruncol] = colstruct
;stop
                  colstruct = 0
                  tmp = 0
                  sptmp = 0

              ENDIF ;; any of these camcols?

          ENDFOR ;; loop over camcols

          ;; Output the file
          mwrfits, outstruct, outfile, /create
          outstruct = 0

      ENDIF ;; any of these runs?
      
  ENDFOR ;; loop over runs

END 
