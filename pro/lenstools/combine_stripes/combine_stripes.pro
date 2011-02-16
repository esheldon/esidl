
PRO combine_stripes, stripes, clr, type, indir=indir, outdir=outdir, hirata=hirata, hudson=hudson

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: combine_stripes, stripes, clr, type, indir=indir, outdir=outdir, /hirata, /hudson'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; this type means what source tags go with these lens catalogs
  CASE type OF
      'specgal': BEGIN 
          typestr = '_nzCuts_specGal'
          columns = ['leafid','clambda','ceta','e1_recorr','e2_recorr',$
                     'e1e1err','e1e2err','e2e2err','photoz_z','photoz_zerr','rsmear']
          nzCuts = 1
      END 
      'maxbcg': BEGIN 
          typestr = '_nzCuts_MaxBCG'
          columns = ['leafid',$
                     'clambda','ceta',$
                     'e1_recorr','e2_recorr',$
                     'e1e1err','e1e2err','e2e2err',$
                     'photoz_z','photoz_zerr',$
                     'rsmear',$
                     'rpetro']
          nzCuts = 1
      END 
      ELSE: message,'Unknown type: '+ntostr(type)
  ENDCASE 

  time = systime(1)
  colors=['u','g','r','i','z']

  ;; parts of the output/input filenames
  nnm = 'srcgal'
  
  nend = '.fit'
  addstr = '.fit'
  addstr2 = '.st'

  IF n_elements(hirata) EQ 0 THEN hirata=1

  IF keyword_set(hirata) THEN BEGIN
      hirstr = '_h' 
      hdr_hir = 'yes'
  ENDIF ELSE BEGIN 
      hirstr = ''
      hdr_hir = 'no'
  ENDELSE 

  IF keyword_set(hudson) THEN addstr='_hudson'+addstr


  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;

  nstripe = n_elements(stripes)
  stripe_str = ntostr(long(stripes))

  primary_bound_multi, stripes, bound

  ;; some arrays
  flam = dblarr(nstripe)
  llam = dblarr(nstripe)
  cmethod = strarr(nstripe)
  baye_ver = strarr(nstripe)
  phtz_ver = strarr(nstripe)
  neach = lonarr(nstripe)

  htmDepth = lonarr(nstripe)

  ;; set up the system variables
  sdssidl_setup
  setup_mystuff

  sdss_shapecorr_dir = sdssidl_config('shapecorr_dir')
  IF n_elements(indir) EQ 0 THEN $
    indir = sdss_shapecorr_dir+'combined/'
  IF n_elements(outdir) EQ 0 THEN outdir = sdss_shapecorr_dir+'combined/'

  hdrstripe_string = ''
  FOR ist=0L, nstripe-1 DO BEGIN 
      hdrstripe_string = hdrstripe_string + stripe_str[ist]+' '
  ENDFOR 

  clr_string = clrarr2string(clr)
  stripe_string = stripearr2string(stripes)


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; setup output file namess
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sname=outdir + 'stripe'+stripe_string+$
    '_srcgal_'+clr_string+hirstr+typestr+addstr
  sname2=outdir + 'stripe'+stripe_string+$
    '_srcgal_'+clr_string+hirstr+typestr+addstr2
  print
  print,'Output file: '+sname
  print,'Output file2: '+sname2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read in runs, find the filled stripe, output
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; read in all the runs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'---------------------------------------------------------'

  ;; first loop through and find numbef of objects in each 
  ;; file

  ntot = 0L
  neach = lonarr(nstripe)
  FOR istripe=0L, nstripe-1 DO BEGIN 

      tfile = SDSS_SHAPECORR_DIR + $
        'combined/stripe'+stripe2string(stripes[istripe])+$
        '_srcgal_'+clr_string+hirstr+addstr

      hdr = headfits(tfile, ext=1)
      neach[istripe] = sxpar(hdr,'naxis2')
      ntot = ntot + neach[istripe]

      flam[istripe] = SXPAR(hdr, "FIRSTLAM")
      llam[istripe] = SXPAR(hdr, "LASTLAM")
      cmethod[istripe] = SXPAR(hdr, "CMETHOD")
      baye_ver[istripe] = SXPAR(hdr, "BAYE_VER")
      phtz_ver[istripe] = SXPAR(hdr, "PHTZ_VER")
      htmDepth[istripe] = SXPAR(hdr, "HTMDEPTH")

  ENDFOR 

  ;; Check if all same cmethod
  print,'Checking CMETHOD'
  wcm = where(cmethod NE cmethod[0], ncm)
  IF ncm NE 0 THEN BEGIN 
      print
      message,'Not all same CMETHOD!!',/inf
      print,cmethod
      message,''
  ENDIF 
  
  ;; check if all same bayesian sg sep
  print,'Checking BAYE_VER'
  wcm = where(baye_ver NE baye_ver[0], ncm)
  IF ncm NE 0 THEN BEGIN 
      print
      message,'Not all same BAYE_VER!!',/inf
      print,baye_ver
      message,''
  ENDIF 

  ;; check if all same bayesian sg sep
;  print,'Checking PHTZ_VER'
;  wcm = where(phtz_ver NE phtz_ver[0], ncm)
;  IF ncm NE 0 THEN BEGIN 
;      print
;      message,'Not all same PHTZ_VER!!',/inf
;      print,phtz_ver
;      message,''
;  ENDIF 

  ;; check if all same HTM depth
  print,'Checking HTMDEPTH'
  whtm = where(htmDepth NE htmDepth[0], nhtm)
  IF nhtm NE 0 THEN BEGIN 
      print
      message,'Not all same HTMDEPTH!!',/inf
      print,htmDepth
      message,''
  ENDIF 

  print,'OK'

  get_scat, stripes, clr, struct, $
    hirata=hirata, hudson=hudson, $
    columns=columns, nzCuts=nzCuts

  print
  colprint,'Stripe    First_lam        Last_lam       Number'
  colprint,stripe_str,flam,llam,neach
  print

  print
  print,'Sorting by HTM leafid'
  s_ind = sort(struct.leafid)
  struct = struct[s_ind]

  print
  print,'Histogramming leafids'
  minid = min(struct.leafid, max=maxid)
  leafHist = histogram(struct.leafid, min=minid, max=maxid, rev=rev)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; print out some info
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Min ceta = ',min(struct.ceta, max=maxceta)
  print,'Max ceta = ',maxceta

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output the file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  print,'Number in stripes: ',ntostr(ntot)

  print
  print,'Writing combined st file: ',sname2
  write_idlstruct, struct, sname2

  print
  print,'Writing combined file: ',sname

  ;; info in both headers
  outhdr = ['END']
  SXADDPAR, outhdr, 'STRIPE', hdrstripe_string, ' Stripe numbers'
  SXADDPAR, outhdr, 'FILTERS', clr_string, ' Bandpasses used'

  SXADDPAR, outhdr, 'CMETHOD', cmethod[0], ' Method used to calculate covariance matrix'
  SXADDPAR, outhdr, 'BAYE_VER', baye_ver[0], ' Version of bayesian s/g separation code used.'
  SXADDPAR, outhdr, 'PHTZ_VER', phtz_ver[0], ' Version of template photoz code used.'
  SXADDPAR, outhdr, 'HIRATA', hdr_hir, ' Did we use the Hirata method?'
  SXADDPAR, outhdr, 'HTMDEPTH', htmDepth[0], ' Depth of HTM leaf ids'

  mwrfits2, struct, sname, outhdr, /create, /destroy

  print
  print,'Appending reverse indices to file: ',sname
  mwrfits, rev, sname

  ptime,systime(1)-time

  return
END 
