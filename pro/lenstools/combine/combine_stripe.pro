
PRO combine_stripe, runs, reruns, clr, indir=indir, outdir=outdir,$
                    nnm=nnm, lcat=lcat,$
                    magcut=magcut, err_cut=err_cut,$
                    nolink=nolink, $
                    seeing_index=seeing_index,$
                    hirata=hirata, hudson=hudson

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    COMBINE_STRIPE
;       
; PURPOSE:
;    Combine shape files from different runs (but the same stripe) into
;    one big catalog.
;
; CALLING SEQUENCE
;    combine_stripe, runs, reruns, clr, indir=indir, outdir=outdir,$
;                   nnm=nnm, /lcat,$
;                   magcut=magcut, err_cut=err_cut,$
;                   nolink=nolink, $
;                   nocheckgaps=nocheckgaps
;
; INPUTS: 
;    runs: an integer vector containing the runs.
;    reruns: an integer vector with reruns.
;    clr: the bandpass to use in integer form (g=1, r=2, i=3)
;
; OPTIONAL INPUTS:
;    indir: override default input dir
;    magcut: the r' magnitude at which to cut the catalog (default is 18)
;    err_cut: cut on error (default is 0.64)
;
; KEYWORD PARAMETERS:
;    /lcat: do lens cat
; OUTPUTS: 
;    One big catalog for whole stripe.
;   
; CALLED ROUTINES:
;    SDSSIDL_SETUP
;    MRDFITS
;    SXADDPAR
;    RMCLOSE_LAMETA
;    MWRFITS
;
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    1-June-2000  Erin Scott Sheldon UofMichigan
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: combine_stripe'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time = systime(1)
  colors=['u','g','r','i','z']

  ;; parts of the output/input filenames
  nnm = 'srcgal'
  
  
  IF keyword_set(hudson) THEN BEGIN 
      magcut = 14.5
      magcutstr = ntostr( rnd(magcut,1), 4)
      fcutstr = '_hudson'
  ENDIF ELSE if n_elements(magcut) EQ 0 THEN BEGIN 
      magcut = 18.
      magcutstr = ntostr( rnd(magcut,1), 4)
      fcutstr = ''
  ENDIF ELSE BEGIN 
      magcutstr = ntostr( rnd(magcut,1), 4)
      IF magcut EQ 18.0 THEN fcutstr = '' $
      ELSE fcutstr = '_magcut'+magcutstr
  ENDELSE 

  IF keyword_set(typecut) THEN BEGIN
      nend = '_typecut.fit'
      addstr = fcutstr+'_typecut.fit' 
  ENDIF ELSE BEGIN
      nend = '.fit'
      addstr = fcutstr+'.fit'
  ENDELSE 

  ;; Will cut on the shape error
  IF n_elements(err_cut) EQ 0 THEN err_cut = 0.4

  IF keyword_set(lcat) THEN BEGIN
      print
      print,'Creating lens catalog: r [16.0:'+magcutstr+']'
      print
      minmag = 16. 
      maxmag = magcut
  ENDIF ELSE BEGIN
      print
      print,'Creating source catalog: r ['+magcutstr+':22.0]'
      print
      minmag = magcut
      maxmag = 22.
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up runs, arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;

  nrun = n_elements(runs)
  nrerun = n_elements(reruns)
  IF nrun NE nrerun THEN BEGIN 
      print,'size(runs) and size(reruns) must be same'
      return
  ENDIF 
  run_str = ntostr(long(runs))
  rerun_str = ntostr(long(reruns))

  ;; some arrays
  flam = dblarr(nrun)
  llam = dblarr(nrun)
  cmethod = strarr(nrun)
  baye_ver = strarr(nrun)
  phtz_ver = strarr(nrun)
  rcut = strarr(nrun)
  neach = lonarr(nrun)

  strip_arr = strarr(nrun)
  stripe_arr = lonarr(nrun)

  ;; hirata now default
  IF n_elements(hirata) EQ 0 THEN hirata=1

  IF keyword_set(hirata) THEN BEGIN
      hirstr = '_h' 
      hdr_hir = 'yes'
  ENDIF ELSE BEGIN 
      hirstr = ''
      hdr_hir = 'no'
  ENDELSE 

  ;; set up the system variables
  sdssidl_setup
  setup_mystuff

  sdss_shapecorr_dir = sdssidl_config('shapecorr_dir')

  IF n_elements(indir) EQ 0 THEN $
    indir = sdss_shapecorr_dir+'corr'+run_str+'/'+rerun_str+'/combined/'
  IF n_elements(outdir) EQ 0 THEN outdir = sdss_shapecorr_dir+'combined/'
  print,'Outputting to file: ',outdir

  ;;;;;;;;;;;;;;;;;;;;;
  ;; Set up header
  ;;;;;;;;;;;;;;;;;;;;;

  run_string = ''
  rerun_string = ''
  FOR irun=0L, nrun-1 DO BEGIN 
      run_string = run_string + run_str[irun]+' '
      rerun_string = rerun_string + rerun_str[irun]+' '
  ENDFOR 

  nclr = n_elements(clr)
  clr_string = ''
  FOR iclr=0L, nclr-1 DO clr_string = clr_string + !colors[clr[iclr]]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read in runs, find the filled stripe, output
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; read in all the runs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  columns = ['run',$
             'rerun',$
             'camcol',$
             'field',$
             'id', $
             'e1',$
             'e2',$
             'e1_recorr',$
             'e2_recorr',$
             'e1e1err',$
             'e1e2err',$
             'e2e2err',$
             'upetro', $
             'gpetro', $
             'rpetro', $
             'ipetro', $
             'zpetro', $
             'grmodel', $
             'rimodel', $
             'rpetrorad', $
             'rsmear', $
             'clambda',$
             'ceta', $
             'objc_prob_psf', $
             'photoz_z',$
             'photoz_zerr', $
             'photoz_quality', $
             'photoz_abscounts', $
             'nrunave', $
             'runcombine_flags']

  print
  print,'---------------------------------------------------------'
  FOR irun = 0L, nrun-1 DO BEGIN 
      
      irunstr = ntostr(irun)
      name = indir[irun] + 'run'+ run_str[irun] +'_'+rerun_str[irun]+'_'+ nnm +'_'+clr_string+hirstr+nend
      
      print,'Reading ',name
      
      ;; read in one of the runs
      comstr = 'rstr'+irunstr+' = mrdfits(name, 1, columns=columns)'
      IF NOT execute(comstr) THEN message,'Could not read run '+irunstr

      comstr='wbad = where(rstr'+irunstr+'.clambda eq -9999., nbad)'
      IF NOT execute(comstr) THEN message,'Could not check for bad lambda'
      IF nbad NE 0 THEN message,'Found bad clambda'

      comstr = 'nrstr'+irunstr+' = n_elements(rstr'+irunstr+')'
      IF NOT execute(comstr) THEN message,'Error checking n_elements'
      comstr = 'neach['+irunstr+'] = nrstr'+irunstr
      IF NOT execute(comstr) THEN message,'Error adding to num array'
      print
      
      ;; get header
      comstr = 'hdr'+irunstr+' = headfits(name)'
      IF NOT execute(comstr) THEN message,'Error reading header'
      
      ;; get first lambda and last lambda
      comstr = 'flam['+irunstr+'] = SXPAR(hdr'+irunstr+', "FIRSTLAM")'
      IF NOT execute(comstr) THEN message,'Error getting FIRSTLAM'
      comstr = 'llam['+irunstr+'] = SXPAR(hdr'+irunstr+', "LASTLAM")'
      IF NOT execute(comstr) THEN message,'Error getting LASTLAM'
      
      ;; get the strip and stripe of this run
      comstr = 'strip_arr['+irunstr+'] = strtrim( SXPAR(hdr'+irunstr+', "STRIP") )'
      IF NOT execute(comstr) THEN message,'Error getting STRIP'
      comstr = 'stripe_arr['+irunstr+'] = SXPAR(hdr'+irunstr+', "STRIPE")'
      IF NOT execute(comstr) THEN message,'Error getting STRIPE'
      
      ;; get the cmethod
      comstr = 'cmethod['+irunstr+'] = SXPAR(hdr'+irunstr+', "CMETHOD")'
      IF NOT execute(comstr) THEN message,'Error getting CMETHOD'
      
      ;; get bayesian version
      comstr = 'baye_ver['+irunstr+'] = SXPAR(hdr'+irunstr+', "BAYE_VER")'
      IF NOT execute(comstr) THEN message,'Error getting BAYE_VER'
      
      ;; get photoz version
      comstr = 'phtz_ver['+irunstr+'] = SXPAR(hdr'+irunstr+', "PHTZ_VER")'
      IF NOT execute(comstr) THEN message,'Error getting PHTZ_VER'
      
      ;; get rcuts
      comstr = 'rcut['+irunstr+'] = SXPAR(hdr'+irunstr+', "RCUT")'
      IF NOT execute(comstr) THEN message,'Error getting RCUT'
      
      ;; cut on uncertainty, Rsmear
      comstr = 'w = where( (rstr'+irunstr+'.e1e1err LT err_cut) AND (rstr'+irunstr+'.e2e2err lt err_cut), ntmp)'
      IF NOT execute(comstr) THEN message,'Error making error cuts.'
      
      IF ntmp NE 0 THEN BEGIN 
          comstr = 'n'+irunstr+' = ntmp'
          IF NOT execute(comstr) THEN message,'Error copying num'
          comstr = 'rstr'+irunstr+' = temporary( rstr'+irunstr+'[w])'
          IF NOT execute(comstr) THEN message,'Error trimming struct rstr'+irunstr
      ENDIF ELSE BEGIN 
          print,'No good'
          return
      ENDELSE 
      
  ENDFOR 
  print,'---------------------------------------------------------'
  print
  
  ;; make sure we got stripe/strip info from headers
  ww=where(stripe_arr EQ 0,nww)
  IF nww NE 0 THEN message,'STRIPE set to zero'
  ww=where(ntostr(strip_arr) EQ '0',nww)
  IF nww NE 0 THEN message,'STRIP set to zero'
  
  ;; Check if all in same stripe
  print,'Checking if same stripe'
  wst=where(stripe_arr NE stripe_arr[0], nst)
  IF nst NE 0 THEN BEGIN 
      print
      message,'Not all same stripe!!',/inf
      message,'Stripes: '+ntostr(stripe_arr)
  ENDIF 
  
  stripestr = stripe2string(stripe_arr[0])
  print,'OK'
  
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

  ;; check if all same bayesian sg sep
  print,'Checking RCUT'
  wcm = where(rcut NE rcut[0], ncm)
  IF ncm NE 0 THEN BEGIN 
      print
      message,'Not all same RCUT!!',/inf
      print,rcut
      message,''
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; setup output file namess
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sname='stripe'+stripestr+$
    '_srcgal_'+clr_string+hirstr+addstr
  sname = outdir+sname

  link_dir  = sdss_shapecorr_dir + 'combined/figures/'
  link_name = sname
  
;  gapdir = outdir + 'figures/'
;  IF NOT exist(gapdir) THEN gapdir = outdir
;  gapfile = gapdir + 'checkgaps_' + repstr(sname, 'fit','ps')


  
  print
  print,'Output file: '+sname
;  print,'Gapfile: '+gapfile
  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find all the strips
  ;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Checking the strips'
  ustrip = strip_arr[ rem_dup(strip_arr) ]
  nstrip = n_elements(ustrip)
  IF nstrip NE 2 THEN BEGIN 
      print
      print,'There is only one strip here'
      return
  ENDIF 
  print,'OK: ',ustrip

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; how many to each strip?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;

  w1 = where(strip_arr EQ ustrip[0], nstrip1)
  w2 = where(strip_arr EQ ustrip[1], nstrip2)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; are there gaps? Probably can
  ;; loosen this criteria done in lambda
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  print
;  print,'Checking for gaps'
;  begplot,name=gapfile,/color
;  names = ntostr(stripe_arr[0])+'-'+ustrip[0]+'-'+run_str[w1]
;  checkgaps, flam[w1], llam[w1], isgap1, names=names

;  names = ntostr(stripe_arr[0])+'-'+ustrip[1]+'-'+run_str[w2]
;  checkgaps, flam[w2], llam[w2], isgap2, names=names
;  endplot 
;  IF (isgap1 OR isgap2) THEN BEGIN 
;      IF NOT keyword_set(nocheckgaps) THEN message,'There are gaps' $
;      ELSE message,'There are gaps: proceeding anyway due to /nocheckgaps',/inf
;  ENDIF 

;  print,'OK'
  print
  s = sort(strip_arr)

  sformat = '(%"%6s %5s %5s %5s %13s %13s %8s")'
  format  = '(%"%6d %5s %5d %5d %13.8f %13.8f %8d")'
  print,format=sformat,$
    "Stripe","Strip","Run","Rerun","FirstLam","LastLam","Number"
  colprint,format=format,$
    stripe_arr[s],strip_arr[s],runs[s],reruns[s],flam[s],llam[s],neach[s]

  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Combine same-strip runs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; look at the first strip
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;

  comstr = 'r1 = temporary(rstr'+ntostr(w1[0])+')'
  IF NOT execute(comstr) THEN message,'Error copying r1'
  IF nstrip1 GT 1 THEN BEGIN    ; more than one run in this strip
                                ; concatenate them
      FOR j=1, nstrip1-1 DO BEGIN 
          irunstr=ntostr(w1[j])
          print,'Concatenating in strip '+ustrip[0]
          print
          comstr='concat_structs, temporary(r1), temporary(rstr'+irunstr+'), temp'

          IF NOT execute(comstr) THEN message,'Error concatenating strip1'
          r1 = temporary(temp)
      ENDFOR 
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; look at the second strip
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;

  comstr = 'r2 = temporary(rstr'+ntostr(w2[0])+')'
  IF NOT execute(comstr) THEN message,'Error copying r2'
  IF nstrip2 GT 1 THEN BEGIN    ; more than one run in this strip
                                ; concatenate them
      FOR j=1, nstrip2-1 DO BEGIN 
          irunstr=ntostr(w2[j])
          print,'Concatenating in strip '+ustrip[1]
          print
          comstr='concat_structs, temporary(r2), temporary(rstr'+irunstr+'), temp'
          IF NOT execute(comstr) THEN message,'Error concatenating strip2'
          r2 = temporary(temp)
      ENDFOR 
  ENDIF 


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; make sure not bad clambda, ceta
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  n1 = n_elements(r1)
  n2 = n_elements(r2)

  w1 = where(r1.clambda NE 0.0d,nw1)
  w2 = where(r2.clambda NE 0.0d,nw2)
  IF (nw1 EQ 0) OR (nw2 EQ 0) THEN message,'All ra/dec are zero!'
  IF nw1 LT n1 THEN BEGIN
      message,'Bad ra/dec in strip '+ustrip[0],/inf
      r1 = temporary(r1[w1])
  ENDIF 
  IF nw2 LT n2 THEN BEGIN
      message,'Bad ra/dec in strip '+ustrip[1],/inf
      r2 = temporary(r2[w1])
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; define overlap region
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  firstlam = max( [min(r1.clambda), min(r2.clambda)] )
  lastlam = min( [max(r1.clambda), max(r2.clambda)] )

  print,'First Lambda: ',firstlam
  print,'Last Lambda: ',lastlam

  print,'Making clambda cuts: ['+ntostr(firstlam)+', '+ntostr(lastlam)+']'

  w1 = where(r1.clambda GE firstlam AND r1.clambda LE lastlam, n1)
  w2 = where(r2.clambda GE firstlam AND r2.clambda LE lastlam, n2)

  IF (n1 EQ 0) OR (n2 EQ 0) THEN BEGIN 
      print,'No good 2'
      help,w1,w2
      return
  ENDIF ELSE BEGIN 
      r1 = temporary( r1[w1] )
      r2 = temporary( r2[w2] )
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; set up new struct
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ntot = n1+n2
  concat_structs, temporary(r1), temporary(r2), struct

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make primary bound cuts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  primary_bound, stripe_arr[0], bound

  wprim = where(struct.ceta GE bound.etamin AND $
                struct.ceta LT bound.etamax, nprim)
  IF nprim EQ 0 THEN message,'None in primary region!'

  nremove = ntot-nprim
  print,'Removed '+ntostr(nremove)+' from non-primary region'
  struct = temporary( struct[wprim] )
  ntot = nprim

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; magnitude cuts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wmag = where(struct.rpetro GE minmag AND struct.rpetro LE maxmag, nmag)
  IF nmag EQ 0 THEN message,'None in primary region!'

  nremove = ntot-nmag
  print,'Removed '+ntostr(nremove)+' in mag cuts'
  struct = temporary( struct[wmag] )
  ntot = nmag

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Get rid of duplicates. 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  s_ind = sort(struct.clambda)
  struct = temporary(struct[s_ind])
  
;  rmclose_lameta_seeing, struct, keep, seeing_index=seeing_index

  rmclose_lameta_meane, struct, keep, out, $
                        new_e1, new_e2, $
                        new_e1_recorr, new_e2_recorr, $
                        new_e1e1err, new_e1e2err, new_e2e2err,$
                        new_rsmear, $
                        photoz_z, photoz_zerr, photoz_quality, $
                        new_cflags, new_nave

  struct = temporary(struct[keep])

  struct.e1 = new_e1
  struct.e2 = new_e2
  struct.e1_recorr = new_e1_recorr
  struct.e2_recorr = new_e2_recorr
  struct.e1e1err = new_e1e1err
  struct.e1e2err = new_e1e2err
  struct.e2e2err = new_e2e2err

  struct.rsmear = new_rsmear

  struct.photoz_z = photoz_z
  struct.photoz_zerr = photoz_zerr
  struct.photoz_quality = photoz_quality

  struct.nrunave = new_nave
  struct.runcombine_flags = new_cflags

  nremove = ntot - n_elements(keep)
  print,ntostr(nremove),' galaxies matched in overlap region'
  ntot = n_elements(keep)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;  Calculate htm leafid
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Looking up HTM leafids'
  depth = 10
  csurvey2eq, struct.clambda, struct.ceta, ra, dec
  htmLookupRadec, ra, dec, depth, leafId

  add_tag, temporary(struct), 'leafid', 0L, struct
  struct.leafId = leafId


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; print out some info
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'First clambda = ',firstlam
  print,'Last clambda = ',lastlam
  print
  print,'Min ceta = ',min(struct.ceta, max=maxceta)
  print,'Max ceta = ',maxceta
  print


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output the file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  print,'Number in stripe: ',ntostr(ntot)
  print
  print,'Writing combined file: ',sname

  strip_string=''
  FOR it=0L,nrun-1 DO strip_string=strip_string+ntostr(strip_arr[it])+' '

  ;; info in both headers
  outhdr = ['END']
  SXADDPAR, outhdr, 'RUN', run_string, ' Imaging run number'
  SXADDPAR, outhdr, 'RERUN', rerun_string, ' Rerun number'
  SXADDPAR, outhdr, 'FILTERS', clr_string, ' Bandpasses used'
  SXADDPAR, outhdr, 'STRIPE', stripe_arr[0], ' Stripe number'
  SXADDPAR, outhdr, 'STRIP', strip_string, ' Strips used N - north, S - south'

  SXADDPAR, outhdr, 'FIRSTLAM', firstlam, ' First lambda in this data'
  SXADDPAR, outhdr, 'LASTLAM', lastlam, ' Last lambda in this data'
  SXADDPAR, outhdr, 'CMETHOD', cmethod[0], ' Method used to calculate covariance matrix'
  SXADDPAR, outhdr, 'BAYE_VER', baye_ver[0], ' Version of bayesian s/g separation code used.'
  SXADDPAR, outhdr, 'PHTZ_VER', phtz_ver[0], ' Version of template photoz code used.'
  SXADDPAR, outhdr, 'ERR_CUT', err_cut, ' Maximum error cut for ellipticity measurement'
  SXADDPAR, outhdr, 'HIRATA', hdr_hir, ' Did we use the Hirata method?'
  SXADDPAR, outhdr, 'HTMDEPTH', depth, ' Depth of HTM leaf ids'

  mwrfits2, struct, sname, outhdr, /create, /destroy

  ptime,systime(1)-time

  return
END 
