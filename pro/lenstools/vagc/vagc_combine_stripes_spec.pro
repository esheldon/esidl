PRO vagc_combine_stripes_spec, stripes, indir=indir, outdir=outdir, lrg=lrg, lss=lss, letter=letter, post=post, sample=sample

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: combine_stripes, stripes, indir=indir, outdir=outdir, /lrg, /lss, letter=letter, post=post, sample=sample'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time = systime(1)

  ;; set up the system variables
  sdssidl_setup
  setup_mystuff

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output file names
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  catname = vagc_catname(lss=lss, letter=letter, post=post, sample=sample)

  IF keyword_set(lrg) THEN BEGIN 
      nnmout = 'lrg_'+catname
  ENDIF ELSE nnmout = catname
  
  nend = '.fit'
  addstr = '.fit'

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;

  nstripe = n_elements(stripes)
  stripe_str = ntostr(long(stripes))

  sdss_shapecorr_dir = sdssidl_config('shapecorr_dir')
  IF n_elements(indir) EQ 0 THEN $
    indir = sdss_shapecorr_dir+'combined/'
  IF n_elements(outdir) EQ 0 THEN outdir = vagc_lensinput_dir(lss=lss)

  hdrstripe_string = ''
  FOR ist=0L, nstripe-1 DO BEGIN 
      hdrstripe_string = hdrstripe_string + stripe_str[ist]+' '
  ENDFOR 

  stripe_string = stripearr2string(stripes)

  files = indir + 'stripe'+stripe2string(stripes)+'_'+nnmout+addstr

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; setup output file namess
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sname=outdir + 'stripe'+stripe_string+'_'+nnmout+addstr

  print
  print,'Output file: '+sname

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; read in all the stripes
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  mrdfits_multi, files, struct

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Get rid of duplicates. 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ntot = n_elements(struct)

  s_ind = sort(struct.clambda)
  struct = temporary(struct[s_ind])
  
  rmclose_lameta, struct.clambda, struct.ceta, keep
  struct = temporary(struct[keep])

  nkeep = n_elements(keep)
  nremove = ntot - nkeep
  print,ntostr(nremove),' galaxies matched in overlap region'
  ntot = nkeep

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output the file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  print,'Number in stripe: ',ntostr(ntot)
  print
  print,'Writing combined file: ',sname

  ;; info in both headers
  outhdr = ['END']
  SXADDPAR, outhdr, 'STRIPE', hdrstripe_string, ' Imaging run number'

  mwrfits2, struct, sname, outhdr, /create, /destroy

  ptime,systime(1)-time

  return
END 
