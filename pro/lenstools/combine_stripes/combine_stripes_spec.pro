
PRO combine_stripes_spec, stripes, indir=indir, outdir=outdir, lrg=lrg

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: combine_stripes, stripes, indir=indir, outdir=outdir, /lrg'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time = systime(1)
  colors=['u','g','r','i','z']

  ;; parts of the output/input filenames
  nnm = 'spec'
  
  nend = '.fit'
  addstr = '.fit'

  IF keyword_set(lrg) THEN nnmout = 'lrg_spec' ELSE nnmout = 'spec'

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;

  nstripe = n_elements(stripes)
  stripe_str = ntostr(long(stripes))

  neach = lonarr(nstripe)

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

  stripe_string = stripearr2string(stripes)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; setup output file namess
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sname=outdir + 'stripe'+stripe_string+'_'+nnmout+addstr

  print
  print,'Output file: '+sname

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read in runs, find the filled stripe, output
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; read in all the runs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'---------------------------------------------------------'
  FOR istripe = 0L, nstripe-1 DO BEGIN 
      
      get_spectra_lcat, stripes[istripe], tmp, /spec, lrg=lrg

      wbad = where(tmp.clambda EQ -9999., nbad)
      IF nbad NE 0 THEN message,'Found bad clambda'

      neach[istripe] = n_elements(tmp)

      IF istripe EQ 0 THEN BEGIN 
          struct = temporary(tmp)
      ENDIF ELSE BEGIN 
          concat_structs, temporary(tmp), temporary(struct), tmp2
          struct = temporary(tmp2)
      ENDELSE 
      
  ENDFOR 
  print,'---------------------------------------------------------'
  print
    
  print,'OK'
  print
  colprint,'Stripe    Number'
  colprint,stripe_str,neach
  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Get rid of duplicates. 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ntot = n_elements(struct)

  s_ind = sort(struct.clambda)
  struct = temporary(struct[s_ind])
  
  rmclose_lameta_matched, struct, keep
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
