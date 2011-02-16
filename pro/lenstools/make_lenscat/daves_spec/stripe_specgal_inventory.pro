PRO stripe_specgal_inventory

  sdss_shapecorr_dir = sdssidl_config('shapecorr_dir')
  indir = SDSS_SHAPECORR_DIR+'combined/'
  outdir = SDSS_SHAPECORR_DIR+'spec_index/'

  lrgstr_arr = ['', 'lrg_']
  lrgstr2_arr = ['', 'LRG']

  cd,indir

  ;; do main and lrg samples
  FOR j=0,1 DO BEGIN 

      lrgstr = lrgstr_arr[j]
      lrgstr2 = lrgstr2_arr[j]

      outfile = outdir+'stripe_'+lrgstr+'specgal_inventory.dat'
  
      files = findfile('stripe[0-9][0-9]_'+lrgstr+'spec.fit')

      nf = n_elements(files)
      stripes = bytarr(nf)
      ngals = lonarr(nf)

      print,'Printing to file: ',outfile
      openw, lun, outfile, /get_lun
      print, '  Stripe   Ngal'
      FOR i=0L, nf-1 DO BEGIN  
          
          hdr = headfits(files[i], ext=1)
          stripes[i] = fix(sxpar(hdr, 'STRIPE'))
          ngals[i] = long(sxpar(hdr, 'NAXIS2'))
          
          print, stripes[i], ngals[i]
          
      ENDFOR 
      
      s=sort(stripes)
      stripes = stripes[s]
      ngals = ngals[s]
      
      printf, lun, 'Ntotal: '+ntostr(long(total(ngals)))+' '+lrgstr2
      printf, lun, '  Stripe   Ngal'
      colprint, stripes, ngals, lun=lun
      
      free_lun, lun
  ENDFOR 

END 
