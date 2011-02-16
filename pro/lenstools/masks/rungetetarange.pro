PRO rungetetarange, stripes, outdir=outdir, indir=indir

  ;; only need to run this on stripes where spectroscopy is
  ;; not rectangle in lambda,eta: 10,82,76,86

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: rungetetarange, stripes, outdir=outdir, indir=indir'
      return
  ENDIF 

  sdssidl_setup,/silent
  setup_mystuff

  IF n_elements(outdir) EQ 0 THEN $
    outdir = sdssidl_config('SHAPECORR_DIR')+'masks/'

  nstripe = n_elements(stripes)

  FOR i=0L, nstripe-1 DO BEGIN 

      stripe=stripes[i]

      ;; spectra etarange
;      get_spectra_lcat, stripe, lcat, indir=indir, /spec
;      IF n_elements(lcat) NE 0 THEN BEGIN 
;          outfile = outdir+'etarange-spectra-stripe'+ntostr(long(stripe))+'.fit'
;          print,outfile
;          getetarange, lcat, binstruct
;          mwrfits, binstruct, outfile, /create
;      ENDIF
;      delvarx,lcat,binstruct

      ;; everything etarange

      get_scat, stripe, [1,2,3], scat, indir=indir
      get_spectra_lcat, stripe, lcat, /spec
      IF n_elements(scat) NE 0 THEN BEGIN 
          outfile = outdir+'etarange-stripe'+ntostr(long(stripe))+'.fit'
          print,outfile
          getetarange, scat, binstruct, oplotstruct=lcat
          mwrfits, binstruct, outfile, /create
      ENDIF
      delvarx,scat,binstruct


  ENDFOR 


END 
