PRO read_etarange, stripe, etarange, spectra=spectra, indir=indir

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: read_etarange, stripe, etarange, spectra=spectra, indir=indir'
      return
  ENDIF 

  IF n_elements(indir) EQ 0 THEN $
    indir = sdssidl_config('SHAPECORR_DIR')+'masks/'

  sdssidl_setup,/silent
  IF NOT keyword_set(spectra) THEN BEGIN 
      etarangefile=indir+'etarange-stripe'+$
        ntostr(long(stripe))+'.fit'
      print
      print,'Reading etarange file: ',etarangefile
      etarange=mrdfits(etarangefile,1,/silent)
  ENDIF ELSE BEGIN 
      etarangefile=indir+'etarange-spectra-stripe'+$
        ntostr(long(stripe))+'.fit'
      print
      print,'Reading etarange_spectra file: ',etarangefile
      etarange=mrdfits(etarangefile,1,/silent)
  ENDELSE 

END 
