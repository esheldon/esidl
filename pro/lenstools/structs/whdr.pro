FUNCTION whdr, wstruct

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: hdr = whdr(wstruct)'
      return,''
  ENDIF 

  nbin = n_elements(wstruct.rmax_act)

  hdr = ['END']
  SXADDPAR, hdr, 'nlenses', wstruct.nlenses
  SXADDPAR, hdr, 'totpairs', wstruct.totpairs
  SXADDPAR, hdr, 'binwidth', wstruct.binsize
  SXADDPAR, hdr, 'rminkpc', wstruct.rmin
  SXADDPAR, hdr, 'rmaxkpc', wstruct.rmax
  SXADDPAR, hdr, 'rmax_act', wstruct.rmax_act[nbin-1]
  SXADDPAR, hdr, 'h', wstruct.h

  return, hdr

END 
