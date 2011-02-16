FUNCTION shhdr, shstruct

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: hdr = sumhdr(shstruct)'
      return,''
  ENDIF 

  nbin = n_elements(shstruct.shear)

  shhdr = ['END']
  SXADDPAR, shhdr, 'nlenses', shstruct.nlenses
  SXADDPAR, shhdr, 'totpairs', shstruct.totpairs
  SXADDPAR, shhdr, 'binwidth', shstruct.binsize
  SXADDPAR, shhdr, 'rmin', shstruct.rmin
  SXADDPAR, shhdr, 'rmax', shstruct.rmax
  SXADDPAR, shhdr, 'rmax_act', shstruct.rmax_act[nbin-1]
  SXADDPAR, shhdr, 'ssh', shstruct.Ssh

  return, shhdr

END 
