FUNCTION sumhdr, sumstruct

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: hdr = sumhdr(sumstruct)'
      return,''
  ENDIF 

  nbin = n_elements(sumstruct.rsum)

  sumhdr = ['END']
  SXADDPAR, sumhdr, 'nlenses', sumstruct.nlenses
  SXADDPAR, sumhdr, 'totpairs', sumstruct.totpairs
  SXADDPAR, sumhdr, 'binwidth', sumstruct.binsize
  SXADDPAR, sumhdr, 'rmin', sumstruct.rmin
  SXADDPAR, sumhdr, 'rmax', sumstruct.rmax
  SXADDPAR, sumhdr, 'rmax_act', sumstruct.rmax_act[nbin-1]
  SXADDPAR, sumhdr, 'sshsum', sumstruct.Sshsum

  return, sumhdr

END 
