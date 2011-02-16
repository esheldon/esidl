FUNCTION zsumhdr, sumstruct

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: hdr = zsumhdr(sumstruct)'
      return,''
  ENDIF 

  nbin = n_elements(sumstruct.rsum)

  sumhdr = ['END']
  sxAddPar, sumhdr, 'nlenses',     sumstruct.nlenses
  sxAddPar, sumhdr, 'totpairs',    sumstruct.totpairs
  sxAddPar, sumhdr, 'binwidth',    sumstruct.binsize
  sxAddPar, sumhdr, 'rminkpc',     sumstruct.rmin
  sxAddPar, sumhdr, 'rmaxkpc',     sumstruct.rmax
  sxAddPar, sumhdr, 'rmax_act',    sumstruct.rmax_act[nbin-1]
  sxAddPar, sumhdr, 'sshsum',      sumstruct.Sshsum
  sxAddPar, sumhdr, 'zsum',        sumstruct.zsum
  sxAddPar, sumhdr, 'scritinvsum', sumstruct.scritinvsum
  sxAddPar, sumhdr, 'h',           sumstruct.h

  IF tag_exist(sumstruct, 'compcut')  THEN sxAddPar, sumhdr, 'compcut',  sumstruct.compcut
  IF tag_exist(sumstruct, 'depth')    THEN sxAddPar, sumhdr, 'depth',    sumstruct.depth
  IF tag_exist(sumstruct, 'comoving') THEN sxAddPar, sumhdr, 'comoving', sumstruct.comoving
  IF tag_exist(sumstruct, 'logbin')   THEN sxAddPar, sumhdr, 'logbin',   sumstruct.logbin

  return, sumhdr

END 
