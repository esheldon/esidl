FUNCTION zshhdr, shstruct

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: hdr = zsumhdr(shstruct)'
      return,''
  ENDIF 

  nbin = n_elements(shstruct.sigma)

  shhdr = ['END']
  sxAddPar, shhdr, 'nlenses',      shstruct.nlenses
  sxAddPar, shhdr, 'totpairs',     shstruct.totpairs
  sxAddPar, shhdr, 'binwidth',     shstruct.binsize
  sxAddPar, shhdr, 'rminkpc',      shstruct.rmin
  sxAddPar, shhdr, 'rmaxkpc',      shstruct.rmax
  sxAddPar, shhdr, 'rmax_act',     shstruct.rmax_act[nbin-1]
  sxAddPar, shhdr, 'ssh',          shstruct.Ssh
  sxAddPar, shhdr, 'zmean',        shstruct.zmean
  sxAddPar, shhdr, 'meanscritinv', shstruct.meanscritinv
  sxAddPar, shhdr, 'h',            shstruct.h

  IF tag_exist(shstruct, 'compcut')  THEN sxAddPar, shhdr, 'compcut',  shstruct.compcut
  IF tag_exist(shstruct, 'depth')    THEN sxAddPar, shhdr, 'depth',    shstruct.depth
  IF tag_exist(shstruct, 'comoving') THEN sxAddPar, shhdr, 'comoving', shstruct.comoving
  IF tag_exist(shstruct, 'logbin')   THEN sxAddPar, shhdr, 'logbin',   shstruct.logbin

  return, shhdr

END 
