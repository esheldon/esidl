PRO zobjshear_define, nbin

  COMMON zobjshear_block, tmpetan,tmperad,tmpetanerr,tmperaderr,$
    tmptansig,tmpradsig,tmptansigerr,tmpradsigerr,$
    tmpwsum,tmprsum,tmpnpsum,tmpwsum_ssh,tmpSsh,$
    rmax_act_tmp,rmin_act_tmp,$
    sumstruct, shstruct, lensumstruct, lensum, groupstruct

  arrval = fltarr(nbin)

  shstruct = zshstruct(arrval)
  sumstruct = zsumstruct(arrval)
  
  lensumstruct = zlensumstruct(arrval)

  ;; temporary
  tmpetan = arrval
  tmperad = arrval
  tmpetanerr = arrval
  tmperaderr = arrval

  tmptansig = arrval
  tmpradsig = arrval
  tmptansigerr = arrval
  tmpradsigerr = arrval

  tmpwsum = arrval
  tmprsum = arrval
  tmpnpsum = arrval
  tmpwsum_ssh = 0.
  tmpSsh = 0.

  rmax_act_tmp = arrval
  rmin_act_tmp = replicate(1.e7, nbin)

END 
