PRO zobjshear_sums, whist, index, ie

  COMMON zobjshear_block, tmpetan,tmperad,tmpetanerr,tmperaderr,$
    tmptansig,tmpradsig,tmptansigerr,tmpradsigerr,$
    tmpwsum,tmprsum,tmpnpsum,tmpwsum_ssh,tmpSsh,$
    rmax_act_tmp,rmin_act_tmp,$
    sumstruct, shstruct, lensumstruct, lensum, groupstruct

  sumstruct.npair[whist] = sumstruct.npair[whist] + tmpnpsum[whist]
  sumstruct.rsum[whist] = sumstruct.rsum[whist] + tmprsum[whist]
  sumstruct.etansum[whist] = sumstruct.etansum[whist] + tmpetan[whist]
  sumstruct.eradsum[whist] = sumstruct.eradsum[whist] + tmperad[whist]
  sumstruct.tansigsum[whist] = sumstruct.tansigsum[whist] + tmptansig[whist]
  sumstruct.radsigsum[whist] = sumstruct.radsigsum[whist] + tmpradsig[whist]
  
  sumstruct.etanerrsum[whist]=sumstruct.etanerrsum[whist]+tmpetanerr[whist]
  sumstruct.eraderrsum[whist]=sumstruct.eraderrsum[whist]+tmperaderr[whist]
  sumstruct.tansigerrsum[whist]=sumstruct.tansigerrsum[whist]+tmptansigerr[whist]
  sumstruct.radsigerrsum[whist]=sumstruct.radsigerrsum[whist]+tmpradsigerr[whist]
  
  sumstruct.Sshsum = sumstruct.Sshsum+tmpSsh
  
  sumstruct.wsum[whist] = sumstruct.wsum[whist] + tmpwsum[whist]
  sumstruct.wsum_ssh = sumstruct.wsum_ssh + tmpwsum_ssh
  
  ;; copy in individual lens stuff
  lensum[index].totpairs = total(tmpnpsum[whist])
  lensum[index].npair[whist] = tmpnpsum[whist]
  lensum[index].ie = ie
  lensum[index].rmax_act[whist] = rmax_act_tmp[whist]
  lensum[index].rmin_act[whist] = rmin_act_tmp[whist]
  lensum[index].rsum[whist] = tmprsum[whist]
  lensum[index].etansum[whist] = tmpetan[whist]
  lensum[index].eradsum[whist] = tmperad[whist]
  lensum[index].tansigsum[whist] = tmptansig[whist]
  lensum[index].radsigsum[whist] =  tmpradsig[whist]
  
  lensum[index].etanerrsum[whist] = tmpetanerr[whist]
  lensum[index].eraderrsum[whist] = tmperaderr[whist]
  lensum[index].tansigerrsum[whist] = tmptansigerr[whist]
  lensum[index].radsigerrsum[whist] = tmpradsigerr[whist]
  
  lensum[index].sshsum = tmpSsh
  
  lensum[index].wsum[whist] = tmpwsum[whist]
  lensum[index].wsum_ssh = tmpwsum_ssh

END 
