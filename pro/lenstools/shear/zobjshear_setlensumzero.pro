PRO zobjshear_setlensumzero, lensum, index

  lensum[index].npair[*] = 0.0
  lensum[index].totpairs = 0.0

  lensum[index].rmax_act[*] = 0.0
  lensum[index].rmin_act[*] = 0.0
  lensum[index].rsum[*] = 0.0

  lensum[index].sigma[*] = 0.0
  lensum[index].sigmaerr[*] = 0.0

  lensum[index].orthosig[*] = 0.0
  lensum[index].orthosigerr[*] = 0.0

  lensum[index].sigerrsum[*] = 0.0
  lensum[index].orthosigerrsum[*] = 0.0

  lensum[index].sshsum = 0.0
  lensum[index].wsum[*] = 0.0
  lensum[index].wsum_ssh = 0.0

  lensum[index].weight = 0.0


END 
