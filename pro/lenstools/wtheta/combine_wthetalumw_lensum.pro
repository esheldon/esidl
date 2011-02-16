PRO combine_wthetalumw_lensum, lensum, binsize, rminkpc, rmaxkpc, hval, sumstruct, wstruct

  ;; this is for weighted npsum, with both lum and sigcrit weighting

  IF n_params() LT 5 THEN BEGIN 
      print,'-Syntax: combine_wtheta_lensum2, lensum, binsize, rminkpc, rmaxkpc, hval, sumstruct, wstruct'
      return
  ENDIF 

  IF tag_exist(lensum, 'lerrsum1') THEN doerrsum=1 ELSE doerrsum=0

  beta = lensum[0].beta
  nbeta = n_elements(beta)
  nlens = n_elements(lensum)
  nbin = n_elements(lensum[0].rsum)
  arrval = fltarr(nbin)
  
  sumstruct = wthetaLumwSumstruct(arrval,beta)

  sumstruct.nlenses = nlens
  sumstruct.binsize = binsize
  sumstruct.rmin = rminkpc
  sumstruct.rmax = rmaxkpc
  sumstruct.h    = hval

  sumstruct.totpairs = total( lensum.totpairs )

  FOR bin=0L, nbin-1 DO BEGIN 

      sumstruct.rmax_act[bin] = max(lensum.rmax_act[bin])
      sumstruct.rmin_act[bin] = min(lensum.rmin_act[bin])
      sumstruct.rsum[bin] = total(lensum.rsum[bin])
      sumstruct.lsum[bin] = total(lensum.lsum[bin])

      IF doerrsum THEN BEGIN 
          sumstruct.lerrsum1[bin] = total(lensum.lerrsum1[bin])
          sumstruct.lerrsum2[bin] = total(lensum.lerrsum2[bin])
          sumstruct.lerrsum3[bin] = total(lensum.lerrsum3[bin])
      ENDIF 

      sumstruct.lwsum[bin] = total(lensum.lwsum[bin])
      FOR ib=0L, nbeta -1 DO BEGIN 
          sumstruct.npsum[bin,ib] = total(lensum.npsum[bin,ib])
      ENDFOR 
      sumstruct.wsum[bin] = total(lensum.wsum[bin])

  ENDFOR 

  combine_wthetalumw_sum, sumstruct, wstruct

  return
END 
