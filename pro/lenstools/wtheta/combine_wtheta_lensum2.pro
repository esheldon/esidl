PRO combine_wtheta_lensum2, lensum, binsize, rminkpc, rmaxkpc, hval, sumstruct, wstruct

  ;; this is for weighted npsum

  IF n_params() LT 5 THEN BEGIN 
      print,'-Syntax: combine_wtheta_lensum2, lensum, binsize, rminkpc, rmaxkpc, hval, sumstruct, wstruct'
      return
  ENDIF 

  wz = getztag(lensum)

  nlens = n_elements(lensum)
  nbin = n_elements(lensum[0].rsum)
  arrval = fltarr(nbin)
  
  wstruct = wthetastruct(arrval)
  sumstruct = wthetasumstruct(arrval)

  sumstruct.nlenses = nlens
  sumstruct.binsize = binsize
  sumstruct.rmin = rminkpc
  sumstruct.rmax = rmaxkpc
  sumstruct.h    = hval

  sumstruct.totpairs = total( lensum.totpairs )

  FOR bin=0L, nbin-1 DO BEGIN 

      sumstruct.rmax_act[bin] = max(lensum.rmax_act[bin])
      sumstruct.rsum[bin] = total(lensum.rsum[bin])
      sumstruct.npsum[bin] = total(lensum.npsum[bin])
      sumstruct.wsum[bin] = total(lensum.wsum[bin])

  ENDFOR 

  wstruct.totpairs = sumstruct.totpairs
  wstruct.nlenses = nlens
  wstruct.binsize = binsize
  wstruct.rmin = rminkpc
  wstruct.rmax = rmaxkpc
  wstruct.h    = hval
  wstruct.rmax_act = sumstruct.rmax_act

  FOR i=0L, nbin-1 DO BEGIN 
      R1 = wstruct.rmin + i*wstruct.binsize
      R2 = wstruct.rmin + (i+1)*wstruct.binsize
  
      wstruct.area[i] = !pi*(R2^2 - R1^2)

      ;; now an average per lens
      wstruct.npair[i] = sumstruct.npsum[i]/sumstruct.wsum[i]

      ;; npair is average per lens, no need to divide by nlenses
      IF wstruct.area[i] NE 0. THEN $
        wstruct.density[i] = wstruct.npair[i]/wstruct.area[i];/wstruct.nlenses
      
      wstruct.meanr[i] = sumstruct.rsum[i]/sumstruct.npsum[i]
      
      ;; Now totals within each rmax_act[i]
      wstruct.tnpair[i] = total(sumstruct.npsum[0:i])/total(sumstruct.wsum[0:i])

  ENDFOR 

  return
END 
