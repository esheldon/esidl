PRO combine_wtheta2, files, outstruct
  
  ;; this is for use with the weighted npsum

  ;; this combines sum files

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: combine_wtheta2, files, outstruct'
      return
  ENDIF 

  nf=n_elements(files)
  FOR i=0, nf-1 DO BEGIN

      sumstruct1 = mrdfits(files[i], 1, hdr, /silent)

      IF i EQ 0 THEN BEGIN

          sumstruct = sumstruct1
          nbin = n_elements(sumstruct1.rsum)
          nbinuse = nbin
          sumstruct1 = 0

      ENDIF ELSE BEGIN 

          nbin2 =  n_elements(sumstruct1.rsum)
          IF nbin2 NE nbin THEN BEGIN 
              message,'Files must contain same number of bins'
          ENDIF 
          IF (  (sumstruct.binsize NE sumstruct1.binsize) OR $
                (sumstruct.rmin NE sumstruct1.rmin) ) THEN BEGIN 
              message,'Either binsize, rmin, or h has changed'
          ENDIF 

          sumstruct.rsum = sumstruct.rsum + sumstruct1.rsum
          sumstruct.npsum = sumstruct.npsum + sumstruct1.npsum
          sumstruct.wsum = sumstruct.wsum + sumstruct1.wsum

          sumstruct.nlenses = sumstruct.nlenses + sumstruct1.nlenses
          sumstruct.totpairs = sumstruct.totpairs + sumstruct1.totpairs

          FOR jj=0L, nbin-1 DO BEGIN
              sumstruct.rmax_act[jj] = max( [sumstruct.rmax_act[jj], sumstruct1.rmax_act[jj]] )
          ENDFOR 

          sumstruct1 = 0
      ENDELSE                   ; check complete bin, also check if already done
      IF (sumstruct.rsum[nbin-1] EQ 0.) AND (nbinuse NE nbin-1) THEN BEGIN 
          ;print,'Last bin is incomplete. Cutting'
          nbinuse=nbin-1
      ENDIF 
  ENDFOR 
 
  arrval = fltarr(nbinuse)
  outstruct = wthetastruct(arrval)

  outstruct.nlenses  = sumstruct.nlenses
  outstruct.totpairs = sumstruct.totpairs
  outstruct.binsize  = sumstruct.binsize
  outstruct.rmin     = sumstruct.rmin
  outstruct.rmax     = sumstruct.rmax
  outstruct.h        = sumstruct.h

  npairold = 0
  FOR i=0L, nbinuse-1 DO BEGIN 

      ;; calculate area and density of background galaxies
      R1 = outstruct.rmin + i*outstruct.binsize
      R2 = outstruct.rmin + (i+1)*outstruct.binsize

      outstruct.area[i] = !pi*(R2^2 - R1^2)

      ;; now an average per lens
      outstruct.npair[i] = sumstruct.npsum[i]/sumstruct.wsum[i]
      ;; npair is average per lens, no need to divide by nlenses
      outstruct.density[i] = outstruct.npair[i]/outstruct.area[i];/outstruct.nlenses

      outstruct.rmax_act[i] = sumstruct.rmax_act[i]

      outstruct.meanr[i]    = sumstruct.rsum[i]/sumstruct.npsum[i]

      ;; Now totals within each rmax_act[i]
      outstruct.tnpair[i] = total(sumstruct.npsum[0:i])/total(sumstruct.wsum[0:i])

  ENDFOR 

  return
END 
