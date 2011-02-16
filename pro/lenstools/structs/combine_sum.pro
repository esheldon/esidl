PRO combine_sum, sumstruct, outstruct

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: combine_sum, sumstruct, outstruct'
      return
  ENDIF 

  nbin = n_elements(sumstruct[0].rsum)

  IF (sumstruct.rsum[nbin-1] EQ 0.) THEN nbinuse=nbin-1 ELSE nbinuse=nbin

  arrval = fltarr(nbinuse)
  outstruct = shstruct(arrval)

  struct_assign, sumstruct, outstruct

  outstruct.Ssh = sumstruct.Sshsum/total(sumstruct.wsum)

  npairold = 0
  FOR i=0L, nbinuse-1 DO BEGIN 

      ;; calculate area and density of background galaxies
      R1 = outstruct.rmin + i*outstruct.binsize
      R2 = outstruct.rmin + (i+1)*outstruct.binsize

      outstruct.area[i] = !pi*(R2^2 - R1^2)

      outstruct.npair[i] = sumstruct.npair[i]
      outstruct.tnpair[i] = total(sumstruct.npair[0:i])

      IF outstruct.area[i] NE 0.0 THEN BEGIN 
          outstruct.density[i] = $
            outstruct.npair[i]/outstruct.area[i]/outstruct.nlenses
      ENDIF 

      outstruct.rmax_act[i] = sumstruct.rmax_act[i]
      outstruct.meanr[i]    = sumstruct.rsum[i]/sumstruct.npair[i]

      outstruct.shear[i]    = sumstruct.etansum[i]/sumstruct.wsum[i]/2.
      outstruct.ortho[i]    = sumstruct.eradsum[i]/sumstruct.wsum[i]/2.
      outstruct.shearerr[i] = sqrt( sumstruct.etanerrsum[i]/sumstruct.wsum[i]^2)/2.
      outstruct.orthoerr[i] = sqrt( sumstruct.eraderrsum[i]/sumstruct.wsum[i]^2)/2.

      ;; Now totals within each rmax_act[i]
      outstruct.tshear[i] = total(sumstruct.etansum[0:i])/total(sumstruct.wsum[0:i])/2.
      outstruct.tshearerr[i] = sqrt( total(sumstruct.etanerrsum[0:i])/total(sumstruct.wsum[0:i])^2 )/2.
      outstruct.tortho[i] = total(sumstruct.eradsum[0:i])/total(sumstruct.wsum[0:i])/2.
      outstruct.torthoerr[i] = sqrt( total(sumstruct.eraderrsum[0:i])/total(sumstruct.wsum[0:i])^2 )/2.

      ;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; halo shape stuff
      ;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; quadrupole moments
      outstruct.tqr[i] = total(sumstruct.qrsum[0:i])/total(sumstruct.wsum[0:i])/2.
      outstruct.tqrerr[i] = sqrt( total(sumstruct.qrerrsum[0:i])/total(sumstruct.wsum[0:i])^2 )/2.
      outstruct.tqi[i] = total(sumstruct.qisum[0:i])/total(sumstruct.wsum[0:i])/2.
      outstruct.tqierr[i] = sqrt( total(sumstruct.qierrsum[0:i])/total(sumstruct.wsum[0:i])^2 )/2.

      outstruct.tqradr[i] = total(sumstruct.qradrsum[0:i])/total(sumstruct.wsum[0:i])/2.
      outstruct.tqradrerr[i] = sqrt( total(sumstruct.qradrerrsum[0:i])/total(sumstruct.wsum[0:i])^2 )/2.
      outstruct.tqradi[i] = total(sumstruct.qradisum[0:i])/total(sumstruct.wsum[0:i])/2.
      outstruct.tqradierr[i] = sqrt( total(sumstruct.qradierrsum[0:i])/total(sumstruct.wsum[0:i])^2 )/2.

      ;; two half spaces
      outstruct.tshear1[i] = total(sumstruct.etan1sum[0:i])/total(sumstruct.w1sum[0:i])/2.
      outstruct.tshear1err[i] = sqrt( total(sumstruct.etan1errsum[0:i])/total(sumstruct.w1sum[0:i])^2  )/2.
      outstruct.tshear2[i] = total(sumstruct.etan2sum[0:i])/total(sumstruct.w2sum[0:i])/2.
      outstruct.tshear2err[i] = sqrt( total(sumstruct.etan2errsum[0:i])/total(sumstruct.w2sum[0:i])^2  )/2.

      outstruct.tortho1[i] = total(sumstruct.erad1sum[0:i])/total(sumstruct.w1sum[0:i])/2.
      outstruct.tortho1err[i] = sqrt( total(sumstruct.erad1errsum[0:i])/total(sumstruct.w1sum[0:i])^2  )/2.
      outstruct.tortho2[i] = total(sumstruct.erad2sum[0:i])/total(sumstruct.w2sum[0:i])/2.
      outstruct.tortho2err[i] = sqrt( total(sumstruct.erad2errsum[0:i])/total(sumstruct.w2sum[0:i])^2  )/2.



  ENDFOR 

  return
END 
