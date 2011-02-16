PRO combine_shear, files, outstruct
 
  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: combine_shear, files, outstruct'
      return
  ENDIF 

  nf=n_elements(files)
  FOR i=0, nf-1 DO BEGIN

      struct1 = mrdfits(files[i], 1, hdr, /silent)

      IF i EQ 0 THEN BEGIN

          struct = struct1
          nbin = n_elements(struct1.rsum)
          nbinuse = nbin
          struct1 = 0

      ENDIF ELSE BEGIN 

          nbin2 =  n_elements(struct1.rsum)
          IF nbin2 NE nbin THEN BEGIN 
              message,'Files must contain same number of bins'
          ENDIF 
          IF (  (struct.binsize NE struct1.binsize) OR $
                (struct.rmin NE struct1.rmin) ) THEN BEGIN 
              message,'Either binsize, rmin, or h has changed'
          ENDIF 

          struct.rsum = struct.rsum + struct1.rsum
          struct.etansum = struct.etansum + struct1.etansum
          struct.etanerrsum = struct.etanerrsum + struct1.etanerrsum
          struct.eradsum = struct.eradsum + struct1.eradsum
          struct.eraderrsum = struct.eraderrsum + struct1.eraderrsum

          struct.qrsum = struct.qrsum + struct1.qrsum
          struct.qrerrsum = struct.qrerrsum + struct1.qrerrsum
          struct.qisum = struct.qisum + struct1.qisum
          struct.qierrsum = struct.qierrsum + struct1.qierrsum

          struct.qradrsum = struct.qradrsum + struct1.qradrsum
          struct.qradrerrsum = struct.qradrerrsum + struct1.qradrerrsum
          struct.qradisum = struct.qradisum + struct1.qradisum
          struct.qradierrsum = struct.qradierrsum + struct1.qradierrsum

          struct.etan1sum = struct.etan1sum + struct1.etan1sum
          struct.etan1errsum = struct.etan1errsum + struct1.etan1errsum
          struct.etan2sum = struct.etan2sum + struct1.etan2sum
          struct.etan2errsum = struct.etan2errsum + struct1.etan2errsum
          struct.erad1sum = struct.erad1sum + struct1.erad1sum
          struct.erad1errsum = struct.erad1errsum + struct1.erad1errsum
          struct.erad2sum = struct.erad2sum + struct1.erad2sum
          struct.erad2errsum = struct.erad2errsum + struct1.erad2errsum

          struct.w1sum = struct.w1sum + struct1.w1sum
          struct.w2sum = struct.w2sum + struct1.w2sum

          struct.wsum = struct.wsum + struct1.wsum
          struct.npair = struct.npair + struct1.npair

          struct.nlenses = struct.nlenses + struct1.nlenses
          struct.totpairs = struct.totpairs + struct1.totpairs
          struct.Sshsum = struct.Sshsum + struct1.Sshsum

          struct.lenswsum = struct.lenswsum + struct1.lenswsum

          FOR jj=0L, nbin-1 DO BEGIN
              struct.rmax_act[jj] = max( [struct.rmax_act[jj], struct1.rmax_act[jj]] )
          ENDFOR 

          struct1 = 0
      ENDELSE                   ; check complete bin, also check if already done
      IF (struct.rsum[nbin-1] EQ 0.) AND (nbinuse NE nbin-1) THEN BEGIN 
          ;print,'Last bin is incomplete. Cutting'
          nbinuse=nbin-1
      ENDIF 
  ENDFOR 

  combine_sum, struct, outstruct

  return
END 
