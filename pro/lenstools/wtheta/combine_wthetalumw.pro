PRO combine_wthetalumw, files, outstruct
  
  ;; this combines sum files
  ;; this is for weighted npsum, with both lum and sigcrit weighting
  ;; this combines sum files

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: combine_wthetalumw, files, outstruct'
      return
  ENDIF 

  nf=n_elements(files)
  FOR i=0, nf-1 DO BEGIN

      print,'Reading file: ',files[i]
      sumstruct1 = mrdfits(files[i], 1, hdr, /silent, status=ok)
      IF ok NE 0 THEN message,'file not found: '+files[i]

      IF i EQ 0 THEN BEGIN

          IF tag_exist(sumstruct1, 'lerrsum1') THEN doerrsum=1 ELSE doerrsum=0
          sumstruct = sumstruct1
          nbin = n_elements(sumstruct1.rsum)
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
          sumstruct.lsum = sumstruct.lsum + sumstruct1.lsum

          IF doerrsum THEN BEGIN 
              sumstruct.lerrsum1 = sumstruct.lerrsum1 + sumstruct1.lerrsum1
              sumstruct.lerrsum2 = sumstruct.lerrsum2 + sumstruct1.lerrsum2
              sumstruct.lerrsum3 = sumstruct.lerrsum3 + sumstruct1.lerrsum3
          ENDIF 

          sumstruct.lwsum = sumstruct.lwsum + sumstruct1.lwsum
          sumstruct.npsum = sumstruct.npsum + sumstruct1.npsum
          sumstruct.wsum = sumstruct.wsum + sumstruct1.wsum

          sumstruct.nlenses = sumstruct.nlenses + sumstruct1.nlenses
          sumstruct.totpairs = sumstruct.totpairs + sumstruct1.totpairs

          FOR jj=0L, nbin-1 DO BEGIN
              sumstruct.rmax_act[jj] = max( [sumstruct.rmax_act[jj], sumstruct1.rmax_act[jj]] )
              sumstruct.rmin_act[jj] = min( [sumstruct.rmin_act[jj], sumstruct1.rmin_act[jj]] )
          ENDFOR 

          sumstruct1 = 0
      ENDELSE                   ; check complete bin, also check if already done

  ENDFOR 
 
  combine_wthetalumw_sum, sumstruct, outstruct

  return
END 
