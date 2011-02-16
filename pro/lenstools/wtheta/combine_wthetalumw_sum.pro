PRO combine_wthetalumw_sum, sumstruct, outstruct

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: combine_wthetalumw_sum, sumstruct, outstruct'
      return
  ENDIF 

  nbin = n_elements(sumstruct.rsum)

  IF tag_exist(sumstruct, 'lerrsum1') THEN doerrsum=1 ELSE doerrsum=0

  nbeta = n_elements(sumstruct.beta)
  arrval = fltarr(nbin)
  outstruct = wthetaLumwStruct(arrval, sumstruct.beta)

  struct_assign, sumstruct, outstruct

  IF sumstruct.binsize EQ -1 THEN logbin=1 ELSE logbin=0

  darea = fltarr(nbin)
  R1 = darea & R2 = R1

  IF keyword_set(logbin) THEN BEGIN 
      rmin = sumstruct.rmin
      logrmin = alog10(rmin)
      logrmax = alog10(sumstruct.rmax)
      binsize_log = ( logrmax - logrmin )/nbin
      r_ratio = 10.^binsize_log
      
      FOR i=0L, nbin-1 DO BEGIN
          R1[i] = rmin*r_ratio^(i)
          R2[i] = rmin*r_ratio^(i+1)
          darea[i] = !pi*(R2[i]^2 - R1[i]^2)
      ENDFOR 
  ENDIF ELSE BEGIN 
      rmin = sumstruct.rmin
      FOR i=0L, nbin-1 DO BEGIN 
          R1[i] = rmin + i*sumstruct.binsize
          R2[i] = rmin + (i+1)*sumstruct.binsize
          darea[i] = !pi*(R2[i]^2 - R1[i]^2)
      ENDFOR 
  ENDELSE 
  outstruct.rmax_act = R2
  outstruct.rmin_act = R1

  FOR i=0L, nbin-1 DO BEGIN 

      outstruct.area[i] = darea[i]

      outstruct.meanr[i] = sumstruct.rsum[i]/sumstruct.npsum[i,0]
      ;; does not include central galaxy
      ml = sumstruct.lsum[i]/sumstruct.lwsum[i]
      outstruct.meanlum[i] = ml

      tarea = !pi*(R2[i]^2 - outstruct.rmin_act[0]^2)
      ia = outstruct.area[0:i]
      ia2=ia^2
      tml = total(sumstruct.lsum[0:i]*ia)/sumstruct.lwsum[i]/tarea

      outstruct.tmeanlum[i] = tml
      IF doerrsum THEN BEGIN 
          outstruct.meanlumerr[i]=sqrt($
                                        (sumstruct.lerrsum1[i] - $
                                         2.*sumstruct.lerrsum2[i]*ml + $
                                         sumstruct.lerrsum3[i]*ml^2 $
                                        ) > 0.  $
                                      )/sumstruct.lwsum[i]

          outstruct.tmeanlumerr[i]=sqrt($
                                      (total(sumstruct.lerrsum1[0:i]*ia2) - $
                                       2.*total(sumstruct.lerrsum2[0:i]*ia2)*tml + $
                                       total(sumstruct.lerrsum3[0:i]*ia2)*tml^2 $
                                      ) > 0.  $
                                       )/sumstruct.lwsum[i]/tarea
      ENDIF 

      FOR ib=0L, nbeta-1 DO BEGIN 
          ;; now an average per lens
          outstruct.npair[i,ib] = sumstruct.npsum[i,ib]/sumstruct.wsum[i]

          ;; Now totals within each rmax_act[i]
          outstruct.tnpair[i,ib]=total(sumstruct.npsum[0:i,ib])/total(sumstruct.wsum[0:i])
          
          outstruct.density[i,ib] = outstruct.npair[i,ib]/outstruct.area[i]
      ENDFOR 


  ENDFOR 

END 
