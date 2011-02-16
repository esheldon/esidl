PRO wtheta_calcmeanerr, lumlensum, meanlum, meanlumerr

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: wtheta_calcmeanerr, lumlensum, meanlum, meanlumerr'
      return
  ENDIF 

  nbin=n_elements(lumlensum[0].rsum)
  
  meanlum = fltarr(nbin)
  meanlumerr=fltarr(nbin)
  
  FOR j=0L, nbin-1 DO BEGIN 
      
      meanlum[j] = total(lumlensum.lsum[j])/total(lumlensum.lwsum[j])
      meanlumerr[j] = sqrt( total( $
                                   (lumlensum.lsum[j] - $
                                    lumlensum.lwsum[j]*meanlum[j])^2 $
                                 )$
                          )/total(lumlensum.lwsum[j])
      
  ENDFOR 
  
END 
