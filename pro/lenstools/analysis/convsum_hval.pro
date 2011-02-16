PRO convsum_hval, sumstruct, hval

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: convsum_hval, shstruct, hval'
      return
  ENDIF 

  IF hval NE sumstruct.h THEN BEGIN 

      hfac = sumstruct.h/hval

      sumstruct.binsize  = sumstruct.binsize*hfac
      sumstruct.rmin     = sumstruct.rmin*hfac
      sumstruct.rmax     = sumstruct.rmax*hfac

      sumstruct.scritinvsum = sumstruct.scritinvsum/hfac
  
      sumstruct.rsum = sumstruct.meanr*hfac
      sumstruct.rmax_act = sumstruct.rmax_act*hfac

      sumstruct.tansigsum = sumstruct.tansigsum/hfac
      sumstruct.radsigsum = sumstruct.radsigsum/hfac
      sumstruct.tansigerrsum = sumstruct.tansigerrsum/hfac
      sumstruct.radsigerrsum = sumstruct.radsigerrsum/hfac

      sumstruct.h = hval
  ENDIF 

  return
END 
