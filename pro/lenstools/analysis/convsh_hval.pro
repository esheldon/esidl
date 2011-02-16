PRO convsh_hval, shstruct, hval

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: convsh_hval, shstruct, hval'
      return
  ENDIF 

  IF hval NE shstruct.h THEN BEGIN 

      hfac = shstruct.h/hval

      shstruct.binsize  = shstruct.binsize*hfac
      shstruct.rmin     = shstruct.rmin*hfac
      shstruct.rmax     = shstruct.rmax*hfac

      shstruct.meanscritinv = shstruct.scritinvsum/hfac
  
      shstruct.meanr = shstruct.meanr*hfac
      shstruct.rmax_act = shstruct.rmax_act*hfac

      shstruct.sigma = shstruct.sigma/hfac
      shstruct.orthosig = shstruct.ortho/hfac
      shstruct.sigmaerr = shstruct.sigmaerr/hfac
      shstruct.orthosigerr = shstruct.orthoerr/hfac

      shstruct.tsigma = shstruct.tsigma/hfac
      shstruct.torthosig = shstruct.tortho/hfac
      shstruct.tsigmaerr = shstruct.tsigmaerr/hfac
      shstruct.torthosigerr = shstruct.torthoerr/hfac

      shstruct.area = shstruct*hfac^2
      shstruct.density = shstruct.density/hfac^2

      shstruct.h = hval
  ENDIF 

  return
END 
