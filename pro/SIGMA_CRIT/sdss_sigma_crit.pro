
;; A convenience function.  Calls sdss_sigma_crit method.
FUNCTION sdss_sigma_crit, omegamat, zlens, $
                          h=h, $
                          scinvstruct=scinvstruct, $
                          npts=npts, $
                          wgood=wgood, cgs=cgs, $
                          photoz=photoz, $
                          limit=limit

  on_error,2

  IF N_params() LT 2 THEN BEGIN 
     print,'-Syntax: result = sdss_sigma_crit(omegmat, zlens, h=, npts=, scinvstruct=, wgood=, /cgs, /photoz)'
     return,-1
  ENDIF 

  sc = obj_new('sdss_sigma_crit')
  scritinv = sc->sigmacritinv(omegamat, zlens, $
                              h=h, $
                              scinvstruct=scinvstruct, $
                              npts=npts, $
                              wgood=wgood, cgs=cgs, $
                              photoz=photoz,$
                              limit=limit)
  obj_destroy, sc
  return, scritinv
  
END 
