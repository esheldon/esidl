PRO fixwthetaconvnan, dstruct

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: fixwthetaconvnan, dstruct'
      return
  ENDIF 

  nn=n_elements(dstruct)

  nr = n_elements(dstruct[0].r)
  FOR i=0L, nn-1 DO BEGIN 

      den_nan = where( finite( dstruct[i].density, /NAN ) , nden)
      
      IF nden NE 0 THEN BEGIN 
          print,"found "+ntostr(nden)+" NAN's at i = "+ntostr(i)+'  Fixing'
          good = lindgen(nr)
          remove, den_nan, good

          ;; interpolate to get correct value
          fixedvals = interpol( dstruct[i].density[good], dstruct[i].r[good], $
                                dstruct[i].r[den_nan] )
          dstruct[i].density[den_nan] = fixedvals

          ;; now re-integrate and set values
          intfunc, dstruct[i].density*2d*!dpi*dstruct[i].r, dstruct[i].r, $
            dens_int, /double
          dens_int = dens_int/!dpi/dstruct[i].r^2

          ;; copy on where we have NAN
          den_int_nan =  where( finite( dstruct[i].density_int, /NAN ) , nden_int)
          dstruct[i].density_int[den_int_nan] = dens_int[den_int_nan]

      ENDIF 

  ENDFOR 


END 
