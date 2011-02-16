PRO apply_primary_bound_cuts, stripes, clambda, ceta, masked, unmasked

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: apply_primary_bound_cuts, stripes/bound_array, clambda, ceta, masked, unmasked'
      return
  ENDIF 

  IF datatype(stripes) EQ 'STC' THEN BEGIN 
      bound_array = stripes 
  ENDIF ELSE BEGIN 
      primary_bound_multi, stripes, bound_array, overall_bound
  ENDELSE 

  nst = n_elements(stripes)
  nobj = n_elements(clambda)

  testarr = intarr(nobj)

  FOR i=0L, nst-1 DO BEGIN 
      good = where( (clambda GE bound_array[i].lammin) AND $
                    (clambda LE bound_array[i].lammax) AND $
                    (ceta GE bound_array[i].etamin) AND $
                    (ceta LE bound_array[i].etamax), $
                    ngood, comp=bad, ncomp=nbad)

      IF ngood NE 0 THEN BEGIN 
          testarr[good] = 1
      ENDIF 

  ENDFOR 
  
  unmasked = where(testarr NE 0, ncomp=ncomp, comp=masked)

END 
