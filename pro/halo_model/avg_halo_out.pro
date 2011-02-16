PRO avg_halo_out, struct, wzstruct, r, deltasig, xigm, $
                  meanz=meanz, physical=physical

  ;; compute the average deltasig and xigm using out weights
  ;; and the dn/dz

  nz = n_elements(struct)

  wz = interpol(wzstruct.weight, wzstruct.z, struct.z)
  num = interpol(wzstruct.num, wzstruct.z, struct.z)

  weights = wz*num
  wsum = total(weights)

  meanz = total( weights*struct.z )/wsum

  IF NOT keyword_set(physical) THEN BEGIN 

      ;; Co-moving (which is the default output of wayne's code)

      nr = n_elements(struct[0].r)
      r = struct[0].r
      deltasig = fltarr(nr)
      xigm = deltasig
      
      FOR i=0L, nr-1 DO BEGIN 
          
          deltasig[i] = total( struct[*].deltasig[i]*weights )/wsum
          xigm[i] = total( struct[*].xigm[i]*weights )/wsum
          
      ENDFOR 

  ENDIF ELSE BEGIN 

      ;; Physical
      nr = n_elements(struct[0].rcommon)
      r = struct[0].rcommon
      deltasig = fltarr(nr)
      xigm = deltasig
      
      FOR i=0L, nr-1 DO BEGIN 
          
          deltasig[i] = total( struct[*].deltasig_physcommon[i]*weights )/wsum
          xigm[i] = total( struct[*].xigmcommon[i]*weights )/wsum
          
      ENDFOR 

  ENDELSE 
  

END 
