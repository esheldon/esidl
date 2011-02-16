PRO zobjshear_get_sources, ra, minra, maxra, LASTRA, w, num, issouth=issouth

  IF (maxra GT 360.0d) OR (minra LT 0.0d) OR $
    (minra GT maxra) THEN BEGIN ; We DO cross [0,360] mark
    ;;  print,'Crossed ra=[180,-180] minra = '+ntostr(minra)+' maxra = '+ntostr(maxra)
      
      IF maxra GT 360. THEN BEGIN 
          diff = maxra - 360.0d
          w=where( (ra GE minra) OR (ra LE diff), num)
          return
      ENDIF ELSE IF minra LT 0.0d THEN BEGIN 
          diff = -minra
          w=where( (ra LE maxra) OR (ra GE 360.0-diff), num)
          return
      ENDIF ELSE BEGIN 
          w=where( (ra GE minra) OR (ra LE maxra), num)
          return
      ENDELSE 
  ENDIF ELSE BEGIN ; We don't cross [0,360] point


      IF NOT keyword_set(issouth) THEN BEGIN 
          ntot=n_elements(ra)
          w=lindgen(ntot)
          binary_search, ra, minra, i1, /round, /edgedefault
          binary_search, ra, maxra, i2, /round, /edgedefault
      ENDIF ELSE BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; need continuous chunk not crossing branch point:
          ;; if south, we won't cross 0.0, so if not crosing 
          ;; the branch point just see if we are greater than
          ;; lastra
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF minra GT LASTRA THEN w=where(ra GT LASTRA,ntot) $
          ELSE w=where(ra LE LASTRA,ntot)
          IF ntot EQ 0 THEN BEGIN 
              w=-1
              num=0
              return
          ENDIF 
          binary_search, ra[w], minra, i1, /round, /edgedefault
          binary_search, ra[w], maxra, i2, /round, /edgedefault
      ENDELSE 
      CASE 1 OF
          (i1 NE -1) AND (i2 NE -1): BEGIN
              w=w[i1:i2]
              num=n_elements(w)
          END 
          (i1 EQ -1) AND (i2 NE -1): BEGIN
              w=w[0:i2]
              num=n_elements(w)
          END 
          (i2 EQ -1) AND (i1 NE -1): BEGIN
              w=w[i1:ntot-1]
              num=n_elements(w)
          END 
          ELSE: BEGIN
              w=-1
              num=0
          END 
      ENDCASE 
  ENDELSE 
  return
END 
