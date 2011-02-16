PRO zobjshear_lambda_get_sources, lambda, minlam, maxlam, w, num, issouth=issouth

  IF (maxlam GT 180.) OR (minlam LT -180.) OR $
    (minlam GT maxlam) THEN BEGIN ; We DO cross -180,180 mark
    ;;  print,'Crossed lambda=[180,-180] minlam = '+ntostr(minlam)+' maxlam = '+ntostr(maxlam)
      
      IF maxlam GT 180. THEN BEGIN 
          diff = maxlam - 180.
          w1 = where(lambda LE 180. AND $
                     lambda GE minlam, n1)
          w2 = where(lambda GE -180. AND $
                     lambda LE (-180.+diff), n2)
      ENDIF ELSE IF minlam LT -180 THEN BEGIN 
          diff = -180. - minlam
          w1 = where(lambda LE maxlam AND $
                     lambda GE -180., n1)
          w2 = where(lambda LE 180. AND $
                     lambda GE (180. - diff), n2)
      ENDIF ELSE BEGIN 
          w=where( (lambda GE minlam) OR $
                   (lambda LE maxlam), num)
          return
      ENDELSE 
      CASE 1 OF
          (n1 NE 0) AND (n2 NE 0): BEGIN 
              w=[w1,w2]
              num=n1+n2
          END 
          (n1 NE 0): BEGIN
              w=w1
              num=n1
          END 
          (n2 NE 0): BEGIN
              w=w2
              num=n2
          END 
          ELSE: BEGIN
              w=-1
              num=0
          END 
      ENDCASE 
  ENDIF ELSE BEGIN ; We don't cross -180,180 mark


      IF NOT keyword_set(issouth) THEN BEGIN 
          ntot=n_elements(lambda)
          w=lindgen(ntot)
          binary_search, lambda, minlam, i1, /round, /edgedefault
          binary_search, lambda, maxlam, i2, /round, /edgedefault
      ENDIF ELSE BEGIN 
          
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; need continuous chunk not crossing branch point:
          ;; if south, we won't cross 0.0, so if not crosing 
          ;; the branch point just see if we are on the 
          ;; positive or negative side
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF maxlam GT 0 THEN w=where(lambda GE 0.,ntot) $
          ELSE w=where(lambda LE 0.,ntot)
          IF ntot EQ 0 THEN BEGIN 
              w=-1
              num=0
              return
          ENDIF 
          binary_search, lambda[w], minlam, i1, /round, /edgedefault
          binary_search, lambda[w], maxlam, i2, /round, /edgedefault
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
