FUNCTION binarea, d_edge, R1, R2, darea=darea, d_edge2=d_edge2, hit=hit

  IF n_elements(d_edge2) NE 0 THEN doboth=1 ELSE doboth=0

  darea = !pi*(R2^2 - R1^2)
  hit=0

  ;; Are we over the edge?
  IF R2 GT d_edge THEN BEGIN 
      hit=hit+1

      adiff2 = R2^2*acos(d_edge/R2) - d_edge*sqrt(R2^2-d_edge^2)
      
      ;; also check inner ring
      IF R1 GT d_edge THEN BEGIN 

          adiff1 =  R1^2*acos(d_edge/R1) - d_edge*sqrt(R1^2-d_edge^2)

      ENDIF ELSE adiff1 = 0.

      adiff = adiff2 - adiff1

      barea = darea-adiff

  ENDIF ELSE BEGIN 

      barea = darea

  ENDELSE 

  ;; may be another edge on opposite side to check
  IF doboth THEN BEGIN 

      IF R2 GT d_edge2 THEN BEGIN 
          
          hit=hit+1

          adiff2 = R2^2*acos(d_edge2/R2) - d_edge2*sqrt(R2^2-d_edge2^2)
          
          ;; also check inner ring
          IF R1 GT d_edge2 THEN BEGIN 
              
              adiff1 =  R1^2*acos(d_edge2/R1) - d_edge2*sqrt(R1^2-d_edge2^2)
              
          ENDIF ELSE adiff1 = 0.
          
          second_adiff = adiff2 - adiff1
                    
      ENDIF ELSE BEGIN 

          second_adiff = 0.
          
      ENDELSE 

  ENDIF 

  IF doboth THEN barea = barea - second_adiff

  return,barea

END 
