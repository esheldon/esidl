PRO des_psf_sampling_getsets, seeing_set, apix_set, seeing, apix, all=all



  ;; The old way, where we tried lots of different combinations

;  IF seeing_set EQ 1 THEN BEGIN 
;      seeing = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75]
;  ENDIF ELSE BEGIN 
;      seeing = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10]
;  ENDELSE 

;  IF apix_set EQ 1 THEN BEGIN 
;      apix = [0.250, 0.260, 0.270, 0.280, 0.290, 0.300, 0.310, 0.320]
;  ENDIF ELSE BEGIN 
;      apix = [0.255, 0.265, 0.275, 0.285, 0.295, 0.305, 0.315, 0.325]
 ; ENDELSE 


  ;; This way we just fix the seeing and vary the pixel scale

  IF keyword_set(all) THEN BEGIN 

      IF n_elements(seeing_set) EQ 0 THEN message,'For /all you still need to input the seeing_set'

      ;; only one seeing set right now

      FOR apix_set=1,4 DO BEGIN 
          des_psf_sampling_getsets, seeing_set, apix_set, $
            seeing, tapix
          
          add_arrval, tapix, apix
      ENDFOR 


      ;; This also sorts
      apix   = apix[rem_dup(apix)]

      ;; rem_dup always returns an array for some stupid reason
      IF n_elements(apix) EQ 1 THEN apix=apix[0]
      
      return
  ENDIF ELSE BEGIN 

      IF seeing_set EQ 1 THEN BEGIN 
          seeing = 0.9
      ENDIF ELSE BEGIN 
          seeing = 0.7
      ENDELSE 
      
      apix_step = 0.005
      napix = 20
      
      CASE apix_set OF
          1: apix = 0.2 + apix_step*findgen(napix)
          2: apix = 0.3 + apix_step*findgen(napix)
          3: apix = 0.4 + apix_step*findgen(napix)
          4: apix = 0.5 + apix_step*findgen(napix)
          ELSE: message,'Unknown apix_set: '+ntostr(apix_set)
      ENDCASE 

  ENDELSE 

END 
