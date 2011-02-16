
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    
;       
; PURPOSE:
;    
;
; CALLING SEQUENCE:
;    
;
; INPUTS: 
;    
;
; OPTIONAL INPUTS:
;    
;
; KEYWORD PARAMETERS:
;    
;       
; OUTPUTS: 
;    
;
; OPTIONAL OUTPUTS:
;    
;
; CALLED ROUTINES:
;    
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO find_radec2, ra, dec, mrun, mrerun, mcamcol, mfield, southern=southern, $
                 silent=silent

  IF N_params() LT 2 THEN BEGIN 
     print,'-Syntax:  find_radec2, ra, dec, mrun, mrerun, mcamcol, mfield, southern=southern, silent=silent'
     print,''
     print,'Use doc_library,""  for more help.'  
     return
  ENDIF 

  eq2survey, ra, dec, lambda, eta
  stripe = eta2stripenum(eta, southern=southern)

  define_stripes

  command = 'run = !stripes.stripe'+ntostr(stripe)+'.run'
  tmp = execute(command)
  IF tmp EQ 0 THEN BEGIN 
      print,'No info on stripe '+ntostr(stripe)
      return
  ENDIF 
  command = 'rerun = !stripes.stripe'+ntostr(stripe)+'.rerun'
  tmp = execute(command)
  IF tmp EQ 0 THEN BEGIN 
      print,'No info on stripe '+ntostr(stripe)
      return
  ENDIF 

  nrun = n_elements(run)
  clr=2                         ;read red astrans file
  mincol = 0L
  minrow = 0L
  maxcol = 2048-1
  maxrow = 1489-1
  FOR i=0L, nrun DO BEGIN 
      FOR camcol=1, 6 DO BEGIN 
          trans=sdss_read('astrans', run[i], camcol, band=clr, rerun=rerun[i],$
              node=node, inc=inc, /silent)
          field = trans.field
          nfield = n_elements(field)
          FOR fi=0L, nfield-1 DO BEGIN 
              
              ;; assume square
              rowcol2munu, trans, field[fi], minrow, mincol, mu1, nu1
              gc2survey, mu1, nu1, node, inc, lambda1, eta1

              rowcol2munu, trans, field[fi], maxrow, maxcol, mu2, nu2
              gc2survey, mu2, nu2, node, inc, lambda2, eta2

              IF (lambda LE lambda2) AND (lambda GE lambda1) AND $
                                 (eta LE eta2) AND (eta GE eta1) THEN BEGIN 
                  mrun = run[i]
                  mrerun = rerun[i]
                  mcamcol = camcol
                  mfield = field[fi]
                  IF NOT keyword_set(silent) THEN BEGIN 
                      print,'Run: ',ntostr(mrun),' Rerun: ',ntostr(mrerun),$
                        ' Camcol: ',ntostr(mcamcol),' Field: ',ntostr(mfield)
                  ENDIF 
                  return
              ENDIF 
          ENDFOR 
      ENDFOR 
  ENDFOR 
  
  print,'RA,DEC not found'

return
END 
