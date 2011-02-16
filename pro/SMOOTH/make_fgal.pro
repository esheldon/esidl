PRO make_fgal, run1, run2

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

  IF n_params() LT 2 THEN BEGIN
      print,'-Syntax: make_fgal, run1, run2'
      print
      print,'Use doc_library,"make_fgal" for more help'
  ENDIF 

  st = create_struct('RA', 0d, 'DEC', 0d)

  r1str = ntostr(long(run1))
  r2str = ntostr(long(run2))

  dir = '/sdss4/data1/esheldon/CORRECTED/'
  clr = [1,2,3]
  colors = ['u','g','r','i','z']
  max_mag = 20.5

  nclr = n_elements(clr)

  FOR i=0, nclr-1 DO BEGIN

      Lg = 0
      Sg = 0
      ostruct = 0

      Lname = dir + 'run' + r1str + '_' + r2str + '_lensgal_' + $
              colors[clr[i]] + '_overlap.fit'
      Sname = dir + 'run' + r1str + '_' + r2str + '_srcgal_' + $
              colors[clr[i]] + '_overlap.fit'

      Fname = dir + 'run' + r1str + '_' + r2str + '_fgal_' + $
              colors[clr[i]] + '.fit'

      IF NOT exist(Lname) OR NOT exist(Sname) THEN BEGIN
          Lname = dir + 'run' + r1str + '_' + r2str + '_lensgal_' + $
                  colors[clr[i]] + '_overlap.fit'
          Sname = dir + 'run' + r1str + '_' + r2str + '_srcgal_' + $
                  colors[clr[i]] + '_overlap.fit'

          Fname = dir + 'run' + r1str + '_' + r2str + '_fgal_' + $
                  colors[clr[i]] + '.fit'
          IF NOT exist(Lname) OR NOT exist(Sname) THEN BEGIN
              print,'No overlap files exist for runs ',r1str,' ',r2str
              return
          ENDIF 
      ENDIF 
      Lg = mrdfits(Lname, 1)
      Sg = mrdfits(Sname, 1)

      nL = n_elements(Lg)

      wS = where(Sg.petrocounts LE max_mag, nS)

      IF nS NE 0 THEN BEGIN 
          ntot = nL + nS
      
          ostruct = replicate(st, ntot)

          ostruct[0:nL-1].ra = Lg.ra
          ostruct[0:nL-1].dec = Lg.dec
          ostruct[nL:nL+nS-1].ra = Sg[wS].ra
          ostruct[nL:nL+nS-1].dec = Sg[wS].dec

          mwrfits, ostruct, Fname, /create
      ENDIF ELSE print,'No gals found in source cat for ',$
        colors[clr[i]],' band'

  ENDFOR 

  return 
END 
