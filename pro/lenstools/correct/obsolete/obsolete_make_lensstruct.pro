PRO make_lensstruct, lensstr
  
  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: make_lensstruct, lensstr'
      return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure out which taglist to use
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  arrval = replicate(1.e10, 5)
  arrval2 = replicate(-1.,5)
  lensstr = create_struct('SEEING', fltarr(5), $
                          'SKY', fltarr(5), $
                          'SKYERR', fltarr(5),$
                          'STARFLAG', 0b, $
                          'E_D_BIT', 0b, $
                          'IXX',fltarr(5), $
                          'IYY',fltarr(5), $
                          'IXY',fltarr(5), $
                          'RHO4',fltarr(5), $
                          'WHYFLAG',bytarr(5), $
                          'PSFIXX', fltarr(5), $
                          'PSFIYY', fltarr(5), $
                          'PSFIXY', fltarr(5), $
                          'PSFRHO4', fltarr(5),$
                          'E1',arrval, $
                          'E2',arrval, $
                          'MOMERR',fltarr(5), $
                          'ROTATION', fltarr(5), $
                          'R',arrval2)

return
END 
