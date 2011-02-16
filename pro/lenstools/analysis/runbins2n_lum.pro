PRO runbins2n_lum, lensum, nbin, cutstruct

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: runbins2n_lum, lensum, nbin, cutstruct'
      return
  ENDIF 

  nrand = 10000L

  arrval = fltarr(nbin-1)
  cutstruct=create_struct('ucuts',arrval,$
                          'gcuts',arrval,$
                          'rcuts',arrval,$
                          'icuts',arrval,$
                          'zcuts',arrval)

  maxlum = [10.0e10, 10.e10, 15.0e10, 30.0e10, 45.0e10]
  
  FOR clr=0,4 DO BEGIN 

      w=where( lensum.lum[clr] GT 0.0 AND $
              lensum.lum[clr] LT maxlum[clr] )

      bins2n_lum, lensum[w], nbin, nrand, 'lum', element=clr,cuts=cuts

      CASE clr OF
          0: cutstruct.ucuts = cuts
          1: cutstruct.gcuts = cuts
          2: cutstruct.rcuts = cuts
          3: cutstruct.icuts = cuts
          ELSE: cutstruct.zcuts = cuts
      ENDCASE 

  ENDFOR 
END 
