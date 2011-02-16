PRO run_make_wthetaconv_func_lumw, sigonly=sigonly, noweight=noweight

  FOR wclr=0,4 DO BEGIN 
      FOR clr=1,3 DO BEGIN 
          make_wthetaconv_func_lumw, clr, wclr, sigonly=sigonly, noweight=noweight
      ENDFOR 
  ENDFOR 

END 
