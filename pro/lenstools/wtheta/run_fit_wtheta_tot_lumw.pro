PRO run_fit_wtheta_tot_lumw, sigonly=sigonly, noweight=noweight

  FOR wclr=0,4 DO BEGIN

      fit_wtheta_tot_lumw,wclr,wtlen,wtrlen, sigonly=sigonly, noweight=noweight,$
        /makeplot

      delvarx, wtlen, wtrlen
  ENDFOR 

END 
