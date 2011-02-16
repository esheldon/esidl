PRO run_rsmear_cuts, run, rerun, purity, clr, struct, $
                     overwrite=overwrite, hirata=hirata

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: run_rsmear_cuts, run, rerun, purity, struct, overwrite=overwrite, hirata=hirata'
      return
  ENDIF 

  ;; First time, reads in all camcols
  rsmear_cuts, run, rerun, purity, clr, struct, meanmag, all_rcuts, $
               overwrite=overwrite, hirata=hirata

  ;; now can do selection on camcol (or if struct wasn't there it would
  ;; just read in that camcol)
  rsmear_cuts, run, rerun, purity, clr, struct, meanmag, rcuts1, camcol=1, $
               overwrite=overwrite, hirata=hirata
  rsmear_cuts, run, rerun, purity, clr, struct, meanmag, rcuts2, camcol=2, $
               overwrite=overwrite, hirata=hirata
  rsmear_cuts, run, rerun, purity, clr, struct, meanmag, rcuts3, camcol=3, $
               overwrite=overwrite, hirata=hirata
  rsmear_cuts, run, rerun, purity, clr, struct, meanmag, rcuts4, camcol=4, $
               overwrite=overwrite, hirata=hirata
  rsmear_cuts, run, rerun, purity, clr, struct, meanmag, rcuts5, camcol=5, $
               overwrite=overwrite, hirata=hirata
  rsmear_cuts, run, rerun, purity, clr, struct, meanmag, rcuts6, camcol=6, $
               overwrite=overwrite, hirata=hirata

  rsmear_cuts_plotall, run, rerun, purity, clr, hirata=hirata

END 
