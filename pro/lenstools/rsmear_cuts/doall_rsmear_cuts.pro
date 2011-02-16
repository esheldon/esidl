PRO doall_rsmear_cuts, purity, overwrite=overwrite
  
    run_status = sdss_runstatus()
    si = sdss_flag_select(run_status.flags, $
            'runstatus', {adatc_exist:'y'}, nsi)

    for i=0l, nsi-1 do begin 

        delvarx, struct
        run = run_status[si[i]].run
        rerun = run_status[si[i]].rerun

        clr=1
        run_rsmear_cuts, run, rerun, purity, clr, struct, overwrite=overwrite

        clr=2
        run_rsmear_cuts, run, rerun, purity, clr, struct, overwrite=overwrite
  
        clr=3
        run_rsmear_cuts, run, rerun, purity, clr, struct, overwrite=overwrite

    endfor

end 
