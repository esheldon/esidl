PRO run_field_range

    run_status = sdss_runstatus()
    w = sdss_flag_select(run_status.flags, 'runstatus', {astrans_exist:'y'}) 

    if w[0] eq -1 then begin 
        print,'No good runs!'
        return
    endif 

    search_runs = run_status[w].run
    search_runs = search_runs[rem_dup(search_runs)]
    nw = n_elements(search_runs)

    colprint,search_runs

    for i=0l, nw-1 do begin 
        run = search_runs[i]
        print
        print,'------------------------------------------------------'
        print,'Processing Run: ',ntostr(run)
        print,'------------------------------------------------------'
        print
        field_range, run, bad=bad
        IF NOT bad THEN add_arrval, run, gsearch_runs
    endfor 

    print
    print,'Good search runs: '
    colprint,gsearch_runs

    return
end 
