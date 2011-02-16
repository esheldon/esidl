pro maxbcg_lensing_runsubs, samples, subtype, proc2run=proc2run

    if n_elements(proc2run) eq 0 then begin
        proc2run=['sub_sample','sub_sample_rand','combine_allrand','corr','jackknife']
    endif
    nsamp=n_elements(samples)

    nrand=24

    for i=0L, nsamp-1 do begin
        sample=samples[i]
        m=obj_new('maxbcg_lensing', sample)

        if in(proc2run,'sub_sample') then m->sub_sample, subtype
        if in(proc2run,'sub_sample_rand') then m->sub_sample, subtype, randnum=lindgen(nrand)

        if in(proc2run,'combine_allrand') then m->combine_allrand, subtype=subtype
        if in(proc2run,'corr') then m->corr, subtype=subtype
        if in(proc2run,'jackknife') then m->jackknife, subtype=subtype

        obj_destroy, m
    endfor

end
