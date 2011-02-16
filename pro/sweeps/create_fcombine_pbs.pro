pro create_fcombine_pbs, run, noclobber=noclobber

    if n_elements(run) eq 0 then message,'usage: create_fcombine_pbs, run, /noclobber'

    ; these are all rerun=301
    window_runlist, runs, minscore=0.101

    ; for some reason these were in Nikhil's list.  They *are* rerun=301.
    runs = [runs,[2856, 4107, 4899]]
    runs = runs[sort(runs)]


    ; putting photoop/idlutils later in case of name space collisions
    setup=[ $
        'setup sdssidl -r ~/exports/sdssidl-work', $
        'setup esidl -r ~/exports/esidl-work', $
        'setup photoop v1_10_15',$
        'setup idlutils v5_4_21', $
        'export PHOTO_CALIB=/clusterfs/riemann/raid006/bosswork/groups/boss/calib/dr8_final', $
        'export PHOTO_RESOLVE=/clusterfs/riemann/raid006/bosswork/groups/boss/resolve/2010-05-23', $
        'export PHOTO_REDUX=/clusterfs/riemann/raid006/dr8/groups/boss/photo/redux']

    pbs_dir=filepath(root='/home/esheldon/pbs', sub='fcombine', run)
    fcombine_dir = filepath(root='/clusterfs/riemann/raid008/bosswork/groups/boss/target',sub='fcombine',run)

    file_mkdir, pbs_dir

    for i=0L, n_elements(runs)-1 do begin
        run=runs[i]

        for camcol=1,6 do begin
            job_name=string(f='("fcomb-",i06,"-",i0)', run,camcol)

            pbs_file = filepath(root=pbs_dir, job_name+'.pbs')
            print,pbs_file
           
            commands=['run='+string(run,f='(i0)'),$
                      'rerun=301',$
                      'camcol='+string(camcol,f='(i0)'), $
                      'outdir="'+fcombine_dir+'"', $
                      'noclobber='+string(keyword_set(noclobber),f='(i0)'), $
                      'fc = obj_new("sdss_fcombine", outdir=outdir, noclobber=noclobber)', $
                      'fc->process_camcol, run, rerun, camcol']
            pbs_riemann_idl, pbs_file, commands, setup=setup, job_name=job_name
        endfor
    endfor
end
