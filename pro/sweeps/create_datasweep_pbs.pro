pro create_datasweep_pbs, sweep_run
    if n_elements(sweep_run) eq 0 then message,'usage: create_fcombine_pbs, sweep_run, /noclobber'

    ; these are all rerun=301
    window_runlist, runs, minscore=0.101

    ; for some reason these were in Nikhil's list.  They *are* rerun=301.
    ;runs = [runs,[2856, 4107, 4899]]
    runs = [[2856, 4107, 4899],runs]
    ;runs = runs[sort(runs)]
    ;for i=0L, n_elements(runs)-1 do begin
    ;    for camcol=1,6 do begin
    ;        gsa = sdss_astrom(runnum, camcol[icol], fields, rerun=301)
    ;    endfor
    ;endfor

    ; putting photoop/idlutils later in case of name space collisions
    ; note explicitly setting lots of directories while the disk failure stuff
    ; gets settled
    setup=[ $
           'setup sdssidl -r ~/exports/sdssidl-work', $
           'setup esidl -r ~/exports/esidl-work', $
           ;'setup photoop v1_10_15',$
           'setup photoop -r ~/svn/photoop-partial',$
           'setup idlutils v5_4_21', $
           'export BOSS_PHOTOOBJ=/clusterfs/riemann/raid006/dr8/groups/boss/photoObj', $
           'export PHOTO_CALIB=/clusterfs/riemann/raid006/bosswork/groups/boss/calib/dr8_final', $
           'export PHOTO_RESOLVE=/clusterfs/riemann/raid006/bosswork/groups/boss/resolve/2010-05-23', $
           'export PHOTO_REDUX=/clusterfs/riemann/raid006/dr8/groups/boss/photo/redux', $
           'export PHOTO_SWEEP=/clusterfs/riemann/raid008/bosswork/groups/boss/target/sweeps/'+sweep_run+'/dr8_final']


    outdir='/home/esheldon/pbs/datasweep/'+sweep_run
    file_mkdir, outdir
    for i=0L, n_elements(runs)-1 do begin
        run=runs[i]
        job_name=string(f='("sweep-",i06)', run)

        pbs_file = filepath(root=outdir, job_name+'.pbs')
        print,pbs_file



        commands=['run='+string(run,f='(i0)'),$
                  'rerun=301',$
                  'datasweep, run, rerun=rerun, /uber, /combine, /clobber, maxmem=5000.0']
        pbs_riemann_idl, pbs_file, commands, setup=setup, job_name=job_name
    endfor

    ; now the indexes and doing checks

    types=['sky', 'star', 'gal']
    for i=0L, n_elements(types)-1 do begin

        ; creating indexes
        job_name='sweep-index-'+types[i]
        pbs_file = filepath(root=outdir, job_name+'.pbs')
        print,pbs_file

        commands=['type="'+types[i]+'"', $
                  'datasweep_index, type=type']

        pbs_riemann_idl, pbs_file, commands, setup=setup, job_name=job_name

        ; creating indexes
        job_name='sweep-check-'+types[i]
        pbs_file = filepath(root=outdir, job_name+'.pbs')
        print,pbs_file

        commands=['type="'+types[i]+'"', $
                  'sweep_fcombine_check, type']

        pbs_riemann_idl, pbs_file, commands, setup=setup, job_name=job_name

    endfor
end
