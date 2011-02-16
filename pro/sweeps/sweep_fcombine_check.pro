pro sweep_fcombine_check, type
    rerun='301'
    old_sweep='/clusterfs/riemann/netapp/dr8/groups/boss/sweeps/dr8_final'
    new_sweep='/clusterfs/riemann/raid008/bosswork/groups/boss/target/sweeps/dr8_final'

    ;; read in the files, sorting them to be sure
    allfiles= file_search(old_sweep+'/'+rerun+'/calibObj-*-*-'+type+'.fits.gz')

    for i=0L, n_elements(allfiles)-1 do begin
        oldfile = allfiles[i]
        bname = file_basename(oldfile)

        newfile = new_sweep+'/'+rerun+'/'+bname
        print,oldfile
        old = mrdfits(oldfile,1)
        new = mrdfits(newfile,1)

        if n_elements(old) ne n_elements(new) then begin
            message,'Different row count'
        endif

        cs = compare_struct(old, new)
        w=where(cs.ndiff gt 0, nw)
        if nw ne 0 then begin
            message,'Difference found'
        endif
    endfor
end
