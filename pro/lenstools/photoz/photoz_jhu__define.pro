function photoz_jhu::init
    return,1
end

function photoz_jhu::dir
    return,sdssidl_config('photoz_jhu_dir')
end 

function photoz_jhu::file, stripes, all=all, count=count, big=big, check=check
    count = 0
    dir = self->dir()

    if keyword_set(all) then begin
        pattern = 's*pz_esheldon.fit'
        files = file_search(dir, pattern, count=count)
        return, files
    endif


    if keyword_set(big) then begin
        tfiles = 'photoz_jhu.st'
    endif else begin
        if n_elements(stripes) eq 0 then begin
            print,'-Syntax: files=pj->file(stripes, /all, /big, /check, count=)'
            on_error, 2
            message,'Halting'
        endif
        tfiles = 's'+ntostr(stripes, format='(I2.2)')+'pz_esheldon.fit'
    endelse

    tfiles = concat_dir(dir, tfiles)

    if keyword_set(check) then begin
        for i=0L, n_elements(tfiles)-1 do begin
            if fexist(tfiles[i]) then begin
                add_arrval, tfiles[i], files
            endif
        endfor
    endif else begin
        files = tfiles
    endelse

    if n_elements(files) eq 0 then begin
        return, -1
    endif else begin
        count = n_elements(files)
        return, files
    endelse
end

function photoz_jhu::read, stripes, all=all, big=big, count=count
    count=0
    files = self->file(stripes, all=all, big=big, count=nf, /check)
    if nf eq 0 then begin
        return, -1
    endif
    if keyword_set(big) then begin
        return, read_idlstruct(files)
    endif else begin
        return, mrdfits_multi(files, count=count) 
    endelse
end

function photoz_jhu::stuffstruct, number
    st = $
        {                       $
            photoid:    0LL,    $
            run:        0,      $
            rerun:      0,      $
            camcol:     0,      $
            field:      0,      $
            id:         0,      $
            stripe:     0,      $
            ra:         0d,     $
            dec:        0d,     $
            z:          0.0,    $
            zerr:       0.0,    $
            t:          0.0,    $
            terr:       0.0,    $
            class:      0,      $
            chisq:      0.0,    $
            quality:    0       $
        }

    if n_elements(number) ne 0 then begin
        st = replicate(st, number)
    endif

    return, st

end


; make one gigantic file.  This will allow us to remove duplicates.
; run with /rmdup to remove the duplicates.  Make sure you run on 
; a high memory machine
pro photoz_jhu::make_big_file, rmdup=rmdup

    outfile = self->file(/big)

    if keyword_set(rmdup) then begin
        print,'Reading file: ',outfile
        t=read_idlstruct(outfile)

        norig = n_elements(t)

        pid = t.photoid
        print,'Getting unique ids'
        uid = uniq( pid, sort(pid) )
        pid = 0
        print,'removing dups'
        t = t[uid]

        nkeep = n_elements(t)   

        print,'Kept '+ntostr(nkeep)+'/'+ntostr(norig)

        print,'Writing file: ',outfile
        write_idlstruct, t, outfile
        return
    endif

    ; Delete the file first
    file_delete, outfile, /quiet

    files = self->file(/all, count=nf)

    for i=0L, nf-1 do begin 
        print,'------------------------------------------------------------'
        print,'Reading file: ',files[i]
        
        st = mrdfits(files[i],1)
        num=n_elements(st)

        outst = self->stuffstruct(num)

        copy_struct, st, outst

        casid_extract, st.objid, run, rerun, camcol, field, id
        outst.run=run
        outst.rerun=rerun
        outst.camcol=camcol
        outst.field=field
        outst.id=id
        outst.photoid = photoid(run,rerun,camcol,field,id)

        write_idlstruct, outst, outfile, /append

    endfor 

end



pro photoz_jhu::stuff

    pg = obj_new('postgres')

    table = 'zphot2'
    conn='user=postgres'
    tmpdir = self->dir()

    struct = self->read(/big)

    pg->struct2table, struct, table, primary_key='photoid', conn=conn, $
                tmpdir=tmpdir, status=status

    if status ne 0 then message,'Stuff failed'

    ;; multi-column index on (run,rerun,camcol,field,id)
    query = $
        'CREATE INDEX '+table+'_rrcfi_index ON '+table+' (run,rerun,camcol,field,id)'
    print,query
    pg->query, query, status=status,conn=conn

    if status ne pg->status_val('no_result') then begin 
        message,'CREATE INDEX failed'
    endif 

    ; Now the rest of the indices
    pg->create_index, table, ['stripe','z','zerr','t','class','chisq','quality'], conn=conn
    
    query = 'ANALYZE '+table
    print,query
    pg->query, query, status=status,conn=conn

    if status ne pg->status_val('no_result') then begin 
        message,'ANALYZE failed'
    endif 

    ; Make sure the sdss user can select
    pg->query, 'GRANT select ON '+table+' TO sdss', conn=conn, status=gstatus
    if gstatus ne pg->status_val('no_result') then message,'Could not grant select to sdss'

    obj_destroy, pg

end 




pro photoz_jhu__define
    struct = { $
        photoz_jhu, $
        photoz_jhu_dummy:0 $
    }
end
