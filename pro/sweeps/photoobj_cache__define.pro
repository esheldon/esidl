function photoobj_cache::init, nmax=nmax, verbose=verbose
    self.verbose = keyword_set(verbose)
    self->init_cache, nmax=nmax
    return, 1
end

pro photoobj_cache::init_cache, nmax=nmax

    heap_free, self.cache

    if n_elements(nmax) eq 0 then nmax=200L
    if nmax eq 0 then message,'Cannot initialize with nmax=0'

    pobj0={run:0, $
           field:0, $
           camcol:0, $
           rerun:' ', $
           fieldid:0ll, $
           obj:ptr_new()}
    cache = replicate(pobj0, nmax)

    self.cache = ptr_new(cache, /no_copy)
    self.nmax = nmax
    self.ip = 0L

end

pro photoobj_cache::add2cache, run, camcol, field, rerun=rerun, readobj=readobj


    ;; check if already cached, if so skip
    fieldid= sdss_fieldid(run, camcol, field, rerun=rerun)
    ic= where(fieldid eq (*self.cache).fieldid, nc)
    if(nc gt 0) then return

    ;; replace the current spot
    (*self.cache)[self.ip].run=run
    (*self.cache)[self.ip].camcol=camcol
    (*self.cache)[self.ip].field=field
    (*self.cache)[self.ip].rerun=rerun
    (*self.cache)[self.ip].fieldid=fieldid
    ptr_free, (*self.cache)[self.ip].obj
    pfile= sdss_name('photoObj', run, camcol, field, rerun=rerun)
    if(file_test(pfile)) then begin
        if self.verbose then begin
            splog,format='("Cacheing photoobj: ",i0," ",i0," ",i0," ",i0)',run,rerun,camcol,field
        endif
        obj= mrdfits(pfile,1,/silent)
    endif else begin
       if(keyword_set(readobj)) then begin

            if self.verbose then begin
                splog,format='("Generating photoobj: ",i0," ",i0," ",i0," ",i0)',run,rerun,camcol,field
            endif
           readobj= sdss_readobj(run, camcol, field, rerun=rerun, $
                                 except='TEXTURE',/silent)
           obj= photoobj_table(readobj, /unsafenan)
           readobj=0
       endif else begin
          message, 'No photoObj file: '+pfile
       endelse
    endelse
    (*self.cache)[self.ip].obj= ptr_new(obj, /no_copy)

    ;; update current spot
    self.ip= (self.ip+1L) MOD self.nmax

end
;
function photoobj_cache::read, run, rerun, camcol, field, id, readobj=readobj


    if(n_elements(run) eq 0 OR $
        n_elements(camcol) eq 0 OR $
        n_elements(rerun) eq 0 OR $
        n_elements(id) eq 0 OR $
        n_elements(field) eq 0) then begin
       message, 'Must specify RUN, CAMCOL, FIELD, ID, RERUN'
    endif
      
    fieldids= sdss_fieldid(run, camcol, field, rerun=rerun)
    isort= sort(fieldids)
    iuniq= uniq(fieldids[isort])
    istart=0L
    for i=0L, n_elements(iuniq)-1L do begin
        iend=iuniq[i]
        icurr= isort[istart:iend]

        tmp_fieldid= fieldids[icurr[0]]
        tmp_run=run[icurr[0]]
        tmp_rerun=rerun[icurr[0]]
        tmp_camcol=camcol[icurr[0]]
        tmp_field=field[icurr[0]]
        
        nc=0L
        if(n_tags((*self.cache)) gt 0) then begin
            ic= where((*self.cache).fieldid eq tmp_fieldid, nc)
            if(nc gt 1) then begin
                message, 'Inconsistency! more than one instance of field in cache!'
            endif
        endif
        if(nc eq 0) then begin
            self->add2cache, tmp_run, tmp_camcol, tmp_field, rerun=tmp_rerun, $
                                readobj=readobj
            ic= where((*self.cache).fieldid eq tmp_fieldid, nc)
            if(nc eq 0) then $
              message, 'Inconsistency! adding to cache failed!'
        endif

        cobj = (*self.cache)[ic].obj

        tmp_obj=(*cobj)[id[icurr]-1L]

        if(n_tags(obj) eq 0) then begin
            obj= replicate(tmp_obj[0], n_elements(run))
        endif

        obj[icurr]= tmp_obj
        istart=iend+1L

    endfor

    return, obj
end

pro photoobj_cache::cleanup
    heap_free, self.cache
end

pro photoobj_cache__define
    struct = {$
        photoobj_cache, $
        verbose:0, $
        nmax: 0L, $
        ip: 0L, $
        cache: ptr_new() $
    }
end
