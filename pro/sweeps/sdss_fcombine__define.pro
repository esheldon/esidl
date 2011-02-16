;
; Notes:
;   I do the "straight" means for nok gt 0 instead of nok gt 1
;
;   I use defaults: -9999 fluxes, 9999 var chi2, 0 for ivar
;   I check calib_status on   filter-by-filter basis.
;
; SDSSIDL dependencies:
;   sdss_photoid
;   combine_ptrlist
;

function sdss_fcombine::init, outdir=outdir, noclobber=noclobber, nmax=nmax
    ; outdir is the basedir for output.  Default is same
    ; as input: BOSS_PHOTOOBJ

    ; using lower than the default cache size
    if n_elements(nmax) eq 0 then nmax=100
    cache_init = self->photoobj_cache::init(nmax=nmax, /verbose)
    if cache_init ne 1 then return, -1

    self->set_outdir_base, outdir=outdir

    self.noclobber = keyword_set(noclobber)

    self.nfilter=5
    self.filters = ['u','g','r','i','z']
    return, 1
end


pro sdss_fcombine::process_camcol, run, rerun, camcol
    rlist = sdss_runlist(run, rerun=rerun) 
    nf = rlist.endfield - rlist.startfield+1
    countarr = lonarr(nf)
    count = 0L
    splog,format='("Run: ",i0," Rerun: ",i0," Camcol: ",i0,"  Nfields: ",i0)',$
        run,rerun,camcol,nf

    for field=rlist.startfield, rlist.endfield do begin
        self->process_field, run, rerun, camcol, field, count=tcount
        count += tcount

        countarr[field-rlist.startfield] = tcount
    endfor

    self->sweep_camcol, run, rerun, camcol
end

pro sdss_fcombine::sweep_camcol, run, rerun, camcol

    rlist = sdss_runlist(run, rerun=rerun) 
    nf = rlist.endfield - rlist.startfield+1
    countarr = lonarr(nf)
    count = 0L
    splog,'Sweeping Camcol'
    splog,format='("    Run: ",i0," Rerun: ",i0," Camcol: ",i0,"  Nfields: ",i0)',$
        run,rerun,camcol,nf

    for field=rlist.startfield, rlist.endfield do begin
        f = self->outfile(run, rerun, camcol, field)
        tcount = self->count_from_header(f)
        count += tcount
        countarr[field-rlist.startfield] = tcount
    endfor

    splog,'    Total objects: ',count
    if count gt 0 then begin
        allcol = self->photoobj_combine_struct(count)

        beg=0L
        for field=rlist.startfield, rlist.endfield do begin
            i=field-rlist.startfield
            if countarr[i] gt 0 then begin
                fname = self->outfile(run, rerun, camcol, field)
                t = mrdfits(fname,1,/silent)
                allcol[beg:beg+countarr[i]-1] = t
                beg += countarr[i]
            endif
        endfor

        self->write_sweep_output, allcol, run, rerun, camcol
        allcol=0

    endif

end


pro sdss_fcombine::process_field, run, rerun, camcol, field, count=count, epochs_count=epochs_count
    if self.noclobber then begin
        outf = self->outfile(run, rerun, camcol, field)
        epochsf = self->epochs_outfile(run, rerun, camcol, field)
        if file_test(outf) and file_test(epochsf) then begin
            splog,'Not clobbering existing files:'
            splog,'    ',outf,format='(a,a)'
            splog,'    ',epochsf,format='(a,a)'
            count = self->count_from_header(outf)
            epochs_count = self->count_from_header(epochsf)
            return
        endif
    endif
    res = self->fcombine_field(run, rerun, camcol, field, epochs=epochs)
    if n_tags(res) eq 0 then begin
        count=0
        epochs_count=0
    endif else begin
        count = n_elements(res)
        epochs_count = n_elements(res)
    endelse
    self->write_output, res, run, rerun, camcol, field
    self->write_epochs_output, epochs, run, rerun, camcol, field
end


function sdss_fcombine::fcombine_field, run, rerun, camcol, field, epochs=epochs

    ; read in the photoobj field data and process the survey_primary
    ; objects
    if (n_elements(run) ne 1) $
            or (n_elements(rerun) ne 1) $
            or (n_elements(camcol) ne 1) $
            or (n_elements(field) ne 1) then begin
        message,'fcombine, run, rerun, camcol, field, obj=obj, epobj=epobj'
    endif

    fieldobj = self->read_field(run, rerun, camcol, field)
    if n_tags(fieldobj) eq 0 then return, 0

    pflag= sdss_flagval('resolve_status', 'survey_primary')
    wprimary= where((fieldobj.resolve_status and pflag) gt 0, nprimary)
    if (nprimary eq 0) then begin
        splog,'No primary objects found' 
        epochs = 0
        return, 0
    endif

    fieldobj = fieldobj[wprimary]

    ;; set up outputs, set up initial values
    obj = self->photoobj_combine_struct(nprimary)
    struct_assign, fieldobj, obj, /nozero

    self->process_singlets, fieldobj, obj, singlet_epochs_ptrs
    self->process_multi, fieldobj, obj, multi_epochs_ptrs

    epochs = self->combine_epochs(singlet_epochs_ptrs, multi_epochs_ptrs)

    if n_tags(obj) ne 0 and n_tags(epochs) ne 0 then begin
        splog,'Verifying....'
        self->check_use_count, obj, epochs
    endif

    return, obj

end

pro sdss_fcombine::process_singlets, fieldobj, obj, epochsp

    if not arg_present(epochsp) then message,'epochsp must be an argument'
    self->clear_epochs, epochsp

    wsingle= where(obj.ndetect eq 1, nsingle)
    for i=0L, nsingle-1L do begin
        cobj= obj[wsingle[i]]

        currobj= self->photoobj_cache::read(cobj.run, $
                                            cobj.rerun, $
                                            cobj.camcol, $
                                            cobj.field, $
                                            cobj.id)
        epochs = self->epoch_struct()


        struct_assign, currobj, epochs, /nozero

        epochs.primary_photoid = sdss_photoid(cobj)

        ; should be smart enough to just copy the primary
        tmp_obj = self->fcombine(epochs)

        struct_assign, fieldobj[wsingle[i]], tmp_obj, /nozero
        obj[wsingle[i]]=tmp_obj

        add_arrval, ptr_new(epochs, /no_copy), epochsp
    endfor
end

pro sdss_fcombine::process_multi, fieldobj, obj, epochsp
    common fcombine_process_multi_block, flist

    if not arg_present(epochsp) then message,'epochsp must be an argument'
    self->clear_epochs, epochsp

    if n_elements(flist) eq 0 then begin
        splog,'Cacheing window flist'
        window_read, flist=flist
    endif
    wmulti= where(obj.ndetect gt 1, nmulti)
    if nmulti eq 0 then return

    dir = getenv('PHOTO_RESOLVE')
    if dir eq '' then message,'PHOTO_RESOLVE is not set'
    indxfile= filepath(root=dir,'thingIndx.fits')
    listfile= filepath(root=dir,'thingList.fits')

    rp = sdss_flagval('resolve_status', 'run_primary')
    for i=0L, nmulti-1L do begin
        cobj= obj[wmulti[i]]

        ;; - read in thingIndx row 
        indx= mrdfits(indxfile, 1, row=cobj.thing_id, /silent)
        
        ;; - read in thingList section
        range= indx.istart+[0L, indx.ndetect-1L]
        list=mrdfits(listfile, 1, range=range, /silent)
        
        ;; - generate list of object ids
        curr_run=flist[list.ifield].run
        curr_rerun=flist[list.ifield].rerun
        curr_camcol=flist[list.ifield].camcol
        curr_field=flist[list.ifield].field
        curr_id=list.id

        ;; - but only keep objects in same rerun as primary
        ikeep= where(curr_rerun eq cobj.rerun, nkeep)
        if(nkeep eq 0) then $
           message, 'Inconsistency! the primary object should be in its own rerun at least!'
        curr_run=curr_run[ikeep]
        curr_rerun=curr_rerun[ikeep]
        curr_camcol=curr_camcol[ikeep]
        curr_field=curr_field[ikeep]
        curr_id=curr_id[ikeep]
        
        ;; - now read in all the appropriate info
        currobj= self->photoobj_cache::read(curr_run, $
                                            curr_rerun, $
                                            curr_camcol, $
                                            curr_field, $
                                            curr_id)
        
        ; fix the one object with bad resolve_status
        self->fix_bad_object, currobj

        ; only send run_primary objects
        w=where( (currobj.resolve_status and rp) ne 0, nw)
        if nw gt 0 then begin
            currobj = currobj[w]

            ; create the epochs struct
            epochs = self->epoch_struct(nw)
            struct_assign, currobj, epochs, /nozero

            epochs.primary_photoid = sdss_photoid(cobj)

            tmp_obj = self->fcombine(epochs)

            struct_assign, fieldobj[wmulti[i]], tmp_obj, /nozero
            obj[wmulti[i]]=tmp_obj

            add_arrval, ptr_new(epochs, /no_copy), epochsp
        endif
    endfor


end

function sdss_fcombine::fcombine, epochs
    ; only send run_primary objects to this function.  At least one of
    ; these should also be survey_primary in the case where none of the
    ; observations are photometric
    ;

    self->fcombine_check, epochs

    outobj0 = self->photoobj_combine_struct()

    ftypes = self->types2sum()

    for iftype=0, n_elements(ftypes)-1 do begin
        ftype = ftypes[iftype]
        self->fcombine_byfilter, epochs, outobj0, ftype
    endfor

    nobj=n_elements(epochs)
    ftypes_clean = self->types2sum(/clean, objc_type=objc_type)
    for iftype=0, n_elements(ftypes_clean)-1 do begin
        ftype = ftypes_clean[iftype]
        w=where(epochs.objc_type eq objc_type[iftype], nw)
        if nw gt 0 then begin
            self->fcombine_byfilter, epochs, outobj0, ftype, index=w, /clean
        endif
    endfor

    return, outobj0
end


pro sdss_fcombine::fcombine_byfilter, epochs, outobj0, ftype, index=index, clean=clean

    ; to be run when there are objects with good calibrations
    ; note for singlets you should just use fcombine_copy_primary

    cflag= sdss_flagval('calib_status', 'photometric')
    for ifilter=0, self.nfilter-1 do begin
        if n_elements(index) eq 0 then begin
            wcalib = where((epochs.calib_status[ifilter] and cflag) ne 0, ncalib)
        endif else begin
            wcalib = where((epochs[index].calib_status[ifilter] and cflag) ne 0, ncalib)
            if ncalib ne 0 then wcalib = index[wcalib]
        endelse
        if ncalib ne 0 then begin
            self->fcombine_filter, epochs, outobj0, ftype, ifilter, index=wcalib, clean=clean
        endif else begin
            self->fcombine_copy_primary_filter, epochs, outobj0, ftype, ifilter, index=index, clean=clean
        endelse

    endfor

end
pro sdss_fcombine::fcombine_filter, epochs, outobj0, ftype, ifilter, index=index, clean=clean

    self->get_tags, ftype, $
        mean_tag, ivar_tag, var_tag, chi2_tag, nuse_tag, used_tag, $
        mjd_maxdiff_tag, mjd_var_tag, $
        clean=clean
    ivar = self->struct_get(epochs, ftype+'flux_ivar', index=index, array_index=ifilter)
    w = where(ivar gt 0, nw)

    if nw gt 0 then begin
        ivar = ivar[w]

        ; now get index back into epochs
        if n_elements(index) ne 0 then w = index[w]

        self->struct_set, outobj0, nuse_tag, nw,  array_index=ifilter

        flux = self->struct_get(epochs, ftype+'flux', array_index=ifilter, index=w)
        weights = self->struct_get(epochs, ftype+'flux_ivar', array_index=ifilter, index=w)
        mjd = self->struct_get(epochs, 'mjd', index=w)

        wsum = total(weights)

        mn = total(weights*flux)/wsum

        self->struct_set, outobj0, mean_tag, mn,   array_index=ifilter
        self->struct_set, outobj0, ivar_tag, wsum, array_index=ifilter

        ; leave these as -9999 if only one object
        if nw gt 1 then begin
            flux_var = total((flux-mn)^2)/(nw-1)
            flux_chi2 = total((flux-mn)^2*weights)

            self->struct_set, outobj0, var_tag,  flux_var,  array_index=ifilter
            self->struct_set, outobj0, chi2_tag, flux_chi2, array_index=ifilter

            mjd_maxdiff = max(mjd)-min(mjd)
            mjd_var = var(mjd)

            self->struct_set, outobj0, mjd_maxdiff_tag, mjd_maxdiff, array_index=ifilter
            self->struct_set, outobj0, mjd_var_tag, mjd_var, array_index=ifilter
        endif

        ; finally, let it be known these objects were used
        self->struct_set, epochs, used_tag, 1, index=w, array_index=ifilter
        self->struct_set, epochs, 'fcombine_used', 1, index=w

    endif
end
pro sdss_fcombine::fcombine_copy_primary_filter, epochs, outobj0, ftype, ifilter, index=index, clean=clean
    ; just copy the data from the survey primary object
    self->get_tags, ftype, $
        mean_tag, ivar_tag, var_tag, chi2_tag, nuse_tag, used_tag, $
        mjd_maxdiff_tag, mjd_var_tag, $
        clean=clean

    ivar = self->struct_get(epochs, ftype+'flux_ivar', index=index, array_index=ifilter)
    rstatus = self->struct_get(epochs, 'resolve_status', index=index)

    pflag= sdss_flagval('resolve_status', 'survey_primary')

    w=where((rstatus and pflag) ne 0 and ivar gt 0, ngood)
    if ngood ne 0 then begin
        if ngood ne 1 then message,'Expected only one primary'

        ivar = ivar[w[0]]

        ; get index back into epochs
        if n_elements(index) ne 0 then w=index[w[0]]

        self->struct_set, outobj0, nuse_tag, 1,  array_index=ifilter

        ; var and chi2 are left at -9999
        flux = self->struct_get(epochs, ftype+'flux', array_index=ifilter, index=w)
        self->struct_set, outobj0, mean_tag, flux, array_index=ifilter
        self->struct_set, outobj0, ivar_tag, ivar, array_index=ifilter

        ; finally, let it be known this object was used
        self->struct_set, epochs, used_tag, 1, index=w, array_index=ifilter
        self->struct_set, epochs, 'fcombine_used', 1, index=w

    endif

end

pro sdss_fcombine::get_tags, ftype, $
        mean_tag, ivar_tag, var_tag, chi2_tag, nuse_tag, used_tag, $
        mjd_maxdiff_tag, mjd_var_tag, $
        clean=clean
    if keyword_set(clean) then begin
        mean_tag = ftype+'flux_clean'
        ivar_tag = ftype+'flux_clean_ivar'
        var_tag =  ftype+'flux_clean_var'
        chi2_tag = ftype+'flux_clean_chi2'
        nuse_tag = ftype+'_clean_nuse'
        used_tag = ftype+'_clean_used'

        mjd_maxdiff_tag = ftype+'_clean_mjd_maxdiff'
        mjd_var_tag     = ftype+'_clean_mjd_var'
    endif else begin
        mean_tag = ftype+'flux_mean'
        ivar_tag = ftype+'flux_mean_ivar'
        var_tag =  ftype+'flux_var'
        chi2_tag = ftype+'flux_chi2'
        nuse_tag = ftype+'_nuse'
        used_tag = ftype+'_used'

        mjd_maxdiff_tag = ftype+'_mjd_maxdiff'
        mjd_var_tag     = ftype+'_mjd_var'
    endelse
end

pro sdss_fcombine::fcombine_check, obj
    nobj = n_elements(obj)
    sp = sdss_flagval('resolve_status', 'survey_primary')
    rp = sdss_flagval('resolve_status', 'run_primary')

    w=where((obj.resolve_status and sp) ne 0, nw)
    if nw eq 0 then message,'At least one object should be survey_primary'

    w=where((obj.resolve_status and rp) ne 0, nw)
    if nw ne nobj then message,'Only send run_primary objects to fcombine'

end

pro sdss_fcombine::fix_bad_object, obj
    w=where(obj.run eq 259 and obj.camcol eq 3 and obj.field eq 567 and obj.id eq 1594, nw)
    if nw ne 0 then begin
        splog,'Fixing bad object: 259 3 567 1594'
        ; fix error where run_primary is not set
        rp = sdss_flagval('resolve_status', 'run_primary')
        if (obj[w].resolve_status and rp) eq 0 then begin
            obj[w].resolve_status += rp
        endif
    endif
end

function sdss_fcombine::struct_get, struct, tagname, index=index, array_index=array_index
    tag_index=where(tag_names(struct) eq strupcase(tagname), nw)
    if nw eq 0 then message,'no such tag: '+string(tagname)

    if n_elements(array_index) ne 0 then begin
        if n_elements(index) ne 0 then begin
            data = struct[index].(tag_index)[array_index]
        endif else begin
            data = struct.(tag_index)[array_index]
        endelse
    endif else begin
        if n_elements(index) ne 0 then begin
            data = struct[index].(tag_index)
        endif else begin
            data = struct.(tag_index)
        endelse
    endelse

    return, data
end
pro sdss_fcombine::struct_set, struct, tagname, data, index=index, array_index=array_index
    tag_index=where(tag_names(struct) eq strupcase(tagname), nw)
    if nw eq 0 then message,'no such tag: '+string(tagname)

    if n_elements(array_index) ne 0 then begin
        if n_elements(index) ne 0 then begin
            struct[index].(tag_index)[array_index] = data
        endif else begin
            struct.(tag_index)[array_index] = data
        endelse
    endif else begin
        if n_elements(index) ne 0 then begin
            struct[index].(tag_index) = data
        endif else begin
            struct.(tag_index) = data
        endelse
    endelse
end


function sdss_fcombine::types2sum, clean=clean, objc_type=objc_type
    if keyword_set(clean) then begin
        objc_type = [6, 3, 3]
        return, [ 'PSF', 'CMODEL', 'MODEL']
    endif else begin
        return, [ 'FIBER', 'FIBER2', 'PSF', 'MODEL', 'CMODEL', 'PETRO', $
                  'DEV', 'EXP']
    endelse
end

function sdss_fcombine::photoobj_combine_struct, n
    common photoobj_combine_struct_block, cstruct

    if n_elements(cst) eq 0 then begin
   
        cstruct= {run:0l, $
                  camcol:0, $
                  field:0l, $
                  id:0l, $
                  rerun:' ',$
                  ra:0.d, $
                  dec:0.d, $
                  thing_id:0l, $
                  resolve_status:0l, $
                  objc_type:0l, $
                  ifield:0l, $
                  ndetect:0l, $
                  nobserve:0l $
                  }

        ftypes= self->types2sum()

        ; use defaults designed to minimize accidental improper use
        arrval = replicate(-9999.0, 5)  ; unlikely to use fluxes < 0
        varval = replicate(-9999.0, 5)  ; var is by def >= 0
        ivarval = fltarr(5)             ; ivar=0 is large noise, zero weight
        mjd_diffval = replicate(-9999L, 5) ; These are both by def >= 0
        mjd_varval = replicate(-9999.0, 5)
        nuseval = intarr(5)

        for i=0L, n_elements(ftypes)-1L do begin
            cstruct= create_struct(cstruct, ftypes[i]+'flux_mean', arrval)
            cstruct= create_struct(cstruct, ftypes[i]+'flux_mean_ivar', ivarval)
            cstruct= create_struct(cstruct, ftypes[i]+'flux_var', varval)
            cstruct= create_struct(cstruct, ftypes[i]+'flux_chi2', varval)
            cstruct= create_struct(cstruct, ftypes[i]+'_nuse', nuseval)
            cstruct= create_struct(cstruct, ftypes[i]+'_mjd_maxdiff', mjd_diffval)
            cstruct= create_struct(cstruct, ftypes[i]+'_mjd_var', mjd_varval)
        endfor
   
        ; clean just means only stars are used for psf, only galaxies are used
        ; for cmodel and model.  probably not the best nomenclature.
        fcleantypes= self->types2sum(/clean)

        for i=0L, n_elements(fcleantypes)-1L do begin
            cstruct= create_struct(cstruct, fcleantypes[i]+'flux_clean', arrval)
            cstruct= create_struct(cstruct, fcleantypes[i]+'flux_clean_ivar', ivarval)
            cstruct= create_struct(cstruct, fcleantypes[i]+'flux_clean_var', varval)
            cstruct= create_struct(cstruct, fcleantypes[i]+'flux_clean_chi2', varval)
            cstruct= create_struct(cstruct, fcleantypes[i]+'_clean_nuse', nuseval)
            cstruct= create_struct(cstruct, fcleantypes[i]+'_clean_mjd_maxdiff', mjd_diffval)
            cstruct= create_struct(cstruct, fcleantypes[i]+'_clean_mjd_var', mjd_varval)
        endfor
    endif

    if n_elements(n) ne 0 then begin
        cstruct = replicate(cstruct, n)
    endif
    return, cstruct
 
end

function sdss_fcombine::epoch_struct, n
    st={primary_photoid: 0LL, $ ; id of primary object we started with
        run:0L,$
        rerun: ' ', $
        camcol: 0, $
        field: 0, $
        id: 0L, $
        mjd: 0L, $
        ra: 0d, $
        dec: 0d, $
        thing_id: 0L, $
        resolve_status:0L, $
        objc_type: 0L, $
        calib_status: lonarr(5),$
        fcombine_used: 0b}

    ftypes= self->types2sum()
    fcleantypes= self->types2sum(/clean)

    arrval = replicate(-9999.0, 5)
    useval = bytarr(5)
    ivarval = fltarr(5)

    for i=0L, n_elements(ftypes)-1L do begin
        ftype = ftypes[i]
        st= create_struct(st, ftype+'flux', arrval)
        st= create_struct(st, ftype+'flux_ivar', ivarval)
        st= create_struct(st, ftype+'_used', useval)
        if in(fcleantypes, ftype) then begin
            st= create_struct(st, ftype+'_clean_used', useval)
        endif
    endfor
   

    if n_elements(n) ne 0 then begin
        st = replicate(st, n)
    endif
    return, st

end

function sdss_fcombine::read_field, run, rerun, camcol, field
    ;; read in data for this field
    splog,format='("Reading field info: ",i0," ",i0," ",i0)',run,camcol,field
    reobj= mrdfits(sdss_name('reObjGlobal', run, camcol, field, rerun=rerun),1,/silent)
    nobj= n_elements(reobj)

    fieldobj= self->photoobj_cache::read(replicate(run,nobj), $
                                         replicate(rerun, nobj), $
                                         replicate(camcol,nobj), $
                                         replicate(field,nobj), $
                                         lindgen(nobj)+1L)

    return, fieldobj
end

function sdss_fcombine::count_from_header, filename
    ; use errmsg to prevent things like 
    ;  % HEADFITS: ERROR - Extension past EOF
    ; going to the stderr

    hdr= headfits(filename, ext=1, errmsg=errmsg)
    if errmsg ne '' and errmsg ne 'Extension past EOF' then begin
        message,errmsg
    endif
    if size(hdr,/tname) eq 'STRING' then begin
        count = long(sxpar(hdr, 'NAXIS2'))
    endif else begin
        count=0L
    endelse
    return, count
end


pro sdss_fcombine::set_outdir_base, outdir=outdir
    if n_elements(outdir) ne 0 then begin
        self.outdir_base=outdir
    endif else begin
        dir=getenv('BOSS_PHOTOOBJ')
        if dir eq '' then message,'BOSS_PHOTOOBJ is not set'
        self.outdir_base = dir
    endelse
end

function sdss_fcombine::outfile, run, rerun, camcol, field
    outfile= sdss_name('photoCombine', run, camcol, field, rerun=rerun)
    outfile = file_basename(outfile)
    outfile = filepath(root=self->outdir(run,rerun,camcol), outfile)
    return, outfile
end
function sdss_fcombine::epochs_outfile, run, rerun, camcol, field
    outfile= sdss_name('photoEpochs', run, camcol, field, rerun=rerun)
    outfile = file_basename(outfile)
    outfile = filepath(root=self->outdir(run,rerun,camcol), outfile)
    return, outfile
end

function sdss_fcombine::sweep_outfile, run, rerun, camcol
    dir=self->sweep_outdir(run, rerun)

    fname='photoCombineCamcol-'+string(run, f='(i6.6)')+'-'+ $
      strtrim(string(camcol),2)+'.fits'
    return,filepath(root=dir, fname)
end


function sdss_fcombine::outdir, run, rerun, camcol
    sdir = self->sweep_outdir(run,rerun)
    cstr = string(camcol,f='(i0)')
    dir = filepath(root=sdir, cstr)
    return, dir
end
function sdss_fcombine::sweep_outdir, run, rerun
    rstr = string(run,f='(i0)')
    rrstr = string(rerun,f='(i0)')
    dir=filepath(root=self.outdir_base, sub=['Combine',rrstr], rstr)
    return, dir
end

pro sdss_fcombine::write_output, res, run, rerun, camcol, field
    dir = self->outdir(run, rerun, camcol)
    file_mkdir, dir
    outf = self->outfile(run, rerun, camcol, field)
    if self.noclobber and file_test(outf) then begin
        splog,'Not overwriting file:',outf,format='(a,a)'
        return
    endif
    splog,'Writing field output: ',outf,format='(a,a)'

    phdr= photoobj_hdr()
    mwrfits, res, outf, phdr, /create
end
pro sdss_fcombine::write_epochs_output, res, run, rerun, camcol, field
    dir = self->outdir(run, rerun, camcol)
    file_mkdir, dir
    outf = self->epochs_outfile(run, rerun, camcol, field)
    if self.noclobber and file_test(outf) then begin
        splog,'Not overwriting file:',outf,format='(a,a)'
        return
    endif
    splog,'Writing field epochs output: ',outf,format='(a,a)'

    phdr= photoobj_hdr()
    mwrfits, res, outf, phdr, /create
end

pro sdss_fcombine::write_sweep_output, res, run, rerun, camcol
    dir = self->sweep_outdir(run, rerun)
    file_mkdir, dir
    outf = self->sweep_outfile(run, rerun, camcol)
    splog,'Writing camcol sweep output: ',outf,format='(a,a)'
    mwrfits, res, outf, /create
end




function sdss_fcombine::read_output, run, rerun, camcol, field, fname=f, silent=silent
    f = self->outfile(run, rerun, camcol, field)
    if not keyword_set(silent) then begin
        splog,'Reading field output: ',f,format='(a,a)'
    endif
    return, mrdfits(f, 1, silent=silent)
end
function sdss_fcombine::read_epochs_output, run, rerun, camcol, field, silent=silent
    f = self->epochs_outfile(run, rerun, camcol, field)
    if not keyword_set(silent) then begin
        splog,'Reading epochs field output: ',f,format='(a,a)'
    endif
    return, mrdfits(f, 1, silent=silent)
end
function sdss_fcombine::read_sweep_output, run, rerun, camcol, field, silent=silent
    f = self->sweep_outfile(run, rerun, camcol, field)
    if not keyword_set(silent) then begin
        splog,'Reading sweep output: ',f,format='(a,a)'
    endif
    return, mrdfits(f, 1, silent=silent)
end




function sdss_fcombine::outdir_base
    return, self.outdir_base
end


pro sdss_fcombine::clear_epochs, epochs
    if n_elements(epochs) ne 0 then begin
        if size(epochs,/tname) eq 'POINTER' then begin
            ptr_free, epochs
        endif
        tmp = temporary(epochs)
    endif
end

function sdss_fcombine::combine_epochs, singlet_epochs_ptrs, multi_epochs_ptrs
    ns = n_elements(singlet_epochs_ptrs)
    nm = n_elements(multi_epochs_ptrs)

    if ns eq 0 and nm eq 0 then begin
        return, 0
    endif else if ns ne 0 and nm ne 0 then begin
        ;sepochs = combine_ptrlist(singlet_epochs_ptrs)
        ;mepochs = combine_ptrlist(multi_epochs_ptrs)
        ;epochs = struct_concat(sepochs, mepochs)
        ;ptrlist = [singlet_epochs_ptrs, multi_epochs_ptrs]
        ;epochs = combine_ptrlist(ptrlist)
        epochs = combine_ptrlist( [singlet_epochs_ptrs, multi_epochs_ptrs] )
    endif else if ns ne 0 then begin
        epochs = combine_ptrlist(singlet_epochs_ptrs)
    endif else begin
        epochs = combine_ptrlist(multi_epochs_ptrs)
    endelse

    self->clear_epochs, singlet_epochs_ptrs
    self->clear_epochs, multi_epochs_ptrs

    wuse = where(epochs.fcombine_used eq 1, nuse)
    if nuse eq 0 then return, 0

    epochs = epochs[wuse]

    return, epochs
end

pro sdss_fcombine::test_camcol, run, rerun, camcol, verbose=verbose
    rlist = sdss_runlist(run, rerun=rerun) 
    nf = rlist.endfield - rlist.startfield+1
    countarr = lonarr(nf)
    count = 0L
    splog,format='("Run: ",i0," Rerun: ",i0," Camcol: ",i0,"  Nfields: ",i0)',$
        run,rerun,camcol,nf

    for field=rlist.startfield, rlist.endfield do begin
        print,format='("run: ",i0," rerun: ",i0," camcol: ",i0," field: ",i0)',$
            run,rerun,camcol,field
        self->test_field, run, rerun, camcol, field, verbose=verbose
    endfor
 
end
pro sdss_fcombine::test_field, run, rerun, camcol, field, verbose=verbose
    ; make sure the epochs file and combine file agree
    ; on things like the number of objects used for averages.
    if keyword_set(verbose) then silent=0 else silent=1
    c = self->read_output(run, rerun, camcol, field, fname=fname, silent=silent)
    e = self->read_epochs_output(run, rerun, camcol, field, silent=silent)
    
    if n_tags(c) eq 0 then begin
        splog,'No combine data in file: '+fname
    endif else begin
        self->check_use_count, c, e, verbose=verbose
    endelse
end
pro sdss_fcombine::check_use_count, c, e, verbose=verbose
    ; c = combined 
    ; e = epochs
    nobj = n_elements(c)

    ftypes = self->types2sum()
    ftypes_clean = self->types2sum(/clean)

    ids = sdss_photoid(c)

    for iobj=0L, nobj-1 do begin

        if keyword_set(verbose) then begin
            print,format='("    run: ",i0," rerun: ",i0," camcol: ",i0," field: ",i0," id: ",i0)',$
                c[iobj].run,c[iobj].rerun,c[iobj].camcol,c[iobj].field,c[iobj].id
        endif
        w=where(e.primary_photoid eq ids[iobj], nw)
        noneuse = 0
        if nw eq 0 then begin
            ; in this case there better not be nuse > 0 for any flux!
            noneuse = 1
        endif

        ; loop over all types
        for i=0L, n_elements(ftypes)-1 do begin
            ftype = ftypes[i]

            ; can pick out specific ones used for clean, only a subset
            for clean=0,1 do begin
                if (clean eq 0) or (clean eq 1 and in(ftypes_clean, ftype)) then begin
                    self->get_tags, ftype, mean_tag, ivar_tag, var_tag, chi2_tag, nuse_tag, used_tag, clean=clean

                    for ifilter=0,4 do begin
                        nuse = self->struct_get(c, nuse_tag, array_index=ifilter, index=iobj)
                        err=''
                        if noneuse then begin
                            if nuse ne 0 then err='        didnt expect any fluxes to be used'
                        endif else begin
                            used = self->struct_get(e, used_tag, array_index=ifilter, index=w)
                            used_total = total(used,/int)

                            if nuse ne used_total then begin
                                err = '        Found nuse: '+strn(nuse)+' but used: '+ntostr(used_total)
                            endif
                        endelse

                        if err ne '' then begin
                                print,'        Photo id: ',ids[iobj]
                                print,'        Clean:    ',clean
                                print,'        Filter:   ',ifilter
                                print,'        nuse:     ',nuse
                                print,err
                                message,'Error, halting'
                        endif

                    endfor ;filter
                endif
            endfor ; clean
        endfor ; ftypes
    endfor
end

pro sdss_fcombine::compare_fruns, frun1, frun2, run, rerun, camcol
    d1=filepath('/clusterfs/riemann/raid008/bosswork/groups/boss/target/fcombine',run1)
    d2=filepath('/clusterfs/riemann/raid008/bosswork/groups/boss/target/fcombine',run2)
end

pro sdss_fcombine__define

	struct = {                           $
		sdss_fcombine,                   $
        noclobber:        0,             $
        outdir_base:      '',            $
        nfilter:          0,             $
        filters:          strarr(5),     $
        inherits          photoobj_cache $
	}
end
