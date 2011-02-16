function regauss_sweeps::init, type
    if n_elements(type) eq 0 then begin
        message,'usage: rg=obj_new("regauss_sweeps", sweep_type)'
    endif
    self.type = type
    self.max_rmag = 22.0
	return,1
end


function regauss_sweeps::process_struct, str, nkeep
	ntot=n_elements(str)
	keep = self->select(str, nkeep, anyresolve=anyresolve, mst=mst)
	if nkeep eq 0 then print,'No objects passed cuts'

	print,'Running on ',nkeep,'/',ntot,' objects',$
		f='(a,i0,a,i0,a)'

	rg=obj_new('regauss')
	res = rg->process_atlas_images(str, indices=keep)

    output = struct_addfields(res, {cmodelmag_dered:fltarr(5), modelmag_dered:fltarr(5)})
    output.modelmag_dered = mst.modelmag_dered
    output.cmodelmag_dered = mst.cmodelmag_dered

	return, output
end



; version for running on a camcol outside of the sweeps__define sort of
; framework.

function regauss_sweeps::process_sweep, run, camcol, rerun=rerun, $
        anyresolve=anyresolve, $
        fieldrange=fieldrange, $
        str=str, $
        nkeep=nkeep

	; NOTE:  In this version, not all objects from the sweep get 
	; returned in the structure, and even -1 could be returned.
	; see run-regauss above for the version that returns a row
	; for each object even if it doesn't get processed by the
	; shape code


    str=self->read_sweep(run,camcol,rerun=rerun, fieldrange=fieldrange)
	ntot=n_elements(str)
	keep = self->select(str, nkeep, anyresolve=anyresolve, mst=mst)
	if nkeep eq 0 then begin
		print,'No objects passed cuts'
		return, -1
	endif
	print,'Running on ',nkeep,'/',ntot,' objects',$
		f='(a,i0,a,i0,a)'

	str=str[keep]
    mst=mst[keep]

    rg = obj_new('regauss')
	res = rg->process_atlas_images(str)

    output = struct_addfields(res, {cmodelmag_dered:fltarr(5), modelmag_dered:fltarr(5)})
    output.modelmag_dered = mst.modelmag_dered
    output.cmodelmag_dered = mst.cmodelmag_dered
	return, output
		           
end

function regauss_sweeps::read_sweep, run, camcol, rerun=rerun, fieldrange=fieldrange
    if n_elements(rerun) eq 0 then rerun=sdss_rerun(run)
    ;str = sdss_read('calibObj.'+self.type,run,camcol,rerun=rerun)
    str = sweep_readobj(run, camcol, rerun=rerun, type=self.type, fieldrange=fieldrange)
    if self.type eq 'star' then begin
        ; convert psf_FWHM to m_rr_cc
        addstr = {m_rr_cc:fltarr(5)}
        newstr = struct_addtags(str, addstr)

        newstr.m_rr_cc = fwhm2mom(newstr.psf_fwhm)
        str=0
        str = temporary(newstr)
    endif
    return, str
end

; this procedural version writes out a file into ~/tmp
pro regauss_sweeps::process_sweep, run, camcol, rerun=rerun, $
		fpobjc=fpobjc, $
		fieldrange=fieldrange

	rg = self->process(run, camcol, rerun=rerun, $
		fpobjc=fpobjc, $
		fieldrange=fieldrange, $
		nkeep=nkeep)
	if nkeep eq 0 then begin
		return
	endif

	outfile = self->output_file(run,camcol, rerun=rerun, $
		fieldrange=fieldrange, $
		fpobjc=fpobjc)
	print,'writing to file: ',outfile
	mwrfits, rg, outfile, /create
		           
end


; we will later use the sweeps class to deal with this
function regauss_sweeps::output_file, run, camcol, rerun=rerun, $
		fieldrange=fieldrange, fpobjc=fpobjc

	if n_elements(rerun) eq 0 then message,'enter a rerun'
	dir=getenv('REGAUSS_DIR')
	file=string(f='("regauss-",a,"-",i06,"-",i0,"-",i0,".fits")',self.type,run,camcol,rerun)

	if keyword_set(fpobjc) then begin
		file=repstr(file,'regauss','regauss-fpobjc')
	endif

	if n_elements(fieldrange) eq 2 then begin
		frstr = string(f='(i04,"-",i04)', fieldrange[0], fieldrange[1])
		file=repstr(file,'.fits', '-'+frstr+'.fits')
	endif

	file=filepath(root=dir, file)
	return, file
end


function regauss_sweeps::field_id, run,rerun,camcol,field
	ten=ulong64(10)
	p1 = 0L
	p2 = 6L
	p3 = 11L
	p4 = 12L
	p5 = 15L

	field_id = ulong64(field)*ten^p2
	field_id = temporary(field_id) + ulong64(camcol)*ten^p3
	field_id = temporary(field_id) + ulong64(rerun)*ten^p4
	field_id = temporary(field_id) + ulong64(run)*ten^p5

	return, field_id
end

function regauss_sweeps::select, objs, nkeep, $
        allfield=allfield, anyresolve=anyresolve, $
        mst=mst
	; select a good subset. Use essentially the same criteria from photoz
	; input files

    nobj = n_elements(objs)
	pz=obj_new('photoz_uchicago')

	mst = pz->derived_struct(objs)

	mag_logic = mst.cmodelmag_dered[2] lt self.max_rmag
	flag_logic = pz->flag_logic(objs)

    if not keyword_set(anyresolve) then begin
        resolve_logic = pz->resolve_logic(objs)
    endif else begin
        resolve_logic = replicate(1L, nobj)
    endelse

    if not keyword_set(allfield) then begin
        goodfield_logic = self->goodfield_logic(objs)
    endif else begin
        goodfield_logic = replicate(1L, nobj)
    endelse
	wmag=where(mag_logic)
	wflag = where(flag_logic)
	wresolve = where(resolve_logic)
	wgoodfield = where(goodfield_logic)
	help,objs,wmag,wflag, wresolve, wgoodfield

	keep= where( $
		goodfield_logic $
		and mag_logic $	
		and flag_logic $
		and resolve_logic, nkeep)

	return, keep

end



pro regauss_sweeps::make_pbs_03, bycamcol=bycamcol, lbl=lbl
    proctype='regauss'
    procrun='03'

    if keyword_set(lbl) then begin
        idlutils_v='v5_4_19'
        photoop_v='v1_10_2'
    endif else begin
        idlutils_v=' -r ~/exports/idlutils-trunk'
        photoop_v=' -r ~/exports/photoop-trunk'
    endelse
	command = $
	  'rs=obj_new("regauss_sweeps","'+self.type+'") & res=rs->process_struct(objs)'

    sw = obj_new('sweeps', self.type)

    ; we are outputting a new place because of bad disks

    setups=$
     ['export SWEEP_REDUCE=~esheldon/oh/sweep-reduce',$
      'setup tree']
   
    ;where_string='str.whyflag_rg[2] eq 0 or str.whyflag_rg[3] eq 0'
    ; since psf never fails, this is just asking if we bothered to run at all
    ; should add a field that says if we ran at all
    where_string='str.whyflag_psf[2] eq 0'

    nper=10
    if not keyword_set(bycamcol) then begin
        sw->create_pbs, proctype, procrun, command, $
            shortname='rg', $
            nper=nper, $
            $
            idlutils_v=idlutils_v, $
            photoop_v=photoop_v, $
            $
            LBL=lbl, $
            /collate, $
            where_string=where_string, $
            extra_setups=setups
    endif else begin
        sw->create_pbs_bycamcol, proctype, procrun, command, $
            shortname='rg', $
            $
            idlutils_v=idlutils_v, $
            photoop_v=photoop_v, $
            $
            LBL=lbl, $
            extra_setups=setups

        ; just for gathering
        sw->create_pbs, proctype, procrun, command, $
            shortname='rg', $
            nper=nper, $
            $
            idlutils_v=idlutils_v, $
            photoop_v=photoop_v, $
            $
            LBL=lbl, $
            where_string=where_string, $
            extra_setups=setups, $
            /collate, $
            /dogather



    endelse

end


pro regauss_sweeps::make_pbs_02
    proctype='regauss'
    procrun='02'

	idlutils_v='v5_4_19'
	photoop_v='v1_10_2'
	command = $
	  'rs=obj_new("regauss_sweeps","'+self.type+'") & res=rs->process_struct(objs)'

    sw = obj_new('sweeps', self.type)

    walltime='48:00:00'
    ; we are outputting a new place because of bad disks

    setups=$
     ['export SWEEP_REDUCE=/clusterfs/riemann/raid008/bosswork/groups/boss/target/esheldon/sweep-reduce']
   
    where_string='str.whyflag_rg[2] eq 0 or str.whyflag_rg[3] eq 0'
    nper=10
    sw->create_pbs, proctype, procrun, command, $
        shortname='rg', $
        nper=nper, $
        $
        idlutils_v=idlutils_v, $
        photoop_v=photoop_v, $
        $
        /LBL, $
        /collate, $
        where_string=where_string, $
        walltime=walltime, $
        extra_setups=setups

end



pro regauss_sweeps::make_pbs_01, dogather=dogather
    proctype='regauss'
    procrun='01'

	idlutils_v='v5_4_19'
	photoop_v='v1_10_2'
	command = $
	  'rs=obj_new("regauss_sweeps","'+self.type+'") & res=rs->process_struct(objs)'

    sw = obj_new('sweeps', self.type)

    if not keyword_set(dogather) then begin
        sw->create_pbs_byrun, proctype, procrun, command, $
            shortname='rg', $
            idlutils_v=idlutils_v, $
            photoop_v=photoop_v, $
            $
            where_string=where_string, $
            /LBL, $
            walltime='48:00:00'
    endif else begin
		where_string='str.whyflag_rg[2] eq 0 or str.whyflag_rg[3] eq 0'
		nper=20
		sw->create_pbs, proctype, procrun, command, $
            shortname='rg', $
			nper=nper, $
			$
			idlutils_v=idlutils_v, $
			photoop_v=photoop_v, $
			$
            /LBL, $
			/collate, $
			where_string=where_string, $
			/dogather

    endelse

end

function regauss_sweeps::goodfield_logic, objs
	field_id = self->field_id(objs.run,objs.rerun,objs.camcol,objs.field)

    bad_field_ids = self->bad_field_ids(nbad)
	goodfield_logic = (field_id ne bad_field_ids[0])
    for i=1L, nbad-1 do begin
	    goodfield_logic = goodfield_logic and (field_id ne bad_field_ids[i])
    endfor
    return, goodfield_logic
end

function regauss_sweeps::bad_field_ids, nbad
	bad_field_id1 =self->field_id(1462,301,3,209)
	bad_field_id2 =self->field_id(6794,301,1,86)
    bad_field_ids = [bad_field_id1, bad_field_id2]
    nbad=n_elements(bad_field_ids)
    return, bad_field_ids
end

pro regauss_sweeps::test_bad_fields, run

    if run eq 1462 then begin
        t=sweep_readobj(1462,3,rerun=301,type='gal')
        w=where(t.field eq 209)
    endif else begin
        t=sweep_readobj(6794,1,rerun=301,type='gal')
        w=where(t.field eq 86)
    endelse
    t=t[w]

    w=self->select(t, nkeep,/allfield)
    t=t[w]
    
    for i=0L, nkeep-1 do begin
        print,'id: ',t[i].id
        read_atlas, t[i], atlas_struct=as
    endfor
end

pro regauss_sweeps__define

  struct = {$
             regauss_sweeps, $
			 type: '', $
             max_rmag: 0.0 $
           }

end 
