;+
; CLASS NAME:
;	sweeps
;
; PURPOSE:
;	Generic parallel processing code for sweeps.  This class
;	begins it's life as a gutting of the target selection code,
;	so may have things left over from that.
;
; CATEGORY:
;	SDSS III/BOSS specific
;
; CALLING SEQUENCE:
;	bt = obj_new('sweeps',pars=)
; METHODS:
;	run IDL> methods, 'sweeps' to see a list of methods.
;
; MODIFICATION HISTORY:
;	Created Erin Sheldon, BNL 2009-??
;
;-
function sweeps::init, type, pars=pars, _extra=_extra

    if n_elements(type) eq 0 then begin
        message,'usage: sw=obj_new("sweeps", sweep_type)'
    endif
    self.type = type

	self->set_default_pars
	self->copy_extra_pars, pars=pars, _extra=_extra
	return, 1
end


function sweeps::pars
	return, self.pars
end

pro sweeps::process_partial_runs, proctype, procrun, command, jobnum, $
		nper=nper, $
		run_index=run_index, $
		addtags=addtags, $
		fpobjc=fpobjc, $
		pars=pars, $
		_extra=_extra

	if n_elements(proctype) eq 0 or n_elements(jobnum) eq 0 then begin
		on_error, 2
		message,'bt->process_partial_runs, proctype, procrun, jobnum, nper=2, run_index=, addtags=, /fpobjc, _extra='
	endif

	if n_elements(nper) eq 0 then begin
		nper=2
	endif
	self->split_runlist, runptrs, rerunptrs, nper=nper
	runs = *runptrs[jobnum]
	reruns = *rerunptrs[jobnum]
	ptr_free, runptrs, rerunptrs

	if n_elements(run_index) ne 0 then begin
		runs=runs[run_index]
		reruns=reruns[run_index]
	endif

	splog,'Using runs: ',runs

	self->process_runs, proctype, procrun, command, $
		runs=runs, fpobjc=fpobjc, addtags=addtags, $
		pars=pars, $
		_extra=_extra

end



;docstart::sweeps::process_runs
;
; NAME:
;	sweeps::process_runs	
;
; PURPOSE:
;	Process a list of runs through a code.  Can process either the
;	datasweeps or fpObjc files.
;
; CALLING SEQUENCE:
;	bt=obj_new('sweeps')
;	bt->process_runs, proctype, procrun, $
;		output_dir=, runs=, logfile=, /fpobjc, _extra=_extra
;
; INPUTS:
;	proctype: 'lrg','qso','std'
;	procrun: Either a string or a number.  Numbers will be converted to
;		a 6-zero padded integer string, e.g. 001732, when including in a
;		file name.
;
; OPTIONAL INPUTS:
;	output_dir: over-ride the default location $SWEEP_REDUCE/procrun	
;	runs: specify a set of runs.  These must be a subset of those returned
;		by a run of bt->runlist
;	addtags=: A set of tags to add to the sweeps from the fpObjc files.
;	/fpobjc: Process fpObjc files instead of sweep files.
;	pars=: A structure with extra parameters to be fed to the underlying 
;		called function.
;	_extra:  Extra keywords for the underlying command
;		See the individual documention.
;	logfile: Specify a file to write log messages
;
;
; OUTPUTS:
;	Writes the status_struct (see ::status_struct) to the first binary
;	table extension.
;
;	Writes the structures returned by the ::select methods to second 
;	extension. 
;
; MODIFICATION HISTORY:
;		* renamed process_datasweeps to process_runs which can now process
;			fpObjc files.
;		* require procrun input
;	2009-06-31:  Write a status structure to the first binary table extension.
;		This is much more powerful than using the header.
;docend::sweeps::process_runs

pro sweeps::process_runs, proctype, procrun, command,  $
		runs=runs, camcols=camcols, $
		addtags=addtags, $
		fpobjc=fpobjc, $
		match_method=match_method, $
		pars=pars, $
		output_dir=output_dir, $
		logfile=logfile, $
		_extra=_extra

	if n_elements(proctype) eq 0 or n_elements(procrun) eq 0 $
			or n_elements(command) eq 0 then begin
		splog,"usage: "
		splog,"  bt=obj_new('sweeps')"
		splog,"  bt->process_runs, proctype, procrun, command, output_dir=, runs=all, camcols=all, logfile=stdout, addtags=, /fpobjc, _extra="
		splog,"  See the ::select method in the individual classes for "
		splog,"  extra keywords that can be sent.  These will be passed along "
		splog,"  with the _extra= technique"
		on_error, 2
		message,'Halting'
	endif

	; begin logging
	splog, filename=logfile

	if n_elements(output_dir) eq 0 then begin
		output_dir = self->output_dir(proctype, procrun)
	endif

	if not file_test(output_dir) then begin
		splog,'Making output_dir: '+output_dir
		;file_mkdir, output_dir
		spawn,'mkdir -pv '+output_dir
	endif

	if not file_test(output_dir,/directory) then begin
		message,'Directory does not exist: '+string(output_dir)
	endif

	; Get some environment variables
	resdir = getenv('PHOTO_RESOLVE')
	if NOT keyword_set(resdir) then begin
		message,'Set PHOTO_RESOLVE'
	endif
	sweepdir = getenv('PHOTO_SWEEP')
	if NOT keyword_set(sweepdir) then begin
		message,'Set PHOTO_SWEEP'
	endif
	calibdir = getenv('PHOTO_CALIB')
	if NOT keyword_set(sweepdir) then begin
		message,'Set PHOTO_CALIB'
	endif

	;targetversion = bosstarget_version()

	; make sure we have our run list
	self->match_runlist, runs, runs2use, reruns2use

	stime0 = systime(1)
   
	nrun = n_elements(runs2use)

	if n_elements(logfile) ne 0 then begin
		splog, 'Log file ' + logfile + ' opened ' + systime()
	endif
	splog, 'Working with Nrun=', nrun
	spawn, 'uname -a', uname
	splog, 'Uname: ' + uname[0]


	if n_elements(camcols) eq 0 then camcols=[1,2,3,4,5,6]

	;---- For loop goes through all the runs

	for irun=0L, nrun-1L do begin
		run = runs2use[irun]
		rerun = reruns2use[irun]

		for icamcol=0, n_elements(camcols)-1 do begin
			camcol=camcols[icamcol]

			status_struct = self->status_struct($
				proctype=proctype,procrun=procrun,$
				run=run,rerun=rerun,camcol=camcol)

			status_struct.date =  systime()
			status_struct.photoop_v = photoop_version()
			status_struct.idlutils_v = idlutils_version()
			status_struct.photo_sweep = file_basename(sweepdir)
			status_struct.photo_resolve = file_basename(resdir)
			status_struct.photo_calib = file_basename(calibdir)

			output_file = self->output_file( $
				proctype, procrun, run, rerun, camcol, $
				output_dir=output_dir, $
				fpobjc=fpobjc, match_method=match_method)

			status_struct.output_file = file_basename(output_file)

			if keyword_set(fpobjc) then begin
				; run on the full fpObjc files
				res = self->process_fpobjc($
					run, rerun, camcol, select_tags, command, $
					pars=pars, $
					nobj=nobj, nres=nres, $
					_extra=_extra)
			endif else if n_elements(match_method) ne 0 then begin
				; take the match_method string as the name of a matching
				; method contained within the class for this type.
				; This should return a structure describing matches between
				; the sweep and some internal catalog
				res = self->match_sweep( $
					run, rerun, camcol, match_method, $
					nobj=nobj, $
					_extra=_extra)
			endif else begin
				; run on the sweeps
				res = self->process_sweep($
					run, rerun, camcol, command, $
					addtags=addtags, $
					pars=pars, $
					nobj=nobj, $
					nres=nres, $
					_extra=_extra, $
					flags=flags)
			endelse

			status_struct.nobj=nobj
			status_struct.nres=nres
			status_struct.flags = flags

			splog,$
				'Writing status structure to file: ',output_file,format='(a,a)'
			mwrfits, status_struct, output_file, /create

			if n_tags(res) ne 0 then begin
                nres=n_elements(res)
				splog, 'Run=',run,' camcol=', camcol,' Nres=',n_elements(res)
				splog,$
					'Appending ',nres,' results to output file: ', output_file,$
					format='(a,i0,a)'
				mwrfits, res, output_file
			endif 
	
		endfor
	endfor

	splog, 'Total time = ', systime(1)-stime0, ' seconds', format='(a,f6.0,a)'
	splog, 'Successful completion at ' + systime()
	splog, /close

	return

end

function sweeps::status_struct, $
		proctype=proctype, procrun=procrun, $
		run=run, rerun=rerun, camcol=camcol, $
		minimal=minimal

	st = { $
		proctype: '', $
		procrun: '', $
		photoop_v: '', $
		;bosstarget_v: '', $
		idlutils_v: '', $
		photo_sweep: '', $
		photo_resolve: '', $
        photo_calib:'' $
	}

	if n_elements(proctype) ne 0 then st.proctype=proctype
	if n_elements(procrun) ne 0 then st.procrun=procrun

	if not keyword_set(minimal) then begin
		st_rest = { $
			run:-9999L,$
			rerun:'-9999L',$
			camcol:-9999L, $
			output_file: '', $
			date: '', $
			nobj: 0L, $
			nres: 0L, $
			flags: 0L}
		st = create_struct(st, st_rest)


		if n_elements(run) ne 0 then st.run=run
		if n_elements(rerun) ne 0 then st.rerun=rerun
		if n_elements(camcol) ne 0 then st.camcol=camcol
	endif

	return, st
end

function sweeps::processing_flag,name
	case strlowcase(name) of
		'read_failed': return, 2L^0
		'noresult': return, 2L^1
		'nofile': return, 2L^2
		'command_failed': return, 2L^3
        'status_read_failed': return, 2L^4
		else: message,'No such flag name: '+string(name)
	endcase
end
function sweeps::processing_flagnames
    return, ['read_failed','noresult','nofile','command_failed','status_read_failed']
end
pro sweeps::print_processing_flags, flags
    fn=self->processing_flagnames()
    for i=0L, n_elements(fn)-1 do begin
        tfn = fn[i]
        if (flags and self->processing_flag(tfn)) ne 0 then begin
            splog,'    Flag set: ',fn
        endif
    endfor
end

function sweeps::procrunstring, procrun
	; if a number is entered, convert it length 6 zero-padded string
	if size(procrun, /tname) ne 'STRING' then begin
		;tstr = string(long(procrun), f='(i06)')
		tstr = strn(long(procrun), len=6,padchar='0')
	endif else begin
		return, procrun
	endelse

end
function sweeps::output_dir, proctype, procrun
	if n_elements(procrun) eq 0 or n_elements(proctype) eq 0 then begin
		on_error,2
		message,'Usage: dir=bt->output_dir(proctype, procrun)'
	endif
	dir = getenv('SWEEP_REDUCE')
	if dir eq '' then message,'SWEEP_REDUCE is not set'

	runstr = self->procrunstring(procrun)
	dir = filepath(root=dir, sub=proctype, runstr)
	return, dir
end


function sweeps::output_file, proctype, procrun, run, rerun, camcol, $
		output_dir=output_dir, $
		fpobjc=fpobjc, $
		extra_name=extra_name, $
		match_method=match_method, $
		gather=gather, $
		collate=collate, $
		verify=verify

	if keyword_set(gather) or keyword_set(collate) or keyword_set(verify) then begin
		gorc=1
		npar = 2
	endif else begin
		gorc=0
		npar = 5
	endelse
	if n_params() lt npar then begin
		on_error,2
		message,'Usage: f=bt->output_file(proctype, procrun, run, rerun, camcol, output_dir=, /fpobjc, /old, /gather, /collate, match_method=)'
	endif

	if n_elements(output_dir) eq 0 then begin
		output_dir=self->output_dir(proctype, procrun)
	endif

	fname='sweep'+self.type
	if keyword_set(fpobjc) then begin
		fname += '-fpobjc'
	endif else if n_elements(match_method) ne 0 then begin
		mstr=repstr(match_method,'_','-')
		fname += '-'+mstr
	endif
	fname+="-"+proctype

	procrunstr = self->procrunstring(procrun)

	if not keyword_set(gorc) then begin
		fname=[fname,strn(long(run),len=6,padchar='0')]
		fname=[fname,strn(camcol)]
		if not keyword_set(old) then begin
			fname=[fname,strn(rerun)]
		endif
	endif 

	if not keyword_set(old) then fname = [fname, procrunstr]

	if keyword_set(collate) then begin
		fname = [fname,'collate']
	endif else if keyword_set(gather) then begin
		fname = [fname,'gather']
	endif else if keyword_set(verify) then begin
		fname = [fname,'verify']
	endif

	if keyword_set(extra_name) then begin
		fname = [fname,extra_name]
	endif
	fname = strjoin(fname, '-') + '.fits'

	fpath = filepath(root=output_dir, fname)
	return, fpath
end



function sweeps::process_sweep, $
		run, rerun, camcol, command, $
		nobj=nobj, nres=nres, $
		pars=pars, $
		addtags=addtags, $
		_extra=_extra, $
		flags=flags

	res = -1
	nobj=0
	nres=0

	flags = 0L

    sf=obj_new('sdss_files')

	if keyword_set(addtags) then begin
		objs = self->augment_sweep(run, camcol, rerun=rerun, type=self.type, $
			addtags=addtags)
	endif else begin
		objs = sf->read('calibobj.'+self.type, run, camcol, rerun=rerun)
	endelse

	if n_tags(objs) ne 0 then begin
		nobj = n_elements(objs)

		if not execute(command) then begin
			splog,'SWEEP_ERROR: Failed to execute command: "',command,"'",format='(a,a,a)'
			splog,'Returning -1'
			flags += self->processing_flag('command_failed')
		endif else begin
			; success
			if n_tags(res) ne 0 then nres=n_elements(res)
		endelse
	endif else begin
		splog,'SWEEP_ERROR: Failed to read file: ',sweep_file
		flags += self->processing_flag('read_failed')
	endelse

	if nres eq 0 and self.pars.require_result then begin
		; only set this if we care
		splog,'SWEEP_ERROR: expected result: ',sweep_file
		flags += self->processing_flag('noresult')
	endif
	
	return, res
end


function sweeps::process_fpobjc, $
		run, rerun, camcol, select_tags, command, $
		nobj=nobj, nres=nres, $
		pars=pars, $
		_extra=_extra, $
		flags=flags

	; we want to encourage selection of a subset!  Too much memory
	; usage otherwise
	if n_elements(select_tags) eq 0 then message,'you must send select_tags!'

	res = -1
	nobj=0
	nres=0

	flags = 0L

	objs = sdss_readobj(run, camcol, rerun=rerun, select_tags=select_tags)

	if n_tags(objs) ne 0 then begin
		nobj = n_elements(objs)

		if not execute(command) then begin
			splog,'Failed to execute command: "',command,"'",form='(a,a,a)'
			splog,'Returning -1'
			flags += self->processing_flag('command_failed')
		endif else begin
			; success
			if n_tags(res) ne 0 then nres=n_elements(res)
		endelse
	endif else begin
		flags += self->processing_flag('read_failed')
		splog,'Failed to read file: ',sweep_file
	endelse

	if nres eq 0 and self.pars.require_result then begin
		; only set this if we care
		flags += self->processing_flag('noresult')
	endif
	
	return, res

end

function sweeps::match_sweep, $
		run, rerun, camcol, match_method, nobj=nobj, _extra=_extra

	if n_params() lt 5 then begin
		on_error,2
		message,'usage: res=bt->match_sweep(run, rerun, camcol, match_method, _extra=)'
	endif
	sweep_file = $
		sdss_name('calibObj.'+self.type,run,camcol,rerun=rerun)

	splog,'Reading file: ',sweep_file,'[.gz]',format='(a,a,a)'
	targ = obj_new('bosstarget_'+proctype)

    sf=obj_new('sdss_files')
	objs = sf->read('calibObj.'+self.type, run, camcol, rerun=rerun)

	if n_tags(objs) ne 0 then begin
		nobj=n_elements(objs)
		command = 'res = targ->'+match_method+'(objs, _extra=_extra)'
		if not execute(command) then begin
			message,'Could not execute command: '+command
		endif
	endif else begin
		nobj=0
		splog,'Failed to read file: ',sweep_file
		res = -1
	endelse

	obj_destroy, targ

	return, res
end




function sweeps::make_bayes_input, struct

	if n_elements(struct) eq 0 then begin
		on_error, 2
		message,'usage: bayes_instruct=bt->make_bayes_input(struct)'
	endif

	arrval=fltarr(5)
	larrval = lonarr(5)
	inst = {$
		psfcounts: arrval, $
		psfcountserr: arrval, $
		counts_dev: arrval, $
		counts_deverr: arrval, $
		counts_exp:arrval, $
		counts_experr:arrval, $
		fracpsf: arrval, $
		m_rr_cc_psf: arrval, $
		flags: larrval,$
		flags2: larrval }	
		
	inst = replicate(inst, n_elements(struct))

	struct_assign, struct, inst, /nozero
	sdss_flux2lups, struct.psfflux, struct.psfflux_ivar, struct.psfcountserr, $
		psfcounts, psfcountserr
	inst.psfcounts = temporary(psfcounts)
	inst.psfcountserr = temporary(psfcountserr)

	sdss_flux2lups, struct.devflux, struct.devflux_ivar, struct.counts_deverr, $
		counts_dev, counts_deverr	
	inst.counts_dev = counts_dev
	inst.counts_deverr = temporary(counts_deverr)

	sdss_flux2lups, struct.expflux, struct.expflux_ivar, struct.counts_experr, $
		counts_exp, counts_experr	
	inst.counts_exp = counts_exp
	inst.counts_experr = temporary(counts_experr)


	return, inst


end

pro sweeps::run_bayes_prob, struct, probgal, probflags

	if n_elements(struct) eq 0 then begin
		on_error, 2
		message,'usage: run_bayes_prob, struct, probgal, probflags'
	endif
	splog,'Making input'
	inst=self->make_bayes_input(struct)
	compute_bayes_prob, inst, probgal, probflags

end


function sweeps::calibobj_type, proctype
	case proctype of
		'pzgal': begin
			return, 'gal'
		end
		'am': begin
			return, 'gal'
		end
		'pzstar': begin
			return, 'star'
		end
        'flstar': begin
            return,'star'
        end
        'flgal': begin
            return,'gal'
        end
		else: message,'Unsupported process type: '+string(proctype)
	endcase

end





function sweeps::read, proctype, procrun, run, rerun, camcol, $
		output_dir=output_dir, $
		pars=send_pars, $
		reselect=reselect, $
		collate=collate, $
		addtags=addtags, $
		columns=columns, $
		where_string=where_string, $
		filename=filename, $
		fpobjc=fpobjc, $
		match_method=match_method, $
		status_struct=status_struct, $
		_extra=_extra


	if n_params() lt 5 then begin
		splog,"usage: "
		splog,"  bt=obj_new('sweeps')"
		splog,"  res = bt->read(proctype, procrun, run, rerun, camcol, output_dir=, /collate, addtags=, columns=, where_string=, filename=, /fpobjc, match_method=, /addlups, status_struct=)"
		splog,"  In the where string, refer the structure 'str'.  e.g. "
		splog,"      str.boss_target1 gt 0 and str.ra gt 336"
		on_error, 2
		message,'Halting'
	endif

	pars=self->pars()

	calibobj_file = $
		sdss_name('calibObj.'+self.type, run, camcol, rerun=rerun)
	output_file = self->output_file( $
		proctype, procrun, run, rerun, camcol, output_dir=output_dir, $
		fpobjc=fpobjc, match_method=match_method)

	if not file_test(output_file) then begin
		splog,'File not found: ',output_file,format='(a,a)'
		return,-1
	endif

	splog,'Reading output file: ',output_file,format='(a,a)'
	if pars.result_ext ne 1 then begin
		status_struct = mrdfits(output_file, pars.status_ext, $
			silent=0, status=rd_status)
		str = mrdfits(output_file, pars.result_ext, $
			silent=0, status=rd_status)
	endif else begin
		str = mrdfits(output_file, pars.result_ext, $
			/silent, status=rd_status)
	endelse
	if rd_status ne 0 then begin
		on_error, 2
		splog,'Error reading file: '+output_file
		return,-1
	endif
	if keyword_set(collate) then begin
		splog,'Collating'

		if keyword_set(fpobjc) then begin
			co = sdss_readobj(run, camcol, rerun=rerun)
		endif else begin
			splog,'Reading calibObj file: ',calibobj_file+'[.gz]',$
				format='(a,a,a)'
			co = self->augment_sweep(run, camcol, rerun=rerun, $
				type=self.type, addtags=addtags)
		endelse


		; make sure the result and sweep match up as expected
		; doesn't make sense for /fpobjc probably
		if self.pars.require_result then begin
			; we require a res for each input
			wbad = where(co.id ne str.id, nbad)
			if nbad ne 0 then begin
				message,"Structures don't line up"	
			endif
		endif else begin
			; just do a match
			sphoto_match, co, str, mco, mstr
			if mco[0] ne -1 then begin
				if n_elements(mstr) eq n_elements(str) then begin
					co = co[mco]
					str = str[mstr]
				endif else begin
					message,'Some objects did not match!'
				endelse
			endif else begin
				message,"no matches found!"
			endelse
		endelse



		; using my own code here since the idlutils code and sdssidl code
		; have namespace conflicts

		; add in extra tags from the output
		tstags = tag_names(str)
		cotags = tag_names(co)
		match, tstags, cotags, mts, mco

		; these will be tags only found in the ts file
		keepind = lindgen(n_elements(tstags))
		remove, mts, tstags
		remove, mts, keepind

		newst = co[0]
		for i=0L, n_elements(tstags)-1 do begin
			newst = create_struct(newst, tstags[i], str[0].(keepind[i]))
		endfor
		newst = replicate(newst, n_elements(str))
		struct_assign, str, newst, /nozero
		struct_assign, co, newst, /nozero

		; rename
		str=0
		str = temporary(newst)

	endif 

	; attempt to re-select based on the input parameters
	; the select function must support the /reselect keyword
	if n_elements(send_pars) ne 0 or keyword_set(reselect) then begin
		splog,'Re-selecting with input pars'
		targ=obj_new('bosstarget_'+proctype)
		res = targ->select(str, pars=send_pars, /reselect, _extra=_extra)

		wdiff=where(res.boss_target1 ne str.boss_target1, ndiff)
		splog,'# differences found in boss_target1: ',ndiff,form='(a,i0)'
		struct_assign, res, str, /nozero
	endif

	; make sure to do the selection *before* we extract
	; tags, to give the user full flexibility
	if n_elements(where_string) ne 0 then begin
		w=self->where_select(str, where_string, nw)
		if nw eq 0 then begin
			return, -1
		endif
		str = str[w]
	endif

	; now extract certain columns(tags)
	if n_elements(columns) ne 0 then begin
		tstr = struct_selecttags(str, select_tags=columns)
		if size(tstr, /tname) ne 'STRUCT' then begin
			message,'None of the requested tags matched: '+$
				'['+strjoin( string(columns), ', ')+']'
		endif
		str=0
		str = temporary(tstr)
	endif

	return, str

end

function sweeps::where_select, str, where_string, nw

	command = 'w = where('+where_string+', nw)'
	splog,'Executing where statement: "',command,'"',$
		form='(a,a,a)'
	if not execute(command) then begin
		message,'Could not execute where_string: '+string(where_string)
	endif
	return, w
end


pro sweeps::cache_runlist, $
		photo_resolve=photo_resolve, $
        minscore=minscore, $
        force=force
	common sweeps_runlist_block, _photo_resolve, _minscore, _flrun, _flrerun

    docache=0

    default_minscore = 0.1 + 0.001
    if n_elements(_photo_resolve) eq 0 then begin
        _photo_resolve=getenv('PHOTO_RESOLVE')
        _minscore=default_minscore
        docache=1
    endif
    if n_elements(photo_resolve) ne 0 then begin
        if photo_resolve ne _photo_resolve then begin
            docache=1
        endif
    endif
    if n_elements(minscore) ne 0 then begin
        if minscore ne _miscore then begin
            docache=1
        endif
    endif

    if keyword_set(docache) or keyword_set(force) then begin

		if n_elements(photo_resolve) ne 0 then begin
            _photo_resolve = photo_resolve
			sold=getenv('PHOTO_RESOLVE')
            splog,'Setting PHOTO_RESOLVE equal to: ',_photo_resolve
			setenv,'PHOTO_RESOLVE='+_photo_resolve
		endif else begin
			_photo_resolve = getenv('PHOTO_RESOLVE')
            splog,'Using PHOTO_RESOLVE: ',_photo_resolve
        endelse

        if n_elements(minscore) ne 0 then begin
            _minscore=minscore
        endif else begin
            _minscore=default_minscore
        endelse

        splog, 'Cacheing window_runlist, ..,  minscore=',_minscore,$
            format='(a,f)'
        window_runlist, _flrun, rerun=_flrerun, minscore=_minscore
        nflrun = n_elements(_flrun)
        splog,'Found: ',nflrun,' runs',format='(a,i0,a)'

		if n_elements(photo_resolve) ne 0 then begin
			setenv,'PHOTO_RESOLVE='+sold
            splog,'Setting PHOTO_RESOLVE back to: ',photo_resolve
		endif

		return

	endif
end

pro sweeps::runlist, flrun, flrerun, $
		minscore=minscore, photo_resolve=photo_resolve, force=force
	common sweeps_runlist_block, _photo_resolve, _minscore, _flrun, _flrerun

	self->cache_runlist, $
			photo_resolve=photo_resolve, force=force, minscore=minscore

	flrun   = _flrun
	flrerun = _flrerun

end

pro sweeps::split_runlist, runs2use, reruns2use, $
		nsplit=nsplit, nper=nper, $
		runs=runs, $
		photo_resolve=photo_resolve, force=force, minscore=minscore, $
        seed=seed

	if n_elements(runs) ne 0 then begin
		self->match_runlist, runs, runs2use, reruns2use, $
			photo_resolve=photo_resolve, force=force, minscore=minscore
	endif else begin
		self->runlist, runs2use, reruns2use, $
			photo_resolve=photo_resolve, force=force, minscore=minscore
	endelse

	if n_elements(nsplit) ne 0 or n_elements(nper) ne 0 then begin
        if n_elements(seed) ne 0 then begin
            ; scramble the runs, since we have multiple large runs
            ; in order 752/756 e.g.
            r=sort(randomu(seed, n_elements(runs2use)))
            runs2use = runs2use[r]
            reruns2use = reruns2use[r]
        endif
		runs2use = self->splitlist(runs2use, nsplit=nsplit, nper=nper)
		reruns2use = self->splitlist(reruns2use, nsplit=nsplit,  nper=nper)
	endif

end

function sweeps::splitlist, list, nsplit=nsplit, nper=nper

	nlist = n_elements(list)
	if nlist eq 0 then begin
		message,'Usage: slist=bt->splitlist(list, nsplit=, nper=)'
	endif

	if n_elements(nper) ne 0 then begin
		nsplit = (nlist/nper) + ((nlist mod nper) gt 0)
	endif else if n_elements(nsplit) ne 0 then begin
		nper = nlist/nsplit
		;nleft = nlist mod nsplit
	endif else begin
		message,'send nsplit= or nper='
	endelse

	split_ptrlist = ptrarr(nsplit)
	
	current = 0LL
	for i=0L, nsplit-1 do begin
		if i eq (nsplit-1) then begin
			; last one, make sure we include the remainder
			split = list[current:nlist-1]
		endif else begin
			split = list[current:current+nper-1]
		endelse

		split_ptrlist[i] = ptr_new(split, /no_copy)
		current = current + nper
	endfor

	return, split_ptrlist

end




pro sweeps::match_runlist, runs, match_runs, match_reruns, $
		photo_resolve=photo_resolve, force=force, minscore=minscore, $
        minput=minput, mflrun=mflrun

	common sweeps_runlist_block, _photo_resolve, _minscore, _flrun, _flrerun
	self->cache_runlist, photo_resolve=photo_resolve, force=force, minscore=minscore

	if n_elements(runs) eq 0 then begin
		splog,'Using all runs'
		match_runs = _flrun
		match_reruns = _flrerun
		return
	endif

	splog,'Matching input runlist'

	rmd = rem_dup(runs)
	match, runs[rmd], _flrun, minput, mflrun
	if mflrun[0] eq -1 then begin
		message,'None of input runs matched window_runlist'
	endif
	if n_elements(mflrun) ne n_elements(rmd) then begin
		message,'Some of the input runs did not match',/inf
	endif

	splog, 'Using ',n_elements(mflrun),'/',n_elements(_flrun),' runs'
	match_runs = _flrun[mflrun]
	match_reruns = _flrerun[mflrun]

end




pro sweeps::gather_partial, proctype, procrun, jobnum, $
		runs=runs, $
		nper=nper, $
		collate=collate, $
		addtags=addtags, $
		pars=send_pars, reselect=reselect, extra_name=extra_name, $
		where_string=where_string, $
		everything=everything, $
		fpobjc=fpobjc, $
		combine=combine,ascii=ascii, $ ; these with regard to combine only
		match_method=match_method, $
		outfile=outfile, $
		noverify=noverify, $
		_extra=_extra

	if n_elements(proctype) eq 0 or (n_elements(jobnum) eq 0 and not keyword_set(combine)) then begin
		on_error, 2
		message,'bt->gather_partial, proctype, procrun, jobnum, nper=2, /collate, addtags=, where_string=, /everything, /fpobjc, /combine, /ascii, extra_name=, match_method=, outfile=, _extra=', /inf
		message,'Send /combine to combine the individual files.  ',/inf
		message,'  /ascii is with regard to the combined file only'
	endif

	pars=self->pars()

	if n_elements(nper) eq 0 then begin
		nper=2
	endif
	
	if not keyword_set(combine) then begin


		self->split_runlist, runptrs, rerunptrs, nper=nper, $
			runs=runs, /force
		njobs=n_elements(runptrs)

		runs = *runptrs[jobnum]
		reruns = *rerunptrs[jobnum]
		ptr_free, runptrs, rerunptrs

		splog, 'Using runs: ',runs


		if n_elements(outfile) eq 0 then begin
			outfile=self->output_file(proctype, procrun, $
				extra_name=extra_name, $
				fpobjc=fpobjc, /gather, collate=collate, $
				match_method=match_method)
			;outfile = $
			;	repstr(outfile,'.fits','-'+string(jobnum,f='(i03)')+'.fits')
			outfile = $
				repstr(outfile,'.fits','-'+strn(jobnum,len=3,padchar='0')+'.fits')
		endif
		splog,'outfile: ',outfile,format='(a,a)'

		self->gather2file, proctype, procrun, $
			pars=send_pars, reselect=reselect, extra_name=extra_name, $
			collate=collate, $
			addtags=addtags, $
			where_string=where_string, $
			everything=everything, $
			runs=runs, $
			columns=columns, /fast,  $
			fpobjc=fpobjc, $
			match_method=match_method, $
			noverify=noverify, $
			outfile=outfile, _extra=_extra

	endif else begin

		splog,'Combining various collated files'
		if n_elements(where_string) ne 0 then begin
			splog,'Combining with where string "'+where_string+'"'
		endif

		if not keyword_set(noverify) then begin
			splog,'First verifying all outputs'

			nbad = self->verify(proctype, procrun, $
				runs=runs, $
				fpobjc=fpobjc, match_method=match_method)
				;extra_name=extra_name)
			if nbad ne 0 then begin
				splog,'Found bad results: ',nbad
				message,'Halting'
			endif
		endif

		self->split_runlist, runptrs, rerunptrs, nper=nper, $
			runs=runs, /force
		njobs=n_elements(runptrs)
		ptr_free, runptrs, rerunptrs

		; combine all the partials
		outfile=self->output_file(proctype, procrun, $
			extra_name=extra_name, $
			fpobjc=fpobjc, match_method=match_method, /gather, $
			collate=collate)


		; Got to roll my own since we don't have sdssidl mrdfits_multi
		flist=strarr(njobs)

		status_ptrlist=ptrarr(njobs)
		res_ptrlist=ptrarr(njobs)
		for job=0L, njobs-1 do begin
			;flist[job] = $
			;	repstr(outfile,'.fits','-'+string(job,f='(i03)')+'.fits')
			flist[job] = $
				repstr(outfile,'.fits','-'+strn(job,len=3,padchar='0')+'.fits')
			splog,'Reading: ',flist[job],form='(a,a)'
			st=mrdfits(flist[job],pars.status_ext, status=stst)

			err=string($
				'Failed to read ',flist[job],$
				' ext=',pars.status_ext,form='(a,a,a,i0)')
			if stst ne 0 then begin
				splog,err
			endif else begin
				status_ptrlist[job] = ptr_new(st, /no_copy)
			endelse
			err=string($
				'Failed to read ',flist[job],$
				' ext=',pars.result_ext,form='(a,a,a,i0)')
			res=mrdfits(flist[job],pars.result_ext, status=tst)
			if tst ne 0 then begin
				splog,err
			endif else begin
				if n_elements(where_string) ne 0 then begin
					nres=n_elements(res)
					w=self->where_select(res,where_string,nw)
					splog,'Keeping ',strn(nw),'/',strn(nres),$
						form='(a,a,a,a)'
					if nw ne 0 then begin
						res=res[w]
					endif else begin
						res=0
					endelse
				endif
				if n_tags(res) ne 0 then begin
					res_ptrlist[job] = ptr_new(res, /no_copy)
				endif
			endelse


		endfor
		tot = combine_ptrlist(res_ptrlist)
		status_tot = combine_ptrlist(status_ptrlist)

		if n_tags(tot) eq 0 then message,'Failed to read'

		; now just keep one,verify above makes sure they are the same
		status_struct = status_tot[0]

		ntot=n_elements(tot)
		if keyword_set(ascii) then begin
			outfile = repstr(outfile, '.fits', '-tab.st')
			splog,'Writing ',ntot,' to file: ',outfile,format='(a,i0,a)'
			; requires sdssidl
			write_idlstruct, tot, outfile, /ascii, hdrstruct=status_struct
		endif else begin
			splog,'Writing status to file: ',outfile,format='(a,a)'
			mwrfits, status_struct, outfile, /create
			splog,'Writing ',ntot,' to file: ',outfile,format='(a,i0,a)'
			mwrfits, tot, outfile
		endelse

	endelse

end








function sweeps::read_gather, proctype, procrun,  $
		collate=collate, $
		output_dir=output_dir, $
		extra_name=extra_name, $
		ascii=ascii, $
		fpobjc=fpobjc, match_method=match_method, $
		status_struct=status_struct, file=file

	if n_elements(proctype) eq 0 or n_elements(procrun) eq 0 then begin
		on_error, 2
		message,'Usage: str=bt->read_gather(proctype, procrun, /collate, output_dir=, /fpobjc, match_method=, status_struct=)', /inf
		message,'Halting'
	endif

	pars=self->pars()
	file=self->output_file(proctype, procrun, $
		output_dir=output_dir, fpobjc=fpobjc, $
		/gather, collate=collate, $
		extra_name=extra_name, $
		match_method=match_method)
	splog,'Reading file: ',file,format='(a,a)'
	if keyword_set(ascii) then begin
		file = repstr(file, '.fits', '-tab.st')
		command=' st=read_idlstruct(file,hdr=status_struct)'
		if not execute(command) then begin
			message,"You probably don't have sdssidl stup"
		endif
	endif else begin
		status_struct = mrdfits(file, pars.status_ext)
		st=mrdfits(file,pars.result_ext)
	endelse
	return, st
end
pro sweeps::gather_usage, proc=proc

	if keyword_set(proc) then begin
		front="bt->gather2file, "
		back=""
	endif else begin
		front="res = bt->gather("
		back=")"
	endelse

	splog,"usage: "
	splog,"  bt=obj_new('sweeps')"
	splog,"  "+front+"proctype, output_dir=, runs=all, logfile=stdout, columns=, where_string=, /everything, /fast"+back+", /noverify"
	splog,"  In the where string, refer the structure 'str'.  e.g. "
	splog,"      str.boss_target1 gt 0 and str.ra gt 336"
	splog,"  /fast uses twice the memory but only does one pass"
end
pro sweeps::gather2file, proctype, procrun, $
		collate=collate, $
		addtags=addtags, $
		output_dir=output_dir, $
		runs=runs, camcols=camcols, nomatch=nomatch, $
		pars=send_pars, reselect=reselect, extra_name=extra_name, $
		logfile=logfile, $
		columns=columns, $
		where_string=where_string, everything=everything, $
		fast=fast, $
		fpobjc=fpobjc, match_method=match_method, $
		noverify=noverify, $
		outfile=outfile, _extra=_extra

	; this procedural version writes to a file

	if n_elements(proctype) eq 0 or n_elements(procrun) eq 0 then begin
		self->gather_usage,/proc
		on_error, 2
		message,'Halting'
	endif

	if n_elements(outfile) eq 0 then begin
		outfile=self->output_file(proctype, procrun, $
			output_dir=output_dir, fpobjc=fpobjc, $
			extra_name=extra_name, $
			match_method=match_method, $
			/gather, collate=collate)
	endif
	splog,'Will write to output file: ',outfile,format='(a,a)'

	res = self->gather($
		proctype, procrun, $
		output_dir=output_dir, $
		runs=runs, camcols=camcols, nomatch=nomatch, $
		pars=send_pars, reselect=reselect, $
		logfile=logfile, columns=columns, $
		where_string=where_string, everything=everything, $
		collate=collate, $
		addtags=addtags, $
		fast=fast, $
		fpobjc=fpobjc,match_method=match_method, $
		noverify=noverify, $
		status_struct=status_struct, _extra=_extra)

	if n_tags(res) eq 0 then begin
		splog,'No results found, not writing file: ',outfile,form='(a,a)'
		return
	endif

	; in order to pass the verify step, the versions must all match so we
	; can just copy in from the first
	st = self->status_struct(/minimal)
	struct_assign, status_struct[0], st, /nozero

	splog,'Writing to output file: ',outfile,format='(a,a)'
	mwrfits, st, outfile, /create
	mwrfits, res, outfile

end

function sweeps::gather, proctype, procrun, $
		output_dir=output_dir, $
		runs=runs, camcols=camcols, nomatch=nomatch, $
		pars=send_pars, reselect=reselect, $
		where_string=where_string, $
		add_where_string=add_where_string, $
		anyresolve=anyresolve, $
		run_primary=run_primary, $
		everything=everything, $
		logfile=logfile, columns=columns, $
		fast=fast, $ ; this is the only option now
		fpobjc=fpobjc, $
		match_method=match_method, $
		noverify=noverify, $
		collate=collate, $
		addtags=addtags, $
		status_struct=status_struct, $
		_extra=_extra



	if n_elements(proctype) eq 0 or n_elements(procrun) eq 0 then begin
		self->gather_usage
		on_error, 2
		message,'Halting'
	endif

	if not keyword_set(noverify) and not keyword_set(nomatch) then begin
		nbad = self->verify(proctype, procrun, $
			output_dir=output_dir, runs=runs, camcols=camcols, $
			logfile=logfile, $
			fpobjc=fpobjc, $
			match_method=match_method)
		if nbad ne 0 then begin
			splog,'Found bad results: ',nbad
			message,'Halting'
		endif
	endif

	if n_elements(add_where_string) ne 0 $
			and n_elements(where_string) ne 0 then begin
		where_string += ' (' + add_where_string + ')'
	endif

	if n_elements(where_string) ne 0 then begin
		splog,'Gathering with where string: "'+where_string+'"'
	endif
	; begin logging
	splog, filename=logfile

	; Get some environment variables
	resdir = getenv('PHOTO_RESOLVE')
	if NOT keyword_set(resdir) then begin
		message,'Set PHOTO_RESOLVE'
	endif

	if not keyword_set(nomatch) then begin
		self->match_runlist, runs, runs2use, reruns2use
	endif else begin
		if n_elements(runs) eq 0 then begin
			message,'send runs if /nomatch is set'
		endif
		runs2use=runs
		splog,'Assuming reruns are 301'
		reruns2use = replicate(301,n_elements(runs))
	endelse
  
	stime0 = systime(1)

	nrun = n_elements(runs2use)

	if n_elements(logfile) ne 0 then begin
		splog, 'Log file ' + logfile + ' opened ' + systime()
	endif
	splog, 'Working with Nrun=', nrun
	spawn, 'uname -a', uname
	splog, 'Uname: ' + uname[0]


	;---- For loop goes through all the runs
	;---- (first count total number)

	ptrlist=ptrarr(nrun*6)
	if arg_present(status_struct) then begin
		sptrlist=ptrarr(nrun*6)
		do_status=1
	endif else begin
		do_status=0
	endelse

	itot=0L
	for irun=0L, nrun-1L do begin
		run = runs2use[irun]
		rerun = reruns2use[irun]
		for camcol=1, 6 do begin

			tstr= self->read(proctype, procrun, run, rerun, camcol, $
				output_dir=output_dir,$
				pars=send_pars, $
				reselect=reselect, $
				where_string=where_string, $
				collate=collate, $
				addtags=addtags, $
				match_method=match_method, $
				columns=columns,$
				fpobjc=fpobjc,  $
				status_struct=tstatus_struct, $
				_extra=_extra)

			if(n_tags(tstr) gt 0) then begin
				splog,'Keeping ',n_elements(tstr),' objects', $
					format='(a,i0,a)'

				ptrlist[itot] = ptr_new(tstr,/no_copy)

				if do_status then begin
					sptrlist[itot] = ptr_new(tstatus_struct,/no_copy)
				endif

				itot = itot+1

			endif else begin
				splog,'Keeping zero objects'
			endelse

		endfor
	endfor

	outst = combine_ptrlist(ptrlist)
	if do_status then begin
		status_struct=combine_ptrlist(sptrlist)
	endif

	if n_tags(outst) ne 0 then begin
		ntot=n_elements(outst)
	endif else begin
		ntot=0
	endelse
	splog,'Keeping total of ',ntot,' objects', format='(a,i0,a)'
	splog,'Total time = ', systime(1)-stime0, ' seconds', format='(a,f6.0,a)'
	splog,'Successful completion at ' + systime()
	splog,/close

	return, outst
end 

function sweeps::verify, proctype, procrun, writebad=writebad, $
		output_dir=output_dir, $
		runs=runs, camcols=camcols, $
		extra_name=extra_name, $
		logfile=logfile, $
		fpobjc=fpobjc, $
		match_method=match_method


	if n_elements(proctype) eq 0 or n_elements(procrun) eq 0 then begin
		splog,'Usage: nbad = bt->verify(proctype, procrun, output_dir=, /fpobjc, match_method=)'
		on_error, 2
		message,'Halting'
	endif


	pars=self->pars()
	if n_elements(camcols) eq 0 then camcols=[1,2,3,4,5,6]

	; begin logging
	splog, filename=logfile

	self->match_runlist, runs, runs2use, reruns2use, /force
  
	stime0 = systime(1)

	nrun = n_elements(runs2use)

	if n_elements(logfile) ne 0 then begin
		splog, 'Log file ' + logfile + ' opened ' + systime()
	endif
	splog, 'Verifying Nrun=', nrun

	ptrlist=ptrarr(nrun*6)
	ii=0LL
	for irun=0L, nrun-1L do begin
		run = runs2use[irun]
		rerun = reruns2use[irun]
		for icamcol=0L, n_elements(camcols)-1 do begin
			camcol=camcols[icamcol]

			output_file = self->output_file( $
				proctype, procrun, run, rerun, camcol, $
				output_dir=output_dir, $
				extra_name=extra_name, $
				fpobjc=fpobjc, match_method=match_method)

			if not file_test(output_file) then begin
				splog,'SWEEP_ERROR: output file missing: ',output_file,form='(a,a)'
				message,'Fatal error: halting'
			endif

            tstatus=mrdfits(output_file, pars.status_ext, /silent, status=ts)
            if ts ne 0 then begin
                splog,'SWEEP_ERROR: Problem reading status from output file: ',output_file,form='(a,a)'
                tstatus = self->status_struct()
                tstatus.flags += self->processing_flag('status_read_failed')
            endif else begin
                if tstatus.flags ne 0 then begin
                    splog,'SWEEP_ERROR: Problem with result in output file: ',output_file,form='(a,a)'
                    self->print_processing_flags, tstatus.flags
                endif
            endelse

			ptrlist[ii] = ptr_new(tstatus, /no_copy)

			ii+=1
		endfor
	endfor

	status_struct = combine_ptrlist(ptrlist)
	help,status_struct,/str

	; version problems are a fatal error, an exception will be thrown
	self->check_version_tags, status_struct

	; only write a file if problems were found
	wbad=where(status_struct.flags ne 0, nbad)
	if nbad ne 0 then begin
		splog,'Found problems with ',nbad,' files'
		if keyword_set(writebad) then begin
			badfile=self->output_file(proctype, procrun, /verify, $
				output_dir=output_dir, fpobjc=fpobjc, $
				extra_name=extra_name, $
				match_method=match_method)
			splog,'Writing verify file: ',badfile,form='(a,a)'
			mwrfits, status_struct, badfile, stname='SWEEPREDUCE_STATUS', $
				/create
		endif
	endif else begin
		splog,'No problems found with status structures'
	endelse

	return, nbad
end

pro sweeps::check_version_tags, status_struct
	; run through all the relevant version tags and make sure they are 
	; equal.  This is fatal and an exception is thrown
	tags=strlowcase(tag_names(status_struct))
    if status_struct[0].proctype eq 'regauss' and status_struct[0].procrun eq '03' then begin
        ; there was a mishap with versioning
        tagcheck = ['procrun', 'proctype', 'photo_sweep', 'photo_resolve']
    endif else begin
        tagcheck = ['photoop_v', 'idlutils_v', $
            'procrun', 'proctype', 'photo_sweep', 'photo_resolve']
    endelse
	for i=0L, n_elements(tagcheck)-1 do begin
		w=where(tags eq tagcheck[i], nw)
		if nw eq 0 then message,'Tag not found: '+tagcheck[i]


		w_unique_tag = rem_dup(status_struct.(w))
        nu = n_elements(w_unique_tag)
		if nu ne 1 then begin
            utag = status_struct[w_unique_tag].(w)
			splog,"Not all '"+tagcheck[i]+"' tags are the same!",/inf
            splog,'Unique values: ',utag
            for j=0L, nu-1 do begin
                help,status_struct[w_unique_tag[j]],/str
            endfor
            message,'halting'
		endif
	endfor
end





; this depends on the idlutils function sdss_readobjlist
function sweeps::augment_sweep, run, camcol, rerun=rerun, $
		addtags=addtags_in, $
		matchstruct=matchstruct, $
		_extra=_extra

	; first read the sweep
    sf=obj_new('sdss_files')
	sweep = sf->read('calibobj.'+self.type, run, camcol, rerun=rerun)

	if n_elements(matchstruct) ne 0 then begin
		; grab a subset
		sphoto_match, sweep, matchstruct, msweep, mmatch
		if msweep[0] eq -1 then message,'no objects matched'
		sweep = sweep[msweep]
	endif

	; get subset of addtags that isn't already in the sweep
	for i=0L, n_elements(addtags_in)-1 do begin
		tag=addtags_in[i]
		if not tag_exist(sweep,tag) then begin
			add_arrval, tag, addtags
		endif
	endfor

	ntags = n_elements(addtags)
	if ntags ne 0 then begin

		select_tags = ['run','rerun','camcol','field','id',addtags]
		;fpobjc=sdss_readobj(run,camcol,rerun=rerun,select_tags=select_tags)
		fpobjc=sdss_readobjlist($
			sweep.run, sweep.camcol, sweep.field, sweep.id,$
			rerun=rerun, $
			select_tags=select_tags)

		wbad=where(sweep.id ne fpobjc.id or sweep.field ne fpobjc.field,nbad)
		if n_elements(fpobjc) ne n_elements(sweep) or nbad ne 0 then begin
			message,'object lists do not match'
		endif

		; now add any tags that matched from the fpobjc files
		st=sweep[0]

		for i=0L, n_elements(addtags)-1 do begin
			tag=addtags[i]
			; add the tag if found in the fpobjc and not already in
			; the sweep
			if tag_exist(fpobjc,tag,index=ind) $
				and not tag_exist(sweep,tag) then begin

				st = create_struct(st, tag, fpobjc[0].(ind))
			endif
		endfor

		; copy in the data: we allow the sweeps to rule here
		st = replicate(st, n_elements(sweep))
		copy_struct, fpobjc, st
		copy_struct, sweep, st

		return, st

	endif else begin
		return, sweep
	endelse

end







pro sweeps::create_pbs_bycamcol, proctype, procrun, command, $
        shortname=shortname, $
		addtags=addtags, $
		esidl_v=esidl_v, $
		sdssidl_v=sdssidl_v, $
		photoop_v=photoop_v, $
		idlutils_v=idlutils_v, $
		photo_resolve=photo_resolve, $
		photo_sweep=photo_sweep, $
        extra_setups=extra_setups, $
		pars=pars, $
		runs=runs, $
		ignore_resolve=ignore_resolve, $
		walltime=walltime, gdl=gdl, $
        $
        LBL=LBL


    if n_elements(shortname) eq 0 then shortname=proctype

    ; default to "current" installs
	if n_elements(photoop_v) eq 0 then photoop_v = ""
	if n_elements(idlutils_v) eq 0 then idlutils_v = ""

    ; my local installs
	if n_elements(esidl_v) eq 0 then esidl_v = "-r ~esheldon/exports/esidl-work"
	if n_elements(sdssidl_v) eq 0 then sdssidl_v = "-r ~esheldon/exports/sdssidl-work"


	add_arrval,"setup photoop "+photoop_v, setups
	add_arrval,"setup idlutils "+idlutils_v, setups
	add_arrval,"setup sdssidl "+sdssidl_v, setups
	add_arrval,"setup esidl "+esidl_v, setups

    if n_elements(extra_setups) ne 0 then begin
        add_arrval, extra_setups, setups
    endif
	
	if n_elements(photo_sweep) eq 0 then begin
		PHOTO_SWEEP=getenv('PHOTO_SWEEP')
	endif
	if n_elements(photo_resolve) eq 0 then begin
		PHOTO_RESOLVE=getenv('PHOTO_RESOLVE')
	endif
	if n_elements(photo_calib) eq 0 then begin
		PHOTO_CALIB=getenv('PHOTO_CALIB')
	endif

	add_arrval, 'export PHOTO_SWEEP='+photo_sweep, setups
	add_arrval, 'export PHOTO_RESOLVE='+photo_resolve, setups
	add_arrval, 'export PHOTO_CALIB='+photo_calib, setups


	;setups = strjoin(setups, ' && ')

	self->cache_runlist, /force, photo_resolve=photo_resolve
	if n_elements(runs) ne 0 then begin
		self->match_runlist, runs, tmpruns, reruns, photo_resolve=photo_resolve
	endif else begin
		self->runlist, runs, reruns, photo_resolve=photo_resolve
	endelse

	print,'Found: ',n_elements(runs),' runs',f='(a,i0,a)'


	home=getenv('HOME')
	pbs_dir=path_join(home,['pbs',proctype,procrun])
	spawn,'mkdir -pv '+pbs_dir

	qsub_file = path_join(pbs_dir, 'submit-'+proctype+'-bycamcol')
	fbase = proctype+'-'+self.type

	qsub_file+='.sh'
	openw, qsub_lun, qsub_file, /get_lun

	nrun=n_elements(runs)
	ntot=nrun*6
	ii=0L
	for i=0L, nrun-1 do begin
		run = runs[i]

		rstr = run2string(run)

		for camcol=1,6 do begin
			cstr=strn(camcol)

		    job_name = shortname+self.type+'-'+rstr+'-'+cstr

			pbs_file = repstr(job_name, proctype, fbase)+'.pbs'
			pbs_file=filepath(root=pbs_dir, pbs_file)


			idl_commands="sw=obj_new('sweeps','"+self.type+"'"
			if keyword_set(ignore_resolve) then begin
				idl_commands += ",/ignore_resolve"
			endif
			idl_commands += ")"


			proc_command = $
				string(f='(%"%s, %s, %s, %s, run=%s, camcol=%s")', $
				"    sw->process_runs", $
				"'"+proctype+"'", "'"+procrun+"'", "'"+command+"'", $
				strn(run), strn(camcol))


			if n_elements(pars) ne 0 then begin
				idl_commands = [idl_commands, 'pars='+tostring(pars)]
				proc_command += ', pars=pars'
			endif 
			if n_elements(addtags) ne 0 then begin
				idl_commands = [idl_commands, 'addtags='+tostring(addtags)]
				proc_command += ', addtags=addtags'
			endif

			idl_commands=[idl_commands,proc_command]

            if keyword_set(LBL) then begin
                pbs_riemann_idl, $
                    pbs_file, idl_commands, setup=setups, job_name=job_name, $
                    walltime=walltime, gdl=gdl
            endif else begin
                pbs_bnl_idl, $
                    pbs_file, idl_commands, setup=setups, job_name=job_name, $
                    walltime=walltime, gdl=gdl
            endelse


			printf, qsub_lun, $
				'echo -n "',ii+1,'/',ntot,' ',pbs_file,' "',$
				format='(a,i0,a,i0,a,a,a)'
			printf, qsub_lun, 'qsub '+pbs_file
			ii=ii+1
		endfor
	endfor


	free_lun, qsub_lun

end


pro sweeps::create_pbs_byrun, proctype, procrun, command, $
        shortname=shortname, $
		addtags=addtags, $
		esidl_v=esidl_v, $
		sdssidl_v=sdssidl_v, $
		photoop_v=photoop_v, $
		idlutils_v=idlutils_v, $
		photo_resolve=photo_resolve, $
		photo_sweep=photo_sweep, $
        extra_setups=extra_setups, $
		pars=pars, $
		runs=runs, $
		ignore_resolve=ignore_resolve, $
		where_string=where_string, gdl=gdl, $
        $
        LBL=LBL, $
        walltime=walltime
	
    if n_elements(shortname) eq 0 then shortname=proctype

    ; default to "current" installs
	if n_elements(photoop_v) eq 0 then photoop_v = ""
	if n_elements(idlutils_v) eq 0 then idlutils_v = ""

    ; my local installs
	if n_elements(esidl_v) eq 0 then esidl_v = "-r ~esheldon/exports/esidl-work"
	if n_elements(sdssidl_v) eq 0 then sdssidl_v = "-r ~esheldon/exports/sdssidl-work"


	add_arrval,"setup photoop "+photoop_v, setups
	add_arrval,"setup idlutils "+idlutils_v, setups
	add_arrval,"setup sdssidl "+sdssidl_v, setups
	add_arrval,"setup esidl "+esidl_v, setups

    if n_elements(extra_setups) ne 0 then begin
        add_arrval, extra_setups, setups
    endif

	if n_elements(photo_sweep) eq 0 then begin
		PHOTO_SWEEP=getenv('PHOTO_SWEEP')
	endif
	if n_elements(photo_resolve) eq 0 then begin
		PHOTO_RESOLVE=getenv('PHOTO_RESOLVE')
	endif
	if n_elements(photo_calib) eq 0 then begin
		PHOTO_CALIB=getenv('PHOTO_CALIB')
	endif

	add_arrval, 'export PHOTO_SWEEP='+photo_sweep, setups
	add_arrval, 'export PHOTO_RESOLVE='+photo_resolve, setups
	add_arrval, 'export PHOTO_CALIB='+photo_calib, setups


	

	;setups = strjoin(setups, ' && ')

	self->cache_runlist, /force
	if n_elements(runs) ne 0 then begin
		self->match_runlist, runs, tmpruns, reruns, photo_resolve=photo_resolve
	endif else begin
		self->runlist, runs, reruns, photo_resolve=photo_resolve
	endelse

	print,'Found: ',n_elements(runs),' runs',f='(a,i0,a)'





	home=getenv('HOME')
	pbs_dir=path_join(home,['pbs',proctype,procrun])
	;file_mkdir, pbs_dir
	spawn,'mkdir -pv '+pbs_dir




	

	qsub_file = path_join(pbs_dir, 'submit-'+proctype+'-byrun')
	fbase = proctype+'-'+self.type



	combine_file = path_join(pbs_dir, 'combine-byrun.sh')
	if n_elements(extra_name) ne 0 then begin
		combine_file = repstr(combine_file, '.sh', '-'+extra_name+'.sh')
		extra_gather = ', extra_name="'+extra_name+'"'
	endif else begin
		extra_gather = ''
	endelse
	openw, lun, combine_file, /get_lun

	printf, lun, 'idl<<EOF'
	printf, lun, "  setenv,'PHOTO_SWEEP="+photo_sweep+"'"
	printf, lun, "  setenv,'PHOTO_RESOLVE="+photo_resolve+"'"
	printf, lun, '  runs='+tostring(runs)
	printf, lun, '  sw=obj_new("sweeps","'+self.type+'")'
	combine_command = string($
		"  sw->gather2file, '",proctype,"','",procrun,"'",+$
		extra_gather, $
		f='(a,a,a,a,a,a)' )
	;if subset then begin
	;	combine_command += ', runs=runs'
	;endif
	if n_elements(where_string) ne 0 then begin
		combine_command += ", where_string='"+where_string+"'"
	endif
	printf, lun, combine_command
	printf, lun, 'EOF'

	free_lun, lun

	if keyword_set(onlycombine) then begin
		print,'Just writing combine script'
		return
	endif


	qsub_file+='.sh'
	openw, qsub_lun, qsub_file, /get_lun


	nrun=n_elements(runs)
	ii=0L
	for i=0L, nrun-1 do begin
		run = runs[i]

		rstr_full = run2string(run)
		rstr = string(run,f='(i0)')

		job_name = shortname+self.type+'-'+rstr
		job_name_full = proctype+'-'+rstr_full

		pbs_file = repstr(job_name_full, proctype, fbase)+'.pbs'
		pbs_file=filepath(root=pbs_dir, pbs_file)


		idl_commands="sw=obj_new('sweeps','"+self.type+"'"
		if keyword_set(ignore_resolve) then begin
			idl_commands += ",/ignore_resolve"
		endif
		idl_commands += ")"


		if n_elements(pars) ne 0 then begin
			idl_commands = [idl_commands, 'pars='+tostring(pars)]
		endif 

		if n_elements(addtags) ne 0 then begin
			idl_commands = [idl_commands, 'addtags='+tostring(addtags)]
		endif



		for camcol=1,6 do begin

			proc_command = $
				string(f='(%"%s, %s, %s, %s, run=%s, camcol=%s")', $
				"    sw->process_runs", $
				"'"+proctype+"'", "'"+procrun+"'", "'"+command+"'", $
				strn(run), strn(camcol))

			if n_elements(pars) ne 0 then begin
				proc_command += ', pars=pars'
			endif
			if n_elements(addtags) ne 0 then begin
				proc_command += ', addtags=addtags'
			endif
			

			idl_commands=[idl_commands,proc_command]
		endfor

        if keyword_set(LBL) then begin
            pbs_riemann_idl, $
                pbs_file, idl_commands, setup=setups, job_name=job_name, gdl=gdl, $
                walltime=walltime
        endif else begin
            pbs_bnl_idl, $
                pbs_file, idl_commands, setup=setups, job_name=job_name, gdl=gdl, $
                walltime=walltime
        endelse


		printf, qsub_lun, $
			'echo -n "',ii+1,'/',nrun,' ',pbs_file,' "',$
			format='(a,i0,a,i0,a,a,a)'
		printf, qsub_lun, 'qsub '+pbs_file
		ii=ii+1
	endfor


	free_lun, qsub_lun

end







pro sweeps::create_pbs, proctype, procrun, command, $
        shortname=shortname, $
		addtags=addtags, $
		runs=runs, $
		esidl_v=esidl_v, $
		sdssidl_v=sdssidl_v, $
		photoop_v=photoop_v, $
		idlutils_v=idlutils_v, $
		photo_resolve=photo_resolve, $
		photo_sweep=photo_sweep, $
        extra_setups=extra_setups, $
		$
		pars=pars, $
		reselect=reselect, extra_name=extra_name, $
		$
		doproc=doproc,$
		dogather=dogather,$
		collate=collate, $
		nper=nper, $
		fpobjc=fpobjc, $
		match_method=match_method, $
		noverify=noverify, $
		where_string=where_string, $
		$
		walltime=walltime, gdl=gdl, $
        $
        LBL=LBL


    if n_elements(shortname) eq 0 then shortname=proctype

	if n_elements(doproc) eq 0 and n_elements(dogather) eq 0 then begin
		doproc=1
		dogather=1
	endif
	if n_elements(nper) eq 0 then begin
		nper=2
	endif

	if n_elements(runs) ne 0 then begin
		runstring=tostring(runs)
	endif

	if n_elements(fpobjc) eq 0 then fpobjc=0




    ; default to "current" installs
	if n_elements(photoop_v) eq 0 then photoop_v = ""
	if n_elements(idlutils_v) eq 0 then idlutils_v = ""

    ; my local installs
	if n_elements(esidl_v) eq 0 then esidl_v = "-r ~esheldon/exports/esidl-work"
	if n_elements(sdssidl_v) eq 0 then sdssidl_v = "-r ~esheldon/exports/sdssidl-work"


	add_arrval,"setup photoop "+photoop_v, setups
	add_arrval,"setup idlutils "+idlutils_v, setups
	add_arrval,"setup sdssidl "+sdssidl_v, setups
	add_arrval,"setup esidl "+esidl_v, setups

    if n_elements(extra_setups) ne 0 then begin
        add_arrval, extra_setups, setups
    endif

	if n_elements(photo_sweep) eq 0 then begin
		PHOTO_SWEEP=getenv('PHOTO_SWEEP')
	endif
	if n_elements(photo_resolve) eq 0 then begin
		PHOTO_RESOLVE=getenv('PHOTO_RESOLVE')
	endif
	if n_elements(photo_calib) eq 0 then begin
		PHOTO_CALIB=getenv('PHOTO_CALIB')
	endif

	add_arrval, 'export PHOTO_SWEEP='+photo_sweep, setups
	add_arrval, 'export PHOTO_RESOLVE='+photo_resolve, setups
	add_arrval, 'export PHOTO_CALIB='+photo_calib, setups

	;setups = strjoin(setups, ' && ')

	sweep_old=getenv('PHOTO_SWEEP')
	resolve_old=getenv('PHOTO_RESOLVE')
	setenv, 'PHOTO_SWEEP='+PHOTO_SWEEP
	setenv, 'PHOTO_RESOLVE='+PHOTO_RESOLVE

	self->split_runlist, runptrs, rerunptrs, nper=nper, /force, $
		runs=runs
	njobs = n_elements(runptrs)
	ptr_free, runptrs, rerunptrs

	setenv, 'PHOTO_SWEEP='+sweep_old
	setenv, 'PHOTO_RESOLVE='+resolve_old



	home=getenv('HOME')
	pbs_dir=path_join(home,['pbs',proctype,procrun])
	spawn,'mkdir -pv '+pbs_dir

	qsub_file = path_join(pbs_dir, 'submit-'+proctype)
	check_file = path_join(pbs_dir, 'check.sh')
	combine_file = path_join(pbs_dir, 'combine.sh')
	if n_elements(extra_name) ne 0 then begin
		combine_file = repstr(combine_file, '.sh', '-'+extra_name+'.sh')
	endif

	fbase = proctype+'-'+self.type

	if fpobjc then begin
		qsub_file+='-fpobjc'
		fbase+='-fpobjc'
	endif

	if n_elements(match_method) ne 0 then begin
		qsub_file += '-match-gather'
		fbase += '-match-gather'
	endif else begin
		if keyword_set(doproc) then begin
			qsub_file += '-proc'
			fbase += '-proc'
		endif
		if keyword_set(dogather) then begin
			qsub_file += '-gather'
			fbase += '-gather'
		endif
	endelse


	if n_elements(extra_name) ne 0 then begin
		qsub_file += '-'+extra_name
		fbase += '-'+extra_name
		extra_gather = ', extra_name="'+extra_name+'"'
	endif else begin
		extra_gather=''
	endelse

	qsub_file+='.sh'
	print,'writing ',qsub_file
	openw, lun, qsub_file, /get_lun


	numstr = strn(njobs-1,len=3,padchar='0')
    printf, lun
    printf, lun, 'for i in `seq -w 0 '+numstr+'`; do' 
    printf, lun, '    file='+fbase+'-${i}.pbs'
    printf, lun, '    echo -n "qsub $file  "'
    printf, lun, '    qsub $file'
    printf, lun, 'done'
	free_lun, lun

	openw, lun, check_file, /get_lun
	printf,lun,'for f in *.pbs; do'
	printf,lun,'    if [ ! -e "$f.log" ]; then'
	printf,lun,'        echo $f'
	printf,lun,'    fi'
	printf,lun,'done'
	free_lun, lun

	print,'writing ',combine_file
	openw, lun, combine_file, /get_lun
	printf, lun, 'idl<<EOF'
	printf, lun, "  setenv,'PHOTO_SWEEP="+photo_sweep+"'"
	printf, lun, "  setenv,'PHOTO_RESOLVE="+photo_resolve+"'"
	printf, lun, '  sw=obj_new("sweeps","'+self.type+'")'
	mess=string($
		"  sw->gather_partial, '",proctype,"','",procrun,"'",+$
		", /combine",extra_gather, $
		f='(a,a,a,a,a,a,a)')
	if n_elements(runs) ne 0 then begin
		printf, lun, '  runs='+runstring
		mess += ', runs=runs'
	endif
	printf, lun, mess
	printf, lun, 'EOF'
	free_lun, lun



	for job=0L, njobs-1 do begin

		jobstr = strn(job,len=3,padchar='0')

		jfbase = fbase + '-'+jobstr
		
		pbs_file=filepath(root=pbs_dir, jfbase+'.pbs')

		proc_command = $
			string(f='(%"%s, %s, %s, %s, %s, nper=%s")', $
			"    sw->process_partial_runs", $
			"'"+proctype+"'", "'"+procrun+"'", "'"+command+"'", $
			strn(job), strn(nper))

		gather_command = string(f='(%"%s, %s, %s, %s, nper=%s")', $
			"    sw->gather_partial", $
			"'"+proctype+"'", "'"+procrun+"'", strn(job), strn(nper))

		if keyword_set(fpobjc) then begin
			proc_command += ", /fpobjc"
			gather_command += ", /fpobjc"
		endif
		if n_elements(match_method) ne 0 then begin
			proc_command += ", match_method='"+match_method+"'"
			gather_command += ", match_method='"+match_method+"'"
		endif 

		idl_commands = "sw=obj_new('sweeps')"
		idl_commands = "sw=obj_new('sweeps','"+self.type+"')"


		if n_elements(runs) ne 0 then begin
			idl_commands = [idl_commands, "runs="+runstring]
		endif


		if keyword_set(reselect) then begin
			gather_command += ", /reselect"
		endif
		if keyword_set(noverify) then begin
			gather_command += ", /noverify"
		endif

		if n_elements(where_string) ne 0 then begin
			gather_command += ", where_string='"+where_string+"'"
		endif


		if n_elements(pars) ne 0 then begin
			idl_commands=[idl_commands, 'pars='+tostring(pars)]
			proc_command += ", pars=pars"
			gather_command += ", pars=pars"
		endif 

		if n_elements(addtags) ne 0 then begin
			idl_commands = [idl_commands, 'addtags='+tostring(addtags)]
			proc_command += ', addtags=addtags'
			gather_command += ', addtags=addtags'
		endif

		if n_elements(extra_name) ne 0 then begin
			gather_command += ", extra_name='"+extra_name+"'"
		endif
			
		if n_elements(runs) ne 0 then begin
			proc_command += ', runs=runs'
			gather_command += ', runs=runs'
		endif

		if keyword_set(collate) then begin
			gather_command += ', /collate'
		endif

		if keyword_set(doproc) then begin
			idl_commands=[idl_commands,proc_command]
		endif
		if keyword_set(dogather) then begin
			idl_commands=[idl_commands,gather_command]
		endif

		job_name = shortname+self.type+"-"+jobstr

        if keyword_set(LBL) then begin
            pbs_riemann_idl, pbs_file, idl_commands, $
                setup=setups, job_name=job_name, $
                walltime=walltime, gdl=gdl
        endif else begin
            pbs_bnl_idl, pbs_file, idl_commands, $
                setup=setups, job_name=job_name, $
                walltime=walltime, gdl=gdl
        endelse
	endfor
end


pro sweeps::create_pbs_photoz_input_pzgal, gdl=gdl

	procrun='m01'

	command = $
	  'pz=obj_new("photoz_uchicago") & res=pz->sweep_input_select(objs)'
	self->create_pbs, 'pzgal', procrun, command, $
		pars=pars, $
		$
		idlutils_v=idlutils_v, $
		photoop_v=photoop_v, $
		$
		photo_sweep=photo_sweep, $
		photo_resolve=photo_resolve, $
		gdl=gdl

	return

end

; this is old, don't use
pro sweeps::create_pbs_admom_sweep, gdl=gdl, dogather=dogather

	procrun='am01'

	command = $
	  'as=obj_new("admom_sweep") & res=as->run_regauss(objs)'

	; on riemann we use the installed versions
	idlutils_v='v5_4_11'
	photoop_v='v1_9_4'

	; we require a row for each result even if not processed
	;pars = {require_result: 1}
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'am', procrun, command, $
			addtags=['m_rr_cc','rowc'], $
			pars=pars, $
			$
			idlutils_v=idlutils_v, $
			photoop_v=photoop_v, $
			$
			photo_sweep=photo_sweep, $
			photo_resolve=photo_resolve, $
			gdl=gdl
	endif else begin

		where_string='str.whyflag_rg[2] eq 0 or str.whyflag_rg[3] eq 0'
		nper=20
		self->create_pbs, 'am', procrun, command, $
			addtags=['m_rr_cc','rowc'], $
			pars=pars, $
			nper=nper, $
			$
			idlutils_v=idlutils_v, $
			photoop_v=photoop_v, $
			$
			photo_sweep=photo_sweep, $
			photo_resolve=photo_resolve, $
			gdl=gdl, $
			/collate, $
			where_string=where_string, $
			/dogather
	endelse

	return

end






function sweeps_default_pars

	pars = { $
		sweeps_parstruct, $
		status_ext:1, $
		result_ext:2, $
		require_result:0 $
	}

	return, pars
end



pro sweeps::set_default_pars
	pars=sweeps_default_pars()
	self.pars = pars
end

pro sweeps::copy_extra_pars, pars=pars, _extra=_extra

	if n_tags(_extra) ne 0 then begin
		tmp = self.pars
		struct_assign, _extra, tmp, /nozero
		self.pars = tmp
	endif
	if n_tags(pars) ne 0 then begin
		tmp=self.pars
		struct_assign, pars, tmp, /nozero
		self.pars=tmp
	endif

end


pro sweeps::test_create, run
    
    setenv,'PHOTO_CALIB=/clusterfs/riemann/raid006/bosswork/groups/boss/calib/dr8_final'
    setenv,'PHOTO_RESOLVE=/clusterfs/riemann/raid006/bosswork/groups/boss/resolve/2010-05-23'
    setenv,'PHOTO_REDUX=/clusterfs/riemann/raid006/dr8/groups/boss/photo/redux'

    setenv,'PHOTO_SWEEP=/clusterfs/riemann/raid008/bosswork/groups/boss/target/sweeps/dr8_final'

    ;readcol,'finalruns_301.dat', run, rerun, format='L, L' 
    ;if (n_elements(imin) GT 0) then run = run[imin:imax]

    ;for ii = 0L, n_elements(run)-1L do begin
        datasweep, run, rerun=301, /uber, /combine, maxmem=5000.0
    ;endfor

end




pro sweeps__define

	defpars = sweeps_default_pars()
	struct = {$
		sweeps, $
        type: '', $
		pars: defpars, $
		dummy: 0 $
	}
end
