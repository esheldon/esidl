function photoz_uchicago::init, type=type

	if n_elements(type) eq 0 then begin 
		;message,'type not initialized.  Setting to "dr7pofz"',/inf
		self.type='dr7pofz'
	endif else begin 
		self.type=type
	endelse 

	return,1

end 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Input files, finding, reading, creating
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function photoz_uchicago::basedir
	basedir = getenv('PHOTOZ_DIR')
	if basedir eq '' then message,'$PHOTOZ_DIR is not set'
	return, basedir
end

function photoz_uchicago::photoz_dir
	basedir = self->basedir()
	return, path_join(basedir, ['uchicago',self.type])
end
function photoz_uchicago::input_dir
	pzdir = self->photoz_dir()
	return, path_join(pzdir, 'input')
end 
function photoz_uchicago::output_dir
	pzdir = self->photoz_dir()
	return, path_join(pzdir, 'output')
end 

function photoz_uchicago::training_dir, original=original
	basedir = self->basedir()
	sub=['uchicago','training']
	if keyword_set(original) then begin
		sub=[sub, 'original']
	endif else begin
		sub=[sub, 'matched']
	endelse
	dir = path_join(basedir, sub)
	return, dir
end 
function photoz_uchicago::training_file, type, original=original
	dir=self->training_dir(original=original)
	if keyword_set(original) then ext='.fits' else ext='.st'
	f=path_join(dir, type+ext)
	return, f
end

function photoz_uchicago::training_read, type, original=original
	f=self->training_file(type, original=original)
	if keyword_set(original) then begin
		t=read_idlstruct(f)
	endif else begin
		t=mrdfits(f,1)
	endelse
	return, t
end



function photoz_uchicago::input_file, runs, reruns=reruns, count=count

	if n_elements(runs) eq 0 then begin 
		print,'-Syntax: files = obj->input_file(runs, reruns=, count=)'
		return,''
	endif 

	sf = obj_new('sdss_files')

	nruns = n_elements(runs)
	nreruns = n_elements(reruns)

	if nreruns ne 0 then begin 
		if nruns ne nreruns then message,'runs/reruns must be same length'
	endif else begin 
		reruns = sf->rerun(runs)
	endelse 

	rstr = sf->run2string(runs)
	rrstr = ntostr(reruns)

	dir = self->input_dir()
	files = 'photoz-input-'+rstr+'-'+rrstr+'.st'
	count = n_elements(files)

	return, concat_dir(dir,files)

end 

function photoz_uchicago::input_structdef_old

  st = create_struct('run', 0L, $
                     'rerun', 0, $
                     'camcol', 0, $
                     'field', 0, $
                     'id', 0L, $
                     'ra', 0d, $
                     'dec', 0d, $
                     'primtarget', 0L, $
                     'objc_type', 0L, $
                     'counts_model', fltarr(5), $
                     'counts_modelerr', fltarr(5), $
                     'objc_prob_psf', 0.0)

  return,st

end 


PRO photoz_uchicago::runs_reruns, runset, runs, reruns

  ;; 5042 just has empty files, how stupid can you get
  w=where(!run_status.tsobj_photo_v ge 5.4 AND $
          !run_status.run NE 5042,nruns)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Split into sets for running on different processors
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  CASE runset OF
      1: BEGIN 
          nruns = 100
          ind = lindgen(100)
      END 
      2: BEGIN 
          nruns = 100
          ind = lindgen(100)+100
      END 
      3: BEGIN 
          nruns = nruns-200
          ind = lindgen(nruns) + 200
      END 
      ELSE: message,'Unknown runset '+ntostr(runset)
  ENDCASE 

  w=w[ind]

  runs = !run_status[w].run
  reruns = !run_status[w].rerun

END 



PRO photoz_uchicago::input_write, runset, overwrite=overwrite

;  IF n_params() EQ 0 THEN BEGIN 
;      print,'-Syntax: obj->input_write, runset, /overwrite'
;      return
;  ENDIF 

  tm=systime(1)

;  self->runs_reruns, runset, runs, reruns

  readcol,'~/tmp/getruns_original.dat',runs,reruns,format='I,I'
  s = sort(runs)
  runs = runs[s]
  reruns = reruns[s]
  nruns = n_elements(runs)
  colprint, runs, reruns

  taglist = ['run','rerun','camcol','field','id',$
             'ra','dec', $
             'primtarget', $
             'objc_type', $
             'counts_model', $
             'counts_modelerr', $
             'reddening', $
             'objc_prob_psf']


  ;; This works off adatc stuff, need to convert,but don't have
  ;; some of these tags in the database!
  wstring = $
    "where(lnew.counts_model[2] NE -9999 AND " + $
    "      lnew.objc_prob_psf GE 0.0 AND " + $
    "      lnew.objc_prob_psf LT (1.-!hardprobcut) AND " + $
    "      ( " + $
    "        (lnew.m_r_h[1] GT 0.0 AND lnew.m_r_h[1] LT !hardrcut) OR " + $
    "        (lnew.m_r_h[2] GT 0.0 AND lnew.m_r_h[2] LT !hardrcut) OR " + $
    "        (lnew.m_r_h[3] GT 0.0 AND lnew.m_r_h[3] LT !hardrcut)    " + $
    "      ) AND "+$
    "      (lnew.cmodel_counts[2]-lnew.reddening[2]) lt 22.0 AND " + $
    "      (lnew.cmodel_counts[2]-lnew.reddening[2]) gt 14.0 AND " + $
    "      (lnew.counts_model[2]-lnew.reddening[2]) lt 22.5 AND " + $
    "      (lnew.counts_model[2]-lnew.reddening[2]) gt 13.0 " + $
    "      )"

  structdef = self->input_structdef()

  outdir = self->input_dir()
  outdir = concat_dir(outdir,'new')

  FOR i=0L, nruns-1 DO BEGIN 

      print
      print,'Processing run #'+ntostr(i+1)+'/'+ntostr(nruns)

      run = runs[i]
      rerun = reruns[i]

      rstr = run2string(run)
      rrstr = ntostr(rerun)
      outFile = 'photoz-input-'+rstr+'-'+rrstr+'.st'
      outFile = concat_dir(outdir, outfile)

      print
      print,'Writing file: ',outFile

      ;; Don't overwrite files unless asked to.  This is so we
      ;; can rerun the script when one run crashes without 
      ;; redoing all runs
      IF keyword_set(overwrite) OR (NOT fexist(outFile)) THEN BEGIN 

          file_delete, outFile, /quiet

          IF fexist(outFile) AND keyword_set(overwrite) THEN BEGIN 
              print,'Overwriting existing file: ',outFile
          ENDIF 
          
          ;; Loop over camcols
          FOR camcol=1,6 DO BEGIN 
              
              read_tsobj, $
                [run,rerun,camcol], str, taglist=taglist, /all, /corr, $
                wstring = wstring, verbose=1
              
              nstruct = n_elements(str)
              
              ;; Anything found?
              IF nstruct NE 0 THEN BEGIN 
                  
                  str.counts_model = str.counts_model - str.reddening
                  
                  outStruct = replicate(structdef, nstruct)
                  copy_struct, str, outStruct

                  write_idlstruct, outStruct, outFile, /ascii, /append

              ENDIF 
              
          ENDFOR 
          

      ENDIF ELSE BEGIN
          print,'File '+outFile+' exists. Skipping'
      ENDELSE 

  ENDFOR 

  print
  ptime,systime(1)-tm

END 


FUNCTION photoz_uchicago::input_read, runs, reruns=reruns, count=count

  IF n_elements(runs) EQ 0 THEN BEGIN 
      print,'-Syntax: files = obj->input_file(runs, reruns=, count=)'
      return,-1
  ENDIF 

  files = self->input_file(runs, reruns=reruns)
  return, read_idlstruct_multi(files, count=count)
  
END 

function photoz_uchicago::calib_logic, objs, flags
	flagval = sdss_flagval('calib_status','photometric')
	calib_logic = $
		(objs.calib_status[0] and flagval) ne 0 $
		and (objs.calib_status[1] and flagval) ne 0 $
		and (objs.calib_status[2] and flagval) ne 0 $
		and (objs.calib_status[2] and flagval) ne 0 $
		and (objs.calib_status[2] and flagval) ne 0
	return,calib_logic
end

function photoz_uchicago::colorflags_logic, flags, objflagtype, flagname, $
		notset=notset
	flagval = sdss_flagval(objflagtype, flagname)
	if not keyword_set(notset) then begin
		check = $
			(flags[0,*] and flagval) ne 0 $
			or (flags[1,*] and flagval) ne 0 $
			or (flags[2,*] and flagval) ne 0 $
			or (flags[3,*] and flagval) ne 0 $
			or (flags[4,*] and flagval) ne 0
	endif else begin
		; checking if it is not set instead if it is set
		check = $
			(flags[0,*] and flagval) eq 0 $
			or (flags[1,*] and flagval) eq 0 $
			or (flags[2,*] and flagval) eq 0 $
			or (flags[3,*] and flagval) eq 0 $
			or (flags[4,*] and flagval) eq 0
	endelse
	return, reform(check)
end


function photoz_uchicago::deredden_flux, flux, extinction
	exponent = 0.4*extinction
	flux_correct = (10.0^exponent)*flux
	return, flux_correct
end

function photoz_uchicago::deredden_ivar, ivar, extinction
	exponent = -0.8*extinction
	ivar_correct = (10.0^exponent)*ivar
	return, ivar_correct
end

function photoz_uchicago::make_cmodelmag, calibobj, cmodelflux=cmodelflux
	
	devflux = calibobj.devflux
	expflux = calibobj.expflux
	fracpsf = calibobj.fracpsf

	cmodelflux = devflux*fracpsf + expflux*(1.0-fracpsf)
	cmodelmag = 22.5-2.5*alog10(cmodelflux > 0.001)

	return, cmodelmag
end



function photoz_uchicago::derived_struct, objs
	st = { $
		modelflux_dered: fltarr(5), $
		modelflux_dered_ivar: fltarr(5), $
		modelmag_dered:fltarr(5), $
		cmodelmag_dered:fltarr(5) $
	}

	st = replicate(st, n_elements(objs))

	; Extract the relevant magnitudes
	; We really should have used cmodelmags, and will later.
	st.modelmag_dered = 22.5 - 2.5*alog10(objs.modelflux > 0.001)
	st.cmodelmag_dered = self->make_cmodelmag(objs)

	; Extinction correct
	st.modelmag_dered = st.modelmag_dered - objs.extinction
	st.cmodelmag_dered = st.cmodelmag_dered - objs.extinction

	st.modelflux_dered = $
		self->deredden_flux(objs.modelflux, objs.extinction)
	st.modelflux_dered_ivar = $
		self->deredden_ivar(objs.modelflux_ivar, objs.extinction)
	return, st

end



function photoz_uchicago::flag_logic, calibobj


	binned1 = sdss_flagval('object1','binned1')
	binned2 = sdss_flagval('object1','binned2')
	binned4 = sdss_flagval('object1','binned4')

	r_binned_logic =  $
		(calibobj.flags[2] and binned1) ne 0 or $
		(calibobj.flags[2] and binned2) ne 0 or $
		(calibobj.flags[2] and binned4) ne 0

	i_binned_logic =  $
		(calibobj.flags[3] and binned1) ne 0 or $
		(calibobj.flags[3] and binned2) ne 0 or $
		(calibobj.flags[3] and binned4) ne 0

	splog,'Cutting to binned in r,i'
	logic = (r_binned_logic and i_binned_logic)
    help,where(logic)


	splog,'Cutting to photometric'
	;notcalib_logic = self->colorflags_logic(calibobj.calib_status,$
	;	    'calib_status','photometric', /notset)
	calib_logic = self->calib_logic(calibobj)
	logic = logic and calib_logic
    help,where(logic)

	; flags that must not be set
	satur = sdss_flagval('object1','satur')
	satur_center = sdss_flagval('object2','satur_center')

	blended = sdss_flagval('object1','blended')
	nodeblend = sdss_flagval('object1','nodeblend')

	oflags = calibobj.objc_flags

	s_logic = (oflags and satur)
    ;sc_logic = (oflags and satur_center)
	;s_or_sc_logic = ( s_logic eq 0 or (s_logic ne 0 and sc_logic eq 0) )

    splog,'cutting satur'
	logic = logic and (s_logic eq 0)
    help,where(logic)


	; famously worded as double negatives
    splog,'cutting bad deblend'
	bdb_logic = $
		((oflags and blended) eq 0) or ((oflags and nodeblend) ne 0)

    logic = logic and bdb_logic
    help,where(logic)

    splog,'Cutting bright'
	bright = sdss_flagval('object1','bright')
	logic = logic and (oflags and bright) eq 0
    help,where(logic)

    splog,'too many peaks'
	too_many_peaks = sdss_flagval('object1','deblend_too_many_peaks')
	logic = logic and (oflags and too_many_peaks) eq 0



	; new stuff to try
    splog,'cutting peakcenter'
	peakcenter = sdss_flagval('object1','peakcenter')
	logic = logic and  (oflags and peakcenter) eq 0
    help,where(logic)

    splog,'cutting notchecked'
	notchecked = sdss_flagval('object1','notchecked')
	logic = logic and (oflags and notchecked) eq 0
    help,where(logic)

    splog,'cutting noprofile'
	noprofile = sdss_flagval('object1','noprofile')
	logic = logic and (oflags and noprofile) eq 0
    help,where(logic)


	; default is cut all saturated since 1/3 of those left over after
	; rdev and tycho cuts were poorly deblended and total was only 1%


	return, logic

end

function photoz_uchicago::resolve_logic, objs
	primary_flag = sdss_flagval('resolve_status','survey_primary')
	return, (objs.resolve_status and primary_flag) ne 0
end

function photoz_uchicago::input_structdef, num

	st = {$
		run:0L, rerun:0,camcol:0,field:0,id:0L, $
		ra:0d, dec:0d, $
		modelflux_dered: fltarr(5), $
		modelflux_dered_ivar: fltarr(5) }

	if n_elements(num) ne 0 then begin
		st=replicate(st, num)
	endif
	return,st

end 



function photoz_uchicago::sweep_input_select, objs, nkeep
	; select objects from datasweep file for inputs to the photoz codes

	mst = self->derived_struct(objs)

	mag_logic = mst.cmodelmag_dered[2] lt 22
	flag_logic = self->flag_logic(objs)
	resolve_logic = self->resolve_logic(objs)

	wmag=where(mag_logic)
	wflag = where(flag_logic)
	wresolve = where(resolve_logic)
	help,objs,wmag,wflag, wresolve

	keep= where( $
		mag_logic $	
		and flag_logic $
		and resolve_logic, nkeep)
	help,keep

	if nkeep eq 0 then begin
		return, -1
	endif else begin
		st = self->input_structdef(nkeep)
		struct_assign, objs[keep], st, /nozero

		st.modelflux_dered = mst[keep].modelflux_dered
		st.modelflux_dered_ivar = mst[keep].modelflux_dered_ivar
		return, st
	endelse

end


function photoz_uchicago::match_struct, base_struct, n

	arr1 = replicate(-9999.0, 5)
	arr2 = fltarr(5)
	arr3 = replicate(9999.0, 5)
	newstruct={$
		matchid: -9999L, $
		modelflux:arr1, $
		modelflux_ivar: arr2, $
		modelmag: arr1, $
		modelmag_err: arr3 $
	}
	outst = create_struct(base_struct, newstruct)

	if n_elements(n) ne 0 then begin
		outst=replicate(outst, n)	
	endif
	return, outst
end

pro photoz_uchicago::match2training, phot=phot, htmid_phot=htmid_phot, revphot=revphot
	types = ['deep2','other','sdssobjids','vvds','zcosmos']
	nt=n_elements(types)
	sw=obj_new('sweeps')

	photfile=file_basename( sw->output_file('pzgal','m01',/gather) )
	if n_elements(phot) eq 0 then begin
		print,'reading photometric catalog'
		phot=sw->read_gather('pzgal','m01')
	endif

	for i=0L, nt-1 do begin
		type=types[i]
		tf= self->training_file(type,/original)
		matchf=self->training_file(type)

		print,'reading training set: ',tf
		t=mrdfits(tf,1)

		outst=self->match_struct(t[0], n_elements(t))
		struct_assign, t, outst, /nozero

		print,'Matching to photometric catalog'
		angle = 2d/3600d*!dpi/180d
		htm_match, t.ra, t.dec, phot.ra, phot.dec, angle, $
			htmid2=htmid_phot, htmrev2=revphot, $
			mt, mphot, maxmatch=1

		help,mphot

		outst[mt].matchid = mphot
		outst[mt].modelflux = phot[mphot].modelflux_dered
		outst[mt].modelflux_ivar = phot[mphot].modelflux_dered_ivar

		outst[mt].modelmag = nmgy2mag(outst[mt].modelflux, $
			                          ivar=outst[mt].modelflux_ivar, $
									  err=err)
		outst[mt].modelmag_err=err

		h=['']
		sxaddpar, h, 'photfile', photfile
		hdr={photfile:photfile}
		print,'Outputting matched catalog: ',matchf
		;mwrfits, outst, matchf, h, /create
		write_idlstruct, outst, matchf, hdr=hdr, /csv
	endfor
end

pro photoz_uchicago::plotmatch, phot=phot
	types = ['deep2','other','sdssobjids','vvds','zcosmos']
	nt=n_elements(types)
	sw=obj_new('sweeps')

	if n_elements(phot) eq 0 then begin
		print,'reading photometric catalog'
		phot=sw->read_gather('pzgal','m01')
	endif

	for i=0L, nt-1 do begin

		type=types[i]
		tf= self->training_file(type)

		print,'reading training set: ',tf
		t=read_idlstruct(tf)

		wm=where(t.matchid ge 0)
		psfile=repstr(tf,'.fits','.eps')

		begplot,psfile,/encap,/color,xsize=11,ysize=8.5

		pplot, t.ra, t.dec, psym=3,/ynozero, $
			xtitle='RA', ytitle='DEC'

		pplot, /over, phot.ra, phot.dec, psym=3
		pplot, /over, t.ra, t.dec, psym=8, symsize=0.5, color='red'
		pplot, /over, t[wm].ra, t[wm].dec, psym=8, symsize=0.5, color='yellow'

		endplot,/trim,/png

	endfor
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output files, finding, reading, creating
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function photoz_uchicago::output_ext, dat=dat
    case self.type of
        'nn': begin
            if keyword_set(dat) then return, '.dat' else return, '.st'
        end
        'dr6cc2': begin
            if keyword_set(dat) then return, '.dat' else return, '.st'
        end
        else: message,'dont yet support type '+self.type
    endcase
end
function photoz_uchicago::output_file, numbers, all=all, dat=dat, big=big, count=count


    ext = self->output_ext(dat=dat)
    case self.type of
        'nn': begin

            dir = self->output_dir()

            if keyword_set(big) then begin
                return, concat_dir(dir, 'DR5_final.st')
            endif else if not keyword_set(all) then begin
                if n_elements(numbers) eq 0 then begin 
                    print,'-Syntax: files = obj->output_file( [numbers, /all, /dat, /big, count=] )'
                    return,''
                endif 

                nstr = ntostr(numbers, format='(I20.3)')
                files = 'DR5_'+nstr+'_final'+ext
                count = n_elements(files)

                return, concat_dir(dir, files)
            endif else begin
                pattern = 'DR5_*_final'+ext
                files = file_search(dir, pattern, count=count)
                return, files
            endelse

        end
        'dr6cc2': begin

            dir = self->output_dir()

            if keyword_set(big) then begin
                return, concat_dir(dir, 'DR6.st')
            endif else if not keyword_set(all) then begin
                if n_elements(numbers) eq 0 then begin 
                    print,'-Syntax: files = obj->output_file( [numbers, /all, /dat, /big, count=] )'
                    return,''
                endif 

                nstr = ntostr(numbers, format='(I20.3)')
                files = 'DR6_'+nstr+ext
                count = n_elements(files)

                return, concat_dir(dir, files)
            endif else begin
                pattern = 'DR6_*'+ext
                files = file_search(dir, pattern, count=count)
                return, files
            endelse


        end
        else: message,'Unsupported type: '+self.type
    endcase

end 

; for reading the .dat files
; #objID  ra   dec   photoz   photoz_err   probgal   probgal_flag
function photoz_uchicago::output_structdef, n, original=original

	case self.type of
		'nn': begin 
			st = {      $ 
				casid:0ll,      $
				ra:0d,            $
				dec:0d,           $
				photoz_z:0.0,     $
				photoz_zerr:0.0,  $
				probgal: 0.0,     $
				probgal_flag: 0   $
				}
		end 
		'dr6cc2': begin 
			if keyword_set(original) then begin
				st = {      $ 
					casid:0ll,      $
					ra:0d,            $
					dec:0d,           $
					zcc2:0.0,     $
					zcc2_err:0.0,  $
					zd1: 0.0, $
					zd1_err: 0.0, $
					u: 0.0, $
					g: 0.0, $
					r: 0.0, $
					i: 0.0, $
					z: 0.0 }
			endif else begin
				st = { $
					casid: 0LL, $
					ra: 0d, $
					dec: 0d, $
					photoz_z: 0.0, $
					photoz_zerr: 0.0 $
				}
			endelse
		end 
		'dr7pofz': begin

			if keyword_set(original) then begin
				;1740 40 1 91 427 0.157377 0.124464 0.115932 0.083024 2 587727225157779883 10.50000507 -10.85634083       0.125581 0.067288 0.060078 0.069887 0.070902 0.088673 0.114152 0.239471 0.309163 0.353196 0.366523 0.295653 0.295341 0.323038 0.311881 0.248626 0.185506 0.099445 0.065554 0.022851 0.011546 0.011546 0 0.011475 0.022963 0.022963 0.022963 0.011488 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

				zvals = arrscl(findgen(100), 0.03, 1.47)
				st = { $
					run:0,$
					rerun:0,$
					camcol:0,$
					field:0,$
					id:0, $
					umg:0.0, $
					gmr:0.0, $
					rmi:0.0, $
					imz:0.0, $
					flags:0, $
					casid: 0LL, $
					ra:0d, $
					dec:0d, $
					pz: fltarr(100) $
				}
			endif else begin
				; might make this <sigma_critinv>(zl)
				st = { $
					ra:0d, $
					dec:0d, $
					flags: 0, $
					pz: fltarr(100) $
				}
			endelse
		end

	endcase 

	if n_elements(n) ne 0 then st=replicate(st,n)
	return,st

END 

function photoz_uchicago::dr7pofz_zvals
	zvals = arrscl(findgen(100), 0.03, 1.47)
	return, zvals
end


function photoz_uchicago::output_read, numbers, all=all, dat=dat, columns=columns, count=count

    count=0L
    files = self->output_file(numbers, all=all, dat=dat, count=fcount)
    if fcount eq 0 then begin
        message,'No files found', /inf
        return, -1
    endif

    if keyword_set(dat) then begin
        return, self->read_datfile(files, count=count)
    endif else begin
        ;return, mrdfits_multi(files, ext=1, count=count)
        return, read_idlstruct_multi(files, columns=columns, count=count)
    endelse
  
end 

function photoz_uchicago::output_skiplines
  case self.type of
      'nn': return, 1
      'dr6cc2': return, 0
	  'dr7pofz': return, 0
      else: message,'Dont support type '+self.type+' yet'
  endcase 
end

FUNCTION photoz_uchicago::read_datfile, files, count=count

    count = 0L
    structdef = self->output_structdef(/original)
    skipline=self->output_skiplines()

    ;; Get number in each file
    nfiles = n_elements(files)
    nlines = lon64arr(nfiles)
    for i=0l, nfiles-1 do begin 
        if not fexist(files[i]) then begin 
            message,'file does not exist, skipping: '+files[i],/inf
        endif else begin 
            spawn,'wc -l '+files[i], result, /sh    
            ; account for header in count, but not in nlines for read_struct
            nlines[i] = long(result[0])
            count = count + nlines[i] - skipline
        endelse 
    endfor 

    if count eq 0 then begin
        message,'no files found', /inf
        return, -1
    endif

    zphot_struct = replicate(structdef, count)
    beg = 0ll
    for i=0l, nfiles-1 do begin 
        if nlines[i] gt 0 then begin 
			print,'Reading file: '+files[i]
            read_struct, files[i], structdef, tmp, $
                skipline=skipline, nlines=nlines[i]

            nobj = n_elements(tmp)
            zphot_struct[beg:beg+nobj-1] = tmp
            tmp=0
        endif 
    endfor 

    return,zphot_struct

end 

pro photoz_uchicago::convert2idlstruct

    files = self->output_file(/all,/dat,count=count)

    for i=0L, count-1 do begin
        file = files[i]
        oext = self->output_ext(/dat)
        next = self->output_ext()
        outfile = repstr(file, oext, next)

        struct = self->read_datfile(file)

		if n_tags(struct) eq 0 then begin
			message,'file empty?: '+file,/inf
		endif else begin
			n = n_elements(struct)
			case self.type of
				'dr6cc2': begin
					ost=temporary(struct)
					struct = self->output_structdef(n)
					struct.casid=ost.casid
					struct.ra=ost.ra
					struct.dec=ost.dec
					struct.photoz_z=ost.zcc2
					struct.photoz_zerr=ost.zcc2_err
				end
				'dr7pofz': begin
					ost=temporary(struct)
					struct = self->output_structdef(n)
					struct.ra=ost.ra
					struct.dec=ost.dec
					struct.flags=ost.flags
					struct.pz=ost.pz
					ost=0
				end
				else:
			endcase

			print
			print,'Writing to file: ',outfile

			write_idlstruct, struct, outfile
		endelse
    endfor

end 









;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Stuffing to database
;; This is a better way than the postgres_zphot class. Growing pains.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; note, probgal_flag not here because we cut on it
function photoz_uchicago::output_stuff_struct, num
 
  struct = $
    {                    $
      photoid: 0LL,      $
      run: 0,            $
      rerun: 0,          $
      camcol: 0,         $
      field: 0,          $
      id: 0,             $
      stripe: 0,         $
      ra: 0d,            $
      dec: 0d,           $
      photoz_z:0.0,      $
      photoz_zerr:0.0,   $
      probgal: 0.0       $
    }

  if n_elements(num) ne 0 then begin 
      struct = replicate(struct, num)
  endif 

  return, struct
end 

function photoz_uchicago::read_runpar
    dir = self->output_dir()
    file = concat_dir(dir, 'run.par')
    yanny_read, file, pdata
    struct = *pdata
    ptr_free, pdata
    return, struct
end
pro photoz_uchicago::getstripe, struct, runpar

    ; histogram from 0 to max run
    h=histogram(struct.run, min=0, rev=rev)
    nh=n_elements(h)
    for i=0L, nh-1 do begin
        if rev[i] ne rev[i+1] then begin
            ; objects in this run
            w=rev[ rev[i]:rev[i+1]-1 ]

            num = n_elements(w)

            run = struct[w[0]].run

            wp = where(runpar.run eq run, nwp)
            if nwp eq 0 then begin
                print,'No matching run: ',run
                struct.stripe = -9999
            endif else begin
                stripe = runpar[wp[0]].stripe
                print,'run = ', run, stripe, num
                struct[w].stripe = stripe
            endelse
        endif
    endfor

end

; pick out galaxies
function photoz_uchicago::stuff_select, struct, nw
    w=where($
        struct.probgal gt 0.8 $
        and struct.probgal le 1 $
        and struct.probgal_flag eq 0,nw)
    return, w
end

; make one gigantic file.  This will allow us to remove duplicates.
; run with /rmdup to remove the duplicates.  Make sure you run on 
; a high memory machine
pro photoz_uchicago::make_big_file, rmdup=rmdup

    outfile = self->output_file(/big)

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

    files = self->output_file(/all, count=nf)

    runpar=self->read_runpar()
    for i=0L, nf-1 do begin 
        print,'------------------------------------------------------------'
        print,'Reading file: ',files[i]
        
        st = read_idlstruct(files[i])

        wkeep = self->stuff_select(st, nkeep)
        if nkeep ne 0 then begin
            print,'kept '+ntostr(nkeep)+'/'+ntostr(n_elements(st))
            st = st[wkeep]
            outst = self->output_stuff_struct(nkeep)

            copy_struct, st, outst

            casid_extract, st.casid, run, rerun, camcol, field, id
            outst.run=run
            outst.rerun=rerun
            outst.camcol=camcol
            outst.field=field
            outst.id=id
            outst.photoid = photoid(run, rerun, camcol, field, id)

            ; add stripe info
            self->getstripe, outst, runpar

            write_idlstruct, outst, outfile, /append
        endif else begin
            message,'No objects passed stuff_select', /inf
        endelse

    endfor 

end



function photoz_uchicago::tablename
    case self.type of
        'nn': return, 'zphot'
        'dr6cc2': return, 'zphotcc2'
        else: message,'dont yet support type: '+self.type
    endcase
end
pro photoz_uchicago::output_stuff

    pg = obj_new('postgres')

    table = self->tablename()

    tmpdir = self->output_dir()
    if self.type eq 'nn' then begin
        ; Read in the big file and stuff it into database
        file = self->output_file(/big)
        print,'Reading file: ',file
        st = read_idlstruct(file)
        pg->struct2table, st, table, $
            primary_key='photoid', conn='user=postgres', $
            tmpdir=tmpdir, status=status

        if status ne 0 then message,'Stuff failed'
        ; Create indices, analyze, and grant read permissions
        print
        print,'Creating indices'
        ;; multi-column index on (run,rerun,camcol,field,id)
        query = $
            'CREATE INDEX '+table+'_rrcfi_index ON '+table+' (run,rerun,camcol,field,id)'
        print,query
        pg->query, query, status=status,conn='user=postgres'

        if status ne pg->status_val('no_result') then begin 
            message,'CREATE INDEX failed'
        endif 

        pg->create_index, table, ['stripe','photoz_z','photoz_zerr','probgal'], conn='user=postgres'

    endif else begin
        files = self->output_file(/all, count=count)
        for i=0L, count-1 do begin
            print,'Reading file: ',files[i]
            st = read_idlstruct(files[i])

            pg->struct2table, st, table, $
                primary_key='casid', conn='user=postgres', $
                tmpdir=tmpdir, status=status
            
        endfor

        pg->create_index, table, ['photoz_z','photoz_zerr'], $
            conn='user=postgres'
    endelse



    query = 'ANALYZE '+table
    print
    print,query
    pg->query, query, status=status,conn='user=postgres'

    if status ne pg->status_val('no_result') then begin 
        message,'ANALYZE failed'
    endif 

    ; Make sure the sdss user can select
    query = 'GRANT select ON '+table+' TO sdss'
    print
    print,query
    pg->query,  query, conn='user=postgres', status=gstatus
    if gstatus ne pg->status_val('no_result') then message,'Could not grant select to sdss'


    obj_destroy, pg


end 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Analysis code
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




FUNCTION photoz_uchicago::specgal_pgsql_get

  ;; Get matches with no other cuts
  query = $
    'SELECT '+$
    'pz.photoid, s.z, s.primtarget, '+$
    'pz.photoz_z, pz.photoz_zerr4 as photoz_zerr, pz.photoz_zwarning '+$
    'FROM specgal as s, zphot as pz ' +$
    'WHERE s.match_photoid != -1 AND s.match_photoid = pz.photoid'

  print,query

  pg=obj_new("postgres")
  st = pg->query(query, status=status)


  return,st
END 

FUNCTION photoz_uchicago::specgal_match_file, fits=fits
  dir = concat_dir(self->output_dir(),'specgal_match')
  IF keyword_set(fits) THEN BEGIN 
      file = concat_dir(dir,'specgal_match_'+self.type+'.fit')
  ENDIF ELSE BEGIN 
      file = concat_dir(dir, 'specgal_match_'+self.type+'.st')
  ENDELSE 

  return,file

END 
FUNCTION photoz_uchicago::specgal_match_get, fits=fits
  file = self->specgal_match_file(fits=fits)
  print
  print,'Reading file: ',file
  IF keyword_set(fits) THEN BEGIN 
      struct = mrdfits(file,1)
  ENDIF ELSE BEGIN 
      struct = read_idlstruct(file)
  ENDELSE 
  return,struct
END 
PRO photoz_uchicago::specgal_match_create

  st = self->specgal_pgsql_get()

  file = self->specgal_match_file(/fits)
  print
  print,'Writing to file: ',file
  mwrfits, st, file, /create

  file = self->specgal_match_file()
  print
  print,'Writing to file: ',file
  write_idlstruct, st, file

END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make plots of the error distribution z-photoz for different samples,
; binned by z and photoz
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO photoz_uchicago::error_dist, st, lrg=lrg

  plot_dir = '~/plots/photoz/uchicago/'

  pg=obj_new('postgres')
  IF n_elements(st) EQ 0 THEN BEGIN 
      query = $
        'SELECT '+$
        's.z, s.primtarget, s.subclass, s.cmodel_counts, '+$
        'pz.photoz_z, pz.photoz_zerr4 as photoz_zerr, pz.photoz_zwarning '+$
        'FROM specgal as s, zphot as pz ' +$
        'WHERE s.match_photoid != -1 AND s.match_photoid = pz.photoid'
      
      st = pg->query(query, status=status)

  ENDIF 


;  w = lindgen(nn)

  IF NOT keyword_set(lrg) THEN BEGIN 
      file_front = 'error_dist_main'
      nperzbin = 10000
      w = $
        where( (st.primtarget AND $
                sdss_flag('target', 'galaxy')) NE 0 AND $
               (st.primtarget AND $
                sdss_flag('target', 'southern_survey')) EQ 0,$
               nw)
  ENDIF ELSE IF 1 THEN BEGIN

      file_front = 'error_dist_lrg'
      nperzbin = 5000
      w = $
        where( (st.primtarget AND $
                sdss_flag('target', 'galaxy_red')) NE 0 AND $
               (st.primtarget AND $
                sdss_flag('target', 'southern_survey')) EQ 0,$
               nw)
  ENDIF 


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot histograms of delta_z in bins of spectroscopic redshift
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  xtitle = 'z - photoz'
  s = sort(st[w].z)
  w = w[s]

  ind = lindgen(nw)

  hist = histogram(ind, bin=nperzbin, reverse_indices=zrev_ind)

  wzh = where(hist EQ nperzbin, nbin)
  print,'nbin = ',nbin

  ny = 6
  nx = nbin/ny
  IF (nbin MOD ny) NE 0 THEN nx = nx+1 
  !p.multi = [0,nx,ny]
  FOR zi=0L, nbin-1 DO BEGIN 

      wzi = zrev_ind[ zrev_ind[zi]:zrev_ind[zi+1] -1 ]
      wzi = w[wzi]

      ;; plot zphot around this z bin

      meanz = mean(st[wzi].z)
      mzstr = '<z> = '+ntostr(meanz, 5, /round)
      
      xrange = [-0.2, 0.2]
      zdiff = st[wzi].z - st[wzi].photoz_z

      plothist, zdiff, bin=0.01, xrange=xrange, title=mzstr, xtitle=xtitle, $
        charsize=1.5
      ;;legend,mzstr,box=0,/right
      ;;key=prompt_kbrd()

  ENDFOR 

  png_file = plot_dir + file_front+'_z_bin.png'
  print,'Writing png file: ',png_file
  write_png, png_file, tvrd(/true)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot histograms of delta_z in bins of spectroscopic redshift
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  key = prompt_kbrd("Hit a key for photoz binning")

  s = sort(st[w].photoz_z)
  w = w[s]

  ind = lindgen(nw)

  hist = histogram(ind, bin=nperzbin, reverse_indices=zrev_ind)

  wzh = where(hist EQ nperzbin, nbin)
  print,'nbin = ',nbin

  ny = 6
  nx = nbin/ny
  IF (nbin MOD ny) NE 0 THEN nx = nx+1 
  !p.multi = [0,nx,ny]
  FOR zi=0L, nbin-1 DO BEGIN 

      wzi = zrev_ind[ zrev_ind[zi]:zrev_ind[zi+1] -1 ]
      wzi = w[wzi]

      ;; plot zphot around this z bin

      meanz = mean(st[wzi].photoz_z)
      mzstr = '<photo_z> = '+ntostr(meanz, 5, /round)
      
      xrange = [-0.2, 0.2]
      zdiff = st[wzi].z - st[wzi].photoz_z

      plothist, zdiff, bin=0.01, xrange=xrange, title=mzstr, xtitle=xtitle, $
        charsize=1.5
      ;;legend,mzstr,box=0,/right
      ;;key=prompt_kbrd()

  ENDFOR 

  png_file = plot_dir + file_front+'_photoz_bin.png'
  print,'Writing png file: ',png_file
  write_png, png_file, tvrd(/true)


  !p.multi=0
END




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Compare the uchicago photoz outputs to old budavari
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO photoz_uchicago::compare_budavari, $
                   stripe, $
                   zphot=zphot, bud_zphot=bud_zphot, $
                   mzphot=mzphot, mbud=mbud, $
                   dops=dops

  IF keyword_set(dops) THEN BEGIN 
      
      dir = '~/plots/photoz/uchicago/'
      file=dir+'compare_zphot_budavari_'+$
        'stripe'+strn(stripe,len=2,padchar='0')+'.ps'
      begplot, name=file, /color
  ENDIF 

;  !p.multi=[0,0,2]

  IF (!d.name EQ 'X' OR !d.name EQ 'Z') THEN BEGIN 
      wzphotclr = !green
      budclr = !cyan
  ENDIF ELSE BEGIN 
      wzphotclr = !blue
      budclr = !red
  ENDELSE 

  w=where(!run_status.stripe EQ stripe)
  
  wruns = rem_dup(!run_status[w].run)
  runs = !run_status[w[wruns]].run
  nruns = n_elements(runs)
  reruns = intarr(nruns)
  FOR i=0L, nruns-1 DO reruns[i] = sdss_rerun(runs[i])

  forprint,runs,reruns

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read the files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(zphot) EQ 0 THEN zphot = self->output_read(runs)
  IF n_elements(bud_zphot) EQ 0 THEN get_scat,stripe,[1,2,3],bud_zphot,/nzcuts

  nzphot = n_elements(zphot)
  nbud = n_elements(bud_zphot)

  w=where(zphot.photoz_zwarning EQ 0,nw)
  print,float(nw)/nzphot,' % objects passed zwarning cut'

  ;; get copies of the redshifts
  zz = zphot.photoz_z
  bzz = bud_zphot.photoz_z
  mm = zphot.counts_model_pogson[2]

  ;; Histograms of all decent measurements (although zwarn is not great)
  binsize = 0.01
  minz = 0.02
  maxz = 1.0
  plothist, zz, min=minz, max=maxz, bin=binsize, /norm, $
    xtitle='photoz', position=aspect(1./!gratio)

  plothist, zz[w], min=minz, max=maxz, bin=binsize, $
    /overplot, color=wzphotclr, /norm, xtitle='zphot'

  plothist, bzz, min=minz, max=maxz, bin=binsize, $
    /overplot, color=budclr, /norm

  legend,['NN','NN zwarn=0','Budavari'], $
    colors=[!p.color, wzphotclr, budclr], line=[0,0,0],/right,box=0, $
    charsize=1

  ;; Match by radec
  IF n_elements(mzphot) EQ 0 AND n_elements(mbud) EQ 0 THEN BEGIN 
      print,'Running sphermatch'
      csurvey2eq, bud_zphot.clambda, bud_zphot.ceta, bra, bdec
      zphotra = zphot.ra
      zphotdec = zphot.dec
      spherematch, $
        zphotra[w], zphotdec[w], bra, bdec, 1.0/3600d, mzphot, mbud, d12, $
        maxmatch=1
  
      mzphot = w[mzphot]
      help,zphot,bud_zphot, mzphot, mbud
  ENDIF 

  key = prompt_kbrd()

  budclr = wzphotclr
  plothist, zz[mzphot], min=minz, max=maxz, bin=binsize, /norm, $
    xtitle='photoz', position=aspect(1./!gratio)

  plothist, bzz[mbud], min=minz, max=maxz, bin=binsize, $
    /overplot, color=budclr, /norm

  legend,['NN matched','budavari matched'], $
    colors=[!p.color, budclr], line=[0,0],/right,box=0, $
    charsize=1

;  outst = {bud_petrocounts:fltarr(5), $
;           bud_photoz_z:0.0, $
;           bud_photoz_zerr: 0.0, $
;           bud_photoz_quality:0.0, $

;  binsize=0.1
;  minm = 15.
;  maxm = 24.
;  plothist, mm, min=minm, max=maxm, bin=binsize, $
;    xtitle='r', /norm
;  plothist, mm[w], min=minm, max=maxm, bin=binsize, $
;    /overplot, color=wzphotclr, /norm
;  plothist, bud_zphot.rpetro,  min=minm, max=maxm, bin=binsize, $
;    /overplot, color=budclr, /norm

;  legend,['zphot','zphot zwarn=0','budavari'], $
;    colors=[!p.color, wzphotclr, budclr], line=[0,0,0],box=0, $
;    charsize=1

;  !p.multi=0

  IF keyword_set(dops) THEN endplot

END 





;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Error plot
;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION photoz_uchicago::read_zdiff
  dir = sdssidl_config('photoz_dir')
  file =  concat_dir(dir, 'photoz_output_nn/zpminzs.dat')
  read_struct, file, {zspec:0.0, zpminzs: 0.0}, struct
  return,struct
END 


;; /color is for ps.  only works for scatter plot
PRO photoz_uchicago::bias_plot, dops=dops, scatter=scatter, color=color

  struct = self->read_zdiff()
  nperbin = 500

  bs = binner(struct.zspec, struct.zpminzs, nperbin=nperbin, method='median')

;  binner_bynum, struct.zspec, struct.zpminzs, nperbin, $
;                xb, yb, yberr, ybinned_sdev = ybsig, $
;                /median

  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(scatter) THEN $
        addstr = 'scatter' ELSE addstr = 'hist2d'
      IF keyword_set(color) THEN $
        addstr = addstr + '_color'
      file = self->plot_dir() + 'bias_plot_'+addstr+'.eps'
      begplot, name=file, /encapsulated, color=color
  ENDIF 

  xtitle = 'z!Dspec!N'
  ytitle = 'z!Dphot!N'+!csym.minus+'z!Dspec!N'
  xrange = [0.0, 0.6]
  yrange = [-0.15, 0.15]

  IF !d.name EQ 'X' THEN BEGIN 
      meanclr = !red
      rangeclr = meanclr
      othick = 2
  ENDIF ELSE BEGIN 
      IF keyword_set(color) THEN BEGIN 
          meanclr = !DarkGreen
          rangeclr = !red
          oline = 0
      ENDIF ELSE BEGIN 
          oline = 2
          IF keyword_set(scatter) THEN BEGIN 
              meanclr = !grey75
              rangeclr = !grey75
          ENDIF 
      ENDELSE 
  ENDELSE 

  nsmooth = 3
  xb = bs.xbinned
  yb = smooth(bs.ybinned, nsmooth)
  ybsig = smooth(bs.ybinned_sdev, nsmooth)



  IF keyword_set(scatter) THEN BEGIN 
      aplot, !gratio, struct.zspec, struct.zpminzs, $
        psym=8, symsize=0.25, $
        xrange=xrange, yrange=yrange, $
        xstyle=3, ystyle=3, $
        xtitle=xtitle, ytitle=ytitle
  ENDIF ELSE BEGIN
      xmold = !x.margin
      !x.margin = [12,3]
      IF keyword_set(color) THEN loadct,0
      ploth, struct.zspec, struct.zpminzs, psym=3, $
             xrange=xrange, yrange=yrange, $
             xstyle=3, ystyle=3, $
             xtitle=xtitle, ytitle=ytitle, /log, range=[0, 2.05]
      !x.margin = xmold
  ENDELSE 
  simpctable
  oplot, bs.xbinned, yb, color=meanclr, thick=othick
  oplot, bs.xbinned, yb+ybsig, color=rangeclr, line=oline, thick=othick
  oplot, bs.xbinned, yb-ybsig, color=rangeclr, line=oline, thick=othick

  IF keyword_set(dops) THEN endplot, /trim_bbox


END 





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Carlos' deconvolutions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION photoz_uchicago::deconv_dir
  return,esheldon_config('lensinput_dir')+'srcgal_dndz/'
END 
FUNCTION photoz_uchicago::deconv_file
  dir = self->photoz_uchicago::deconv_dir()
  return,dir + 'dndzspec.erin.0.040.dat'
END 
FUNCTION photoz_uchicago::deconv_read, $
  original=original, $
  hybrid=hybrid, $
  photoz=photoz, $
  deconv=deconv

  ;; Default is hybrid

; From carlos:
; 1: centroid of redshift bin
; 2: reconstructed distribution
; 3: photo-z distribution
; 4: distribution of training set
; 5: 30-bootstrap-averaged reconstruction. 
; I would use the 2nd column up to the bin centered at 0.62 
; and the photometric redshift distribution afterwards. 

; I hard-coded this into the reader below by indexing the arrays


  file = self->deconv_file()
  struct = {z: 0.0, $
            deconv_dndz: 0.0, $
            photoz_dndz: 0.0, $
            training_dndz: 0.0, $
            boot_dndz:0.0}

  read_struct, file, struct, st
  nst = n_elements(st)

  IF keyword_set(original) THEN return,st

  outstruct = {z: 0.0, $
               dndz: 0.0}

  outstruct = replicate(outstruct, nst)
  outstruct.z = st.z

  IF keyword_set(photoz) THEN BEGIN 
      print,'Returning photoz dndz'
      outstruct.dndz = st.photoz_dndz
  ENDIF ELSE IF keyword_set(deconv) THEN BEGIN 
      print,'Returning deconvolved dndz'
      outstruct.dndz = st.deconv_dndz
  ENDIF ELSE BEGIN 
      ;; I'm following Carlos' advice above.
      print,'Returning hybrid dndz'
      outstruct[0:15].dndz = st[0:15].deconv_dndz
      outstruct[16:nst-1].dndz = st[16:nst-1].photoz_dndz
  ENDELSE 
  return,outstruct
            
END 

FUNCTION photoz_uchicago::plot_dir
  return,'~/plots/photoz/uchicago/'
END 
PRO photoz_uchicago::deconv_plot, dops=dops

  IF keyword_set(dops) THEN BEGIN 
      deconv_file = self->plot_dir()+'deconv_plot.eps'
      deconv_allfile = self->plot_dir()+'deconv_allplot.eps'
  ENDIF  

  ost = self->deconv_read(/original)
  st = self->deconv_read()

  IF keyword_set(dops) THEN begplot,name=deconv_file, /encap
  


  xrange = [0,1.0]
  xtitle = 'z'
  ytitle = 'dn/dz'
  aplot, !gratio, st.z, st.dndz, psym=10, $
        xtitle=xtitle, ytitle=ytitle, xrange=xrange

  key = prompt_kbrd('hit a key')
  IF keyword_set(dops) THEN endplot, /trim_bbox



  IF keyword_set(dops) THEN BEGIN 
      begplot,name=deconv_allfile, /color, /encap
;      dclr = !blue
      pclr = !grey50

      pthick = !p.thick
      uthick = 2*!p.thick

  ENDIF ELSE BEGIN 
;      dclr = !green
      pclr = !green

      pthick = !p.thick
      uthick = 2*!p.thick
  ENDELSE

;  dline = 0
  pline = 0
  uline = 0


  aplot, !gratio, st.z, st.dndz, psym=10, $
        xtitle=xtitle, ytitle=ytitle, xrange=xrange

;  oplot, ost.z, ost.deconv_dndz, psym=10, color=dclr, line=dline
  oplot, ost.z, ost.photoz_dndz, psym=10, color=pclr, line=pline
  oplot, st.z, st.dndz, psym=10, line=uline

  legend, ['dn/dz deconv','dn/dz photoz'],$
          line=[uline,pline], color=[!p.color, pclr], $
          box=0, /right, charsize=1

  IF keyword_set(dops) THEN endplot, /trim_bbox

END 



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Marcos' matching deconvolutions.  He assigns weights to the spec galaxies
; in such a way that the weighted hist of the spec match the normal hist of
; the photometric sample.  The the idea is that the  weighted specz hist
; might be close to that of the photometric sample
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function photoz_uchicago::wroot
    return,sdssidl_config('photoz_dir')
end
function photoz_uchicago::psample
    return,'princeton6'
end
function photoz_uchicago::wdir, dtype, subtype=subtype, number=number, createdir=createdir
    if n_elements(dtype) eq 0 or n_elements(sample) eq 0 then begin
        print,'-Syntax: dir = pz->wdir(sample, dtype, subtype=, /createdir)'
        print,' dtype = input|output'
    endif
    f=self->wfile(dtype, subtype=subtype, number=number, createdir=createdir, dir=dir)
    return, dir 
end
function photoz_uchicago::wfile, pathel, number=number, dir=dir, createdir=createdir

    if n_elements(pathel) eq 0 then begin
        print,'-Syntax: files = pz->wfile(pathel, number=, dir=, /createdir)'
        print,' dtype = input|output'
    endif

    ext='.dat'
    root=self->wroot()
    sample=self->wsample()

    file = datafile(pathel, root=root, project='weight', sample=sample,  $
                    number=number, ext=ext, dir=dir)
    return, file
end



; Write the photometric input file
pro photoz_uchicago::winput_write, sample, number=number
    file = self->wfile(sample, 'input', sub='photo', number=number)
end

function photoz_uchicago::wread, nneigh, catnum
    f=self->wfile(nneigh, catnum)
    sdef = {zphot:0.0, zspec:0.0, w:0.0, cm:fltarr(5)}
    read_struct,f, sdef, struct
    return, struct
end






function photoz_uchicago::wplot_dir
    pd = self->plot_dir()
    pd = filepath(root=pd, 'weight')
    file_mkdir, pd  
    return, pd
end
pro photoz_uchicago::wcompare_legend, colors
    legend, ['p zphot', 'w spec zphot', 'w spec zspec'], $
        line=0, color=colors, $
        /right, box=0, charsize=1
end
pro photoz_uchicago::wcompare, nneigh, catnum, ws, ps, dops=dops, $
        nboot=nboot, zbin=zbin, zmin=zmin, zmax=zmax

	pg=obj_new('postgres')
    pdir = self->wplot_dir()
    psfile = 'zphot-compare-weighed.'+ntostr(nneigh)+'.'+ntostr(catnum)+'.ps'
    psfile=filepath(root=pdir, psfile)

    if n_elements(ws) eq 0 then begin
        ws = self->wread(nneigh, catnum)
    endif
    if n_elements(ps) eq 0 then begin
        q = 'select photoz_z as zphot, counts_model as cm from scat_princeton6 limit 1000000'
        print,q
        ps=pg->query(q)
    endif

    if keyword_set(dops) then begplot,psfile,/color

    ; Create index array for bootstrap samples
    ndata = n_elements(ws)
    ;print,'Getting bootstrap sample indices'
    ;barr = boot_indarray(ndata, nboot)

    ; For z binning
    if n_elements(zbin) eq 0 then zbin = 0.03
    if n_elements(zmin) eq 0 then zmin = 0.0
    if n_elements(zmax) eq 0 then zmax = 1.2

    if !d.name eq 'PS' then zpcolor=!blue else zpcolor=!green
    zpcolor=!darkgreen
    zscolor=!red

    if 1 then begin
        ; First plot all
        print,'Doing histograms for all galaxies'
        plothist, ps.zphot, /norm, $
            bin=zbin, min=zmin, max=zmax, $
            xtitle='z', ytitle='P(z)', aspect=!gratio, xrange=[0,zmax], /fill
        plothist, ws.zphot, weights=ws.w, /overplot, /norm, nboot=nboot, $
            bin=zbin, min=zmin, max=zmax, bhist=bhist, /verbose, color=zpcolor, $
            line=0
        plothist, ws.zspec, weights=ws.w, /overplot, /norm, nboot=nboot, $
            bin=zbin, min=zmin, max=zmax, bhist=bhist, /verbose, color=zscolor

        self->wcompare_legend, [!p.color, zpcolor, zscolor]

        key=prompt_kbrd('hit a key')

        !p.charsize=1
        erase & multiplot, [3,2], /square 

        mbin = 0.2
        mmin = 16. & mmax=28.
        xrange=[mmin,mmax]
        yrange=[0,1.2]

        for i=0,4 do begin
            if i gt 2 then xt='mag' else xt=''
            if i eq 0 or i eq 3 then yt='dN/dmag' else yt=''
            plothist, ps.cm[i], bin=mbin, min=mmin, max=mmax, peak=1, $
                xrange=xrange, xsty=1, yrange=yrange, ystyle=1, xtitle=xt, ytitle=yt
            plothist, ws.cm[i], bin=mbin, min=mmin, max=mmax, peak=1, $
                weights=ws.w, /overplot, color=zscolor
              
            legend, !colors[i], box=0, charsize=1
            ;if i eq 0 then begin
            ;    legend,['photometric','spec weighted'], $
            ;        line=0, color=[!p.color, zscolor], /right, box=0, $
            ;        charsize=0.5, number=0.5
            ;endif
            if i ne 4 then multiplot
        endfor
        multiplot, /default

        key=prompt_kbrd('hit a key')

        cbin = [0.1,0.1,0.05,0.05]
        erase & multiplot, [2,2], /square
        for i=0,3 do begin
            if i gt 1 then xt='color' else xt=''
            if i eq 0 or i eq 2 then yt='dN/color' else yt=''

            color = ps.cm[i] - ps.cm[i+1]
            cmin=-2. & cmax=4.
            xrange=[cmin,cmax]
            yrange=[0,1.2]
            plothist, color, bin=cbin[i], min=cmin, max=cmax, peak=1, $
                xrange=xrange, xsty=1, yrange=yrange, ystyle=1, xtitle=xt, ytitle=yt
            color = ws.cm[i] - ws.cm[i+1]
            plothist, color, bin=cbin[i], min=cmin, max=cmax, peak=1, $
                weights=ws.w, /overplot, color=zscolor

            cstr = !colors[i]+' - '+!colors[i+1]
            legend, cstr, box=0, charsize=1
            if i ne 3 then multiplot
        endfor
        multiplot, /default
        key=prompt_kbrd('hit a key')
    endif

    if keyword_set(dops) then endplot
    return
    
    zpbin = 0.015
    nbin = 16
    zpminarr = 0.1 + 0.05*findgen(nbin)
    zpmaxarr = 0.1 + 0.05*findgen(nbin) + 0.05

    erase & multiplot, [4,4], /square
    !p.charsize=0.75

    yrange = [0,10]
    for i=0L, n_elements(zpminarr)-1 do begin

        zpmin=zpminarr[i]
        zpmax=zpmaxarr[i]
        
        if i eq 0 or i eq 4 or i eq 8 or i eq 12 then yt='P(z)' else yt=''
        if i ge 12 then xt='z' else xt=''

;        pi=where(ps.zphot ge zpmin and ps.zphot le zpmax)
;        plothist, ps[pi].zphot, /norm, $
;            bin=zbin, min=zmin, max=zmax, $
;            xtitle=xt, ytitle=yt, xrange=[0,zmax], yrange=yrange


        si=where(ws.zphot ge zpmin and ws.zphot le zpmax)
;        plothist, ws[si].zphot, weights=ws[si].w, /overplot, /norm, nboot=nboot, $
;            bin=zbin, min=zmin, max=zmax, bhist=bhist, /verbose, color=zpcolor

        plothist, ws[si].zspec, weights=ws[si].w, /norm, nboot=nboot, $
            bin=zbin, min=zmin, max=zmax, bhist=bhist, /verbose, $
            xtitle=xt, ytitle=yt, xrange=[0,zmax], yrange=yrange
        ;self->wcompare_legend, [!p.color, zpcolor, zscolor]

        zpmean = (zpmin+zpmax)/2
        pplot, [zpmean,zpmean], [0,100], color=zpcolor, /over

        legend,string(zpmin,f='(f4.2)')+' < z < '+string(zpmax,f='(f4.2)'),$
            /left,box=0,charsize=0.75
        if i ne n_elements(zpminarr)-1 then multiplot
    endfor

    multiplot,/default

    !p.multi=0
    return


end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Princeton
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 

function photoz_uchicago::princeton_input_dir
    dir=sdssidl_config('photoz_dir')
    dir=concat_dir(dir,'photoz_input_princeton')
    return,dir
end 

function photoz_uchicago::princeton_input_file, stripe
    dir = self->princeton_input_dir()
    return,concat_dir(dir,'photoz-input-'+stripe2string(stripe)+'.st')
end 


function photoz_uchicago::princeton_input_structdef
    return,{photoid:0ULL, $
            ra: 0d, $
            dec: 0d, $
            counts_model: fltarr(5) }
end 
function photoz_uchicago::princeton_input_read, stripe, status=status
    file = self->princeton_input_file(stripe)
    print,'Reading file: ',file
    t = read_idlstruct(file, status=status)
    return,t
end 
pro photoz_uchicago::princeton_input_write, stripes=stripes


    ms = obj_new('make_scat')
    minflux = 0.001               ; mag of 30
    maxflux = 1.e5 

    ;; not all will be found, but no big deal
    if n_elements(stripes) eq 0 then stripes = 1+lindgen(86)
    nst = n_elements(stripes)
    for i=0l, nst-1 do begin 
        srcfile = ms->princeton_bystripe_infile(stripes[i])
        if fexist(srcfile) then begin 

            print,'-----------------------------------------------------------'
            print,'Reading princeton srcfile: ',srcfile
            t = mrdfits(srcfile,1)

            nt = n_elements(t)
            outstruct = replicate(self->princeton_input_structdef(), nt)

            outstruct.photoid = $
                photoid(t.run,t.rerun,t.camcol,t.field,t.id)
            outstruct.ra = t.ra
            outstruct.dec = t.dec
            ; no bounds checking needed because asinh are well behaved.
            outstruct.counts_model = k_maggies2lups(t.modelflux*1.e-9)

            outfile = self->princeton_input_file(stripes[i])
            print
            print,'Writing input file: ',outfile
            write_idlstruct, outstruct, outfile, /ascii
        endif 
    endfor 

    obj_destroy,ms

end 

function photoz_uchicago::princeton_output_dir
  dir = concat_dir(sdssidl_config('photoz_dir'),'photoz_output_princeton_nn')
  return,dir
end 

function photoz_uchicago::princeton_output_filelist, count=count
  dir = self->princeton_output_dir()
  pattern = concat_dir(dir,'photoz-input-[0-9][0-9].wzflag.etbl.gz')
  files=file_search(pattern, count=count)
  return,files
end

function photoz_uchicago::princeton_output_file, stripe
  if n_elements(stripe) eq 0 then begin 
      print,'file = pz->princeton_output_file(stripe)'
      on_error, 2
      message,'Halting'
  endif 
  dir = self->princeton_output_dir()
  files = concat_dir(dir,'photoz-input-'+stripe2string(stripe)+'.wzflag.etbl.gz')
  return,files
end 
function photoz_uchicago::princeton_output_structdef
  struct = {photoz_z:0.0, $
            photoz_zerr4: 0.0, $
            photoid: 0ULL, $
            counts_model_pogson: fltarr(5),$
            photoz_zwarning:0b}

  return,struct
            
end 
function photoz_uchicago::princeton_output_read, stripe_or_file, count=count
    if size(stripe_or_file, /tname) eq 'STRING' then begin
        files = stripe_or_file
    endif else begin
        stripe = stripe_or_file
        files = self->princeton_output_file(stripe)
    endelse

    structdef = self->princeton_output_structdef()
    return, self->read_multi(files, structdef, count=count)
end 

function photoz_uchicago::princeton_output_stuff_struct, num

  struct = $
    {                    $
      photoid: 0LL,      $
      run: 0,            $
      rerun: 0,          $
      camcol: 0,         $
      field: 0,          $
      id: 0,             $
      stripe: 0,         $
      photoz_z:0.0,      $
      photoz_zerr:0.0,   $
      photoz_zwarning:0b $
    }

  IF n_elements(num) NE 0 THEN BEGIN 
      struct = replicate(struct, num)
  ENDIF 

  return, struct
END 


pro photoz_uchicago::princeton_output_stuff

    pg = obj_new('postgres')

    table = 'zphot_princeton'
    tmpdir = self->output_dir()


    for stripe=1L, 86 do begin 
        print,'------------------------------------------------------------'
        st = self->princeton_output_read(stripe, count=count)

        if count gt 0 then begin

            outst = self->princeton_output_stuff_struct(count)

            photoid_extract, st.photoid, run, rerun, camcol, field, id

            outst.photoid=st.photoid
            outst.run=run
            outst.rerun=rerun
            outst.camcol=camcol
            outst.field=field
            outst.id=id

            outst.stripe=stripe

            outst.photoz_z=st.photoz_z
            outst.photoz_zerr=st.photoz_zerr4
            outst.photoz_zwarning=st.photoz_zwarning

            pg->struct2table, $
                outst, table, primary_key='photoid', conn='user=postgres', $
                tmpdir=tmpdir, status=status

            if status ne 0 then message,'Stuff failed'
        endif
    endfor 


    ;; multi-column index on (run,rerun,camcol,field,id)
    query = $
        'CREATE INDEX '+table+'_rrcfi_index ON '+table+' (run,rerun,camcol,field,id)'
    print,query
    pg->query, query, status=status,conn='user=postgres'

    ;; Note this may already exist
    if status ne pg->status_val('no_result') then begin 
        message,'CREATE INDEX failed'
    endif 

    ; Now the rest of the indices
    pg->create_index, table, ['stripe','photoz_z','photoz_zerr','photoz_zwarning'], conn='user=postgres'
    
    query = 'ANALYZE '+table
    print,query
    pg->query, query, status=status,conn='user=postgres'

    if status ne pg->status_val('no_result') then begin 
        message,'ANALYZE failed'
    endif 

    ; Make sure the sdss user can select
    pg->query, 'GRANT select ON '+table+' TO sdss', conn='user=postgres', status=gstatus
    if gstatus ne pg->status_val('no_result') then message,'Could not grant select to sdss'

    obj_destroy, pg

end 


pro photoz_uchicago::test_sigmacritinv_calc, str
	zvals = self->dr7pofz_zvals()

end


function photoz_uchicago::cleanup
  return,1
end 

pro photoz_uchicago__define

  struct = {$
             photoz_uchicago, $
             type:'' $
           }

end 
