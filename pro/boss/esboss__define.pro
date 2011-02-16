;+
; If you want to do a run:
;   create_pbs or create_pbs_bycamcol (for qsos usually)
;	if you ran bycamcol you will need a separate /dogather 
;	For lrgs we also have the noknown gather.
;
; For examples, see 
;	esboss::create_pbs_lrg_sweep20090928_local
;	esboss::create_pbs_qso_sweep20090928_local
;
;
;
;	Recent tuning on qsos with new likelihood, chunk2:
;		created loose cut on likelihood run '2009-12-10-newlike2-loose'
;		Tuned like to 35/sq degree with ::tune_chunk2_like
;		create new gather '2009-12-10-newlike2-norank' to verify and for
;			adam to use (see below).
;		
;		Create a rank/inindow file from norank for adam to create rank 
;		cuts on a grid using bosstarget_qso_recoverprobs_v2
;
;		Figure out which cuts give 20 and 60/sq degree using 
;			::explore_orcuts to find the right "original density goal" that
;				gives the actual density we want.
;			::orcut_plot to interpolate the 3 ranked cuts based on the
;				above density
;		
;		Create new fial runs with these cuts.  Run names to be determined.
;
;
;
; Then some analysis:
;	::compare_densities
;	::plothist_tile_density
;	::plot_targets_and_tiles
;
; You can tune a 1-d parameter using
;	::tune_chunk_density1d
;		convenience functions:
;			::tune_chunk2_nn  ::tune_chunk2_like
;
; I tuned 2-d with things like:
;	::contour_chunk2_bonus
;
; And x2_star with this on the ngc
;	::tune_chunk2_qso_density_bonus_by_ra
;  !!!!! make sure your run with bosstarget_qso is doing the simple
;  x2_star cut only!
;-





function esboss::init
	return, 1
end

; for commissioning 2 I ran
;	For lrg test
;	::create_pbs_lrg_sweep20090928_local
;	::create_pbs_lrg_sweep20090928_local, /noknown
;   then ran the pbs jobs
;	then ran bosstarget::gather_partial, 'lrg', '2009-10-01-test', /combine
;
;   For qso test, first is by camcol, then do the gathers
;   ::create_pbs_qso_sweep20090928_local
;	ran pbs
;   ::create_pbs_qso_sweep20090928_local, /dogather
;	bosstarget::gather_partial, 'qso', '2009-10-01-test', /combine

function esboss::orflags, flags, types
	; sdss_flagval sums an array of flags
	all_flagfals = sdss_flagval('boss_target1',types)
	; for "or" we just check it is not zero
	logic = (flags and all_flagvals) ne 0
	return, logic
end

function esboss::andflags, flags, types, notset=notset
	; sdss_flagval sums an array of flags
	flagtot = sdss_flagval('boss_target1',types)

	; for and we must check that it actually *equals* the sum
	if not keyword_set(notset) then begin
		logic = (flags and flagtot) eq flagtot
	endif else begin
		logic = (flags and flagtot) ne flagtot
	endelse
	return, logic
end



function esboss::get_field_stats, struct

	field_area = (2048d*1361d)*(0.396d/60d/60d)^2

	h=histogram(struct.field, $
		min=min(struct.field), max=max(struct.field), rev=rev)

	nh = n_elements(h)

	stdef = {$
		field: 0L, $
		mra:0d, mdec:0d, ml:0d, mb:0d, $
		mseeingi:-9999., $
		count:0L, $
		density:0d}

	stats = replicate(stdef, nh)

	field = min(struct.field)

	for i=0L, nh-1 do begin
		stats[i].field = field
		if rev[i] ne rev[i+1] then begin
			wf = rev[ rev[i]:rev[i+1]-1 ]

			stats[i].mra = median_check(struct[wf].ra)
			stats[i].mdec = median_check(struct[wf].dec)
			stats[i].mseeingi = median_check(struct[wf].psf_fwhm[3])
			stats[i].count = h[i]
			stats[i].density = stats[i].count/field_area
		endif
		field += 1
	endfor


	glactc, stats.mra, stats.mdec, 2000.0, ml, mb, 1, /degree
	stats.ml = ml
	stats.mb = mb

	return, stats

end

function esboss::btlogic, str, target_flags_or_type, $
		all=all, $
		only=only, $
		fixup=fixup, $
		objc_type=objc_type, $
		primary=primary, $
		run_primary=run_primary, $
		silent=silent

	if n_elements(str) eq 0 or n_elements(target_flags_or_type) eq 0 then begin
		on_error,2
		print,'Usage: logic=eb->btlogic(struct, target_flags_or_type, objc_type=, /primary, /run_primary, /all, /only, /fixup, /silent, count=)'
		print,'Send /all to require all the required flags are set'
		print,'Send /only to require all the required flags are set, and only those are set'
		message,'Halting'
	endif

	; these are flags that must *not* be set
	notflags = 0L
	logic = replicate(1, n_elements(str))

	if strmatch(target_flags_or_type[0], 'qso*') then begin

		target_type = 'qso'
		if target_flags_or_type[0] eq 'qso' then begin
			; the generic 'qso' type was sent
			;target_flags=$
			;	['qso_nn',$
			;	 'qso_like',$
			;	 'qso_kde', $
			;	 'qso_core_main', $
			;	 'qso_bonus_main', $
			;	 'qso_known_midz',$
			;	 'qso_first_boss']
			 target_flags=$
				 ['qso_core_main', $
				  'qso_bonus_main', $
				  'qso_known_midz',$
				  'qso_first_boss']
		endif else begin
			target_flags = target_flags_or_type
		endelse


		; never target qso_known_lohiz
		notflags += sdss_flagval('boss_target1','qso_known_lohiz')
	endif else if ($
			strmatch(target_flags_or_type[0], 'lrg*') $
			or strmatch(target_flags_or_type[0], 'gal*') ) then begin

		target_type = 'lrg'
		if target_flags_or_type[0] eq 'lrg' then begin
			; generic type 'lrg' given
			target_flags = ['gal_loz','gal_cmass','gal_cmass_sparse','gal_cmass_all']
		endif else begin
			target_flags = target_flags_or_type
		endelse

	endif else if strmatch(target_flags_or_type[0],'std*') then begin

		target_type='std'
		if target_flags_or_type[0] eq 'std' then begin
			; generic 'std' type sent
			target_flags = ['std_fstar','std_wd','std_qso']
		endif else begin
			target_flags = target_flags_or_type
		endelse

	endif

	if not keyword_set(silent) then begin
		splog,'Selecting target type "',target_flags_or_type,'"';, form='(a,a,a)'
	endif

	; select primary
	if keyword_set(primary) then begin
		if not keyword_set(silent) then begin
			splog,'  Only primary'
		endif
		primary = sdss_flagval('resolve_status','survey_primary')
		logic = logic and (str.resolve_status and primary) ne 0
	endif else if keyword_set(run_primary) then begin
		if not keyword_set(silent) then begin
			splog,'  Only run_primary'
		endif
		run_primary = sdss_flagval('resolve_status','run_primary')
		logic = logic and (str.resolve_status and run_primary) ne 0
	endif

	; only select objc_type if requested
	if n_elements(objc_type) ne 0 then begin
		if not keyword_set(silent) then begin
			splog,'  Only objc_type = ',objc_type,form='(a,i0)'
		endif
		logic = logic and str.objc_type eq objc_type
	endif


	orflags = sdss_flagval('boss_target1',target_flags)
	logic = (str.boss_target1 and notflags) eq 0

	if keyword_set(only) then begin
		; all the flags are set and no others are set
		logic = logic and (str.boss_target1 eq orflags)
	endif else if keyword_set(all) then begin
		; all the flags are set
		logic = logic and (str.boss_target1 and orflags) eq orflags
	endif else begin
		; any of the flags are set
		logic = logic and (str.boss_target1 and orflags) ne 0
	endelse

	if keyword_set(fixup) then begin
		; add some more restrictive cuts
		if n_elements(target_flags) ne 1 then message,'/fixup is only for single flag'
		if target_flags[0] eq 'qso_bonus' then begin
			faint=str.gfaint ne 0
			bright=str.gfaint eq 0

			newlogic = $
				(str.gfaint ne 0 $
				and alog10(str.kde_qsodens_faint) ge -0.57 $
				and alog10(str.kde_stardens_faint) lt -0.52) $
				or $
				(str.gfaint eq 0 $
				and alog10(str.kde_qsodens_bright) ge -0.57 $
				and alog10(str.kde_stardens_bright) lt -0.095)
			logic = logic and newlogic
		endif else if target_flags[0] eq 'qso_like' then begin
			newlogic = str.like_ratio gt 0.255
			logic = logic and newlogic
		endif else if target_flags[0] eq 'qso_nn' then begin
			newlogic = str.nn_xnn gt 0.781
			logic = logic and newlogic
		endif
	endif


	return, logic

end

function esboss::btselect, str, target_flags_or_type, $
		objc_type=objc_type, $
		all=all, $
		only=only, $
		primary=primary, $
		count=count, $
		complement=complement, $
		ncomplement=ncomplement, $
		fixup=fixup, $
		silent=silent

	if n_elements(str) eq 0 or n_elements(target_flags_or_type) eq 0 then begin
		on_error,2
		print,'Usage: w=eb->btselect(struct, target_flags_or_type, objc_type=, /primary, /all, /only, /silent, count=)'
		print,'Send /all to require all the required flags are set'
		print,'Send /only to require all the required flags are set, and only those are set'
		message,'Halting'
	endif

	if keyword_set(fixup) then begin
		if keyword_set(all) or keyword_set(only) then begin
			message,'fixup only works for ored logic currently'
		endif
		logic = replicate(1L, n_elements(str))
		for i=0L, n_elements(target_flags_or_type)-1 do begin
			tlogic = self->btlogic(str, target_flags_or_type[i], $
				all=all, only=only, $
				objc_type=objc_type, $
				primary=primary, fixup=fixup, silent=silent)

			if i eq 0 then begin
				logic=tlogic
			endif else begin
				logic = tlogic or logic
			endelse
		endfor
	endif else begin
		logic = self->btlogic(str, target_flags_or_type, $
			all=all, only=only, $
			objc_type=objc_type, $
			silent=silent, $
			primary=primary, fixup=fixup)

	endelse

	wtargets = where(logic, count, comp=complement, ncomp=ncomplement)
	tf = '['+strjoin(string(target_flags_or_type), ', ')+']'
	if not keyword_set(silent) then begin
		splog,'Found:',count,tf, form='(a," ",i0," ",a)'
	endif
	return, wtargets

end


function esboss::observed_platedir
	dir=getenv('BOSS_TARGET')
	return, path_join(dir, 'observed-plates')
end
function esboss::observed_platefile, chunk=chunk, orig=orig
	dir=self->observed_platedir()
	if keyword_set(orig) then begin
		fname='plate-data-orig.dat'
	endif else begin
		if n_elements(chunk) eq 0 then message,'Enter chunk(s)'
		fname=string('plate-data-chunk',chunk,'.par',f='(a,i0,a)')
	endelse
	return, path_join(dir, fname)
end
function esboss::read_observed_plates, chunk=chunk, orig=orig
	file=self->observed_platefile(orig=orig, chunk=chunk)
	if keyword_set(orig) then begin
		read_struct, file, {plateid:0L, ra:0d, dec:0d}, outstruct
	endif else begin
		outstruct=yanny_readone(file)
	endelse
	return, outstruct
end


pro esboss::observed_plates_addarea, chunk
	data = self->read_observed_plates(/orig)
	poly=self->bosstile_poly_read(chunk)
	areas = self->circle_poly_area(data.ra, data.dec, poly)
	colprint,areas,areas/(!pi*1.49^2)

	add_tags, data, ['area'], ['0d'], newstruct
	newstruct.area = areas
	file=self->observed_platefile(chunk=chunk)
	print,'Writing to file: ',file

	stname=string('CHUNK',chunk,'OBSERVED',f='(a,i0,a)')
	yanny_write, file, ptr_new(newstruct), stnames=stname
end

function esboss::mmt_plateinfo
	mmtdir = getenv("BOSSTARGET_DIR")
	mmtdir=path_join(mmtdir, "data")

	mmtname = path_join(mmtdir, "mmt-plate-data.dat")

	stdef = {ra:0d, dec:0d, rad:0d}
	read_struct, mmtname, stdef, tplateinfo, skiplines=1

	add_tags, $
		tplateinfo, $
		['clambda','ceta','area','platenum'], $
		['0d','0d','0d','0L'], plateinfo
	plateinfo.area = !dpi*plateinfo.rad^2
	plateinfo.platenum = lindgen(n_elements(plateinfo))
	eq2csurvey, tplateinfo.ra, tplateinfo.dec, clambda, ceta
	plateinfo.clambda=clambda
	plateinfo.ceta=ceta

	return, plateinfo
end
function esboss::mmt_match_dir, type, all=all
	boss_target = getenv("BOSS_TARGET")
	dir = path_join(boss_target, subdir=["esheldon","mmt"])
	dir = path_join(dir, 'sweep-match-'+type)

	if keyword_set(all) then begin
		dir = dir +'-all'
	endif
	return, dir
end
function esboss::mmt_match_psdir, type, all=all
	dir = '~/public_html/bosstarget/density/mmt'
	return, dir
end
function esboss::mmt_match_file, type, index, all=all, targets=targets

	dir = self->mmt_match_dir(type, all=all)
	if keyword_set(all) then addstr='all-' else addstr=''

	if not keyword_set(targets) then begin
		name = string($
			'mmt-sweep-match-',addstr,type,index,'.fits', $
			format='(a,a,a,i02,a)')	
	endif else begin
		name='mmt-sweep-match-'+addstr+type+'-target.fits'
	endelse

	name = path_join(dir, name)
	return, name
end

pro esboss::mmt_combine_and_target, type, str, all=all

	outfile = self->mmt_match_file(type, 0, all=all)
	outfile=repstr(outfile,'00.fits','-target.fits')
	print,'Will write to file: ',outfile
	
	plateinfo = self->mmt_plateinfo()

	; read the files
	nplates = n_elements(plateinfo)
	
	print,'Reading the data'
	if n_elements(str) eq 0 then begin
		for i=0L, nplates-1 do begin
			f=self->mmt_match_file(type, i, all=all)
			print,'Attepting to read file: ',f
			t=mrdfits(f,1,status=status)
			if status eq 0 then begin
				add_tags, t, 'platenum', '0L', newt
				newt.platenum = i
				add_arrval, ptr_new(newt,/no_copy), ptrlist
			endif 
		endfor

		str=combine_ptrlist(ptrlist)
	endif

	print,'Running qso target selection'

	bq=obj_new('bosstarget_qso')

	res=bq->select(str, /struct)

	add_tags, str, $
		['boss_target1','known_qso_matchflags','known_qso_id'],$
		['0L','0L','0L'], $
		newstr
	newstr.boss_target1 = res.boss_target1
	newstr.known_qso_matchflags = res.known_qso_matchflags 
	newstr.known_qso_id = res.known_qso_id

	print,'Writing to file: ',outfile
	mwrfits, newstr, outfile, /create
end

pro esboss::mmt_density, target_flags, str=str, bin=bin

	substring = strjoin(target_flags, ' or ')
	subfilestring = strjoin(target_flags, '-')
	subfilestring=repstr(subfilestring,'_','-')


	f=self->mmt_match_file('star', /target)	
	if n_elements(str) eq 0 then begin
		print,'reading: ',f
		str=mrdfits(f,1)
	endif

	plateinfo=self->mmt_plateinfo()
	
	dir=self->mmt_match_dir('star')
	psdir=self->mmt_match_psdir('star')
	file_mkdir, psdir

	foot_psfile=repstr(f, '.fits', '-footprint-'+subfilestring+'.eps')
	hist_psfile=repstr(f, '.fits', '-densities-'+subfilestring+'.eps')
	foot_psfile = repstr(foot_psfile, dir, psdir)
	hist_psfile = repstr(hist_psfile, dir, psdir)

	wtarg=self->btselect(str, target_flags, /primary)
	h=histogram(str[wtarg].platenum, min=0, max=max(plateinfo.platenum))

	w=where(h gt 0)
	densities = h[w]/plateinfo[w].area

	print,'Densities for ',substring
	print,densities

	med_density = median(densities)
	mean_density = mean(densities)
	sdev_density = stdev(densities)

	; plots

	begplot, foot_psfile, /encap, xsize=18, ysize=3, /color
	eq2csurvey, str[wtarg].ra, str[wtarg].dec, clambda, ceta
	plot, clambda, ceta, psym=3, /ynozero, $
		xtitle=textoidl('\lambda_c'), ytitle=textoidl('\eta_c'), $
		iso=1, charsize=1, yrange=[140,155], ystyle=3

	pplot, plateinfo.clambda, plateinfo.ceta, psym=7, symsize=0.5, thick=2,$
		color='red',/over
	for i=0L, n_elements(plateinfo)-1 do begin
		tvcircle, plateinfo[i].rad, plateinfo[i].clambda, plateinfo[i].ceta, $
			/data, color=c2i('red'), thick=2
	endfor
	legend, substring



	endplot

	if n_elements(bin) eq 0 then bin=7
	begplot, hist_psfile, /encap, xsize=8.25, ysize=6.375
	plothist, densities, bin=bin, $
		xtitle='Density (#/square degree)'

	legend, substring

	legend,$
		['median density: '+strn(round(med_density)), $
		'mean density: '+strn(round(mean_density)), $
		'sdev density: '+strn(round(sdev_density))], $
		/right



	endplot



end







pro esboss::make_run_html_run, reload=reload

	common plot_run_density_run_blk, t745, t1345

	dir='/mount/early1/bosstarget'

	if n_elements(t745) eq 0 or keyword_set(reload) then begin
		bt=obj_new('bosstarget')
		t745 = bt->read('lrg', 745, 1, target_dir=dir, /collate)
		t1345 = bt->read('lrg', 1345, 1, target_dir=dir, /collate)
	
		obj_destroy, bt
	endif

	self->make_run_html, t1345, ['gal_loz']
	self->make_run_html, t1345, ['gal_cmass']
	self->make_run_html, t1345, ['gal_cmass_notred']

	self->make_run_html, t745, ['gal_loz']
	self->make_run_html, t745, ['gal_cmass']
	self->make_run_html, t745, ['gal_cmass_notred']
end


function esboss::orflags, flags, types
	logic = (flags and sdss_flagval('boss_target1',types[0])) ne 0
	for i=1L, n_elements(types)-1 do begin
		logic = logic or (flags and sdss_flagval('boss_target1',types[i])) ne 0
	endfor

	return, logic
end

function esboss::andflags, flags, types, notset=notset
	flagtot = sdss_flagval('boss_target1',types)

	if not keyword_set(notset) then begin
		logic = (flags and flagtot) eq flagtot
	endif else begin
		logic = (flags and flagtot) ne flagtot
	endelse
	return, logic
end



function esboss::get_field_stats, struct

	field_area = (2048d*1361d)*(0.396d/60d/60d)^2

	h=histogram(struct.field, $
		min=min(struct.field), max=max(struct.field), rev=rev)

	nh = n_elements(h)

	stdef = {$
		field: 0L, $
		mra:0d, mdec:0d, ml:0d, mb:0d, $
		mseeingi:-9999., $
		count:0L, $
		density:0d}

	stats = replicate(stdef, nh)

	field = min(struct.field)

	for i=0L, nh-1 do begin
		stats[i].field = field
		if rev[i] ne rev[i+1] then begin
			wf = rev[ rev[i]:rev[i+1]-1 ]

			stats[i].mra = median_check(struct[wf].ra)
			stats[i].mdec = median_check(struct[wf].dec)
			stats[i].mseeingi = median_check(struct[wf].psf_fwhm[3])
			stats[i].count = h[i]
			stats[i].density = stats[i].count/field_area
		endif
		field += 1
	endfor


	glactc, stats.mra, stats.mdec, 2000.0, ml, mb, 1, /degree
	stats.ml = ml
	stats.mb = mb

	return, stats

end

function esboss::bosstile_radius
	return, 1.49d
end



;+
; NAME:
;   circle_poly
; PURPOSE:
;   Create the polygons needed to describe the input circle(s).
; CALLING SEQUENCE:
;   poly=circle_poly(ra,dec,tilerad)
; INPUTS:
;   ra, dec - Center of the circle. (J2000 deg) Can be arrays.
;   tilerad - radius of circle in degrees.
; OUTPUTS:
;   Polygons describing the circle. If ra/dec are arrays so this this.
; REVISION HISTORY:
;   5-Oct-2009 MRB, NYU
;-

function esboss::circle_poly, ra, dec, tilerad

	pi=!DPI
	ncaps=1L

	cap1=replicate(construct_cap(), ncaps)
	polystr1= construct_polygon(ncaps=1,nelem=1)
	poly1=create_struct($
		polystr1, $
		'locationid',-1L, $
		'areaname', ' ', $
		'ra',0.D, $
		'dec',0.D)

	poly=replicate(poly1, n_elements(ra))
	for i=0L, n_elements(ra)-1L do $
		poly[i].caps=ptr_new(cap1)
	poly.ra= ra
	poly.dec= dec
	poly.ncaps=1
	poly.weight=1.
	set_use_caps,poly,[0]

	for i=0L, n_elements(ra)-1L do begin
		(*poly[i].caps)[0]=circle_cap(ra=ra[i], dec=dec[i], tilerad)
		poly[i].str= garea(poly[i])
	endfor

	return, poly

end

;+
; NAME:
;   circle_poly_area
; PURPOSE:
;   Calculate the area of one or more tiles within some geometry.  Only
;		the fraction inside the goemetry is counted.  If the geometry
;		is *not* sent then just the plain areas are returned.
; CALLING SEQUENCE:
;   area= circle_poly_area(ra, dec [, geom, tilerad=])
; INPUTS:
;   ra, dec - [Ntiles] center of tile(s) (J2000 deg)
; OPTIONAL INPUTS:
;   geom - [Ngeom] polygons describing area (full tile area given
;          otherwise)
;   tilerad - radius of tiles to assume (default 1.49 deg)
; OUTPUTS:
;   area - [Ntiles] area in square deg
; REVISION HISTORY:
;   5-Oct-2009 MRB, NYU
;-


function esboss::circle_poly_area, ra, dec, geom, tilerad=tilerad;, total_area=total_area

	; find the area of a circle thatis contained in the window
	; geometry.  If ra/dec are arrays, then an array of areas
	; is returned
	; if geom is not given, then just add up the area of all
	; the tiles

	if(n_elements(ra) eq 0) then $
		message, 'Must input RA and DEC'
	if(n_elements(ra) ne n_elements(dec)) then $
		message, 'RA and DEC must have same size'
	if(n_elements(tilerad) gt 1) then $
		message, 'TILERAD input must be a scalar'
	if(n_elements(tilerad) eq 0) then $
		tilerad=self->bosstile_radius()

	; this doesn't work
	;if arg_present(total_area) then begin
	;	circpoly= self->circle_poly(ra, dec, tilerad)
	;	where_polygons_overlap, $
	;		circpoly, geom, geom_match_ind, nmatch, $
	;		areamatch=areamatch
	;	total_area=total(areamatch, /double)*(180.D/!DPI)^2
	;endif


	area= dblarr(n_elements(ra))
	for i=0L, n_elements(area)-1L do begin
		circpoly= self->circle_poly(ra[i], dec[i], tilerad)
		if(keyword_set(geom)) then begin
			where_polygons_overlap, $
				circpoly, geom, geom_match_ind, nmatch, $
				areamatch=areamatch
			area[i]=total(areamatch, /double)*(180.D/!DPI)^2
		endif else begin
			area[i]=circpoly.str*(180.D/!DPI)^2
		endelse
		destruct_polygon, circpoly

	endfor

	return, area

end

function esboss::tile_area_old, ra, dec, geom, tilerad=tilerad

	; this Jeremy's version
	if(n_elements(ra) eq 0) then $
		message, 'Must input RA and DEC'
	if(n_elements(ra) ne n_elements(dec)) then $
		message, 'RA and DEC must have same size'
	if(n_elements(tilerad) gt 1) then $
		message, 'TILERAD input must be a scalar'
	if(n_elements(tilerad) eq 0) then $
		tilerad=self->bosstile_radius()

	xx = vmid(geom)
	x_to_angles, xx, phi, theta
	ra_geom = transpose(phi)
	dec_geom = (90.D)-transpose(theta)

	ntiles=n_elements(ra)
	area = dblarr(ntiles)

	for i=0L, ntiles-1 do begin
		;; make the circle polygon for each tile here, call it "loc"

		tilepoly = self->circle_poly(ra[i], dec[i], tilerad)

		; this only works because we're at the center of the circle,
		; it is not general
		ii = is_in_window(tilepoly, ra=ra_geom, dec=dec_geom)
		area[i] = total(geom[where(ii gt 0)].str)*(180d/!dpi)^2
		
		destruct_polygon, tilepoly
	endfor

	return, area
end


pro esboss::_random_radec_circle, ra, dec, R, rand_ra, rand_dec
	; R in degrees

	rrad = R*!dpi/180d

	n=n_elements(ra)

	; generate uniformly in R^2 (radians)
	rand_r = randomu(seed, n)

	rand_r = sqrt(rand_r)*rrad

	; generate theta uniformly from 0 to 2*pi
	rand_psi = 2d*!dpi*randomu(seed, n)


	cospsi = cos( rand_psi )

	

end

function esboss::bosstile_read, chunk

	dir=getenv('BOSSTILELIST_DIR')
	if dir eq '' then message,'bosstilelist is not setup'

	chunkstr=string(chunk,f='(i0)')
	bossname='boss'+chunkstr

	fname = 'tiles-'+bossname+'.par'
	file = path_join([dir,'outputs',bossname,fname])

	;print,'Reading tile file: ',file
	tiles = yanny_readone(file)
	if n_tags(tiles) eq 0 then message,'Could not read file: '+file

	tiles = rename_tags(tiles, ['racen','deccen'], ['ra','dec'])

	eq2csurvey, tiles.ra, tiles.dec, clambda, ceta
	glactc, tiles.ra, tiles.dec, 2000.0, l, b, 1, /degree


	tiles = struct_addtags($
		tiles, $
		['clambda','ceta','l','b'],$
		['0d','0d','0d','0d'])
	tiles.clambda=clambda
	tiles.ceta=ceta
	tiles.l=l
	tiles.b=b


	nchunks = n_elements(chunks)
	if nchunks ne 0 then begin
		match, tiles.chunk, [chunks], mtiles, mchunks, /sort
		if mtiles[0] eq -1 then message,'No matches found to input chunks'
		tiles=tiles[mtiles]
	endif

	return, tiles


end


pro esboss::bosstile_match_poly, $
		ra, dec, chunks, $
		allmatch_obj, allmatch_tiles, unique_obj, unique_tiles, $
		tiles=tiles, tilepoly=tilepoly

	if n_elements(chunks) eq 0 or n_elements(ra) eq 0 $
			or n_elements(dec) eq 0 then begin
		on_error, 2
		print,'Usage: '
		print,'  eb->bosstile_match_poly,'
		print,'       ra, dec, chunks, '
		print,'       allmatch_obj, allmatch_tiles, unique_obj, unique_tiles, '
		print,'       tiles=, tilepoly=, '
		message,'Halting'
	endif

	; first trim objects to the window, then match to the
	; tile circle centers to get tile matches


	tiles=self->bosstile_read(chunks)

	; polygon descriptions of the tiles and overall region
	tilepoly=self->bosstile_poly_read(chunks)


	inwindow=is_in_window(tilepoly, ra=ra, dec=dec)
	inwindow = where(inwindow, nwindow)

	allmatch_obj=-1
	allmatch_tiles=-1
	unique_obj=-1
	unique_tiles=-1
	if nwindow ne 0 then begin

		; the match lists will contain duplicats, objects can be in more
		; than one tile.
		spherematch, $
			ra[inwindow], dec[inwindow], $
			tiles.ra, tiles.dec,  $
			self->bosstile_radius(), $
			allmatch_obj, allmatch_tiles, distances, maxmatch=0

		if allmatch_tiles[0] ne -1 then begin
			allmatch_obj = inwindow[allmatch_obj]
			unique_obj = allmatch_obj[rem_dup(allmatch_obj)]

			unique_tiles=allmatch_tiles[rem_dup(allmatch_tiles)]
		endif
	endif

	if not arg_present(tilepoly) then destruct_polygon, tilepoly

end


pro esboss::make_inwindow, target_type, target_run, chunk, $
		extra_name=extra_name, bounds=bounds, trim=trim
	tilepoly=self->bosstile_poly_read(chunk, bounds=bounds)

	bt=obj_new('bosstarget')
	t=bt->read_collated(target_type, target_run, $
		file=file, status_struct=status_struct, extra_name=extra_name)

	add_tags, t, 'inwindow', '0', newt
	newt.inwindow = is_in_window(tilepoly, ra=t.ra, dec=t.dec)
	print,'Area = ',total(tilepoly.str)*(180d/!dpi)^2

	chstr=strjoin( string(chunk,f='(i0)'), '-')
	newfile=repstr(file, '.fits', '-inchunk'+chstr+'.fits')
	if keyword_set(trim) then begin
		newfile=repstr(newfile,'.fits','-trim.fits')
		w=where(newt.inwindow,nw)
		newt = newt[w]
	endif
	if newfile eq file then message,'error making file name'

	print,'Writing new file: ',newfile
	mwrfits, status_struct, newfile, /create
	mwrfits, newt, newfile
end

pro esboss::make_fullcache_and_inwindow, target_type, target_run, chunk, $
		extra_name=extra_name, restrict=restrict

	tilepoly=self->bosstile_poly_read(chunk)

	bt=obj_new('bosstarget')
	struct=bt->read_collated(target_type, target_run, $
		file=file, status_struct=status_struct, extra_name=extra_name)

	photo_sweep = status_struct.photo_sweep
	print,'photo_sweep = ',photo_sweep
	sold=getenv('PHOTO_SWEEP')
	setenv,'PHOTO_SWEEP='+photo_sweep


	print,'Adding full like cache'
	bq=obj_new('bosstarget_qso')
	bq->get_unique_run_rerun_camcol, $
		struct.run, struct.rerun, struct.camcol, $
		urun, urerun, ucamcol

	nrun=n_elements(urun)
	for i=0L, nrun-1 do begin
		if bq->cache_exists('like',urun[i],urerun[i],ucamcol[i]) then begin
			ls=bq->cache_read('like',urun[i],urerun[i],ucamcol[i])

			if n_elements(newstruct) eq 0 then begin
				ignoretags=['run','rerun','camcol','field','id']
				ls0 = ls[0]
				for j=0,n_tags(ls0)-1 do begin
					ls0.(j) = -9999
				endfor
				tmp = remove_tags(ls[0], ignoretags)
				newstruct = create_struct(struct[0], tmp)
				newstruct = create_struct(newstruct, 'inwindow',0)
				newstruct=replicate(newstruct, n_elements(struct))
				struct_assign, struct, newstruct, /nozero
			endif

			sphoto_match, newstruct, ls, mstruct, mls
			if mstruct[0] ne -1 then begin
				tmpstruct = newstruct[mstruct]
				struct_assign, ls[mls], tmpstruct, /nozero
				newstruct[mstruct] = tmpstruct
			endif
		endif
	endfor

	setenv,'PHOTO_SWEEP='+sold

	print,'Adding inwindow'
	eb=obj_new('esboss')
	newstruct.inwindow = $
		self->bosstile_is_in_window(struct.ra,struct.dec,chunk)

	if keyword_set(restrict) then begin
		newstruct = newstruct[where(newstruct.inwindow)]
	endif

	newfile=repstr(file, '.fits', '-fullcache.fits')
	if newfile eq file then message,'error making file name'

	print,'Writing new file: ',newfile
	mwrfits, status_struct, newfile, /create
	mwrfits, newstruct, newfile

end




pro esboss::bosstile_match, struct, mstruct, mtile, distances, $
		chunks=chunks, tiles=tiles, tile_subset=tile_subset
	if n_elements(struct) eq 0 then begin
		on_error, 2
		message,'Usage: eb->bosstile_match, struct, mstruct, mtiles, dist, chunks=, tiles=]'
	endif

	tiles=self->bosstile_read(chunks=chunks)

	if n_elements(tile_subset) ne 0 then begin
		;splog,'Using a subset of the tiles'
		tiles=tiles[tile_subset]
	endif
	spherematch, struct.ra, struct.dec, tiles.ra, tiles.dec,  $
		self->bosstile_radius(), $
		mstruct, mtile, distances, maxmatch=0

end

function esboss::bosstile_chunk_bounds_polygon_file, chunk
	; describes the overall boundary
	dir=getenv('BOSSTILELIST_DIR')
	if dir eq '' then message,'bosstilelist is not setup'
	
	bosschunk=string('boss',chunk,f='(a,i0)')
	fname='trim-'+bosschunk+'.ply'
	;fname='trim-'+bosschunk+'.fits'

	file=path_join(dir, ['inputs',bosschunk, fname])
	return, file
end

function esboss::bosstile_chunk_polygon_file, chunk

	bt=obj_new('bosstarget')
	return, bt->chunk_polygon_file(chunk)

end

function esboss::bosstile_poly_read, chunk, bounds=bounds, verbose=verbose
	
	bt=obj_new('bosstarget')
	return, bt->chunk_polygon_read(chunks,bounds=bounds,verbose=verbose)

end



function esboss::bosstile_is_in_window, ra, dec, chunk, bounds=bounds
	if (n_elements(ra) eq 0 $
		  or n_elements(dec) eq 0 $
		  or n_elements(chunk) eq 0) then begin
		message,'usage: inwindow=eb->'+$
			'bosstile_is_in_window(ra,dec,chunk,/bounds)'
	endif

	tilepoly=self->bosstile_poly_read(chunk)
	inwindow=is_in_window(tilepoly, ra=ra, dec=dec, bounds=bounds)
	destruct_polygon, tilepoly
	return, inwindow
end


function esboss::bosstile_inwindow, ra, dec, chunk, bounds=bounds, count=count
	if (n_elements(ra) eq 0 $
		  or n_elements(dec) eq 0 $
		  or n_elements(chunk) eq 0) then begin
		message,'usage: inwindow=eb->'+$
			'bosstile_inwindow(ra,dec,chunk,/bounds,count=)'
	endif

	tilepoly=self->bosstile_poly_read(chunk,bounds=bounds)
	inwindow=where(is_in_window(tilepoly, ra=ra, dec=dec), count)
	destruct_polygon, tilepoly
	return, inwindow
end

function esboss::bosstile_mean_density, ra, dec, chunk, bounds=bounds
	if (n_elements(ra) eq 0 $
		  or n_elements(dec) eq 0 $
		  or n_elements(chunk) eq 0) then begin
		message,'usage: density=eb->bosstile_mean_density(ra,dec,chunk, bounds=bounds)'
	endif

	tilepoly=self->bosstile_poly_read(chunk,bounds=bounds)
	inwindow=self->bosstile_inwindow(ra, dec, chunk, count=count)

	total_area=total(tilepoly.str)*(180d/!dpi)^2
	density = count/total_area
	destruct_polygon, tilepoly
	return, density
end


function esboss::bosstile_density, ra, dec, chunks, $
		getmean=getmean, $
		$
		tiles=tiles, $
		matched_tiles=matched_tiles, $
		$
		matched_objects=matched_objects, $
		nmatched_objects=nmatched_objects, $
		$
		total_area=total_area, $
		mean_density=mean_density, $
		$
		rasort=rasort, $
		csort=csort, $
		$
		tilepoly=tilepoly

	if n_params() lt 3 then begin
		on_error, 2
		message,'Usage: density_struct=eb->bosstile_density(struct, chunks=, tiles=, matched_tiles=, matched_objects=, tilepoly=, total_area=)'
	endif

	if n_elements(chunks) eq 0 then begin
		message,'chunks= is now a required keyword'
	endif
	if n_elements(chunks) ne 1 then message,'chunks must be length 1 for now'


	self->bosstile_match_poly, $
		ra, dec, chunks, $
		allmatch_obj, allmatch_tiles, $ ; contain dups
		matched_objects, matched_tiles, $ ;these are the unique ones
		tiles=tiles, tilepoly=tilepoly

	total_area=total(tilepoly.str)*(180d/!dpi)^2

	nmatched_objects=0L
	if matched_objects[0] ne -1 then begin
		nmatched_objects=n_elements(matched_objects)
	endif

	mean_density=nmatched_objects/total_area



	; just return the average density
	if keyword_set(getmean) then begin
		if not arg_present(tilepoly) then destruct_polygon, tilepoly
		return, mean_density
	endif





	; output structure
	; add the area calculated from the polygons for this chunk
	; this will help with fractional polygon coverage
	tile_areas = self->circle_poly_area(tiles.ra, tiles.dec, tilepoly)


	density_struct = struct_addtags($
		tiles, $
		['area','count','density'],$
		['0d','0LL','0d'])
	density_struct.area = tile_areas



	if nmatched_objects ne 0 then begin
		ntile=n_elements(tiles)
		bs=binner(allmatch_tiles, rev=rev, min=0, max=ntile-1)

		density_struct.count = bs.hist
		density_struct.density = bs.hist/density_struct.area
	endif



	; clean up memory
	if not arg_present(tilepoly) then destruct_polygon, tilepoly


	if keyword_set(csort) then begin
		s=sort(density_struct.clambda)
		density_struct = density_struct[s]
	endif else if keyword_set(rasort) then begin
		s=sort(density_struct.ra)
		density_struct = density_struct[s]
	endif
	return, density_struct

end



; tune the qso_core_main + qso_bonus_main to be 40/sq degree in 
; each of the subregions.  If the density in core is already >= 40
; then zero bonus are added.

pro esboss::subreg_density_tune, $
		target_run, chunks, maxres, $
		force=force, $
		str=str, status_struct=status_struct


	bt=obj_new('bosstarget')
	tfile=bt->target_file('qso',target_run,/collate)
	new_target_file=repstr(tfile,'.fits','-maskngc.fits')
	print,'will write to new target file: ',new_target_file
	if n_elements(str) eq 0 or n_elements(status_struct) eq 0 then begin
		str=bt->read_collated('qso',target_run, status_struct=status_struct)
	endif

	target_density = 40.0

	; memory managed internally by the cache mechanism, do not destroy
	self->bosstile_sub_regions, chunks, maxres, $
		polygons, regions, geometry

	fname = self->density_tune_file(target_run, chunks, maxres)


	nreg = n_elements(regions)


	coreflag = sdss_flagval('boss_target1','qso_core_main')
	bonusflag = sdss_flagval('boss_target1','qso_bonus_main')

	extra_flagnames=['qso_known_midz','qso_first_boss']
	extra_flags = $
		sdss_flagval('boss_target1',extra_flagnames)

	print,'getting core,bonus, and extra logic'
	; the core sample
	core_logic = (str.boss_target1 and coreflag) ne 0
	; the bonus sample not including the core sample
	bonus_logic = $
		(str.boss_target1 and bonusflag) ne 0 $
		and (str.boss_target1 and coreflag) eq 0
	extra_logic = (str.boss_target1 and extra_flags) ne 0


	print,'Checking against window'
	inwin = is_in_window( polygons, ra=str.ra, dec=str.dec )
	win_core = where(inwin and core_logic, ncore)
	win_bonus = where(inwin and bonus_logic, ncore)
	win_extra = where(inwin and extra_logic, nextra)

	help, win_core, win_bonus, win_extra


	vl='-9999L'
	vd='-9999d'
	; highest threshold
	dthresh='0.97d'
	add_tags, regions, $
		['count','density','ncore','core_density','nbonus','bonus_density','thresh','ra','dec','clambda','ceta','l','b'], $
		[     vl,       vd,     vl,            vd,      vl,             vd, dthresh,  vd,   vd,       vd,   vd,  vd, vd],$
		density_struct

	total_count = 0L
	total_area = 0d

	keepflag = intarr(n_elements(str))
	
	; we always keep these flag types that are within window
	keepflag[win_extra] = 1

	for i=0L,nreg-1 do begin
		print,'   ',i+1,'/',nreg,form='(a,i0,a,i0)'

		winthis_core = where( is_in_window( $
			polygons[ regions[i].pstart:regions[i].pend], $
			ra=str[win_core].ra, dec=str[win_core].dec), nthis_core)

		densthis_core = nthis_core/regions[i].area

		density_struct[i].ncore = nthis_core
		density_struct[i].core_density = nthis_core/regions[i].area

		if nthis_core gt 0 then begin
			win_this = win_core[winthis_core]
		endif

		density_struct[i].nbonus = 0
		density_struct[i].bonus_density = 0
		if density_struct[i].core_density lt target_density then begin

			winthis_bonus = where( is_in_window( $
				polygons[ regions[i].pstart:regions[i].pend], $
				ra=str[win_bonus].ra, $
				dec=str[win_bonus].dec), nthis_bonus)

			if nthis_bonus ne 0 then begin
				winthis_bonus = win_bonus[winthis_bonus]

				nn_value = str[winthis_bonus].nn_value
				; reverse sort, will take the top n_need
				s=reverse( sort(nn_value) )

				dens_need =  $
					(target_density - density_struct[i].core_density)
				n_need = dens_need*regions[i].area
				n_need = long(n_need) < nthis_bonus > 1

				winthis_bonus = winthis_bonus[ s[0:n_need-1] ]

				density_struct[i].thresh = min(str[winthis_bonus].nn_value)

				add_arrval, winthis_bonus, win_this

				density_struct[i].nbonus = n_need
				density_struct[i].bonus_density = n_need/regions[i].area

			endif


		endif


		density_struct[i].count = $
			density_struct[i].ncore + density_struct[i].nbonus
		density_struct[i].density = $
			density_struct[i].count/regions[i].area

		total_count += density_struct[i].count
		total_area += regions[i].area

		if density_struct[i].count gt 0 then begin
			if density_struct[i].count eq 1 then begin
				mra = str[win_this].ra
				mdec = str[win_this].dec
			endif else begin
				mra=median(str[win_this].ra)
				mdec=median(str[win_this].dec)
			endelse
			eq2csurvey, mra, mdec, clambda,ceta
			glactc, mra, mdec, 2000.0, l, b, 1, /degree

			density_struct[i].ra=mra
			density_struct[i].dec=mdec
			density_struct[i].clambda=clambda
			density_struct[i].ceta=ceta
			density_struct[i].l = l
			density_struct[i].b = b

			keepflag[win_this] = 1
		endif

		delvarx, win_this
	endfor ; over regions

	print,'Writing cache file: ',fname
	hdr={total_count: total_count, $
		total_area:total_area, $
		total_density:total_count/total_area}
	;write_idlstruct, density_struct, fname, hdr=hdr

	fhdr=['']
	sxaddpar, fhdr, 'tcount',total_count
	sxaddpar, fhdr, 'tarea',total_area
	sxaddpar, fhdr, 'tdens',total_count/total_area
	;mwrfits, density_struct, repstr(fname,'.st','.fits'), fhdr,/create


	keep = where(keepflag, nkeep)
	print,'keeping total of ',nkeep,'/',n_elements(str),'targets', $
		format='(a,i0,a,i0,a)'
	keepstruct = str[keep]
	print,'writing new target file: ',new_target_file
	mwrfits, status_struct, new_target_file, /create
	mwrfits, keepstruct, new_target_file
end




function esboss::density_tune_file, target_run, chunks, maxres
	boss_target = getenv('BOSS_TARGET')
	dir=filepath(root=boss_target,sub='esheldon','density-tune')
	if not file_test(dir) then begin
		file_mkdir, dir
	endif
	chstr = [target_run,'subreg-density-tune']
	for i=0L, n_elements(chunks)-1 do begin
		chstr_tmp=strtrim(string(chunks[i]),2)

		chstr = [chstr, chstr_tmp]
	endfor

	chstr = [chstr, 'maxres'+strtrim(string(maxres),2)]

	fname = strjoin(chstr, '-')+'.st'
	fname=filepath(root=dir,fname)
	return, fname
end

function esboss::density_tune_read, target_run, chunks, maxres, hdr=hdr
	file=self->density_tune_file(target_run, chunks, maxres)
	if not file_test(file) then begin
		self->subreg_density_tune, target_run, chunks, maxres
	endif
	return, read_idlstruct(file, hdr=hdr)
end


function esboss::density_cache_file, chunks, maxres, nra
	boss_target = getenv('BOSS_TARGET')
	dir=filepath(root=boss_target,sub='esheldon','density-cache')
	if not file_test(dir) then begin
		file_mkdir, dir
	endif
	chstr = ['subreg-density-cache']
	for i=0L, n_elements(chunks)-1 do begin
		chstr_tmp=strtrim(string(chunks[i]),2)

		chstr = [chstr, chstr_tmp]
	endfor

	chstr = [chstr, 'maxres'+strtrim(string(maxres),2)]
	chstr = [chstr, strtrim(string(nra),2)]

	fname = strjoin(chstr, '-')+'.st'
	fname=filepath(root=dir,fname)
	return, fname
end

pro esboss::density_cache_write, $
		chunks, maxres, nra, density_struct, total_area, total_density

	fname = self->density_cache_file(chunks, maxres, nra)
	hdr = {total_area:total_area, total_density:total_density}
	print,'Writing density cache: ',fname
	dir=file_dirname(fname)
	if not file_test(dir) then begin
		file_mkdir, dir
	endif
	write_idlstruct, density_struct, fname, hdr=hdr

end

function esboss::density_cache_read, chunks, maxres, nra

	fname = self->density_cache_file(chunks, maxres, nra)
	print,'Reading density cache: ',fname
	density_struct = read_idlstruct(fname,hdr=hdr)
	outstruct = { $
		total_area: hdr.total_area, $
		total_density: hdr.total_density, $
		dstruct:density_struct}

	return, outstruct
end


pro esboss::cache_bosstile_subreg_density, $
		ra, dec, chunks, maxres, $
		force=force, $
		polygons=polygons

	; memory managed internally by the cache mechanism, do not destroy
	self->bosstile_sub_regions, chunks, maxres, $
		polygons, regions, geometry

	nra = n_elements(ra)
	fname = self->density_cache_file(chunks, maxres, nra)

	if keyword_set(force) or not file_test(fname) then begin

		print,'Getting overall density'	
		total_area = total(geometry.str)*(180d/!dpi)^2
		w=where(is_in_window(geometry, ra=ra, dec=dec),nw)
	
		total_density = nw/total_area

		nreg = n_elements(regions)

		vl='-9999L'
		vd='-9999d'
		add_tags, regions, $
			['count','density','ra','dec','clambda','ceta','l','b'], $
			[     vl,       vd,  vd,   vd,       vd,   vd,  vd, vd],$
			density_struct

		if nw gt 0 then begin
			print,'getting subreg densities'
			raw = ra[w]
			decw = dec[w]

			for i=0L,nreg-1 do begin
				print,'   ',i+1,'/',nreg,form='(a,i0,a,i0)'
				
				inwin = is_in_window( $
					polygons[ regions[i].pstart:regions[i].pend], $
					ra=raw, dec=decw)
				w2=where(inwin, nw2)

				density_struct[i].count = nw2
				density_struct[i].density = nw2/regions[i].area

				if nw2 gt 0 then begin
					if nw2 eq 1 then begin
						mra = raw[w2]
						mdec = decw[w2]
					endif else begin
						mra=median(raw[w2])
						mdec=median(decw[w2])
					endelse
					eq2csurvey, mra, mdec, clambda,ceta
					glactc, mra, mdec, 2000.0, l, b, 1, /degree

					density_struct[i].ra=mra
					density_struct[i].dec=mdec
					density_struct[i].clambda=clambda
					density_struct[i].ceta=ceta
					density_struct[i].l = l
					density_struct[i].b = b
				endif
			endfor
		endif

		print,'Writing cache file: ',fname
		hdr={total_area:total_area, $
			 total_density:total_density}
		write_idlstruct, density_struct, fname, hdr=hdr
	endif	
end

function esboss::bosstile_subreg_density, ra, dec, chunks, $
		maxres=maxres, $
		polygons=polygons, $
		force=force

	if n_elements(maxres) eq 0 then maxres=5
	self->cache_bosstile_subreg_density, $
		ra, dec, chunks, maxres, $
		polygons=polygons, force=force

	nra=n_elements(ra)
	outstruct = self->density_cache_read(chunks, maxres, nra)
	return, outstruct
end



pro esboss::bosstile_denscut_area, method, str=str, dens=dens


	common dpsdd60_block, oplot_chunk, isinwin

	goal=40

	plot_chunk = ['sgc','ngc']
	target_type='qso'
	; pre-tuned each method to ~60/sq degree on chunk 5. Note
	; the result for kde ended up 50 on the test chunk
	target_run='2010-03-03d60'


	if n_elements(oplot_chunk) eq 0 then begin
		oplot_chunk = 'none'
		isinwin=0
	endif



	bt=obj_new('bosstarget')
	if n_elements(str) eq 0 then begin
		; this is the all-sky version
		str=bt->read_collated(target_type,target_run,extra_name='all')
	endif
	ntot = n_elements(str)


	pcstr = strjoin( strtrim(string(plot_chunk),2) ,'-')
	opcstr = strjoin( strtrim(string(oplot_chunk),2) ,'-')
	if pcstr ne opcstr then begin
		ply=self->bosstile_poly_read(plot_chunk,/verbose)
		print,'finding inwindow'
		isinwin=is_in_window(ply, ra=str.ra, dec=str.dec )
		destruct_polygon, ply
		oplot_chunk=plot_chunk
	endif




	print,'trimming to thresh and chunk'
	thresh_struct = self->get_ngc_large_thresholds()




	wtype=where(thresh_struct.method eq method $
				and thresh_struct.density eq goal, ntype)
	if ntype ne 1 then begin
		message,string('Bad method/goal: ',method,goal,f='(a,a,"/",f0.3)')
	endif

	target_flags = thresh_struct[wtype].flagname

	if not tag_exist(str, thresh_struct[wtype].tagname, index=tagind) then begin
		message,'bad tag name: ',thresh_struct[wtype].tagname
	endif

	wuse = where( str.(tagind) gt thresh_struct[wtype].thresh $
		          and isinwin, nuse)
	if nuse eq 0 then begin
		message,'No objects passed cuts on thresh/intest_region'
	endif

	use_str = str[wuse]


	wstr_type = self->btselect(use_str, target_flags,count=nuse)
	use_str = use_str[wstr_type]


	print,nuse,ntot,f='(i0,"/",i0," passed thresh/intest_region")'

	dens = self->bosstile_subreg_density(use_str.ra,use_str.dec,plot_chunk)


	nthresh = 20
	dthresh = arrscl(findgen(nthresh), 40.0, 80.0)
	area = fltarr(nthresh)

	for i=0L, nthresh-1 do begin
		w=where(dens.dstruct.density lt dthresh[i])
		area[i] = total(dens.dstruct[w].area)
	endfor



	
	; area vs threshold
	psdir=filepath( $
		root=getenv('BOSS_TARGET'), subdir=[target_run[0]],'plots')
	mstr=repstr(method,'_','-')
	chstr = strjoin(plot_chunk, '-')
	psfile=[target_run,target_type,chstr,mstr,'areathresh']
	psfile=strjoin(psfile,'-')+'.eps'
	psfile=path_join(psdir, psfile)

	print,'psfile: ',psfile

	begplot,psfile,/color,ysize=8,xsize=8.5
	!p.charsize=!p.charsize*1.3

	plot, dthresh, area, $
		xtitle='density threshold', ytitle=textoidl('area [deg^2]')
	plegend,method

	endplot,/trim,/png



	; contour plot in the ngc



	psfile=[target_run,target_type,chstr,mstr,'ngc-contour']
	psfile=strjoin(psfile,'-')+'.eps'
	psfile=path_join(psdir, psfile)

	print,'psfile: ',psfile

	begplot,psfile,/color,xsize=11,ysize=11
	!p.charsize=!p.charsize*1.3




	w=where(dens.dstruct.ra gt 100 and dens.dstruct.ra lt 270)


	ply_ngc = self->bosstile_poly_read('ngc')
	self->plot_poly, ply_ngc, color=c2i('black'),xstyle=3,ystle=3
	destruct_polygon, ply_ngc


	levels=[50,60,70,80,90]
	nlevel=n_elements(levels)
	c_color = (make_rainbow(nlevel+1))[1:nlevel]
	contour, $
		/irreg, /over, $
		dens.dstruct[w].density, dens.dstruct[w].ra, dens.dstruct[w].dec, $
		levels=levels, c_color=c_color


	plegend,method
	plegend,string(levels,f='(i0)'),line=0, color=c_color,/right

	endplot,/trim,/png


	; contour plot in the sgc
	psfile=[target_run,target_type,chstr,mstr,'sgc-contour']
	psfile=strjoin(psfile,'-')+'.eps'
	psfile=path_join(psdir, psfile)

	print,'psfile: ',psfile

	begplot,psfile,/color,xsize=11,ysize=11
	!p.charsize=!p.charsize*1.3

	sra= shiftra(dens.dstruct.ra, 180)
	w=where(sra gt 130 and sra lt 230)


	ply_sgc = self->bosstile_poly_read('sgc')
	self->plot_poly, ply_sgc, color=c2i('black'),xstyle=3,ystle=3, offset=180
	destruct_polygon, ply_sgc


	levels=[50,60,70,80,90]
	nlevel=n_elements(levels)
	c_color = (make_rainbow(nlevel+1))[1:nlevel]
	contour, $
		/irreg, /over, $
		dens.dstruct[w].density, sra[w], dens.dstruct[w].dec, $
		levels=levels, c_color=c_color


	plegend,method
	plegend,string(levels,f='(i0)'),line=0, color=c_color,/right

	endplot,/trim,/png





	return



end

pro esboss::plot_subreg_density, ra, dec, chunks, psfile=psfile, ind=ind, $
		maxres=maxres, $
		xsize=xsize, ysize=ysize, xrange=xrange, yrange=yrange, $
		offset=offset, $
		hammer=hammer, $
		density_struct=ds, $
		over=over, dmin=dmin, dmax=dmax

	; make 2-d a density map with sub regions as the areas

	if not keyword_set(over) then begin
		loadct, 13 ; rainbow
	endif

	if n_elements(chunks) eq 1 then begin

		if strn(chunks[0]) eq 'ngc' or strn(chunks[0]) eq 'ngcgood' then begin
			if keyword_set(hammer) then begin
				scale_x = 19000
				scale_y = 4500
				scale_ysize=9500

				plat=0d
				plon=180d
				;lim = [-10d,-90d,75d,90d]
				lim=[-10d,110d,90d,-80d]
				map_set, plat, plon, /hammer, /isotropic, lim=lim

				latlab=180d
				over=1

			endif else begin
				scale_x = 18000
				scale_y = 3300
				scale_ysize=12000
			endelse
			scale_step=5
			scale_charsize=1.0
		endif else if strn(chunks[0]) eq 'sgc' then begin
			if keyword_set(hammer) then begin
				scale_x = 19000
				scale_y = 5000
				scale_ysize=9500

				plat=0d
				plon=0d
				;lim = [-10d,-90d,75d,90d]
				;lim=[-10d,110d,90d,-80d]
				lim=[-15,315,38,60]
				map_set, plat, plon, /hammer, /isotropic, lim=lim
				over=1

			endif else begin
				offset=90
				scale_x = 6000
				scale_y = 3300

				scale_ysize=15000

			endelse
			scale_step=5
			scale_charsize=1.0

		endif else if chunks[0] eq 1 then begin
			scale_x = 34000
			scale_y = 3300
			scale_ysize=5300
			scale_step=10
			scale_charsize=0.65
			if n_elements(offset) eq 0 then offset=90
			if n_elements(maxres) eq 0 then maxres=4
		endif else if chunks[0] eq 2 then begin
			scale_x = 18000
			scale_y = 3300
			scale_ysize=11000
			scale_step=5
			scale_charsize=0.8
		endif else if chunks[0] eq 3 then begin
			scale_x = 18000
			scale_y = 3300
			scale_ysize=8000
			scale_step=5
			scale_charsize=0.8
		endif else if chunks[0] eq 4 then begin

			scale_x = 6000
			scale_y = 3300

			scale_ysize=15000

			scale_step=5
			scale_charsize=1.0
			;xrange=[110.0,220.0]
			;!x.style=3
			;if n_elements(maxres) eq 0 then maxres=4
		endif else if chunks[0] eq 5 then begin
			xrange=[180,250]
			xstyle=1
			iso=1
			scale_x = 30000
			scale_y = 3300

			scale_ysize=7500

			scale_step=5
			scale_charsize=1.0

		endif else begin
			scale_x = 6000
			scale_y = 3300

			scale_ysize=15000

			scale_step=5
			scale_charsize=1.0
		endelse

		if strn(chunks[0]) eq '1' then !p.charsize=1
	endif else begin
		;message,'multiple chunks not supported'
		if n_elements(dmin) eq 0 or n_elements(dmax) eq 0 then begin
			message,'you must send dmin= dmax= for multiple chunks'
		endif
		scale_x = 18000
		scale_y = 3300
		scale_ysize=11000
		scale_step=5
		scale_charsize=0.8
	endelse

	ds=self->bosstile_subreg_density(ra, dec, chunks, polygons=polygons, $
		maxres=maxres)



	; now scale lowest density to 0 and highest density to 255
	; this only works if there is one polygon per region!!

	if n_elements(polygons) ne n_elements(ds.dstruct) then begin
		;message,'cannot use this 1-1 coloration when regions contain '+$
		;	'multiple polygons'
	endif

	w=where(ds.dstruct.density gt 0 and ds.dstruct.area gt 0.2, nw)
	colors = replicate(0, n_elements(ds.dstruct))

	if n_elements(dmin) eq 0 then dmin_use = min(ds.dstruct[w].density) else dmin_use=dmin
	if n_elements(dmax) eq 0 then dmax_use = max(ds.dstruct[w].density) else dmax_use=dmax
	colors[w] = bytscl(ds.dstruct[w].density, min=dmin_use, max=dmax_use)

	!p.charsize=2

	; now associate these colors with the polygons: there may be more
	; than one polygon in a region
	for i=0L, nw-1 do begin
		reg = ds.dstruct[w[i]].region
		wmatch = where(polygons.region eq reg,nmatch)
		if nmatch ne 0 then begin
			color = colors[w[i]]
			add_arrval, replicate(color, nmatch), pcolors
			add_arrval, polygons[wmatch], ppolygons
		endif
	endfor

	;self->plot_poly, polygons[w], /fill, color=colors[w], $
	;	xrange=xrange, yrange=yrange, offset=offset, $
	;	over=over, iso=iso, xstyle=xstyle, ystyle=ystyle
	self->plot_poly, ppolygons, /fill, color=pcolors, $
		xrange=xrange, yrange=yrange, offset=offset, $
		over=over, iso=iso, xstyle=xstyle, ystyle=ystyle


	if n_elements(chunks) eq 1 then begin
		if strn(chunks) eq '2' then begin
			nover=30
			line_ceta = arrscl( findgen(nover), 30,42 )
			line_clambda = replicate(-41, nover)
			csurvey2eq, line_clambda, line_ceta, line_ra, line_dec
			pplot, line_ra,line_dec,/overplot, color=0
			print,line_ra,line_dec
		endif
	endif

	; over plot b=30
	if strn(chunks[0]) ne '1' and strn(chunks[0]) ne 'sgc' then begin
		if keyword_set(hammer) then bcolor=255
		glactc, ds.dstruct.ra, ds.dstruct.dec, 2000.0, l, b, 1, /degree
		maxl = max(l, min=minl)
		if minl gt 0 then minl = 0.9*minl else minl=1.1*minl
		if maxl gt 0 then maxl = 1.1*maxl else maxl=0.9*maxl

		bplot_vals = [20d,25d,30d]
		lines = [0,1,2]
		for j=0L,n_elements(bplot_vals)-1 do begin
			bval=bplot_vals[j]
			line=lines[j]

			nover = 100
			lvals = arrscl(findgen(nover), minl, maxl)
			bvals = replicate(bval, nover)
			glactc, ravals, decvals, 2000.0, lvals, bvals, 2, /degree
			pplot, ravals, decvals, /over, line=line, color=bcolor
		endfor
		mess=string(bplot_vals,f='(i0)')
		plegend,'b='+mess,/left,line=lines, charsize=1.0, color=bcolor
	endif
	
	if keyword_set(hammer) then begin
		lons=-180d + 15*findgen(25)
		; use ra [0,360] type names instead of the -180,180 longitude values
		names = strarr(25)
		w=where(lons ge 0)
		names[w] = string(lons[w],f='(i0)')
		w=where(lons lt 0)
		names[w] = string(lons[w] + 360,f='(i0)')
		;map_grid, latdel=15, londel=15,  /label, $
		;	charsize=1.0,lonalign=1.0, latalign=1.0, latlab=180d,lonlab=0d
		map_grid, /label, $
			lons=lons, lonnames=names, latdel=15, $
			charsize=1.0,lonalign=1.0, latalign=1.0, latlab=latlab,lonlab=0d
	endif


	if not keyword_set(over) or keyword_set(hammer) then begin
		put_clr_scl, scale_x, scale_y, [dmin_use,dmax_use], scale_step, $
			ysize=scale_ysize, charsize=scale_charsize
	endif

end


function esboss::get_ngc_large_thresholds

	; now will apply stricter cuts to get to 40/sq degree
	thresh_struct = { $
		method: '',   $
		tagname: '',  $
		flagname: '', $
		density: 0,   $ ;note integer
		thresh: 0.0   $
	}

	t = replicate(thresh_struct, 2*5)

	t[0].method = 'like'
	t[0].tagname = 'like_ratio'
	t[0].flagname = 'qso_like'
	t[0].density = 20
	t[0].thresh = 0.543

	t[1].method = 'mclike'
	t[1].tagname = 'like_ratio_pat'
	t[1].flagname = 'qso_like'
	t[1].density = 20
	t[1].thresh = 0.262

	t[2].method = 'kde'
	t[2].tagname = 'kde_prob'
	t[2].flagname = 'qso_kde'
	t[2].density = 20
	t[2].thresh = 0.903

	t[3].method = 'nn'
	t[3].tagname = 'nn_xnn'
	t[3].flagname = 'qso_nn'
	t[3].density = 20
	t[3].thresh = 0.852

	t[4].method = 'nn_value'
	t[4].tagname = 'nn_value'
	t[4].flagname = 'qso_bonus_main'
	t[4].density = 20
	t[4].thresh = 0.853



	t[0+5].method = 'like'
	t[0+5].tagname = 'like_ratio'
	t[0+5].flagname = 'qso_like'
	t[0+5].density = 40
	t[0+5].thresh = 0.234

	t[1+5].method = 'mclike'
	t[1+5].tagname = 'like_ratio_pat'
	t[1+5].flagname = 'qso_like'
	t[1+5].density = 40
	t[1+5].thresh = 0.108

	t[2+5].method = 'kde'
	t[2+5].tagname = 'kde_prob'
	t[2+5].flagname = 'qso_kde'
	t[2+5].density = 40
	t[2+5].thresh = 0.599

	t[3+5].method = 'nn'
	t[3+5].tagname = 'nn_xnn'
	t[3+5].flagname = 'qso_nn'
	t[3+5].density = 40
	t[3+5].thresh = 0.563

	t[4+5].method = 'nn_value'
	t[4+5].tagname = 'nn_value'
	t[4+5].flagname = 'qso_bonus_main'
	t[4+5].density = 40
	t[4+5].thresh = 0.573


	return, t


end

pro esboss::do_plot_subreg_density_2010_03_03d60, method, goal, str=str, $
		chunk=chunk, maxres=maxres

	common dpsdd60_block, oplot_chunk, isinwin


	if n_elements(chunk) ne 0 then plot_chunk=chunk else plot_chunk = 'ngc'
	target_type='qso'
	; pre-tuned each method to ~60/sq degree on chunk 5. Note
	; the result for kde ended up 50 on the test chunk
	target_run='2010-03-03d60'


	if n_elements(oplot_chunk) eq 0 then begin
		oplot_chunk = 'none'
		isinwin=0
	endif



	bt=obj_new('bosstarget')
	if n_elements(str) eq 0 then begin
		; this is the all-sky version
		str=bt->read_collated(target_type,target_run,extra_name='all')
	endif
	ntot = n_elements(str)


	if plot_chunk[0] ne oplot_chunk[0] then begin
		ply=self->bosstile_poly_read(plot_chunk)
		print,'finding inwindow'
		isinwin=is_in_window(ply, ra=str.ra, dec=str.dec )
		destruct_polygon, ply
		oplot_chunk=plot_chunk
	endif




	print,'trimming to thresh and chunk'
	thresh_struct = self->get_ngc_large_thresholds()

	wtype=where(thresh_struct.method eq method $
				and thresh_struct.density eq goal, ntype)
	if ntype ne 1 then begin
		message,string('Bad method/goal: ',method,goal,f='(a,a,"/",f0.3)')
	endif


	if not tag_exist(str, thresh_struct[wtype].tagname, index=tagind) then begin
		message,'bad tag name: ',thresh_struct[wtype].tagname
	endif

	wuse = where( str.(tagind) gt thresh_struct[wtype].thresh $
		          and isinwin, nuse)
	if nuse eq 0 then begin
		message,'No objects passed cuts on thresh/intest_region'
	endif

	use_str = str[wuse]
	print,nuse,ntot,f='(i0,"/",i0," passed thresh/intest_region")'



	extra_name = string(goal,f='("alld",i0)')
	if method eq 'mclike' then begin
		extra_name += '-pat'
	endif
	print,'extra_name: ',extra_name


	if goal eq 40 then begin
		dmin=30
		dmax=70
	endif else if goal eq 20 then begin
		dmin = 10
		dmax = 50
	endif
	self->plot_subreg_density_byrun, $
		target_type, target_run, plot_chunk, $
		target_flags = thresh_struct[wtype].flagname, $
		str=use_str, $
		dmin=dmin, dmax=dmax, $
		extra_name=extra_name,/hammer, $
		/testregoplot, $
		maxres=maxres

end

pro esboss::plot_subreg_density_byrun, target_type, target_run, chunks, $
		target_flags=target_flags, $
		where_string=where_string, $ ; in addition to target_flags
		goal=goal, $
		dmin=dmin, dmax=dmax, $
		extra_name=extra_name, $
		xsize=xsize, ysize=ysize, xrange=xrange, yrange=yrange, $
		maxres=maxres, $
		hammer=hammer, $
		offset=offset, $
		$
		testregoplot=testregoplot, $
		psfile=psfile, $
		str=str

	if n_elements(target_flags) eq 0 then begin
		case target_type of 
			'qso': begin
				target_flags=$
					['qso_bonus_main','qso_first_boss','qso_known_midz']
			end
			'lrg': begin
				target_flags = $
					['gal_loz','gal_cmass','gal_cmass_sparse']
			end
			else: message,'qso or lrg only'
		endcase
	endif


	if size(chunks[0],/tname) eq 'STRING' then begin
		chstr=strjoin(chunks,'-')
	endif else begin
		chstr=strjoin(string(chunks,f='(i0)'), '-')
	endelse
	chstr='chunk'+chstr
	fmid = chstr
	if n_elements(extra_name) ne 0 then fmid=[extra_name,fmid]
	if n_elements(target_flags) ne 0 then begin
		tf=strjoin(target_flags, '-')
		tf=repstr(tf,'_','-')
		fmid=[fmid,tf]
	endif
	tr = strjoin(target_run, '-')
	if n_elements(psfile) eq 0 then begin
		psdir=filepath( $
			root=getenv('BOSS_TARGET'), subdir=[target_run[0]],'plots')
		psfile=[tr,target_type,fmid,'densmap']
		if n_elements(maxres) ne 0 then begin
			psfile=[psfile,'maxres'+ntostr(long(maxres))]
		endif
		if keyword_set(hammer) then begin
			psfile=[psfile,'hammer']
		endif
		psfile=strjoin(psfile,'-')+'.eps'
		psfile=path_join(psdir, psfile)
	endif else begin
		psdir = file_dirname(psfile)
	endelse

	if not file_test(psdir) then file_mkdir, psdir


	nchunk = n_elements(chunks)
	;if nchunk gt 1 then begin
	;	if n_elements(dmin) eq 0 or n_elements(dmax) eq 0 then begin
	;		message,'you must send dmin= dmax= for multiple chunks'
	;	endif
	;endif

	if n_elements(chunks) eq 1 then begin
		case chunks[0] of
			1: begin
				if n_elements(ysize) eq 0 then ysize=4.0
				if n_elements(xsize) eq 0 then xsize=15.0
			end
			4: begin
				if n_elements(xsize) eq 0 then xsize=15.0
			end
			5: begin
				xsize=15
				ysize=5
			end
			'possible5': begin
				xsize = 15.0
			end
			'ngc': begin
				densmap_dpi=130
			end
			'ngcgood': begin
				densmap_dpi=130
			end
			'sgc': begin
				densmap_dpi=130
			end
			else:  begin
				densmap_dpi=85
			end
		endcase
	endif

	if n_elements(xsize) eq 0 then xsize=8.5
	if n_elements(ysize) eq 0 then ysize=8.0

	bt=obj_new('bosstarget')
	
	begplot,psfile,/color,/encap,xsize=xsize,ysize=ysize
	for i=0L, nchunk-1 do begin
		if i gt 0 then begin
			over=1
		endif
		chunk=chunks[i]	

		if n_elements(target_run) eq nchunk then begin
			run=target_run[i]
		endif else begin
			run=target_run
		endelse
		if n_elements(str) eq 0 or nchunk gt 1 then begin
			str=bt->read_collated(target_type,run,extra_name=extra_name)
		endif
		w=self->btselect(str, target_flags)

		if n_elements(where_string) ne 0 then begin
			ww=struct_select(str[w], where_string, nww,/verbose)
			
			w=w[ww]
		endif

		self->plot_subreg_density, str[w].ra, str[w].dec, chunk, $
			xsize=xsize, ysize=ysize, xrange=xrange, yrange=yrange, $
			dmin=dmin, dmax=dmax, $
			offset=offset, over=over, $
			hammer=hammer, $
			maxres=maxres, $
			density_struct=ds
	endfor


	if keyword_set(testregoplot) then begin
		testply = self->bosstile_poly_read('ngc-large')
		plot_poly, testply, color=c2i('white'), /over,outline_thick=2
		destruct_polygon, testply
	endif

	endplot,/png,dpi=densmap_dpi,/trim

	if n_elements(target_run) gt 1 then begin
		print,'not doing density plots for multiple runs'
		return
	endif


	; versus galactic latitude b
	psfile=[tr,target_type,fmid,'densb']
	psfile=strjoin(psfile,'-')+'.eps'
	psfile=path_join(psdir, psfile)

	begplot,psfile,/color,/encap,xsize=8.5,ysize=8.0
	!p.charsize=2
	b=ds.dstruct.b
	dens=ds.dstruct.density
	w=where(b gt -9000, nw)
	b=b[w]
	dens=dens[w]

	s=sort(b)
	pplot, b[s], dens[s], psym=-8, xtitle='b', ytitle=textoidl('#/deg^2')
	if n_elements(goal) ne 0 then begin
		pplot, !x.crange, [goal,goal], line=2, color='red', $
			/over
	endif
	endplot,/png,dpi=55



	; versus RA
	psfile=[tr,target_type,fmid,'densra']
	psfile=strjoin(psfile,'-')+'.eps'
	psfile=path_join(psdir, psfile)

	begplot,psfile,/color,/encap,xsize=8.5,ysize=8.0
	!p.charsize=2
	ra=ds.dstruct.ra
	dens=ds.dstruct.density
	w=where(ra gt -9000, nw)
	ra=ra[w]
	dens=dens[w]

	if chunks[0] eq 1 then begin
		ra=shiftra(ra,/wrap)
	endif

	s=sort(ra)
	pplot, ra[s], dens[s], psym=-8, xtitle='RA', ytitle=textoidl('#/deg^2')
	if n_elements(goal) ne 0 then begin
		pplot, !x.crange, [goal,goal], line=2, color='red', $
			/over
	endif

	endplot,/png, dpi=55



end

pro esboss::do_heat_map, types

	nt=n_elements(types)
	for i=0L, nt-1 do begin

		delvarx,xrange,yrange,extra_name
		type=strlowcase(types[i])

		if strmatch(type, '*qso*') then begin
			target_type = 'qso'
			goal=60
			dmin=45
			dmax=150
		endif else begin
			target_type = 'lrg'
			goal=130
			dmin=75
			dmax=200
			extra_name='noknown'
		endelse

		case type of
			'lrg1': begin
				target_runs='2009-12-10-newlrg1'
				chunks=2
			end
			'lrg2': begin
				target_runs='2009-12-10-newlrg2'
				chunks=2
			end
			'lrg3': begin
				target_runs='main001'
				chunks=3
			end
			'lrg4': begin
				target_runs='main001'
				chunks=4
			end

			'qso1': begin
				target_runs='2009-12-10-newlike1'
				chunks=2
			end
			'qso2': begin
				target_runs='2009-12-10-newlike2'
				chunks=2
			end
			'qso3': begin
				target_runs='main001'
				chunks=3
			end
			'qso4': begin
				target_runs='main001'
				chunks=4
			end

			'lrg2-3': begin
				target_runs = $
					['2009-12-10-newlrg2',$
					 'main001']
				chunks=[2,3]
				xrange=[100,140]
				yrange=[20,60]
			end
			'qso2-3': begin
				target_runs = $
					['2009-12-10-newlike2',$
					 'main001']
				chunks=[2,3]
				xrange=[100,140]
				yrange=[20,60]
			end


			'lrg2-3-4': begin
				target_runs = $
					['2009-12-10-newlrg2',$
					 'main001',$
					 'main001']
				chunks=[2,3,4]
				xrange=[100,220]
				yrange=[-10,60]
			end
			'qso2-3-4': begin
				target_runs = $
					['2009-12-10-newlike2',$
					 'main001',$
					 'main001']
				chunks=[2,3,4]
				xrange=[100,220]
				yrange=[-10,60]
			end
		endcase

		self->plot_subreg_density_byrun,target_type, target_runs, chunks, $
			goal=goal, dmin=dmin, dmax=dmax, xrange=xrange, yrange=yrange, $
			extra_name=extra_name
	endfor
end

function esboss::sloane_tiles
	ntile = 9722
	tile_radius = 1.49
	tile_area = !dpi*tile_radius^2
	sloane, ntile, theta=theta, phi=phi

	theta=reform(theta)
	phi=reform(phi)

	ra = phi
	dec = 90d - theta

	eq2csurvey, ra, dec, clambda, ceta
	glactc, ra, dec, 2000.0, l, b, 1, /degree

	st={ntile:ntile, $
		tile_radius: tile_radius, $
		tile_area: tile_area, $
		tileid: lindgen(ntile), $
		theta:theta, phi:phi, $
		ra:ra, dec:dec, $
		clambda:clambda, ceta:ceta, $
		l:l, b:b}

	return, st
end
pro esboss::sloane_match, struct, mstruct, msloane, distances
	if n_elements(struct) eq 0 then begin
		message,'Usage: eb->sloane_match, struct [, mstruct, msloane, dist]'
	endif
	sloane_struct=self->sloane_tiles()

	spherematch, struct.ra, struct.dec, sloane_struct.ra, sloane_struct.dec,  $
		sloane_struct.tile_radius, $
		mstruct, msloane, distances, maxmatch=0
end




pro esboss::run_plothist_tile_density, target_type, target_run, qso=qso, lrg=lrg, anyresolve=anyresolve, fpobjc=fpobjc, noradcheck=noradcheck

	if target_type eq 'qso' then begin
		self->plothist_tile_density, 'qso', target_run, struct=qso, $
			['qso_core','qso_bonus'], fpobjc=fpobjc, noradcheck=noradcheck
		self->plothist_tile_density, 'qso', target_run, struct=qso, $
			'qso_core', fpobjc=fpobjc, noradcheck=noradcheck
		self->plothist_tile_density, 'qso', target_run, struct=qso, $
			'qso_bonus', fpobjc=fpobjc, noradcheck=noradcheck

		self->plothist_tile_density, 'qso', target_run, struct=qso, $
			['qso_core','qso_bonus','qso_known_midz'], fpobjc=fpobjc, $
			noradcheck=noradcheck
		self->plothist_tile_density, 'qso', target_run, struct=qso, $
			['qso_core','qso_known_midz'], fpobjc=fpobjc, $
			noradcheck=noradcheck
		self->plothist_tile_density, 'qso', target_run, struct=qso, $
			['qso_bonus','qso_known_midz'], fpobjc=fpobjc, $
			noradcheck=noradcheck

	endif

	if target_type eq 'lrg' then begin
		self->plothist_tile_density, 'lrg', target_run, struct=lrg, $
			['gal_loz','gal_cmass'], fpobjc=fpobjc,noradcheck=noradcheck
		self->plothist_tile_density, 'lrg', target_run, struct=lrg, $
			'gal_loz', fpobjc=fpobjc, noradcheck=noradcheck
		self->plothist_tile_density, 'lrg', target_run, struct=lrg, $
			'gal_cmass', fpobjc=fpobjc, noradcheck=noradcheck
	endif

end


; this the chunks, or a requested chunk
pro esboss::plothist_tile_density, $
		chunks, target_type, target_run, $
		btflags=btflags, $
		where_string=where_string, $
		extra_name=extra_name, $
		struct=struct, reselect=reselect, outdir=outdir

	if n_params() lt 3 then begin
		on_error,2
		print,'Usage: '
		print,'  eb->plothist_tile_density, chunks, target_type, target_run, '
		print,'     btflags=, '
		print,'     where_string=, '
		print,'     extra_name=, '
		print,'     /reselect, '
		print,'     outdir=, '
		print
		message,'Halting'
	endif


	if n_elements(chunks) ne 1 then begin
		message,'chunks= is now required and must be length 1'
	endif
	bt=obj_new('bosstarget')
	if n_elements(struct) eq 0 or keyword_set(reselect) then begin
		struct=bt->read_collated(target_type, target_run, fpobjc=fpobjc,$
			extra_name=extra_name)
	endif
	nstruct=n_elements(struct)


	if keyword_set(reselect) then begin
		if chunks eq 2 then comm2=1
		if chunks eq 1 then commissioning=1
		btypeobj=obj_new('bosstarget_'+target_type)
		struct.boss_target1 = btypeobj->select(struct, $
			comm2=comm2, commissioning=commissioning, /reselect, pars=pars)
	endif


	if n_elements(btflags) eq 0 then begin
		case target_type of
			'qso': begin
				btflags= [$
					'qso_core','qso_bonus','qso_known_midz','qso_known_lohiz',$
					'qso_nn','qso_ukidss','qso_kde_coadd','qso_like',$
					'qso_first_boss']
			end
			'std': begin
				btflags = ['std_fstar']
			end
			'lrg': begin
				btflags = ['gal_loz', 'gal_cmass','gal_cmass_sparse']
			end
			else: message,"target_type must be 'qso','lrg','std'"
		endcase
		btflags_string='any'
		btflags_string_leg='any'
	endif else begin
		btflags_string = strjoin(btflags, '-')
		btflags_string_leg = strjoin(btflags, ' or ')
	endelse


	keep = self->btselect(struct, btflags)

	extra=repstr(btflags_string,'_','-')
	mess1 = btflags_string_leg+' targets'

	help,struct,keep

	if n_elements(where_string) ne 0 then begin
		mess1=where_string
		comm='keep2 = where('+where_string+', nkeep2)'
		if not execute(comm) then message,'command: '+comm+' failed'
		keep=keep[keep2]
	endif
	help,struct,keep



	; match up to the tiles and get density
	tstruct = struct[keep]
	density_struct=self->bosstile_density(tstruct.ra, tstruct.dec, chunks,$
		tiles=tiles, matched_objects=mstruct)


	if n_elements(chunks) eq 0 then begin
		chunks=tiles[rem_dup(tiles.chunk)].chunk
	endif
	chunkstr=strjoin(string(chunks,f='(i03)'), '-')
	if n_elements(chunks) eq 1 then chunkleg='chunk: ' else chunkleg='chunks: '
	chunkleg=chunkleg+strjoin(string(chunks,f='(i03)'), ',')
	extra=[extra,chunkstr]

	if n_elements(outdir) eq 0 then begin
		outdir = bt->target_dir(target_run)
		outdir = path_join(outdir, 'plots')
	endif
	if not file_test(outdir) then file_mkdir, outdir

	if n_elements(extra_name) ne 0 then extra=[extra,extra_name]

	if keyword_set(reselect) then extra=[extra,'reselect']
	extra = strjoin(extra,'-')






	; positions in ra/dec
	psfile = path_join(outdir, target_run+'-'+target_type+'-tile-positions-eq-'+extra+'.eps')
	begplot, psfile, /encap, /color, xsize=11, ysize=8.5

	;yrange=[120,160]
	;xrange=[-50, 50]

	if target_type eq 'std' then begin
		all_psym=8
		intile_psym=8
		all_symsize=0.25
		intile_symsize=0.75 
	endif else begin
		all_psym=3
		intile_psym=8
		intile_symsize=0.1
	endelse

	plot, struct[keep].ra, struct[keep].dec, $
		psym=all_psym, symsize=all_symsize, $$
		xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, $
		xtitle=textoidl('RA'), ytitle=textoidl('DEC')

	pplot, /over, $
		struct[keep[mstruct]].ra, struct[keep[mstruct]].dec, $
		psym=intile_psym, symsize=intile_symsize, color='blue'


	;self->bosstile_polygons_read, chunks, poly, polyids
	tilepoly=self->bosstile_poly_read(chunks)
	plot_poly, tilepoly, /over,outline_thick=2, color=!p.color

	if 0 and n_elements(density_struct.ra) lt 50 then begin
		for ii=0L, n_elements(density_struct.ra)-1 do begin
			xyouts, density_struct[ii].ra, density_struct[ii].dec, string(density_struct[ii].density, f='(f0.1)')
		endfor
	endif

	legend, [mess1, 'within tile'],$
		psym=8,color=[!p.color,c2i('blue')], /left
	legend, ['',chunkleg], /right

	pplot, /over, tiles.ra, tiles.dec, psym=7, color='red'
	wlow=where(density_struct.area lt self->bosstile_radius()^2*!pi*0.99, nlow)
	if nlow ne 0 then begin
		pplot, /over, density_struct[wlow].ra, tiles[wlow].dec, $
			psym=7, color='orange'
	endif
	legend, $
		['filled tiles','partially filled'], psym=7, color=['red','orange'],$
		/right,/bottom
	n=100
	eq_ra = arrscl( findgen(n), 0, 360)
	eq_dec = replicate(0, n)
	pplot, eq_ra, eq_dec,/overplot, thick=1
	endplot


	if 0 then begin

		; positions in clambda-ceta
		psfile = path_join(outdir, target_run+'-'+target_type+'-tile-positions-'+extra+'.eps')
		begplot, psfile, /encap, /color, xsize=11, ysize=8.5

		eq2csurvey, struct.ra, struct.dec, clambda, ceta
		;yrange=[120,160]
		;xrange=[-50, 50]
		plot, clambda[keep], ceta[keep], psym=3, $
			xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, $
			xtitle=textoidl('\lambda_c'), ytitle=textoidl('\eta_c')

		color='red'
		tvcircle, $
			self->bosstile_radius(), $
			tiles.clambda, tiles.ceta, $
			color=c2i(color), /data, thick=3

		if n_elements(density_struct.clambda) lt 50 then begin
			for ii=0L, n_elements(density_struct.clambda)-1 do begin
				;xyouts, density_struct[ii].clambda, density_struct[ii].ceta, strn(density_struct[ii].count)
				xyouts, density_struct[ii].clambda, density_struct[ii].ceta, string(density_struct[ii].density, f='(f0.1)')
			endfor
		endif
		if n_elements(mess1) ne 0 then legend, mess1, /left
		legend, ['',chunkleg], /right

		n=100
		eq_ra = arrscl( findgen(n), 0, 360)
		eq_dec = replicate(0, n)
		eq2csurvey, eq_ra, eq_dec, eq_clambda, eq_ceta
		pplot, eq_clambda, eq_ceta,/overplot, thick=1
		endplot



		; positions in l-b
		psfile = path_join(outdir, target_run+'-'+target_type+'-tile-positions-lb-'+extra+'.eps')
		begplot, psfile, /encap, /color, xsize=11, ysize=8.5


		;brange=[-75, -20]
		;lrange=[30,210]
		glactc, struct.ra, struct.dec, 2000.0, l, b, 1, /degree
		plot, l[keep], b[keep], psym=3, $
			yrange=brange, ystyle=3, xrange=lrange, xstyle=3, $
			xtitle='l', ytitle='b'

		color='red'
		tvcircle, $
			self->bosstile_radius(), $
			tiles.l, tiles.b, $
			color=c2i(color), /data, thick=3

		if n_elements(mess1) ne 0 then legend, mess1, /left
		legend, ['',chunkleg], /right

		n=100
		eq_ra = arrscl( findgen(n), 0, 360)
		eq_dec = replicate(0, n)
		glactc, eq_ra, eq_dec, 2000.0, eq_l, eq_b, 1, /degree
		pplot, eq_l, eq_b,/overplot, thick=1
		endplot

	endif



	psfile = path_join(outdir, target_run+'-'+target_type+'-tile-density-hist-'+extra+'.eps')
	begplot, psfile, /encap, xsize=11, ysize=8.5

	case target_type of
		'qso': begin
			binsize=5
		end
		'lrg': begin
			binsize=7
		end
		'std': begin
			binsize=1
		end
		else: message,"target_type must be 'qso','lrg','std'"
	endcase

	plothist, density_struct.density, $
		bin=binsize, xhist, yhist, /noplot;min=densmin, /noplot

	plothist, density_struct.density, $
		bin=binsize, xtitle='density #/square degree', $;/xlog, $
		/ylog, yrange=[0.9, max(yhist)*1.5], ystyle=3;, $
		;xrange=[densmin, max(xhist)*1.5], xstyle=3;, min=densmin

	med_density = median(density_struct.density)
	mean_density = mean(density_struct.density)
	sdev_density = stdev(density_struct.density)
	mess = $
		['median density: '+strn(round(med_density)), $
		'mean density: '+strn(round(mean_density)), $
		'sdev density: '+strn(round(sdev_density))]

	if n_elements(mess1) ne 0 then begin
		mess=[btflags_string_leg+' tiled targets: '+strn(n_elements(mstruct)),$
			mess]
	endif
	legend,mess, /right
	legend, chunkleg, /left
	endplot



	; plot density versus galactic latitude

	psfile = path_join(outdir, target_run+'-'+target_type+'-tile-density-vs-b-'+extra+'.eps')

	begplot, psfile, /encap, xsize=11, ysize=8.5, /color
	plot, density_struct.b, density_struct.density, psym=8, $
		xtitle='galactic latitude', ytitle='#/square degree',/ynozero

	if n_elements(mess1) ne 0 then legend, mess1, /left
	legend, chunkleg, /right
	endplot

	; vs ra
	psfile = path_join(outdir, target_run+'-'+target_type+'-tile-density-vs-ra-'+extra+'.eps')

	begplot, psfile, /encap, xsize=11, ysize=8.5, /color
	plot, shiftra(density_struct.ra,/wrap), density_struct.density, psym=8, $
		xtitle='RA', ytitle='#/square degree', /ynozero

	if n_elements(mess1) ne 0 then legend, mess1, /left
	legend, ['',chunkleg], /right
	endplot




	fitsfile = $
		path_join(outdir, target_run+'-'+target_type+'-tile-density-'+extra+'.fits')
	print,'Writing to file: ',fitsfile
	mwrfits, density_struct, fitsfile, /create


end



; this does the whole sgc
pro esboss::plothist_sgc_tile_density, target_type, target_run, btflags, $
		struct=struct, fpobjc=fpobjc, noradcheck=noradcheck, $
		psdir=psdir

	bt=obj_new('bosstarget')
	if n_elements(struct) eq 0 then begin
		struct=bt->read_collated(target_type, target_run, fpobjc=fpobjc)
	endif

	btflags_string = strjoin(btflags, '-')
	extra=repstr(btflags_string,'_','-')
	mess1 = strjoin(btflags, ' or ')+' targets'

	if keyword_set(noradcheck) then begin
		extra=extra+'-noradcheck'
		mess1 += '-noradcheck'
	endif


	keep = self->btselect(struct, btflags)
	help,struct,keep

	; order important since the tiles go over the pole
	tstruct = struct[keep]

	isplit = 21.3
	if target_type eq 'qso' then begin
		bq=obj_new('bosstarget_qso')
		tlups=bq->get_lups(tstruct)
		lowi=where(tlups[3,*] lt isplit, comp=highi)
		density_struct_lowi=self->sloane_tile_density($
			tstruct[lowi], sloane=sloane_struct)
		density_struct_highi=self->sloane_tile_density($
			tstruct[highi], sloane=sloane_struct)

		wlowi = where(density_struct_lowi.keep)
		whighi = where(density_struct_highi.keep)
	endif




	density_struct=self->sloane_tile_density($
		tstruct, sloane=sloane_struct)

	w=where(density_struct.keep)

	if n_elements(psdir) eq 0 then begin
		psdir = bt->target_dir(target_run)
		psdir = path_join(psdir, 'plots')
	endif
	if not file_test(psdir) then file_mkdir, psdir




	; positions in clambda-ceta
	psfile = path_join(psdir, target_run+'-tile-positions-'+extra+'.eps')
	begplot, psfile, /encap, /color, xsize=11, ysize=8.5

	eq2csurvey, struct.ra, struct.dec, clambda, ceta
	plot, clambda[keep], ceta[keep], psym=3, $
		xrange=[-50, 50], yrange=[120, 160], xstyle=3, ystyle=3, $
		xtitle=textoidl('\lambda_c'), ytitle=textoidl('\eta_c')

		;color='yellow'
	color='red'
	tvcircle, $
		sloane_struct.tile_radius, $
		sloane_struct.clambda[w], sloane_struct.ceta[w], $
		color=c2i(color), /data, thick=3

	legend, mess1, /left

	n=100
	eq_ra = arrscl( findgen(n), 0, 360)
	eq_dec = replicate(0, n)
	eq2csurvey, eq_ra, eq_dec, eq_clambda, eq_ceta
	pplot, eq_clambda, eq_ceta,/overplot, thick=1
	endplot



	; positions in l-b
	psfile = path_join(psdir, target_run+'-tile-positions-lb-'+extra+'.eps')
	begplot, psfile, /encap, /color, xsize=11, ysize=8.5


	glactc, struct.ra, struct.dec, 2000.0, l, b, 1, /degree
	plot, l[keep], b[keep], psym=3, $
		yrange=[-75, -20], ystyle=3, xrange=[30,210], xstyle=3, $
		xtitle='l', ytitle='b'

	color='red'
	tvcircle, $
		sloane_struct.tile_radius, $
		sloane_struct.l[w], sloane_struct.b[w], $
		color=c2i(color), /data, thick=3

	legend, mess1, /left

	n=100
	eq_ra = arrscl( findgen(n), 0, 360)
	eq_dec = replicate(0, n)
	glactc, eq_ra, eq_dec, 2000.0, eq_l, eq_b, 1, /degree
	pplot, eq_l, eq_b,/overplot, thick=1
	endplot






	psfile = path_join(psdir, target_run+'-tile-density-hist-'+extra+'.eps')
	begplot, psfile, /encap, xsize=11, ysize=8.5

	if target_type eq 'qso' then begin
		if in(btflags,'qso_bonus') then begin
			binsize=7
		endif else begin
			binsize=1
		endelse
	endif else begin
		binsize=7
	endelse

	plothist, density_struct.density[w], $
		bin=binsize, xhist, yhist, /noplot;min=densmin, /noplot

	plothist, density_struct.density[w], $
		bin=binsize, xtitle='density #/square degree', $;/xlog, $
		/ylog, yrange=[0.9, max(yhist)*1.5], ystyle=3;, $
		;xrange=[densmin, max(xhist)*1.5], xstyle=3;, min=densmin

	med_density = median(density_struct.density[w])
	mean_density = mean(density_struct.density[w])
	sdev_density = stdev(density_struct.density[w])
	legend,$
		[mess1, $
		'median density: '+strn(round(med_density)), $
		'mean density: '+strn(round(mean_density)), $
		'sdev density: '+strn(round(sdev_density))], $
		/right
	endplot



	; plot density versus galactic latitude
	; also split by g=21.3 

	psfile = path_join(psdir, target_run+'-tile-density-vs-b-'+extra+'.eps')

	begplot, psfile, /encap, xsize=11, ysize=8.5, /color
	plot, density_struct.b[w], density_struct.density[w], psym=8, $
		xtitle='galactic latitude', ytitle='#/square degree', $
		xrange=[-80,-20], xstyle=3
		
	err=replicate(1.0, n_elements(w))
	fitlin, density_struct.b[w], density_struct.density[w], err, $
		intercept, sig_int, slope, sig_slope

	allb = density_struct.b[w]
	allb = allb[sort(allb)]
	oplot, allb, intercept + slope*allb

	afit = string(intercept, format='(%"a=%0.2f")')
	bfit = string(slope,format='(%"b=%0.2f")')

	legend, mess1, /left
	legend, [afit,bfit],/right
	endplot


	if target_type eq 'qso' then begin
		; low i < 21.3
		psfile = path_join(psdir, $
			target_run+'-tile-density-vs-b-'+extra+'-lowi.eps')

		begplot, psfile, /encap, xsize=11, ysize=8.5, /color
		plot, density_struct_lowi.b[wlowi], $
			density_struct_lowi.density[wlowi], $
			psym=8, $
			xtitle='galactic latitude', ytitle='#/square degree', $
			xrange=[-80,-20], xstyle=3
			
		legend, [mess1,'i < '+string(isplit,f='(f0.1)')], /left
		endplot

		; high i > 21.3
		psfile = path_join(psdir, $
			target_run+'-tile-density-vs-b-'+extra+'-highi.eps')

		begplot, psfile, /encap, xsize=11, ysize=8.5, /color
		plot, $
			density_struct_highi.b[whighi], $
			density_struct_highi.density[whighi], $
			psym=8, $
			xtitle='galactic latitude', ytitle='#/square degree', $
			xrange=[-80,-20], xstyle=3
			
		legend, [mess1,'i > '+string(isplit,f='(f0.1)')], /left
		endplot
	endif










	if n_elements(dperp0) ne 0 then begin
		add_tags, density_struct, 'dperp0', '0d', tmp
		tmp.dperp0 = dperp0
		density_struct = tmp
	endif
	if n_elements(ihi) ne 0 then begin
		add_tags, density_struct, 'ihi', '0d', tmp
		tmp.ihi = ihi
		density_struct = tmp
	endif

	fitsfile = $
		path_join(psdir, target_run+'-tile-density-'+extra+'.fits')
	print,'Writing to file: ',fitsfile
	mwrfits, density_struct, fitsfile, /create


end

function esboss::sgc_lrg_tile_density_pars
	nihi=21
	nrlo=16
	ntri=21
	ndperp0=10
	return, {$
		nihi:nihi, $
		ihi:arrscl(findgen(nihi), 19.9, 20.0), $
		$
		nrlo:nrlo, $
		rlo:13.25+0.01*findgen(nrlo), $
		$
		ntri:ntri, $
		trifrac:arrscl(findgen(ntri), 0.0 ,1.0), $
		$
		ndperp0:ndperp0, $
		dperp0:arrscl(findgen(ndperp0), 0.45, 0.55) }

end

function esboss::sgc_qso_tile_density_pars, chi2_pars=chi2_pars, old=old

	x2star_permissive = 3.0
	x2star_restrictive = 6.0
	chi2_pars = { $
		x2star_permissive: x2star_permissive, $
		x2star_restrictive: x2star_restrictive}

	nbonus_bright = 20	
	nbonus_faint = 20	

	ncore_bright = 20	
	ncore_faint = 20	


	if not keyword_set(old) then begin

		; defaults in main code are -=0.70, -0.80
		bright_restrictive_def = -1.0
		faint_restrictive_def  = -1.0

		;bright_permissive_def = 0.20
		;faint_permissive_def = -0.60
		bright_permissive_def = -1.0
		faint_permissive_def = -1.0

		;core_maxoffset=2.0
		core_maxoffset=0.5
		bonus_maxoffset=2.0
		; this returns combined pars
		return, {$
			$ ; core sample
			ncore_bright: ncore_bright, $
			ncore_faint: ncore_faint, $
			$
			logsdmax_bright_restr_default:bright_restrictive_def, $
			logsdmax_faint_restr_default: faint_restrictive_def, $
			$
			logsdmax_bright_restr_offsets: $
				arrscl(findgen(ncore_bright), 0.0, core_maxoffset), $
			logsdmax_faint_restr_offsets: $
				arrscl(findgen(ncore_faint), 0.0, core_maxoffset), $
			$
			$ ; bonus sample
			nbonus_bright: nbonus_bright, $
			nbonus_faint: nbonus_faint, $
			$
			logsdmax_bright_perm_default:bright_permissive_def, $
			logsdmax_faint_perm_default: faint_permissive_def, $
			$
			logsdmax_bright_perm_offsets: $
				arrscl(findgen(nbonus_bright), 0.0, bonus_maxoffset), $
			logsdmax_faint_perm_offsets: $
				arrscl(findgen(nbonus_faint), 0.0, bonus_maxoffset) $
		}

	endif else begin
		logstardensmax_bright_permissive = 0.20
		logstardensmax_faint_permissive = -0.30

		logstardensmax_bright_restrictive = -0.60
		logstardensmax_faint_restrictive  = -1.10

		maxoffset=1.45
		; this returns combined pars
		return, {$
			nbonus_bright: nbonus_bright, $
			nbonus_faint: nbonus_faint, $
			logsdmax_bright_perm_default: 0.20, $
			logsdmax_faint_perm_default: -0.30, $
			logsdmax_bright_perm_offsets: $
				arrscl(findgen(nbonus_bright), 0.0, maxoffset), $
			logsdmax_faint_perm_offsets: $
				arrscl(findgen(nbonus_faint), 0.0, maxoffset), $
			$
			ncore_bright: ncore_bright, $
			ncore_faint: ncore_faint, $
			logsdmax_bright_restr_default: -0.60, $
			logsdmax_faint_restr_default: -1.10, $
			logsdmax_bright_restr_new: -0.55, $
			logsdmax_faint_restr_new: -0.30, $
			logsdmax_bright_restr_offsets: $
				arrscl(findgen(ncore_bright), 0.0, maxoffset), $
			logsdmax_faint_restr_offsets: $
				arrscl(findgen(ncore_faint), 0.0, maxoffset), $
			x2star_permissive: x2star_permissive, $
			x2star_restrictive: x2star_restrictive}
	endelse

		 
end

function esboss::qso_tile_density_file, target_run, qso_types, pars, extra_name


	if n_params() lt 4 then begin
		on_error, 2
		message,'usage: f=eb->qso_tile_density-file(target_run, qso_type, pars, extra_name)'
	endif

	bt=obj_new('bosstarget')
	tmp_name=bt->target_file('qso', target_run, /collate, extra_name=extra_name)

	tmp_name=file_basename(tmp_name)

	outdir = bt->target_dir(target_run)
	outdir = path_join(outdir, 'density-'+extra_name)

	file_mkdir, outdir

	qso_types_str = strjoin(qso_types, '-')
	extra=repstr(qso_types_str,'_','-')

	if n_elements(pars) ne 0 then begin
		if in(qso_types, 'qso_core') and in(qso_types,'qso_bonus') then begin
			extra += $
				'-sdbp-'+string(pars.logstardensmax_bright_permissive,$
				f='(f0.3)')
			extra += $
				'-sdfp-'+string(pars.logstardensmax_faint_permissive,$
				f='(f0.3)')
		endif else if qso_types eq 'qso_core' or qso_types eq 'qso_kde_coadd' then begin

			if tag_exist(pars,'logstardensmax_bright_restrictive') then begin
				extra += $
					'-sdbr-'+string(pars.logstardensmax_bright_restrictive,$
					f='(f0.3)')
			endif
			if tag_exist(pars,'logstardensmax_faint_restrictive') then begin
				extra += $
					'-sdfr-'+string(pars.logstardensmax_faint_restrictive,$
					f='(f0.3)')
			endif

			if tag_exist(pars,'x2star_restrictive_ra_hidens') then begin
				extra += $
					'-x2star-core-hidens'+$
					string(pars.x2star_restrictive_ra_hidens,$
					       f='(f0.3)')
			endif

			if tag_exist(pars,'x2star_restrictive_region1') then begin
				extra += $
					'-x2star-core-region1'+$
					string(pars.x2star_restrictive_region1,$
					       f='(f0.3)')
			endif
			if tag_exist(pars,'x2star_restrictive_region2') then begin
				extra += $
					'-x2star-core-region2'+$
					string(pars.x2star_restrictive_region2,$
					       f='(f0.3)')
			endif
			if tag_exist(pars,'x2star_restrictive_region3') then begin
				extra += $
					'-x2star-core-region3'+$
					string(pars.x2star_restrictive_region3,$
					       f='(f0.3)')
			endif

			if tag_exist(pars,'x2star_kde_coadd_region1') then begin
				extra += $
					'-x2star-kde-coadd-region1-'+$
					string(pars.x2star_kde_coadd_region1,$
					       f='(f0.3)')
			endif
			if tag_exist(pars,'x2star_kde_coadd_region2') then begin
				extra += $
					'-x2star-kde-coadd-region2-'+$
					string(pars.x2star_kde_coadd_region2,$
					       f='(f0.3)')
			endif
			if tag_exist(pars,'x2star_kde_coadd_region3') then begin
				extra += $
					'-x2star-kde-coadd-region3-'+$
					string(pars.x2star_kde_coadd_region3,$
					       f='(f0.3)')
			endif



		endif
	endif else begin
		extra+='-density'
	endelse
	fitsfile = repstr(tmp_name, '.fits', '-'+extra+'.fits')

	fitsfile=path_join(outdir, fitsfile)

	return, fitsfile
end


pro esboss::select_by_ra, str, minra, maxra, chunk, indices, tile_subset
	sra = shiftra(str.ra, /wrap)

	indices=where(sra ge minra and sra lt maxra, nw)


	tiles=self->bosstile_read(chunk=chunk)
	sracen = shiftra(tiles.ra,/wrap)
	rad = self->bosstile_radius()

	tile_subset=where( (sracen - rad) gt minra and (sracen+rad) lt maxra)
end

function esboss::_select_by_name, struct, cutname, count=count, $
		tile_subset=tile_subset

	sra = shiftra(struct.ra, /wrap)
	bq=obj_new('bosstarget_qso')
	pars=bq->pars()

	rabounds=pars.ra_dens_bounds

	tiles=self->bosstile_read(chunk=1)
	sracen = shiftra(tiles.ra,/wrap)
	rad = self->bosstile_radius()

	case cutname of
		'lodens': begin
			w=where(sra gt -5. and sra lt 15.)
			tile_subset=where( $
				(sracen - rad) gt -5. and (sracen+rad) lt 15.)
		end
		'region1': begin
			w = where(sra lt rabounds[0], count)
			tile_subset=where( (sracen+rad) lt rabounds[0] )
		end
		'region2': begin
			w=where(sra gt rabounds[0] and sra lt rabounds[1], count)
			tile_subset=where( $
				(sracen - rad) gt rabounds[0] and (sracen+rad) lt rabounds[1])
		end
		'region3': begin
			w=where(sra gt rabounds[1], count)
			tile_subset=where( (sracen-rad) gt rabounds[1] )
		end
		else: begin
			splog,'Unknown cut type: ',cutname,': skipping',form='(a,a,a)'
			count=n_elements(struct)
			w=lindgen(count)
		end
	endcase

	return, w
end


pro esboss::get_kde_coadd_matches, target_run, struct, kde, tile_subset, $
		epars=epars, region=region, notloose=notloose

	bt=obj_new('bosstarget')
	bq=obj_new('bosstarget_qso')

	; read the loosely chosen targets
	if not keyword_set(notloose) then begin
		struct=bt->read_collated('qso', target_run, extra_name='loosecoadd')
	endif else begin
		struct=bt->read_collated('qso', target_run)
	endelse
	; read the kde_coadd catalog
	kde=bq->kde_coadd_read()


	; get the objects that matched by position
	wmatch=where(struct.kde_coadd_id ge 0, nmach)
	help,wmatch
	; these now line up
	struct = struct[wmatch]
	kde = kde[struct.kde_coadd_id]


	wtrim = bq->kde_coadd_trim(kde, ntrim, epars=epars)
	help,wtrim
	struct=struct[wtrim]
	kde=kde[wtrim]


	; down to just the region
	if n_elements(region) ne 0 then begin
		wreg = self->_select_by_name(struct,region,tile_subset=tile_subset)
		help,wreg
		struct = struct[wreg]
		kde = kde[wreg]
	endif


end


; 1-d tuning on sparse region ramod
pro esboss::tune_chunk2_ramod, outdir=outdir, goal=goal
	if n_elements(goal) eq 0 then goal=5.0
	chunk = 2
	target_run = '2009-12-10-newlrg2-loose-sparse'

	bt=obj_new('bosstarget')
	str = bt->read_collated('lrg',target_run, extra_name='noknown')

	print,'selecting window and gal_cmass_sparse'
	target_flags = ['gal_cmass_sparse']

	tilepoly=self->bosstile_poly_read(chunk)
	inwindow=where(is_in_window(tilepoly, ra=str.ra, dec=str.dec), count)

	str=str[inwindow]
	total_area = total(tilepoly.str)*(180d/!dpi)^2

	num=100
	ramod_vals = indgen(num)

	bl=obj_new('bosstarget_lrg')

	densities = fltarr(num)

	for i=0L, num-1 do begin
		pars = {ramod_thresh: ramod_vals[i]}	
		str.boss_target1 = bl->select(str, pars=pars)
		w=self->btselect(str, target_flags, count=count)

		densities[i] = count/total_area
		print,ramod_vals[i], densities[i]
	endfor



	if n_elements(outdir) eq 0 then begin
		outdir=getenv('BOSS_TARGET')
		outdir=path_join(outdir, [target_run,'plots'])
	endif
	file_mkdir, outdir
	psname=string(target_run+'-tune-sparse-'+strn(goal)+'.eps')
	psfile=path_join(outdir, psname)


	begplot,psfile,/encap,xsize=8.5,ysize=8,/color
	!p.charsize=2

	if n_elements(xtitle) eq 0 then xtitle='ramod threshold'
	pplot, ramod_vals, densities, $
		xtitle=xtitle, ytitle='#/square degree', /ynozero

	cutval=interpol(ramod_vals, densities, goal)
	pplot, /over, [cutval,cutval], [0,goal], line=2, color='red'
	pplot, /over, [0,cutval],[goal,goal],line=2,color='red'


	mess=[string('goal=',goal,'/sq deg',f='(a,f0.1,a)'),$
		  string('thresh=',strn(cutval),f='(a,a)')]
	plegend, mess, /right, charsize=1.5, spacing=2.0


	endplot

	colprint, ramod_vals, densities

	return
end



; 1-d tuning on ihi_cmass
pro esboss::tune_chunk2_ihi_cmass, outdir=outdir, goal=goal, num=num
	if n_elements(goal) eq 0 then goal=130
	if n_elements(num) eq 0 then num=20
	chunk = 2
	target_run = '2009-12-09-newlrg2-loose'

	bt=obj_new('bosstarget')
	str = bt->read_collated('lrg',target_run, extra_name='noknown')

	print,'selecting window and lrg'
	target_flags = ['gal_loz','gal_cmass']

	tilepoly=self->bosstile_poly_read(chunk)
	inwindow=where(is_in_window(tilepoly, ra=str.ra, dec=str.dec), count)

	str=str[inwindow]
	total_area = total(tilepoly.str)*(180d/!dpi)^2



	ihi_cmass_vals = arrscl(findgen(num), 19.8, 20.0)


	bl=obj_new('bosstarget_lrg')

	densities = fltarr(num)

	for i=0L, num-1 do begin
		pars = {ihi_cmass: ihi_cmass_vals[i]}	
		str.boss_target1 = bl->select(str, pars=pars)
		w=self->btselect(str, target_flags, count=count)

		densities[i] = count/total_area
	endfor



	if n_elements(outdir) eq 0 then begin
		outdir=getenv('BOSS_TARGET')
		outdir=path_join(outdir, [target_run,'plots'])
	endif
	file_mkdir, outdir
	psname=string(target_run+'-tune-ihi-cmass-'+strn(goal)+'.eps')
	psfile=path_join(outdir, psname)


	begplot,psfile,/encap,xsize=8.5,ysize=8,/color
	!p.charsize=2

	if n_elements(xtitle) eq 0 then xtitle='ihi_cmass threshold'
	pplot, ihi_cmass_vals, densities, $
		xtitle=xtitle, ytitle='#/square degree', /ynozero

	cutval=interpol(ihi_cmass_vals, densities, goal)
	pplot, /over, [cutval,cutval], [0,goal], line=2, color='red'
	pplot, /over, [0,cutval],[goal,goal],line=2,color='red'


	mess=[string('goal=',goal,'/sq deg',f='(a,f0.1,a)'),$
		  string('thresh=',strn(cutval),f='(a,a)')]
	plegend, mess, /right, charsize=1.5, spacing=2.0


	endplot



	return
end



pro esboss::contour_ihi_hiz_vs_ihi_cmass, target_run, chunk, n, $
		struct=struct, doscatter=doscatter

	if n_elements(n) eq 0 then n=10
	n=long64(n)

	bt=obj_new('bosstarget')
	bl=obj_new('bosstarget_lrg', /commiss)
	limits = bl->limits()
	target_flags = ['gal_loz','gal_cmass']
	orflags = sdss_flagval('boss_target1',target_flags)



	; the output files
	outdir = bt->target_dir(target_run)
	outdir = path_join(outdir, 'plots')
	densfile = string($
		target_run+'-density-ihi-hiz-ihi-cmass-',n,'x',n,'.fits',$
		f='(a,i02,a,i02,a)')
	densfile=path_join(outdir, densfile)
	psfile=repstr(densfile,'.fits','.eps')
	scatterfile=repstr(densfile,'.fits','-dperp-mag.eps')


	if not fexist(densfile) then begin

		print,'creating file: ',densfile


		if n_elements(struct) eq 0 then begin
			struct=bt->read_collated('lrg',target_run,extra_name='noknown')
		endif


		; begin with limits from commissioning and go brighter
		help,limits,/str

		ihi_hiz_low = 19.8
		ihi_hiz = arrscl( findgen(n), ihi_hiz_low, limits.ihi_hiz  )
		ihi_cmass_low = 19.8
		ihi_cmass = arrscl( findgen(n), ihi_cmass_low, limits.ihi_cmass  )

		; windows for calculating densities
		tilepoly=self->bosstile_poly_read(chunk)
		total_area=total(tilepoly.str)*(180d/!dpi)^2

		ntot=n*n
		dens = fltarr(n,n)
		itot=0LL
		for i=0LL,n-1 do begin
			limits.ihi_hiz = ihi_hiz[i]
			for j=0L, n-1 do begin
				limits.ihi_cmass = ihi_cmass[j]

				boss_target1 = bl->select(struct, pars=limits,/commissioning)

				w=where( (boss_target1 and orflags) ne 0,nw )

				if nw ne 0 then begin	

					inwindow=where( $
						is_in_window(tilepoly, $
							ra=struct[w].ra, dec=struct[w].dec), nwin)

					dens[i,j] = nwin/total_area
				endif
				print,(itot+1),'/',ntot,'  density: ',dens[i,j],$
					f='(i0,a,i0,a,g0)'
				itot+=1
			endfor
		endfor


		destruct_polygon, tilepoly
		destruct_polygon, boundspoly

		dens_struct = {ihi_hiz: ihi_hiz, ihi_cmass:ihi_cmass, dens:dens}
		print,'Writing file: ',densfile
		mwrfits, dens_struct, densfile, /create
		print,'Run again to make the contour plot'
		return
	endif else begin
		print,'Reading file: ',densfile
		dens_struct = mrdfits(densfile, 1)
	endelse

	begplot, psfile, /encap, xsize=8.5, ysize=8
	nlevel=50
	levels=100+2*lindgen(nlevel)
	c_label = replicate(1,nlevel)
	contour, dens_struct.dens, dens_struct.ihi_hiz, dens_struct.ihi_cmass,$
		levels=levels, c_label=c_label, $
		xtitle='ihi_hiz', ytitle='ihi_cmass', $
		xstyle=3, ystyle=3

	endplot

	if keyword_set(doscatter) then begin
		if n_elements(struct) eq 0 then begin
			struct=bt->read_collated('lrg',target_run,extra_name='noknown')
		endif
		begplot,scatterfile,/encap, xsize=8.5, ysize=8, /color

		mst=bl->magstruct(struct)
		plot, mst.cmodelmag[3], mst.dperp, psym=3, $
			xtitle='cmodelmag[3]', ytitle='dperp', $
			xrange=[17.5, 20.5], yrange=[0.25,1.2], $
			xstyle=3, ystyle=3

		limits.ihi_hiz = 19.90
		limits.ihi_cmass= 19.86
		;limits.ihi_cmass= 19.5
		boss_target1 = bl->select(struct, pars=limits,/commissioning)
		w=where( (boss_target1 and orflags) ne 0,nw, comp=comp )
		pplot, mst[comp].cmodelmag[3], mst[comp].dperp, psym=3, $
			color='red', /overplot

		legend, ['all','removed'], psym=8, $
			color=[!p.color, c2i('red')], /left
		legend, $
			['ihi_hiz='+string(limits.ihi_hiz,f='(f0.2)'), $
			 'ihi_cmass='+string(limits.ihi_cmass,f='(f0.2)')], $
			/right,/bottom

		endplot
	endif

end



pro esboss::contour_qso_chunk1_tile_density_lodens, target_run, $
		nvals=nvals, $
		door=door, $
		commissioning=commissioning, $
		density=density, bfratio=bfratio, generate=generate, struct=bigstruct


	if n_params() lt 1 then begin
		on_error, 2
		message,'usage: '
	endif
	
	qso_type='qso_kde_coadd'

	if n_elements(nvals) eq 0 then nvals = 20
	min_bright = -0.1
	max_bright = 0.3
	min_faint = -0.1
	max_faint = 0.3
	sdbr = arrscl(findgen(nvals), min_bright, max_bright)
	sdfr = arrscl(findgen(nvals), min_faint, max_faint)


	bt=obj_new('bosstarget')
	bq =obj_new('bosstarget_qso')
	pars=bq->pars()
	print,pars.kde_coadd_x2_star_quadfit
	if n_elements(density) eq 0 or n_elements(bfratio) eq 0 then begin



		; read the loosely chosen targets
		if n_elements(bigstruct) eq 0 then begin
			bigstruct=$
				bt->read_collated('qso', target_run, extra_name='loosecoadd')
		endif
		; read the kde_coadd catalog
		kde=bq->kde_coadd_read()


		; get the objects that matched by position
		wmatch=where(bigstruct.kde_coadd_id ge 0, nmach)

		; these now line up
		struct = bigstruct[wmatch]
		kde = kde[struct.kde_coadd_id]



		; some basic trimming
		epars={kde_coadd_logstardens_bright: 5.0, $
			kde_coadd_logstardens_faint: 5.0,$
			kde_coadd_x2_star_quadfit: [3.0,0.0,0.0]}

		wtrim = bq->kde_coadd_trim(kde, ntrim, epars=epars)
		struct=struct[wtrim]
		kde=kde[wtrim]



		; down to just the lodens
		wreg = self->_select_by_name(struct,'lodens',tile_subset=tile_subset)
		struct = struct[wreg]
		kde = kde[wreg]


		; we'll need these
		lups=bq->get_lups(struct, /deredden)

		ntot=nvals*nvals
		itot=0

		density=dblarr(nvals,nvals)
		bfratio = density
		for i=0L, nvals-1 do begin
			stardens_bright = 10d^sdbr[i]
			for j=0L, nvals-1 do begin
				stardens_faint = 10d^sdfr[j]

				;print,stardens_bright,stardens_faint

				blogic = kde.kde_stardens_bright lt stardens_bright
				flogic = kde.kde_stardens_faint lt stardens_faint
				if keyword_set(door) then begin
					logic=blogic or flogic 
				endif else begin
					logic=blogic and flogic
				endelse
				w=where(logic, nw)

				chunk=1
				; this is not broken since we removed tile_subset keyword
				ds=self->bosstile_density(kde[w].ra, kde[w].dec, chunk, $
						tile_subset=tile_subset)

				density[i,j] = median(ds.density)

				wfaint = where(lups[1,w] gt pars.gsplit, nfaint)
				wbright = where(lups[1,w] lt pars.gsplit, nbright)


				ra = kde[w[wfaint]].ra
				dec = kde[w[wfaint]].dec

				ds_faint=self->bosstile_density(ra,dec,chunk,$
					tile_subset=tile_subset)
				ra = kde[w[wbright]].ra
				dec = kde[w[wbright]].dec
				ds_bright=self->bosstile_density(ra, dec, chunk, $
					tile_subset=tile_subset)

				bfratio[i,j] = $
					median(ds_bright.density)/median(ds_faint.density)

				itot += 1
				print,itot,'/',ntot,density[i,j],bfratio[i,j],$
					form='(i0,a,i0," ",f0.2," ",f0.2)'
			endfor
		endfor
	endif

	bt=obj_new('bosstarget')
	outdir = bt->target_dir(target_run)
	outdir = path_join(outdir, 'density')

	if keyword_set(door) then begin
		short_type='kc-or'
	endif else begin
		short_type='kc'
	endelse
	file=self->qso_tile_density_file(target_run,qso_type,none,'loosecoadd')
	psfile=repstr(file,'.fits','-'+short_type+'-contours-2par.eps')

	print,psfile
	begplot,psfile,/encap,xsize=8.5,ysize=8, /color
	




	levels = 30 + 1*lindgen(60)
	;yrange=[0.0,0.6]
	yrange=[min_faint-0.1, max_faint]

	nlevel=n_elements(levels)
	c_label = replicate(1,nlevel)

	xtitle='log(star dens max bright)'
	ytitle='log(star dens max faint)'


	;nlevel=10
	contour, density, sdbr, sdfr, $
		levels=levels, $
		c_label=c_label,$
		xtitle=xtitle, ytitle=ytitle, $
		xstyle=3, $
		yrange=yrange, ystyle=3


	;if in(qso_type, 'qso_bonus') then step=0.10 else step=0.20
	ratstep=0.02
	levels = 0.1+ratstep*findgen(100)
	nlevel = n_elements(levels)
	c_label = replicate(1,nlevel)

	ocolor=c2i('blue')
	oline=0
	contour, bfratio, sdbr, sdfr, $
		levels=levels, c_label=c_label, /overplot, c_color=ocolor, $
		c_line=oline

	legend,/bottom,$
		['qso',target_run,$
		strjoin(qso_type,' or '),$
		'x2_star=3']
	legend,/right,/bottom,$
		['density #/sqdeg','bright/faint ratio'], $
		line=[0,oline], $
		color=[!p.color,ocolor]

	;pplot, [-0.02], [-0.02], psym=7, color='red', /overplot
	pplot, [0.065], [0.065], psym=7, color='red', /overplot

	endplot


end

pro esboss::compare_chunks_tune1d, $
		chunks, target_type, target_run, goal, target_flags, $
		extra_names=extra_names, range=range

	if n_elements(outdir) eq 0 then begin
		outdir=getenv('BOSS_TARGET')
		outdir=path_join(outdir, [target_run,'plots'])
	endif

	file_mkdir, outdir
	tf=strjoin(target_flags, '-')

	tstr = target_run
	if n_elements(extra_name) ne 0 then begin
		tstr += '-'+extra_name
	endif

	chstr=strjoin(string(chunks,f='(i0)'),'-')
	goalstr=string(goal,f='(i0)')

	psname=string($
		tstr, chstr, goal, repstr(tf, '_', '-'), $
		f='("compare-",a,"-chunk",a,"-tune",i0,"-",a,".eps")')

	psfile=path_join(outdir, psname)
	print,psfile

	begplot,psfile,/encap,xsize=8.5,ysize=8,/color
	!p.charsize=2

	nextra=n_elements(extra_names)

	nchunk=n_elements(chunks)
	colors=make_rainbow(nchunk)

	ptrlist=ptrarr(nchunk)

	lines=[0,2,3,4]
	lines=lines[0:nchunk-1]
	for i=0L,nchunk-1 do begin
		if n_elements(extra_names) gt 0 then begin
			if nextra eq 1 then begin
				extra_name=extra_names 
			endif else begin
				extra_name=extra_names[i]
			endelse
		endif
		case target_flags of 
			'qso_nn': begin
				self->tune_nn1d,chunks[i],target_run,goal,$
					denstruct=denstruct, extra_name=extra_name, /noplot, $
					range=range
				xtitle='xnn2 threshold'
			end
			'qso_like': begin
				self->tune_like1d,chunks[i],target_run,goal,$
					denstruct=denstruct, extra_name=extra_name, /noplot, $
					range=range

				xtitle='likelihood threshold'
			end
			'qso_kde': begin
				self->tune_kde1d,chunks[i],target_run,goal,$
					denstruct=denstruct, extra_name=extra_name, /noplot, $
					range=range

				xtitle='kdeprob threshold'
			end
			else:message,'bad target_flags: '+string(target_flags)
		endcase

		if i eq 0 then begin
			pplot, denstruct.vals, denstruct.density, $
				xtitle=xtitle, ytitle='#/square degree', /ynozero, $
				line=lines[i]
		endif

		pplot, /overplot, denstruct.vals, denstruct.density, color=colors[i], $
			line=lines[i]

		pplot, /over, $
			[denstruct.bestval,denstruct.bestval], [0,goal], line=1, $
			color=colors[i]
		pplot, /over, [0,denstruct.bestval],[goal,goal],line=1, $
			color=colors[i]


		add_arrval, string(chunks[i],denstruct.bestval, $
			f='("chunk",i0,": ",g0)'), mess
		delvarx, denstruct
	endfor

	plegend, mess, line=lines, colors=colors,/right
	endplot, /png, dpi=75, /trim


end



pro esboss::tune_all1d, chunk, target_run, goal, str=str, extra_name=extra_name


	self->tune_nn1d, chunk, target_run, goal, $
		str=str, extra_name=extra_name

	self->tune_kde1d, chunk, target_run, goal, $
		str=str, extra_name=extra_name

	self->tune_like1d, chunk, target_run, goal, $
		str=str, extra_name=extra_name

	self->tune_mclike1d, chunk, target_run, goal, $
		str=str, extra_name=extra_name


	self->tune_nn_value1d, chunk, target_run, goal, $
		str=str, extra_name=extra_name


end


pro esboss::tune_nn_value1d, chunk, target_run, goal, range=range, $
		str=str, denstruct=denstruct, number=number, $
		reselect=reselect, outdir=outdir, extra_name=extra_name, $
		noplot=noplot


	if n_elements(range) eq 0 then range=[0.3,1.0]

	threshtype='gt'
	target_type='qso'
	; these should be very loose
	target_flags=['qso_like','qso_nn','qso_kde']
	;target_flags='qso_bonus_main'
	tagname='nn_value'

	if n_elements(str) eq 0 or keyword_set(reselect) then begin
		bt=obj_new('bosstarget')
		str=bt->read_collated(target_type, target_run,extra_name=extra_name)
	endif

	if n_elements(extra_name) ne 0 then begin
		; to distinguish from the normal likelihood
		ename = extra_name + '-nnval'
	endif else begin
		ename = 'nnval'
	endelse



	self->tune_chunk_density1d, $
		threshtype, chunk, target_type, target_flags, target_run, tagname, $
		range, goal, $
		extra_name=ename, $
		str=str, denstruct=denstruct, number=number, $
		reselect=reselect, $
		outdir=outdir, $
		xtitle='NN value', $
		noplot=noplot

	return
end

pro esboss::tune_mclike1d, chunk, target_run, goal, range=range, $
		str=str, denstruct=denstruct, number=number, $
		reselect=reselect, outdir=outdir, extra_name=extra_name, $
		noplot=noplot


	if n_elements(range) eq 0 then range=[0.03,0.8]

	; WARNING: this only works because we have been doing such loose
	; cuts on like_ratio > 0, which also implies like_ratio_pat > 0
	threshtype='gt'
	target_type='qso'
	target_flags='qso_like'
	tagname='like_ratio_pat'


	if n_elements(str) eq 0 or keyword_set(reselect) then begin
		bt=obj_new('bosstarget')
		str=bt->read_collated(target_type, target_run,extra_name=extra_name)
	endif

	if n_elements(extra_name) ne 0 then begin
		; to distinguish from the normal likelihood
		ename = extra_name + '-pat'
	endif else begin
		ename = 'pat'
	endelse
	self->tune_chunk_density1d, $
		threshtype, chunk, target_type, target_flags, target_run, tagname, $
		range, goal, $
		extra_name=ename, $
		str=str, denstruct=denstruct, number=number, $
		reselect=reselect, $
		outdir=outdir, $
		xtitle='weighted likelihood threshold', $
		noplot=noplot

	return
end



pro esboss::tune_like1d, chunk, target_run, goal, range=range, $
		str=str, denstruct=denstruct, number=number, $
		reselect=reselect, outdir=outdir, extra_name=extra_name, $
		noplot=noplot


	if n_elements(range) eq 0 then range=[0.03,0.8]

	threshtype='gt'
	target_type='qso'
	target_flags='qso_like'
	tagname='like_ratio'
	self->tune_chunk_density1d, $
		threshtype, chunk, target_type, target_flags, target_run, tagname, $
		range, goal, $
		extra_name=extra_name, $
		str=str, denstruct=denstruct, number=number, $
		reselect=reselect, $
		outdir=outdir, $
		xtitle='likelihood threshold', $
		noplot=noplot

	return
end

pro esboss::tune_nn1d, chunk, target_run, goal, range=range, $
		str=str, denstruct=denstruct, number=number, $
		reselect=reselect, outdir=outdir, extra_name=extra_name, $
		noplot=noplot


	if n_elements(range) eq 0 then range=[0.2,1.0]

	threshtype='gt'
	target_type='qso'
	target_flags='qso_nn'
	;tagname='nn_xnn2'
	tagname='nn_xnn'
	self->tune_chunk_density1d, $
		threshtype, chunk, target_type, target_flags, target_run, tagname, $
		range, goal, $
		extra_name=extra_name, $
		str=str, denstruct=denstruct, number=number, $
		reselect=reselect, $
		outdir=outdir, $
		xtitle='xnn threshold', $
		noplot=noplot

	return
end

pro esboss::tune_kde1d, chunk, target_run, goal, range=range, $
		str=str, denstruct=denstruct, number=number, $
		reselect=reselect, outdir=outdir, extra_name=extra_name, $
		noplot=noplot

	if n_elements(range) eq 0 then range=[0.0,1.0]

	threshtype='gt'
	target_type='qso'
	target_flags='qso_kde'
	tagname='kde_prob'
	self->tune_chunk_density1d, $
		threshtype, chunk, target_type, target_flags, target_run, tagname, $
		range, goal, $
		extra_name=extra_name, $
		str=str, denstruct=denstruct, number=number, $
		reselect=reselect, $
		outdir=outdir, $
		xtitle='kdeprob threshold', $
		noplot=noplot

	return
end





pro esboss::tune_chunk2_like, target_run, goal, str=str, denstruct=denstruct, n=n, $
		reselect=reselect, outdir=outdir

	; note it's not really contour here since it is one dimensional
	; this is "commissioning 2, or comm2"
	; will look for about 20-25/square degree for both bright and faint?
	chunk = 2
	threshtype='gt'
	goal=35.
	tagname='like_ratio'
	range=[0.1,1.0]
	self->tune_chunk_density1d, $
		'gt', chunk, 'qso', 'qso_like', target_run, tagname, range, goal, $
		str=str, denstruct=denstruct, n=n, $
		reselect=reselect, $
		outdir=outdir, $
		xtitle='likelihood threshold'
	return
end

pro esboss::tune_chunk2_nn, target_run, str=str, denstruct=denstruct, n=n, $
		reselect=reselect, outdir=outdir
	; note it's not really contour here since it is one dimensional
	; this is "commissioning 2, or comm2"
	; will look for about 20-25/square degree for both bright and faint?
	chunk = 2
	threshtype='gt'
	goal=20.0
	tagname='nn_xnn'
	range=[0.65,0.9]
	self->tune_chunk_density1d, $
		'gt', chunk, 'qso', 'qso_nn', target_run, tagname, range, goal, $
		str=str, denstruct=denstruct, n=n, $
		reselect=reselect, $
		outdir=outdir, $
		xtitle='nn_xnn threshold'
	return
end



pro esboss::tune_chunk_density1d, $
		threshtype, $
		chunk, $
		target_type, $
		target_flags, $
		target_run, $
		tagname, $
		range, $
		goal, $
		noplot=noplot, $
		extra_name=extra_name, $
		reselect=reselect, $
		str=str, $
		denstruct=denstruct, $
		number=n, $
		outdir=outdir, $
		xtitle=xtitle

	if n_params() lt 8 then begin
		print,'Usage: '
		print,'  eb->tune_chunk_density1d, '
		print,'      threshtype, chunk, target_type, target_flags, target_run, tagname, range, goal, '
		print,'      str=, /reselect,dens=, n=, outdir=, xtitle='
		print
		print,"  threshtype is 'gt' or 'lt'"
		print,"  target_type should be 'lrg', 'qso', etc."
		print,"  target_flags should be 'qso_like' etc."
		print,"  tagname is the tag name to tune on"
		print,"  goal is the desired density goal"
		on_error, 2
		message,'Halting'
	endif

	bt=obj_new('bosstarget')
	if n_elements(str) eq 0 or keyword_set(reselect) then begin
		str=bt->read_collated(target_type, target_run,extra_name=extra_name)
	endif

	;tilepoly=self->bosstile_poly_read(chunk)
	tilepoly=self->bosstile_poly_read(chunk,/bounds)
	total_area=total(tilepoly.str)*(180d/!dpi)^2


	if keyword_set(reselect) then begin
		if chunk eq 2 then comm2=1
		if chunk eq 1 then commissioning=1

		btypeobj=obj_new('bosstarget_'+target_type)
		res = btypeobj->select(str, $
			comm2=comm2, commissioning=commissioning, /reselect)
		str.boss_target1 = res.boss_target1
	endif

	if not tag_exist(str, tagname, ind=tagind) then begin
		message,"Tagname '",tagname,"' does not exists",f='(a,a,a)'
	endif

	; initial selection to make sure *other* selections are still in place
	wtype=self->btselect(str, target_flags)

	; select on the window first
	ra=str.ra
	dec=str.dec
	isin=is_in_window(tilepoly, ra=ra, dec=dec)

	; now targetrs in the window
	wtype2 = where( isin[wtype] )
	wtype = wtype[wtype2]

	if n_elements(n) eq 0 then n=100L


	;testvals=arrscl( findgen(n), range[0], range[1] )

	if n_elements(denstruct) eq 0 then begin
		print,'Calculating densities'

		denstruct = {$
			goal: goal, $
			bestval: 0.0, $
			vals: fltarr(n), $
			density: fltarr(n) $
		}
		denstruct.vals = arrscl( findgen(n), range[0], range[1] )

		for i=0L, n-1L do begin
			if threshtype eq 'gt' then begin
				w=where( str[wtype].(tagind) gt denstruct.vals[i], nw)	
			endif else begin
				w=where( str[wtype].(tagind) lt denstruct.vals[i], nw)	
			endelse

			denstruct.density[i] = nw/total_area

			if ((i MOD 10) EQ 0) then print,i,n,string(13b),format='(i,i,a,$)'
		endfor
	endif

	denstruct.bestval = interpol(denstruct.vals, denstruct.density,  goal)

	if keyword_set(noplot) then return


	if n_elements(outdir) eq 0 then begin
		outdir=getenv('BOSS_TARGET')
		outdir=path_join(outdir, [target_run,'plots'])
	endif
	file_mkdir, outdir
	tf=strjoin(target_flags, '-')

	tstr = target_run
	if n_elements(extra_name) ne 0 then begin
		tstr += '-'+extra_name
	endif

	if size(chunk,/tname) eq 'STRING' then begin
		chstr = strjoin(chunk,'-')
	endif else begin
		chstr=strjoin(string(chunk,f='(i0)'),'-')
	endelse

	goalstr=string(goal,f='(i0)')

	psname=string($
		tstr, chstr, goal, repstr(tf, '_', '-'), $
		f='(a,"-chunk",a,"-tune",i0,"-",a,".eps")')

	psfile=path_join(outdir, psname)

	begplot,psfile,/encap,xsize=8.5,ysize=8,/color
	!p.charsize=2

	if n_elements(xtitle) eq 0 then xtitle='threshold'
	pplot, denstruct.vals, denstruct.density, $
		xtitle=xtitle, ytitle='#/square degree', /ynozero

	pplot, /over, $
		[denstruct.bestval,denstruct.bestval], [0,goal], line=2, color='red'
	pplot, /over, [0,denstruct.bestval],[goal,goal],line=2,color='red'


	mess=[string('goal=',goal,'/sq deg',f='(a,f0.1,a)'),$
		  string('thresh=',strn(denstruct.bestval),f='(a,a)')]
	plegend, mess, /right, charsize=1.5, spacing=2.0


	endplot, /png, dpi=75, /trim

end


pro esboss::tune_qso_postcomm_loose, chunk, target_run, $
		str=str, inwindow=inwindow
	; tune to the loose ~100/sq degree we want for the ranking method

	target_density = 100.0 ;per square degree
	;target_density = 80.0 ;per square degree

	if n_elements(outdir) eq 0 then begin
		outdir=getenv('BOSS_TARGET')
		outdir=path_join(outdir, [target_run,'plots'])
	endif
	if not file_test(outdir) then begin
		file_mkdir, outdir
	endif


	tilepoly=self->bosstile_poly_read(chunk)
	total_area=total(tilepoly.str)*(180d/!dpi)^2


	if n_elements(str) eq 0 then begin
		bt=obj_new('bosstarget')
		str=bt->read_collated('qso',target_run)
	endif

	if n_elements(inwindow) eq 0 then begin
		print,'Checking against window'
		inwindow=where(is_in_window(tilepoly, ra=str.ra, dec=str.dec), nwin)
	endif

	;
	; kde core+bonus
	;
	if 0 then begin
		; we already tuned this separately
		target_flags = ['qso_core','qso_bonus']
		w=self->btselect(str[inwindow], target_flags)
		w=inwindow[w]

		nvals=100
		x2star_thresh = arrscl(findgen(nvals), 0.1, 2)

		dens=fltarr(nvals)
		print,'Calculating kde densities'
		for i=0L, nvals-1 do begin
			; note, this won't affect core since it is already x2_star > 7
			wthresh=where(str[w].x2_star gt x2star_thresh[i], nthresh)
			dens[i] = nthresh/total_area
		endfor

		thresh_best=interpol(x2star_thresh, dens, target_density)

		psfile=path_join(outdir, $
			target_run+'-qso-tune-loose-kde.eps')
		begplot,psfile,xsize=8.5,ysize=8.0,/color,/encap

		pplot, x2star_thresh, dens, $
			xtitle='x2_star permissive threshold', ytitle='density #/sq deg', $
			/ynozero

		pplot, [0,thresh_best],[target_density,target_density], $
			/overplot,color='red'	
		pplot, [thresh_best,thresh_best],[0,target_density], /overplot,color='red'
		plegend,['['+strjoin(target_flags,',')+'] target density: '+$
						  string(target_density,f='(i0)'),$
				'threshold: '+string(thresh_best,f='(f0.2)')],/right	

		endplot
	endif


	;
	; NN
	;
	target_flags='qso_nn'
	w=self->btselect(str[inwindow], target_flags)
	w=inwindow[w]

	nvals=100
	xnn_thresh = arrscl(findgen(nvals), 0.4, 0.8)

	dens=fltarr(nvals)
	print,'Calculating NN densities'
	for i=0L, nvals-1 do begin
		wthresh=where(str[w].nn_xnn gt xnn_thresh[i], nthresh)
		dens[i] = nthresh/total_area
	endfor

	thresh_best=interpol(xnn_thresh, dens, target_density)

	psfile=path_join(outdir, $
		target_run+'-qso-tune-loose-nn.eps')
	begplot,psfile,xsize=8.5,ysize=8.0,/color,/encap

	pplot, xnn_thresh, dens, $
		xtitle='XNN threshold', ytitle='density #/sq deg'

	pplot, [0,thresh_best],[target_density,target_density], $
		/overplot,color='red'	
	pplot, [thresh_best,thresh_best],[0,target_density], /overplot,color='red'

	plegend,['['+strjoin(target_flags,',')+'] target density: '+$
		              string(target_density,f='(i0)'),$
			'threshold: '+string(thresh_best,f='(f0.2)')],/right	
	endplot

	;
	; likelihood
	;
	target_flags='qso_like'
	w=self->btselect(str[inwindow], target_flags)
	w=inwindow[w]

	nvals=100
	like_thresh = arrscl(findgen(nvals), 0.02, 0.2)

	dens=fltarr(nvals)
	print,'Calculating like densities'
	for i=0L, nvals-1 do begin
		wthresh=where(str[w].like_ratio gt like_thresh[i], nthresh)
		dens[i] = nthresh/total_area
	endfor

	thresh_best=interpol(like_thresh, dens, target_density)

	psfile=path_join(outdir, $
		target_run+'-qso-tune-loose-likelihood.eps')
	begplot,psfile,xsize=8.5,ysize=8.0,/color,/encap

	pplot, like_thresh, dens,$
		xtitle='like ratio threshold', ytitle='density #/sq deg'
	pplot, [0,thresh_best],[target_density,target_density], $
		/overplot,color='red'

	pplot, [thresh_best,thresh_best],[0,target_density], /overplot,color='red'
	plegend,['target density: '+string(target_density,f='(i0)'),$
			'threshold: '+string(thresh_best,f='(f0.2)')],/right	
	plegend,['['+strjoin(target_flags,',')+'] target density: '+$
		              string(target_density,f='(i0)'),$
			'threshold: '+string(thresh_best,f='(f0.2)')],/right	

	endplot




	destruct_polygon, tilepoly

end



pro esboss::contour_chunk_bonus, chunk, target_run, num, $
		denstruct=denstruct


	flag='qso_bonus'
	if keyword_set(addcore) then flag=[flag,'qso_core']
	flagval=sdss_flagval('boss_target1',flag)


	; values to check
	;log_stardensmax_bright_vals = arrscl(findgen(n), -0.5, 0.32)
	;log_stardensmax_faint_vals = arrscl(findgen(n), -1.5, -0.21)
	;log_stardensmax_bright_vals = arrscl(findgen(n), -0.1-0.2, -0.1+0.2)
	;log_stardensmax_faint_vals = arrscl(findgen(n), -.5-0.2, -0.5+0.2)

	log_stardensmax_bright_vals = arrscl(findgen(num), 0.2,0.4)
	log_stardensmax_faint_vals = arrscl(findgen(num), -0.2,0.3)



	if n_elements(denstruct) eq 0 then begin

		bt=obj_new('bosstarget')
		str=bt->read_collated('qso', target_run)


		; First get the ones in the window.  This will greatly speed up things
		; later
		print,'Getting in-window objects'
		tilepoly=self->bosstile_poly_read(chunk)
		total_area=total(tilepoly.str)*(180d/!dpi)^2
		ra =str.ra
		dec=str.dec
		inwindow=where( is_in_window(tilepoly, ra=ra, dec=dec), nwin)
		destruct_polygon,tilepoly

		faint  = where($
			(str[inwindow].boss_target1 and flagval) ne 0 $
			and str[inwindow].gfaint ne 0 $
			and alog10(str[inwindow].kde_qsodens_faint) ge -0.57 )
		bright  = where($
			(str[inwindow].boss_target1 and flagval) ne 0 $
			and str[inwindow].gfaint eq 0 $
			and alog10(str[inwindow].kde_qsodens_bright) ge -0.57)

		faint = inwindow[faint]
		bright = inwindow[bright]

		log_stardens_bright = alog10(str[bright].kde_stardens_bright)
		log_stardens_faint = alog10(str[faint].kde_stardens_faint)


		print,'Calculating densities'
		base=fltarr(num,num)
		denstruct={ $
			dens:base, $
			bdens:base, $
			fdens:base, $
			bfratio:base }

		ntot=num*num
		i=0L
		for bi=0L, num-1 do begin
			wb=where( $
				log_stardens_bright lt log_stardensmax_bright_vals[bi], nb )
			bdensity = nb/total_area
			for fi=0L, num-1 do begin

				wf=where( $
					log_stardens_faint lt log_stardensmax_faint_vals[fi], nf )

				fdensity = nf/total_area

				denstruct.bdens[bi, fi] = bdensity
				denstruct.fdens[bi, fi] = fdensity
				denstruct.dens[bi, fi] = bdensity + fdensity

				if fdensity gt 0 then begin
					denstruct.bfratio[bi, fi] = $
						denstruct.bdens[bi, fi]/denstruct.fdens[bi, fi]
				endif


				if ((i MOD 10) EQ 0) then begin
					print,i,'/',ntot,' ',format='(i,a,i,a,$)'
				endif
				i+=1
			endfor
		endfor
	endif
	print


	outdir=getenv('BOSS_TARGET')
	outdir=path_join(outdir, [target_run,'plots'])
	file_mkdir, outdir
	psfile=path_join(outdir, target_run+'-contour-bonus-bright-faint.eps')
	if keyword_set(addcore) then psfile=repstr(psfile,'bonus','core-bonus')

	begplot,psfile,/encap,xsize=8.5,ysize=8,/color
	!p.charsize=2

	nlevel=50
	levels = 30 + 2*lindgen(nlevel)
	c_label = replicate(1,nlevel)
	contour, denstruct.dens, log_stardensmax_bright_vals, log_stardensmax_faint_vals, $
		levels=levels, $
		c_label=c_label, $
		xtitle='log(stardens_bright)', ytitle='log(stardens_faint)'



	; add the bright/faint ratio
	ratstep=0.05
	nlevel=200
	levels = 0.1+ratstep*findgen(nlevel)
	c_label = replicate(1,nlevel)

	ocolor=c2i('blue')
	oline=0
	contour, $
		denstruct.bfratio, log_stardensmax_bright_vals, log_stardensmax_faint_vals, $
		levels=levels, c_label=c_label, /overplot, c_color=ocolor, $
		c_line=oline

	qdb=0.325
	qdf=0.107
	pplot, [qdb], [qdf], psym=1, color='red', symsize=2, /overplot
	mess=string('sdb=',qdb,' sdf=',qdf,f='(a,f0.3,a,f0.3)')
	plegend, psym=1, color='red', mess, /right, $
		charsize=1.5


	endplot
end



pro esboss::contour_chunk2_bonus_qsodens, target_run, str=str, $
		denstruct=denstruct, n=n

	; this is not used any more

	; this is "commissioning 2, or comm2"
	; will look for about 20-25/square degree for both bright and faint?
	chunk = 2

	bt=obj_new('bosstarget')
	if n_elements(str) eq 0 then begin
		str=bt->read_collated('qso', target_run)
	endif

	pstr = extract_tags(str, ['ra', 'dec'])

	if n_elements(n) eq 0 then n=10L


	;log_qsodens_bright_vals = arrscl(findgen(n), -1.4, 1.0)
	;log_qsodens_faint_vals = arrscl(findgen(n), -5.5, 1.0)
	log_qsodens_bright_vals = arrscl(findgen(n), 0.035-0.1, 0.035+0.1)
	log_qsodens_faint_vals = arrscl(findgen(n), -0.35-0.2, -0.35+0.2)

	wbonus = self->btselect(str, 'qso_bonus')
	faint=where(str[wbonus].gfaint)
	bright=where(str[wbonus].gfaint eq 0)
	faint=wbonus[faint]
	bright=wbonus[bright]

	log_kde_qsodens_bright = alog10(str[bright].kde_qsodens_bright)
	log_kde_qsodens_faint = alog10(str[faint].kde_qsodens_faint)


	if n_elements(denstruct) eq 0 then begin
		print,'Calculating densities'
		base=fltarr(n,n)
		denstruct={ $
			dens:base, $
			bdens:base, $
			fdens:base, $
			bfratio:base }

		ntot=n*n
		i=0L
		for bi=0L, n-1 do begin
			; note log_kde_qsodens_bright is already subscripted with bright
			wb=where( $
				log_kde_qsodens_bright gt log_qsodens_bright_vals[bi], nb )
			
			ra=pstr[bright[wb]].ra
			dec=pstr[bright[wb]].dec
			dsb=self->bosstile_density(ra,dec,chunk)
			for fi=0L, n-1 do begin

				; note log_kde_qsodens_faint is already subscripted with faint
				wf=where( $
					log_kde_qsodens_faint gt log_qsodens_faint_vals[fi], nf )
				ra=pstr[faint[wf]].ra
				dec=pstr[faint[wf]].dec
				dsf=self->bosstile_density(ra,dec,chunk)

				denstruct.bdens[bi, fi] = median( dsb.density )
				denstruct.fdens[bi, fi] = median( dsf.density )
				denstruct.dens[bi, fi] = $
					denstruct.bdens[bi, fi] + denstruct.fdens[bi, fi]

				if denstruct.fdens[bi, fi] gt 0 then begin
					denstruct.bfratio[bi, fi] = $
						denstruct.bdens[bi, fi]/denstruct.fdens[bi, fi]
				endif


				if ((i MOD 10) EQ 0) then print,i,ntot,string(13b),format='(i,i,a,$)'
				i+=1
			endfor
		endfor
	endif


	outdir=getenv('BOSS_TARGET')
	outdir=path_join(outdir, [target_run,'plots'])
	file_mkdir, outdir
	pngfile=path_join(outdir, 'contour-bonus-bright-faint.png')
	psfile=path_join(outdir, 'contour-bonus-bright-faint.eps')

	begplot,psfile,/encap,xsize=8.5,ysize=8,/color
	;devold=!d.name
	;setupplot, 'Z'
	;device, set_resolution=[800,800]
	;!p.color=c2i('black', rct=rct, gct=gct, bct=bct)
	;!p.background=c2i('white')
	;!p.charsize=2
	;print,!p.color,!p.background
	;help,rct,gct,bct

	nlevel=50
	levels = 15 + 1*lindgen(nlevel)
	c_label = replicate(1,nlevel)
	contour, denstruct.dens, log_qsodens_bright_vals, log_qsodens_faint_vals, $
		levels=levels, $
		c_label=c_label, $
		xtitle='log(qsodens_bright)', ytitle='log(qsodens_faint)'



	; add the bright/faint ratio
	ratstep=0.1
	nlevel=50
	levels = 0.1+ratstep*findgen(nlevel)
	c_label = replicate(1,nlevel)

	ocolor=c2i('blue')
	oline=0
	contour, $
		denstruct.bfratio, log_qsodens_bright_vals, log_qsodens_faint_vals, $
		levels=levels, c_label=c_label, /overplot, c_color=ocolor, $
		c_line=oline

	qdb=0.035
	qdf=-0.370
	pplot, [qdb], [qdf], psym=7, color='red', symsize=0.8, /overplot
	mess=string('qdb=',qdb,' qdf=',qdf,f='(a,f0.3,a,f0.2)')
	legend, psym=7, color='red', mess, /right


	;write_png, pngfile, tvrd(), rct, gct, bct
	;setupplot, devold
	endplot
end



pro esboss::plot_tuned_densities, str
	bq=obj_new('bosstarget_qso')
	res=bq->select(str, /reselect, /struct, /commiss)
	tstr = str
	tstr.boss_target1 = res.boss_target1

	wk=self->btselect(tstr, 'qso_kde_coadd')
	chunk=1
	ds = self->bosstile_density(tstr[wk].ra,tstr[wk].dec,chunk)

	print,'Median density: ',median(ds.density)

	begplot,'~/tmp/test.ps', /landscape
	pplot, shiftra(ds.ra,/wrap), ds.density, psym=8, $
		xtitle='RA', ytitle='#/square degree'
	endplot,/landfix
end


pro esboss::tune_chunk1_qso_tile_density_kde_coadd_by_ra, $
		target_run, nvals=nvals, regenerate=regenerate, $
		struct=struct, kde=kde

	; the density we are aiming for in each region
	density_goal = 50.
	chunk=1


	; define the ra regions
	nregions = 15	
	ravals = arrscl( findgen(nregions+1), -40, 40 )



	; number of values in x2_star to check
	if n_elements(nvals) eq 0 then nvals=20




	; if we re-do this after we will have to instead put in a 
	; flat interpolation
	if n_elements(struct) eq 0 or n_elements(kde) eq 0 then begin
		; loose flat cut
		epars={kde_coadd_x2_star_quadfit: [2.0,0.0,0.0]}
		self->get_kde_coadd_matches, target_run, struct, kde, tile_subset, $
			epars=epars
	endif



	; range of x2star to explore
	minx2star=2.
	maxx2star=20.

	; values of x2star to check
	x2star_vals = $
		arrscl(findgen(nvals), minx2star, maxx2star)


	x2star_best = fltarr(nregions)

	; loop over the regions and get the x2star that gives us the required
	; density goal

	bt=obj_new('bosstarget')
	outdir = bt->target_dir(target_run)
	dirpart='tune-chunk1-density-vs-ra'
	outdir = path_join(outdir, dirpart)
	file_mkdir, outdir
	psfile=target_run+'-'+dirpart+'-each.eps'
	psfile=path_join(outdir, psfile)


	begplot,psfile,/encap,xsize=8.5,ysize=8, /color
	colors=make_rainbow(nregions)
	pplot, [0],/nodata,xrange=[minx2star,maxx2star], yrange=[0,100], $
		xtitle='x2star', $
		ytitle='#/square degree'
	

	for ireg=0L, nregions-1 do begin
		minra = ravals[ireg]
		maxra = ravals[ireg+1]

		midra=[minra+maxra]/2.
		add_arrval, midra, middle_ravals

		self->select_by_ra, struct, minra, maxra, chunk, indices, tile_subset

		density=fltarr(nvals)
		for i=0L, nvals-1 do begin

			x2star=x2star_vals[i]

			w=where(kde[indices].x2_star gt x2star)
			
			ra=struct[indices[w]].ra
			dec=struct[indices[w]].dec
			chunk=1
			ds=self->bosstile_density(ra,dec,chunk, tile_subset=tile_subset)
			density[i] = median(ds.density)
		endfor

		; interpolate to get our 50
		x2star_best[ireg]=interpol(x2star_vals, density, density_goal) > 3

		print,ireg+1,nregions,minra,maxra,x2star_best[ireg],$
			form='(i0,"/",i0,"  [",f0.2,",",f0.2,"]  best x2star: ",f0.2)'

		pplot, x2star_vals, density, $
			color=colors[ireg],/overplot

	endfor



	pplot, [0,100], [density_goal,density_goal], $
		line=2, /over

	endplot

	psfile=target_run+'-'+dirpart+'-best.eps'
	psfile=path_join(outdir, psfile)

	begplot,psfile,/encap,xsize=8.5,ysize=8, /color

	pplot, middle_ravals, x2star_best, $
		xtitle='RA', $
		ytitle='Best x2_star'

	sra=shiftra(struct.ra,/wrap)

	quadfit=poly_fit(middle_ravals,x2star_best,2)
	yfit = poly(middle_ravals, quadfit)

	pplot, middle_ravals, yfit, color='blue', /over
	
	endplot


	outfile='~/tmp/kde-coadd-chunk1-density-tune-vs-ra.fits'
	print,'Writing to file: ',outfile
	ost=replicate({ra:0d, x2_star_best:0d}, nregions)
	ost={$
		ra:middle_ravals, $
		x2_star_best: x2star_best, $
		quadfit:quadfit[*], $
		yfit: yfit}
	mwrfits, ost, outfile, /create

	print,'kde_coadd_x2_star_ravals:       ',tostring(middle_ravals), form='(a,a)'
	print,'kde_coadd_x2_star_best: ',tostring(x2star_best), form='(a,a)'
	print,'quadfit: ',quadfit

end


pro esboss::tune_chunk2_qso_density_bonus_by_ra, target_run, $
		nvals=nvals, regenerate=regenerate,  $
		outdir=outdir

	; tune bonus objects so that we meet the *overall* density goal 
	; in each region

	; the density we are aiming for in each region
	density_goal = 80.

	chunk=2




	; number of values in x2_star to check
	if n_elements(nvals) eq 0 then nvals=20

	bt=obj_new('bosstarget')
	struct=bt->read_collated('qso', target_run)
	pstr=extract_tags(struct, ['ra','dec'])


	bq=obj_new('bosstarget_qso')
	struct.boss_target1=bq->select(struct, /comm2, /reselect)

	types = [$
		'qso_core','qso_bonus','qso_known_midz','qso_known_lohiz',$
		'qso_nn','qso_ukidss','qso_kde_coadd','qso_like','qso_first_boss']

	; make sure you aren't selecting on x2_star using anything fancy!
	wtarg=self->btselect(struct, types, /fixup)
	bonusflag=sdss_flagval('boss_target1','qso_bonus')

	; range of x2star to explore
	minx2star=3
	maxx2star=16.
	;maxx2star=10.

	; values of x2star to check
	x2star_vals = arrscl(findgen(nvals), minx2star, maxx2star)
	print,'x2star_vals: ',x2star_vals



	; loop over the regions and get the x2star that gives us the required
	; density goal

	bt=obj_new('bosstarget')
	dirpart='tune-chunk2-density-vs-ra'
	if n_elements(outdir) eq 0 then begin
		outdir = bt->target_dir(target_run)
		outdir = path_join(outdir, dirpart)
	endif
	file_mkdir, outdir
	psfile=target_run+'-'+dirpart+'-each.eps'
	psfile=path_join(outdir, psfile)




	; match positions to the tiles
	self->bosstile_match_poly, $
		pstr[wtarg].ra, pstr[wtarg].dec, chunk, $
		allmatch_obj, allmatch_tiles, $
		matched_objects, matched_tiles, $
		tiles=alltiles

	wtarg = wtarg[matched_objects]

	; sort the tiles by ra
	;eq2csurvey, alltiles.ra, alltiles.dec, clam,ceta
	;s=sort(clam[matched_tiles])

	s=sort(alltiles[matched_tiles].ra)

	; now this should be a sorted list back into the entire tile list
	matched_tiles = matched_tiles[s]

	; now many tiles matched?
	nmatch=n_elements(matched_tiles)
	;x2star_best = fltarr(nmatch)

	tile_window=3


	begplot,psfile,/encap,xsize=8.5,ysize=8, /color
	colors=make_rainbow(nmatch)
	pplot, [0],/nodata,xrange=[minx2star,maxx2star], yrange=[30,100], $
		xtitle='x2star', $
		ytitle='#/square degree'
	
	for ireg=0L, nmatch-1 do begin

		regwindow = [ireg-tile_window, ireg+tile_window] > 0 < (nmatch-1)
		;print,'regwindow = ',regwindow

		tiles2use = matched_tiles[regwindow[0]:regwindow[1]]
		;print,tiles2use

		;if n_elements(tiles2use) ge 2*tile_window then begin
		if n_elements(tiles2use) gt 0 then begin

			tiles_in_reg = alltiles[tiles2use]
			midra=mean(tiles_in_reg.ra)

			add_arrval, midra, middle_ravals

			density=fltarr(nvals)
			for i=0L, nvals-1 do begin

				x2star=x2star_vals[i]

				w=where(struct[wtarg].x2_star gt x2star)

				bonusonly=1
				if bonusonly then begin
					; this where selects objects outside of the region of 
					; interest, we will take the complement. Note we will 
					; remove bonus-only objects only.
					wbad=where($
						(struct[wtarg].boss_target1 eq bonusflag) $
						and struct[wtarg].x2_star lt x2star, $
						comp=keep, ncomp=nkeep)
				endif else begin
					wbad=where($
						(struct[wtarg].boss_target1 and bonusflag) ne 0 $
						and struct[wtarg].x2_star lt x2star, $
						comp=keep, ncomp=nkeep)
				endelse

				if nkeep gt 0 then begin
					keep=wtarg[keep]
					dstruct = $
						self->bosstile_density(pstr[keep].ra,pstr[keep].dec, $
						                       chunk, $
											   tile_subset=tiles2use, $
											   wincheck=0)
					;density[i] = median(dstruct.density)
					density[i] = mean(dstruct.density)
					;density[i] = total(dstruct.count)/total(dstruct.area)
					;print,$
					;	density[i],$
					;	total(dstruct.count)/total(dstruct.area), $
					;	median(dstruct.density)
				endif
			endfor

			;print,density
			; interpolate to get our 50
			;x2star_best[ireg]=interpol(x2star_vals, density, density_goal) > 3

			wok = where(density ge density_goal, nok)
			if nok eq 0 then begin
				best = minx2star
			endif else begin
				best=interpol(x2star_vals, density, density_goal)
			endelse
			add_arrval, best, x2star_best

			print,ireg+1,nmatch,midra,best,$
				form='(i0,"/",i0,"  [",f0.2,"]  best x2star: ",f0.3)'

			pplot, x2star_vals, density, $
				color=colors[ireg],/overplot
		endif
	endfor



	pplot, [0,100], [density_goal,density_goal], $
		line=2, /over

	endplot

	psfile=target_run+'-'+dirpart+'-best.eps'
	psfile=path_join(outdir, psfile)

	begplot,psfile,/encap,xsize=8.5,ysize=8, /color
	!p.charsize=2

	pplot, middle_ravals, x2star_best, $
		xtitle='RA', $
		ytitle='Best x2_star'

	sra=shiftra(struct.ra,/wrap)

	order=2
	fitpars2=poly_fit(middle_ravals,x2star_best,2)
	yfit2 = poly(middle_ravals, fitpars2)


	order=3
	fitpars3=poly_fit(middle_ravals,x2star_best,3)
	yfit3 = poly(middle_ravals, fitpars3)

	order=4
	fitpars4=poly_fit(middle_ravals,x2star_best,4)
	yfit4 = poly(middle_ravals, fitpars4)

	pplot, middle_ravals, yfit2, color='blue', /over
	pplot, middle_ravals, yfit3, color='red', /over
	pplot, middle_ravals, yfit4, color='darkgreen', /over
	plegend,['order=2','order=3','order=4'],color=['blue','red','darkgreen'], $
		line=0, $
		/right
	
	endplot


	outfile='tune-chunk2-bonus-density-vs-ra.fits'
	outfile=path_join(outdir, outfile)
	print,'Writing to file: ',outfile
	ost=replicate({ra:0d, x2_star_best:0d}, nmatch)
	ost={$
		ra:middle_ravals, $
		x2_star_best: x2star_best, $
		$
		fitpars2:fitpars2[*], $
		yfit2: yfit2, $
		$
		fitpars3:fitpars3[*], $
		yfit3: yfit3, $
		$
		fitpars4:fitpars4[*], $
		yfit4: yfit4 $
		}
	mwrfits, ost, outfile, /create

	;print,'middle_ravals:       ',tostring(middle_ravals), form='(a,a)'
	;print,'x2star_best: ',tostring(x2star_best), form='(a,a)'
	;print,'fitpars2: ',tostring(fitpars[*])

end



pro esboss::tune_chunk1_tile_density_qso_kde_coadd_region, $
		target_run, region, nvals=nvals, regenerate=regenerate

	qso_type='qso_kde_coadd'
	extra_cutname=region


	bt=obj_new('bosstarget')
	bq=obj_new('bosstarget_qso')

	if 0 then begin
		struct=bt->read_collated('qso',target_run)
		kde=bq->kde_coadd_read()

		kcflag=sdss_flagval('boss_target1', 'qso_kde_coadd')
		wkde_match=where(struct.kde_coadd_id ge 0 and $
			(struct.boss_target1 and kcflag) ne 0)
		help,wkde_match
		struct=struct[wkde_match]
		kde = kde[struct.kde_coadd_id]
	endif

	epars={ $
		kde_coadd_x2star_region1:2.0, $; 7.42
		kde_coadd_x2star_region2:2.0, $; 3.0
		kde_coadd_x2star_region3:2.0 $ ; 9.07
		}
	self->get_kde_coadd_matches, target_run, struct, kde, tile_subset, $
		epars=epars, region=region



	wreg=self->_select_by_name(kde, region, tile_subset=tile_subset)
	help,wreg
	struct=struct[wreg]
	kde=kde[wreg]



	extra_name = region

	if n_elements(nvals) eq 0 then nvals=20

	case region of
		'region1': begin
			minx2star=2.
			maxx2star=12.0
		end
		'region2': begin
			minx2star=2.0
			maxx2star=8.0
		end
		'region3': begin
			minx2star=2.
			maxx2star=12.0
		end
		else: message,'bad region: '+region
	endcase


	x2star_vals = $
		arrscl(findgen(nvals), minx2star, maxx2star)

	density=fltarr(nvals)
	for i=0L, nvals-1 do begin

		x2star=x2star_vals[i]

		w=where(kde.x2_star gt x2star)

		chunk=1
		ds=self->bosstile_density(kde[w].ra,kde[w].dec, chunk, $
			tile_subset=tile_subset)
		density[i] = median(ds.density)
	endfor

	bt=obj_new('bosstarget')
	outdir = bt->target_dir(target_run)
	outdir = path_join(outdir, 'density')

	file=self->qso_tile_density_file(target_run,qso_type,none,extra_name)
	short_type='qc'
	psfile=repstr(file,'.fits','-'+short_type+'-x2starcheck.eps')

	begplot,psfile,/encap,xsize=8.5,ysize=8, /color
	

	pplot, x2star_vals, density, $
		xtitle='x2star', $
		ytitle='#/square degree'


	goal=50
	pplot, [0,100], [goal,goal], color='darkgreen', line=2, /over

	xval=interpol(x2star_vals, density, goal)
	pplot, [xval], [goal], psym=7, color='red', /over


	legend, [region], /left
	legend,[textoidl('x2star \sim '+string(xval,f='(f0.2)'))], /right,/bottom
	endplot


end

pro esboss::contour_chunk1_tile_density_qso_region, $
		region, qso_type, nvals=nvals, regenerate=regenerate


	target_run = '2009-08-02c'
	;qso_type='qso_core'
	extra_name='looseqso'
	extra_cutname=region

	_extra_name = [extra_name, extra_cutname]

	if n_elements(nvals) eq 0 then nvals=20

	case region of
		'region1': begin
			minx2star=8
			maxx2star=11.0
		end
		'region2': begin
			minx2star=5.0
			maxx2star=8.0
		end
		'region3': begin
			minx2star=7.0
			maxx2star=10.0
		end
		else: message,'bad region: '+region
	endcase



	x2star_restrictive_region_vals = $
		arrscl(findgen(nvals), minx2star, maxx2star)

	density=fltarr(nvals)
	bfratio=density
	for i=0L, nvals-1 do begin

		val=x2star_restrictive_region_vals[i]

		case region of
			'region1': begin
				pars={x2star_restrictive_region1:val, gsplit:21.0}
			end
			'region2': begin
				pars={x2star_restrictive_region2:val, gsplit:21.0}
			end
			'region3': begin
				pars={x2star_restrictive_region3:val, gsplit:21.0}
			end
			else: message,'bad region: '+region
		endcase


		fitsfile=$
			self->qso_tile_density_file(target_run,qso_type,pars,_extra_name)

		if not file_test(fitsfile) or keyword_set(regenerate) then begin
			self->calculate_chunk1_tile_density_qso, $
				target_run, qso_type, pars, extra_name, $
				struct=struct, $
				extra_cutname=extra_cutname,/comm
		endif


		splog,'file: ',fitsfile
		tall=mrdfits(fitsfile,1)
		tbright=mrdfits(fitsfile,2)
		tfaint=mrdfits(fitsfile,3)

		all_dens = median(tall.density)
		bright_dens = median(tbright.density)
		faint_dens = median(tfaint.density)
		if faint_dens gt 0 then begin
			bfratio[i] = bright_dens/faint_dens
		endif


		density[i] = all_dens
	endfor

	bt=obj_new('bosstarget')
	outdir = bt->target_dir(target_run)
	outdir = path_join(outdir, 'density')

	file=self->qso_tile_density_file(target_run,qso_type,none,_extra_name)
	short_type='qc'
	psfile=repstr(file,'.fits','-'+short_type+'-contours.eps')

	begplot,psfile,/encap,xsize=8.5,ysize=8, /color
	

	pplot, x2star_restrictive_region_vals, density, $
		xtitle='x2star_restrictive (ra-hidens)', $
		ytitle='#/square degree'

	bfcolor='blue'
	pplot, x2star_restrictive_region_vals, bfratio, color=bfcolor,$
		/overplot

	goal=20
	pplot, [0,100], [goal,goal], color='darkgreen', line=2, /over

	xval=interpol(x2star_restrictive_region_vals, density, goal)
	pplot, [xval], [goal], psym=7, color='red', /over


	legend, [region], /left
	legend, ['density','bright/faint ratio'],color=[!p.color,c2i(bfcolor)], $
		/right, line=0
	legend,[textoidl('x2star \sim '+string(xval,f='(f0.2)'))], /right,/bottom
	endplot
end

pro esboss::calculate_chunk1_tile_density_qso, $
		target_run, qso_types, pars, extra_name, $
		extra_cutname=extra_cutname, $
		struct=struct, commissioning=commissioning

	if n_params() lt 4 then begin
		on_error, 2
		message,'usage: eb->calculate_chunk1_tile_density_qso, target_run, qso_type, extra_name, struct=, /usewhere, commissioning='
	endif


	bt=obj_new('bosstarget')
	bq=obj_new('bosstarget_qso')

	chunk=1

	if n_elements(extra_cutname) ne 0 then begin
		_extra_name = [extra_name, extra_cutname]
	endif else begin
		_extra_name = extra_name
	endelse
	fitsfile=$
		self->qso_tile_density_file(target_run, qso_types, pars, _extra_name)
	splog,'Will write to file: ',fitsfile,format='(a,a)'


	if n_tags(struct) eq 0 then begin
		struct=bt->read_collated('qso', target_run, extra_name=extra_name)
	endif

	lups = bq->get_lups(struct, /deredden)
	w=self->qso_select_bypar(struct, qso_types, pars, $
		count=count, usewhere=usewhere, res=res, $
		commissioning=commissioning)

	if count gt 0 and n_elements(extra_cutname) ne 0 then begin
		w2=self->_select_by_name(struct[w], extra_cutname, count=count, $
			tile_subset=tile_subset)
		if count gt 0 then w=w[w2]
	endif
	;endif
	if count gt 0 then begin
		ra=struct[w].ra
		dec=struct[w].dec
		density_struct=self->bosstile_density(ra,dec,chunk, $
			tile_subset=tile_subset)
		splog,'Writing to file: ',fitsfile,format='(a,a)'
		mwrfits, density_struct, fitsfile, /create

		wb = where( lups[1,w] lt pars.gsplit, nb)
		wf = where( lups[1,w] ge pars.gsplit, nf)
		if nb ne 0 or nf ne 0 then begin
			if nb ne 0 then begin
				wb=w[wb]
				splog,'Found ',nb,' bright',format='(a,i0,a)'
				ra=struct[wb].ra
				dec=struct[wb].dec
				density_struct=self->bosstile_density(ra,dec, chunk,$
					tile_subset=tile_subset)
			endif else begin
				struct_assign, {tmp:0}, density_struct
			endelse
			splog,'Appending bright to file: ',fitsfile,format='(a,a)'
			mwrfits, density_struct, fitsfile

			if nf ne 0 then begin
				wf=w[wf]
				splog,'Found ',nf,' faint',format='(a,i0,a)'
				ra=struct[wf].ra
				dec=struct[wf].dec
				density_struct=self->bosstile_density(ra,dec,chunk,$
					tile_subset=tile_subset)
			endif else begin
				struct_assign, {tmp:0}, density_struct
			endelse
			splog,'Appending faint to file: ',fitsfile,format='(a,a)'
			mwrfits, density_struct, fitsfile

		endif

	endif else begin
		splog,'No objects passed cuts'
	endelse


end


pro esboss::calculate_chunk1_tile_density_qso_pbs, $
		target_run, qso_type, extra_name, $
		bosstarget_v=bosstarget_v, $
		photoop_v=photoop_v, $
		idlutils_v=idlutils_v, $
		photo_calib=photo_calib, $
		photo_resolve=photo_resolve, $
		photo_sweep=photo_sweep, $
		commissioning=commissioning, $
		extra_cutname=extra_cutname, $
		walltime=walltime

	if n_params() lt 3 then begin
		on_error, 2
		message,'usage: eb->calculate_sgc_tile_density_qso_pbs, target_run, qso_type, extra_name, walltime=, extra_cutname='
	endif




	if n_elements(bosstarget_v) eq 0 then begin
		bosstarget_v = "-r /home/esheldon/exports/bosstarget-work"
	endif
	if n_elements(photoop_v) eq 0 then photoop_v = "v1_9_4"
	if n_elements(idlutils_v) eq 0 then begin
		idlutils_v = "-r /home/esheldon/exports/idlutils-work"
	endif

	if n_elements(commissioning) eq 0 then commissioning=0

	;setups = "source ~esheldon/.dotfiles/bash/riemann/idl_setup.sh"
	setups = "setup sas bosswork"
	setups=[setups,"setup photoop "+photoop_v]
	setups=[setups,"setup idlutils "+idlutils_v] 
	setups=[setups,"setup bosstarget "+bosstarget_v] 
	setups=[setups,"setup sdssidl -r /home/esheldon/svn/sdssidl"] 
	setups=[setups,"setup esidl -r /home/esheldon/exports/esidl-work"] 

	if n_elements(photo_sweep) ne 0 then begin
		setups=[setups,'export PHOTO_SWEEP='+photo_sweep]
	endif
	if n_elements(photo_resolve) ne 0 then begin
		setups=[setups,'export PHOTO_RESOLVE='+photo_resolve]
	endif
	if n_elements(photo_calib) ne 0 then begin
		setups=[setups,'export PHOTO_CALIB='+photo_calib]
	endif

	setups = strjoin(setups, ' && ')






	; bonus means permissive or bonus
	select_pars=self->sgc_qso_tile_density_pars()

	pbs_dir='/home/esheldon/pbs/testqsopars'
	pbs_dir=path_join(pbs_dir,target_run)

	file_mkdir, pbs_dir

	if tag_exist(select_pars, 'logsdmax_bright_perm_default') then begin
		sdbp_def = select_pars.logsdmax_bright_perm_default
		sdfp_def = select_pars.logsdmax_faint_perm_default
	endif

	if tag_exist(select_pars, 'logsdmax_bright_restr_default') then begin
		sdbr_def = select_pars.logsdmax_bright_restr_default
		sdfr_def = select_pars.logsdmax_faint_restr_default
	endif

	; for the pars I've set up qso_bonus is really qso_bonus or qso_core
	if in(qso_type,'qso_bonus') then begin
		nbright = select_pars.nbonus_bright
		nfaint = select_pars.nbonus_faint
	endif else if qso_type eq 'qso_core' then begin
		nbright = select_pars.ncore_bright
		nfaint = select_pars.ncore_faint
	endif else  message,'bad qso type'

	qso_type_str = "'"+strjoin(qso_type, "', '")+"'"
	if n_elements(qso_type) gt 1 then begin
		qso_type_str = "["+qso_type_str+"]"
	endif

	for i=0L, nbright-1 do begin
		for j=0L, nfaint-1 do begin

			if in(qso_type,'qso_bonus') and in(qso_type,'qso_core') then begin
				sdbp_off=select_pars.logsdmax_bright_perm_offsets[i]
				sdfp_off=select_pars.logsdmax_faint_perm_offsets[j]

				bval=(sdbp_def + sdbp_off)
				fval=(sdfp_def + sdfp_off)


				pars={logstardensmax_bright_permissive: bval, $
					  logstardensmax_faint_permissive: fval, $
					  gsplit:21.0}

				commands='pars='+tostring(pars)
				extra_commands=''


				short_type='qb'
				nstr = string(bval,'-',fval,f='(f0.3,a,f0.3)')
				name = short_type+'-'+nstr
				pbs_file='core-bonus-'+nstr+'.pbs'
			endif else if qso_type eq 'qso_core' then begin
				sdbr_off=select_pars.logsdmax_bright_restr_offsets[i]
				sdfr_off=select_pars.logsdmax_faint_restr_offsets[j]

				bval=(sdbr_def + sdbr_off)
				fval=(sdfr_def + sdfr_off)

				pars={logstardensmax_bright_restrictive: bval, $
					  logstardensmax_faint_restrictive: fval, $
					  gsplit:21.0}

				commands='pars='+tostring(pars)
				extra_commands=''

				short_type='qc'
				nstr = string(bval,'-',fval,f='(f0.3,a,f0.3)')
				name = short_type+'-'+nstr
				pbs_file='core-'+nstr+'.pbs'
			endif else message,'bad qso type'


			pbs_file=path_join(pbs_dir, pbs_file)

			main_command = $
				"eb->calculate_chunk1_tile_density_qso, target_run, qso_types, pars, extra_name"


			if keyword_set(commissioning) then begin
				pbs_file=repstr(pbs_file, '.pbs', '-comm.pbs')
				main_command += ', /commissioning'
			endif
			if n_elements(extra_cutname) ne 0 then begin
				main_command += ", extra_cutname='"+extra_cutname+"'"
				pbs_file=repstr(pbs_file, '.pbs', '-'+extra_cutname+'.pbs')
			endif

			print,pbs_file
			idl_commands = [ $
				"eb=obj_new('esboss')", $
				"bq=obj_new('bosstarget_qso')",$
				"target_run='"+target_run+"'", $
				commands,$
				extra_commands, $
				"qso_types="+qso_type_str, $
				"extra_name='"+extra_name+"'", $
				main_command]

			pbs_riemann_idl, pbs_file, idl_commands, job_name=name, $
				setup=setups, $
				walltime=walltime

		endfor
	endfor


end


pro esboss::contour_chunk1_tile_density_qso, $
		target_run, qso_type, extra_name, $
		extra_cutname=extra_cutname, $
		commissioning=commissioning, $
		density=density, bfratio=bfratio, generate=generate, struct=struct


	if n_params() lt 3 then begin
		on_error, 2
		message,'usage: eb->contour_sgc_tile_density_qso, target_run, qso_type, extra_name, struct=, /usewhere'
	endif
	
	; bonus means permissive or bonus
	select_pars=self->sgc_qso_tile_density_pars()

	if tag_exist(select_pars, 'logsdmax_bright_perm_default') then begin
		sdbp_def = select_pars.logsdmax_bright_perm_default
		sdfp_def = select_pars.logsdmax_faint_perm_default
	endif

	if tag_exist(select_pars, 'logsdmax_bright_restr_default') then begin
		sdbr_def = select_pars.logsdmax_bright_restr_default
		sdfr_def = select_pars.logsdmax_faint_restr_default
	endif


	if in(qso_type, 'qso_bonus') then begin
		nbright = select_pars.nbonus_bright
		nfaint = select_pars.nbonus_faint
		short_type='qb'
		pars={logstardensmax_bright_permissive:0.0, $
			  logstardensmax_faint_permissive:0.0,$
			  gsplit:21.0}
	endif else if qso_type eq 'qso_core' then begin
		nbright = select_pars.ncore_bright
		nfaint = select_pars.ncore_faint
		short_type='qc'

		pars={logstardensmax_bright_restrictive:0.0, $
			  logstardensmax_faint_restrictive:0.0,$
			  gsplit:21.0}
	endif else  message,'bad qso type'
	bq=obj_new('bosstarget_qso')
	;pars=bq->pars()


	if n_elements(extra_cutname) ne 0 then begin
		_extra_name = [extra_name, extra_cutname]
	endif else begin
		_extra_name = extra_name
	endelse


	if n_elements(density) eq 0 or n_elements(bfratio) eq 0 then begin
		ntot=nbright*nfaint
		itot=0

		density=dblarr(nbright,nfaint)
		bfratio = density
		for i=0L, nbright-1 do begin
			for j=0L, nfaint-1 do begin
				itot+=1

				if in(qso_type,'qso_bonus') then begin
					sdbp_off=select_pars.logsdmax_bright_perm_offsets[i]
					sdfp_off=select_pars.logsdmax_faint_perm_offsets[j]

					bval=(sdbp_def + sdbp_off)
					fval=(sdfp_def + sdfp_off)

					pars.logstardensmax_bright_permissive = bval
					pars.logstardensmax_faint_permissive = fval

					nstr = string(bval,'-',fval,f='(f0.3,a,f0.3)')
					name = short_type+'-'+nstr
				endif else if qso_type eq 'qso_core' then begin
					sdbr_off=select_pars.logsdmax_bright_restr_offsets[i]
					sdfr_off=select_pars.logsdmax_faint_restr_offsets[j]

					bval=(sdbr_def + sdbr_off)
					fval=(sdfr_def + sdfr_off)

					pars.logstardensmax_bright_restrictive = bval
					pars.logstardensmax_faint_restrictive = fval

					nstr = string(bval,'-',fval,f='(f0.3,a,f0.3)')
					name = short_type+'-'+nstr
				endif else message,'bad qso type'

				fitsfile=$
					self->qso_tile_density_file(target_run,qso_type,pars,_extra_name)
				if keyword_set(generate) and not file_test(fitsfile) then begin
					self->calculate_chunk1_tile_density_qso, $
						target_run, qso_type, pars, extra_name, $
						extra_cutname=extra_cutname, $
						commissioning=commissioning, $
						struct=struct
				endif
				print
				print,'Reading ',itot,'/',ntot,': ',fitsfile, $
					f='(a,i0,a,i0,a,a)'
				tall=mrdfits(fitsfile,1)
				tbright=mrdfits(fitsfile,2)
				tfaint=mrdfits(fitsfile,3)

				all_dens = median(tall.density)
				bright_dens = median(tbright.density)
				faint_dens = median(tfaint.density)
				if faint_dens gt 0 then begin
					bfratio[i,j] = bright_dens/faint_dens
				endif

				density[i,j] = all_dens

				print,all_dens, bright_dens, faint_dens

			endfor
		endfor
	endif

	bt=obj_new('bosstarget')
	outdir = bt->target_dir(target_run)
	outdir = path_join(outdir, 'density')

	file=self->qso_tile_density_file(target_run,qso_type,none,_extra_name)
	psfile=repstr(file,'.fits','-'+short_type+'-contours-2par.eps')

	begplot,psfile,/encap,xsize=8.5,ysize=8, /color
	



	if in(qso_type, 'qso_bonus') then begin
		sdbp_off=select_pars.logsdmax_bright_perm_offsets
		sdfp_off=select_pars.logsdmax_faint_perm_offsets
		bval=(sdbp_def + sdbp_off)
		logstardensmax_bright_permissive = bval
		fval=(sdfp_def + sdfp_off)
		logstardensmax_faint_permissive = fval

		levels = 30 + 5*lindgen(30)
		yrange=[-1.2,1.0]
		ratstep=0.2
	endif else if qso_type eq 'qso_core' then begin
		sdbr_off=select_pars.logsdmax_bright_restr_offsets
		sdfr_off=select_pars.logsdmax_faint_restr_offsets

		bval=(sdbr_def + sdbr_off)
		logstardensmax_bright_restrictive = bval
		fval=(sdfr_def + sdfr_off)
		logstardensmax_faint_restrictive = fval

		levels = 1*lindgen(100)
		yrange=[-1.1,-0.5]

		ratstep=0.1
	endif else message,'bad qso type'

	nlevel=n_elements(levels)
	c_label = replicate(1,nlevel)


	xtitle='log(star dens max bright)'
	ytitle='log(star dens max faint)'


	print,yrange
	print,bval
	print
	print,fval
	;nlevel=10
	contour, density, bval, fval, $
		levels=levels, $
		c_label=c_label,$
		xtitle=xtitle, ytitle=ytitle, $
		xstyle=3, $
		yrange=yrange, ystyle=3


	;if in(qso_type, 'qso_bonus') then step=0.10 else step=0.20
	levels = 0.8+ratstep*findgen(40)
	nlevel = n_elements(levels)
	c_label = replicate(1,nlevel)

	ocolor=c2i('blue')
	oline=0
	contour, bfratio, bval, fval, $
		levels=levels, c_label=c_label, /overplot, c_color=ocolor, $
		c_line=oline

	legend,/bottom,$
		['qso',target_run,strjoin(qso_type,' or ')]
	legend,/right,/bottom,$
		['density #/sqdeg','bright/faint ratio'], $
		line=[0,oline], $
		color=[!p.color,ocolor]

	if qso_type eq 'qso_core' then begin
		pplot, [-0.63], [-0.65], psym=7, color='red', /overplot
	endif

	endplot


end








pro esboss::calculate_sgc_tile_density_qso, $
		target_run, qso_types, pars, extra_name, $
		extra_cutname=extra_cutname, $
		struct=struct, usewhere=usewhere

	if n_params() lt 4 then begin
		on_error, 2
		message,'usage: eb->calculate_sgc_tile_density_qso, target_run, qso_type, extra_name, struct=, /usewhere'
	endif


	bt=obj_new('bosstarget')
	bq=obj_new('bosstarget_qso')


	if n_elements(extra_cutname) ne 0 then begin
		_extra_name = [extra_name, extra_cutname]
	endif else begin
		_extra_name = extra_name
	endelse
	fitsfile=$
		self->qso_tile_density_file(target_run, qso_types, pars, _extra_name)
	splog,'Will write to file: ',fitsfile,format='(a,a)'


	if n_tags(struct) eq 0 then begin
		struct=bt->read_collated('qso', target_run, extra_name=extra_name)
	endif

	lups = bq->get_lups(struct, /deredden)
	;w18=where(lups[2,*] gt 18, count)
	;if count gt 0 then begin
		;w=self->qso_select_bypar(struct[w18], qso_types, pars, $
		;	count=count, usewhere=usewhere, res=res)
	w=self->qso_select_bypar(struct, qso_types, pars, $
		count=count, usewhere=usewhere, res=res)

	if count gt 0 and n_elements(extra_cutname) ne 0 then begin
		w2=self->_select_by_name(struct[w], extra_cutname, count=count)
		if count gt 0 then w=w[w2]
	endif
	;endif
	if count gt 0 then begin
		density_struct=self->sloane_tile_density(struct[w])
		splog,'Writing to file: ',fitsfile,format='(a,a)'
		mwrfits, density_struct, fitsfile, /create

		wb = where( lups[1,w] lt pars.gsplit, nb)
		wf = where( lups[1,w] ge pars.gsplit, nf)
		if nb ne 0 or nf ne 0 then begin
			if nb ne 0 then begin
				wb=w[wb]
				splog,'Found ',nb,' bright',format='(a,i0,a)'
				density_struct=self->sloane_tile_density(struct[wb])
			endif else begin
				struct_assign, {tmp:0}, density_struct
			endelse
			splog,'Appending bright to file: ',fitsfile,format='(a,a)'
			mwrfits, density_struct, fitsfile

			if nf ne 0 then begin
				wf=w[wf]
				splog,'Found ',nf,' faint',format='(a,i0,a)'
				density_struct=self->sloane_tile_density(struct[wf])
			endif else begin
				struct_assign, {tmp:0}, density_struct
			endelse
			splog,'Appending faint to file: ',fitsfile,format='(a,a)'
			mwrfits, density_struct, fitsfile

		endif

	endif else begin
		splog,'No objects passed cuts'
	endelse


end



pro esboss::calculate_sgc_tile_density_qso_pbs, $
		target_run, qso_type, extra_name, $
		extra_cutname=extra_cutname, $
		walltime=walltime

	if n_params() lt 3 then begin
		on_error, 2
		message,'usage: eb->calculate_sgc_tile_density_qso_pbs, target_run, qso_type, extra_name, walltime=, extra_cutname='
	endif

	; bonus means permissive or bonus
	select_pars=self->sgc_qso_tile_density_pars()

	pbs_dir='/home/esheldon/pbs/testqsopars'
	pbs_dir=path_join(pbs_dir,target_run)

	file_mkdir, pbs_dir

	if tag_exist(select_pars, 'logsdmax_bright_perm_default') then begin
		sdbp_def = select_pars.logsdmax_bright_perm_default
		sdfp_def = select_pars.logsdmax_faint_perm_default
	endif

	if tag_exist(select_pars, 'logsdmax_bright_restr_default') then begin
		sdbr_def = select_pars.logsdmax_bright_restr_default
		sdfr_def = select_pars.logsdmax_faint_restr_default
	endif

	; for the pars I've set up qso_bonus is really qso_bonus or qso_core
	if in(qso_type,'qso_bonus') then begin
		nbright = select_pars.nbonus_bright
		nfaint = select_pars.nbonus_faint
	endif else if qso_type eq 'qso_core' then begin
		nbright = select_pars.ncore_bright
		nfaint = select_pars.ncore_faint
	endif else  message,'bad qso type'

	qso_type_str = "'"+strjoin(qso_type, "', '")+"'"
	if n_elements(qso_type) gt 1 then begin
		qso_type_str = "["+qso_type_str+"]"
	endif

	for i=0L, nbright-1 do begin
		for j=0L, nfaint-1 do begin

			if in(qso_type,'qso_bonus') and in(qso_type,'qso_core') then begin
				sdbp_off=select_pars.logsdmax_bright_perm_offsets[i]
				sdfp_off=select_pars.logsdmax_faint_perm_offsets[j]

				bval=(sdbp_def + sdbp_off)
				fval=(sdfp_def + sdfp_off)


				pars={logstardensmax_bright_permissive: bval, $
					  logstardensmax_faint_permissive: fval, $
					  gsplit:21.0}

				commands='pars='+tostring(pars)
				extra_commands=''


				short_type='qb'
				nstr = string(bval,'-',fval,f='(f0.3,a,f0.3)')
				name = short_type+'-'+nstr
				pbs_file='core-bonus-'+nstr+'.pbs'
			endif else if qso_type eq 'qso_core' then begin
				sdbr_off=select_pars.logsdmax_bright_restr_offsets[i]
				sdfr_off=select_pars.logsdmax_faint_restr_offsets[j]

				bval=(sdbr_def + sdbr_off)
				fval=(sdfr_def + sdfr_off)

				pars={logstardensmax_bright_restrictive: bval, $
					  logstardensmax_faint_restrictive: fval, $
					  gsplit:21.0}

				commands='pars='+tostring(pars)
				extra_commands=''

				short_type='qc'
				nstr = string(bval,'-',fval,f='(f0.3,a,f0.3)')
				name = short_type+'-'+nstr
				pbs_file='core-'+nstr+'.pbs'
			endif else message,'bad qso type'


			pbs_file=path_join(pbs_dir, pbs_file)

			main_command = $
				"eb->calculate_sgc_tile_density_qso, target_run, qso_types, pars, extra_name"
			if n_elements(extra_cutname) ne 0 then begin
				main_command += ", extra_cutname='"+extra_cutname+"'"
				pbs_file=repstr(pbs_file, '.pbs', '-'+extra_cutname+'.pbs')
			endif

			print,pbs_file
			idl_commands = [ $
				"eb=obj_new('esboss')", $
				"bq=obj_new('bosstarget_qso')",$
				"target_run='"+target_run+"'", $
				commands,$
				extra_commands, $
				"qso_types="+qso_type_str, $
				"extra_name='"+extra_name+"'", $
				main_command]

			pbs_riemann_idl, pbs_file, idl_commands, job_name=name, $
				walltime=walltime

		endfor
	endfor


end

pro esboss::contour_sgc_tile_density_qso, target_run, qso_type, extra_name, $
		density=density, bfratio=bfratio, generate=generate, struct=struct, $
		kludge_name=kludge_name, usewhere=usewhere


	if n_params() lt 3 then begin
		on_error, 2
		message,'usage: eb->contour_sgc_tile_density_qso, target_run, qso_type, extra_name, struct=, /usewhere'
	endif
	
	; bonus means permissive or bonus
	select_pars=self->sgc_qso_tile_density_pars()

	if tag_exist(select_pars, 'logsdmax_bright_perm_default') then begin
		sdbp_def = select_pars.logsdmax_bright_perm_default
		sdfp_def = select_pars.logsdmax_faint_perm_default
	endif

	if tag_exist(select_pars, 'logsdmax_bright_restr_default') then begin
		sdbr_def = select_pars.logsdmax_bright_restr_default
		sdfr_def = select_pars.logsdmax_faint_restr_default
	endif


	if in(qso_type, 'qso_bonus') then begin
		nbright = select_pars.nbonus_bright
		nfaint = select_pars.nbonus_faint
		short_type='qb'
	endif else if qso_type eq 'qso_core' then begin
		nbright = select_pars.ncore_bright
		nfaint = select_pars.ncore_faint
		short_type='qc'
	endif else  message,'bad qso type'
	bq=obj_new('bosstarget_qso')
	pars=bq->pars()

	if n_elements(density) eq 0 or n_elements(bfratio) eq 0 then begin
		ntot=nbright*nfaint
		itot=0

		density=dblarr(nbright,nfaint)
		bfratio = density
		for i=0L, nbright-1 do begin
			for j=0L, nfaint-1 do begin
				itot+=1

				if in(qso_type,'qso_bonus') then begin
					sdbp_off=select_pars.logsdmax_bright_perm_offsets[i]
					sdfp_off=select_pars.logsdmax_faint_perm_offsets[j]

					bval=(sdbp_def + sdbp_off)
					fval=(sdfp_def + sdfp_off)

					pars.logstardensmax_bright_permissive = bval
					pars.logstardensmax_faint_permissive = fval

					nstr = string(bval,'-',fval,f='(f0.3,a,f0.3)')
					name = short_type+'-'+nstr
				endif else if qso_type eq 'qso_core' then begin
					sdbr_off=select_pars.logsdmax_bright_restr_offsets[i]
					sdfr_off=select_pars.logsdmax_faint_restr_offsets[j]

					bval=(sdbr_def + sdbr_off)
					fval=(sdfr_def + sdfr_off)

					pars.logstardensmax_bright_restrictive = bval
					pars.logstardensmax_faint_restrictive = fval

					nstr = string(bval,'-',fval,f='(f0.3,a,f0.3)')
					name = short_type+'-'+nstr
				endif else message,'bad qso type'

				fitsfile=$
					self->qso_tile_density_file(target_run,qso_type,pars,extra_name)
				if keyword_set(generate) and not file_test(fitsfile) then begin
					self->calculate_sgc_tile_density_qso, $
						target_run, qso_type, pars, extra_name, struct=struct, $
						usewhere=usewhere
				endif
				print
				print,'Reading ',itot,'/',ntot,': ',fitsfile, $
					f='(a,i0,a,i0,a,a)'
				t=mrdfits(fitsfile,1)
				tbright=mrdfits(fitsfile,2)
				tfaint=mrdfits(fitsfile,3)

				if tfaint.median gt 0 then begin
					bfratio[i,j] = tbright.median/tfaint.median
				endif

				density[i,j] = t.median

				print,t.median,tbright.median,tfaint.median

			endfor
		endfor
	endif

	bt=obj_new('bosstarget')
	outdir = bt->target_dir(target_run)
	outdir = path_join(outdir, 'density')

	file=self->qso_tile_density_file(target_run,qso_type,none,extra_name)
	psfile=repstr(file,'.fits','-'+short_type+'-contours.eps')

	begplot,psfile,/encap,xsize=8.5,ysize=8, /color
	;levels = 5+5*lindgen(30)
	;levels = 5+lindgen(200)
	



	if in(qso_type, 'qso_bonus') then begin
		sdbp_off=select_pars.logsdmax_bright_perm_offsets
		sdfp_off=select_pars.logsdmax_faint_perm_offsets
		bval=(sdbp_def + sdbp_off)
		logstardensmax_bright_permissive = bval
		fval=(sdfp_def + sdfp_off)
		logstardensmax_faint_permissive = fval

		levels = 30 + 5*lindgen(30)
	endif else if qso_type eq 'qso_core' then begin
		sdbr_off=select_pars.logsdmax_bright_restr_offsets
		sdfr_off=select_pars.logsdmax_faint_restr_offsets

		bval=(sdbr_def + sdbr_off)
		logstardensmax_bright_restrictive = bval
		fval=(sdfr_def + sdfr_off)
		logstardensmax_faint_restrictive = fval

		levels = 2*lindgen(100)
	endif else message,'bad qso type'

	nlevel=n_elements(levels)
	c_label = replicate(1,nlevel)


	xtitle='log(star dens max bright)'
	ytitle='log(star dens max faint)'


	;nlevel=10
	contour, density, bval, fval, $
		levels=levels, $
		c_label=c_label,$
		xtitle=xtitle, ytitle=ytitle, $
		xstyle=3, $
		yrange=[-1.2,1.0], ystyle=3


	;if in(qso_type, 'qso_bonus') then step=0.10 else step=0.20
	step=0.20
	levels = 0.8+step*findgen(40)
	nlevel = n_elements(levels)
	c_label = replicate(1,nlevel)

	ocolor=c2i('blue')
	oline=0
	contour, bfratio, bval, fval, $
		levels=levels, c_label=c_label, /overplot, c_color=ocolor, $
		c_line=oline

	legend,/bottom,$
		['qso',target_run,strjoin(qso_type,' or ')]
	legend,/right,/bottom,$
		['density #/sqdeg','bright/faint ratio'], $
		line=[0,oline], $
		color=[!p.color,ocolor]

	endplot


end



pro esboss::calculate_sgc_tile_density_dperp0_ihi_pbs, target_run

	target_type='lrg'
	naxis=10
	dperp0=arrscl(findgen(naxis), 0.45, 0.55)
	ihi=arrscl(findgen(naxis), 19.9, 20.0)

	tmpdir='/home/esheldon/pbs/test_commiss'

	;outdir=path_join(tmpdir,'data')
	pbs_dir=path_join(tmpdir,'pbs')

	;file_mkdir, outdir
	file_mkdir, pbs_dir

	job=0L
	for i=0L, naxis-1 do begin
		for j=0L, naxis-1 do begin

			dperp0_str = string(dperp0[i], f='(f0.2)')
			ihi_str = string(ihi[j], f='(f0.2)')

			pbs_file=path_join(pbs_dir, 'grid-'+dperp0_str+'-'+ihi_str+'.pbs')
			name = target_type+"-"+dperp0_str+"-"+ihi_str
			idl_commands = [ $
				"eb=obj_new('esboss')", $
				"target_type='"+target_type+"'", $
				"target_run='"+target_run+"'", $
				"btflags=['gal_loz','gal_cmass']", $
				"dperp0="+dperp0_str, $
				"ihi="+ihi_str, $
				"eb->calculate_sgc_tile_density, target_type, target_run, btflags, dperp0=dperp0, ihi=ihi"]

			pbs_riemann_idl, pbs_file, idl_commands, job_name=name

		endfor
	endfor

end

pro esboss::calculate_sgc_tile_density_ihi_rlo_pbs, target_run, $
		noclobber=noclobber

	pars=self->sgc_lrg_tile_density_pars()

	target_type='lrg'
	dperp0=0.55
	dperp0_str=string(dperp0,f='(f0.3)')
	trifrac=0.0
	tri_str=string(trifrac,f='(f0.3)')

	pbs_dir='/home/esheldon/pbs/test_rlo_ihi'

	file_mkdir, pbs_dir

	for i=0L, pars.nihi-1 do begin
		ihi = pars.ihi[i]
		ihi_str = string(ihi, f='(f0.3)')

		for j=0L, pars.nrlo-1 do begin
			rlo = pars.rlo[j]
			rlo_str = string(rlo, f='(f0.3)')

			if keyword_set(noclobber) then begin
				fitsfile=self->tile_density_file(target_type,target_run,$
					['gal_loz','gal_cmass','gal_triangle'], $
					ihi=ihi, rlo=rlo, triangle_frac=trifrac, dperp0=dperp0)
				if file_test(fitsfile) then begin
					print,'Skipping because fitsfile exists:',fitsfile
				endif
				continue
			endif

			pbs_file=path_join(pbs_dir, 'grid-'+rlo_str+'-'+ihi_str+'.pbs')
			print,pbs_file
			name = target_type+"-"+rlo_str+"-"+ihi_str
			idl_commands = [ $
				"eb=obj_new('esboss')", $
				"target_type='"+target_type+"'", $
				"target_run='"+target_run+"'", $
				"btflags=['gal_loz','gal_cmass','gal_triangle']", $
				"ihi="+ihi_str, $
				"rlo="+rlo_str, $
				"trifrac="+tri_str, $
				"dperp0="+dperp0_str, $
				"eb->calculate_sgc_tile_density, target_type, target_run, btflags, ihi=ihi, rlo=rlo, triangle_frac=trifrac,dperp0=dperp0"]

			pbs_riemann_idl, pbs_file, idl_commands, job_name=name

		endfor
	endfor

end


pro esboss::calculate_sgc_tile_density_triangle_pbs, target_run, $
		noclobber=noclobber

	target_type='lrg'
	nihi=21
	ntri=21

	ihi_arr=arrscl(findgen(nihi), 19.9, 20.0)
	triangle_frac_arr=arrscl(findgen(ntri), 0.0 ,1.0)
	dperp0=0.55
	dperp0_str=string(dperp0,f='(f0.3)')

	pbs_dir='/home/esheldon/pbs/test_triangle'

	file_mkdir, pbs_dir

	for i=0L, nihi-1 do begin
		ihi = ihi_arr[i]
		ihi_str = string(ihi, f='(f0.3)')

		for j=0L, ntri-1 do begin
			trifrac = triangle_frac_arr[j]
			tri_str = string(trifrac, f='(f0.3)')

			if keyword_set(noclobber) then begin
				fitsfile=self->tile_density_file(target_type,target_run,$
					['gal_loz','gal_cmass','gal_triangle'], $
					ihi=ihi, triangle_frac=trifrac, dperp0=dperp0)
				if file_test(fitsfile) then begin
					print,'Skipping because fitsfile exists:',fitsfile
				endif
				continue
			endif

			pbs_file=path_join(pbs_dir, 'grid-'+tri_str+'-'+ihi_str+'.pbs')
			print,pbs_file
			name = target_type+"-"+tri_str+"-"+ihi_str
			idl_commands = [ $
				"eb=obj_new('esboss')", $
				"target_type='"+target_type+"'", $
				"target_run='"+target_run+"'", $
				"btflags=['gal_loz','gal_cmass','gal_triangle']", $
				"ihi="+ihi_str, $
				"trifrac="+tri_str, $
				"dperp0="+dperp0_str, $
				"eb->calculate_sgc_tile_density, target_type, target_run, btflags, ihi=ihi, triangle_frac=trifrac,dperp0=dperp0"]

			pbs_riemann_idl, pbs_file, idl_commands, job_name=name

		endfor
	endfor

end



pro esboss::contour_sgc_tile_densities_ihi_dperp0, target_run, $
		density=density

	pars=self->sgc_lrg_tile_density_pars()

	target_type='lrg'
	btflags=['gal_loz','gal_cmass']

	if n_elements(density) eq 0 then begin
		density = dblarr(naxis,naxis)

		for i=0L, pars.ndperp0-1 do begin
			for j=0L, pars.nihi-1 do begin
				fitsfile=self->tile_density_file($
					target_type, target_run, btflags, $
					dperp0=pars.dperp0[i], ihi=pars.ihi[j])

				print,'reading: ',fitsfile
				t=mrdfits(fitsfile,1)

				density[i,j] = t.median
			endfor
		endfor
	endif

	bt=obj_new('bosstarget')
	outdir = bt->target_dir(target_run)
	outdir = path_join(outdir, 'density')

	file=self->tile_density_file(target_type,target_run,btflags)
	psfile=repstr(file,'density.fits','density-dperp0-ihi-contours.eps')

	begplot,psfile,/encap,xsize=8.5,ysize=8
	levels = 160+5*lindgen(20)
	
	nlevel=n_elements(levels)
	c_label = replicate(1,nlevel)
	;nlevel=10
	contour, density, dperp0, ihi, $
		levels=levels, $
		c_label=c_label,$
		xtitle='dperp0', ytitle='ihi'
	legend,/bottom,$
		[target_type,target_run,strjoin(btflags,' or ')]

	endplot
end

pro esboss::plot_sgc_tile_densities_cmass_ihi, struct=struct, generate=generate

	target_type='lrg'
	target_run ='2009-07-01t'
	target_subtypes=['gal_loz','gal_cmass']

	bt=obj_new('bosstarget')
	bl=obj_new('bosstarget_lrg')

	if n_elements(struct) eq 0 then begin
		struct=bt->read_collated(target_type, target_run)
	endif

	file='~/tmp/cmass-ihi-density-test-'+strjoin(target_subtypes,'-')+'.fits'
	if not file_test(file) or keyword_set(generate) then begin

		splog,'Will write to file: ',file,form='(a,a)'
		ncmass_ihi = 100
		cmass_ihi_min = 20.0
		cmass_ihi_max = 20.5
		cmass_ihi_array = $
			arrscl(findgen(ncmass_ihi), cmass_ihi_min, cmass_ihi_max)

		density=fltarr(ncmass_ihi)

		for i=0L, ncmass_ihi-1 do begin

			cmhstr=string(cmass_ihi_array[i], form='(f0.3)')
			psfile=repstr(file,'.fits','-'+cmhstr+'.eps')

			splog,'ihi_cmass: ',cmhstr
			pars={ihi_cmass: cmass_ihi_array[i]}
			res = bl->select(struct, pars=pars, /commiss, /struct)

			wall=self->btselect(res, ['gal_loz','gal_cmass','gal_triangle'])
			w=self->btselect(res, target_subtypes)


			tstruct=struct[w]
			density_struct=self->sloane_tile_density($
				tstruct, sloane=sloane_struct)
			density[i] = density_struct.median

			splog,'Found density of: ',density[i],form='(a,g0)'


			; plot of the space

			mst=bl->magstruct(struct)
			begplot,psfile,/color,/land
			plotrand, mst[wall].cmodelmag[3], mst[wall].dperp, psym=3, $
				xrange=[17.5,20.5], yrange=[0.50, 1.0], $
				xstyle=1, ystyle=1, $
				xtitle='cmodelmag[3]', ytitle='dperp'
			plotrand, mst[w].cmodelmag[3], mst[w].dperp, psym=3, $
				/over, color='blue'
			legend, string(density[i],form='(f0.1)')+'/sq deg',/right, $
				charsize=1
			endplot,/landfix



		endfor

		density_struct={cmass_ihi: cmass_ihi_array, density:density}
		mwrfits, density_struct, file,/create
	endif else begin
		density_struct=mrdfits(file,1)
	endelse

	psfile=repstr(file,'.fits','.eps')	
	begplot,psfile,/encap,xsize=8.5,ysize=8.5,/color
	pplot, density_struct.cmass_ihi, density_struct.density, $
		xtitle='cmass_ihi', ytitle='#/sq deg', /ynozero, $
		xrange=[20.0,20.35],xstyle=3, aspect=1
	xval=interpol(density_struct.cmass_ihi, density_struct.density, 230d)
	pplot, [xval], [230d], psym=8, color='red',/over
	pplot, [0,xval],[230d,230d],line=2,/over
	pplot, [xval,xval],[0,230d],line=2,/over

	mess=[strjoin(target_subtypes,' or '), $
		textoidl('density \sim 230 at cmass ihi \sim ')+string(xval,form='(f0.2)')]
	legend, mess
	endplot
end

pro esboss::contour_sgc_tile_densities_triangle, target_run, $
		density=density

	target_type='lrg'
	btflags=['gal_loz','gal_cmass','gal_triangle']

	nihi=21
	ntri=21

	ihi=arrscl(findgen(nihi), 19.9, 20.0)
	triangle_frac=arrscl(findgen(ntri), 0.0 ,1.0)
	dperp0=0.55

	if n_elements(density) eq 0 then begin
		density = dblarr(ntri, nihi)

		for i=0L, ntri-1 do begin
			for j=0L, nihi-1 do begin
				fitsfile=self->tile_density_file($
					target_type, target_run, btflags, $
					triangle_frac=triangle_frac[i], ihi=ihi[j], dperp0=dperp0)

				print,'reading: ',fitsfile
				t=mrdfits(fitsfile,1)

				density[i,j] = t.median
			endfor
		endfor
	endif

	bt=obj_new('bosstarget')
	outdir = bt->target_dir(target_run)
	outdir = path_join(outdir, 'density')

	file=self->tile_density_file(target_type,target_run,btflags)
	psfile=repstr(file,'density.fits','density-ihi-trifac-contours.eps')

	begplot,psfile,/encap,xsize=8.5,ysize=8
	levels = 160+5*findgen(50)
	
	nlevel=n_elements(levels)
	c_label = replicate(1,nlevel)
	;nlevel=10
	contour, density, triangle_frac, ihi, $
		levels=levels, $
		c_label=c_label,$
		ytitle='ihi', xtitle='triangle fraction'
	legend,/bottom,$
		[target_run,$
		string('dperp0=',dperp0,f='(a,f0.3)'),$
		strjoin(btflags,' or ')]
	;legend,/right,/bottom,$
	;	textoidl('density [#/deg^2'), line=0

	endplot
end





pro esboss::contour_sgc_tile_densities_ihi_vs_rlo, target_run, $
		density=density

	pars=self->sgc_lrg_tile_density_pars()
	target_type='lrg'
	btflags=['gal_loz','gal_cmass','gal_triangle']


	dperp0=0.55
	triangle_frac=0.0

	if n_elements(density) eq 0 then begin
		density = dblarr(pars.nrlo, pars.nihi)

		for i=0L, pars.nrlo-1 do begin
			for j=0L, pars.nihi-1 do begin
				fitsfile=self->tile_density_file($
					target_type, target_run, btflags, $
					triangle_frac=triangle_frac, $
					dperp0=dperp0, $
					ihi=pars.ihi[j], rlo=pars.rlo[i])

				print,'reading: ',fitsfile
				t=mrdfits(fitsfile,1)

				density[i,j] = t.median
			endfor
		endfor
	endif

	bt=obj_new('bosstarget')
	outdir = bt->target_dir(target_run)
	outdir = path_join(outdir, 'density')

	file=self->tile_density_file(target_type,target_run,btflags)
	psfile=repstr(file,'density.fits','density-ihi-rlo-contours.eps')

	begplot,psfile,/encap,xsize=8.5,ysize=8
	levels = 160+5*findgen(50)
	
	nlevel=n_elements(levels)
	c_label = replicate(1,nlevel)
	;nlevel=10
	contour, density, pars.rlo, pars.ihi, $
		levels=levels, $
		c_label=c_label,$
		ytitle='ihi', xtitle='rlo'
	legend,/bottom,$
		[target_run,$
		string('dperp0=',dperp0,f='(a,f0.3)'),$
		string('trifrac=',triangle_frac,f='(a,f0.3)'), $
		strjoin(btflags,' or ')]

	endplot
end






pro esboss::calculate_sgc_tile_density, target_type, target_run, btflags, $
		struct=struct, fpobjc=fpobjc, noradcheck=noradcheck, $
		triangle_frac=triangle_frac, seed=seed, $
		dperp0=dperp0, ihi=ihi, $
		rlo=rlo, $
		outdir=outdir

	bt=obj_new('bosstarget')
	if n_elements(struct) eq 0 then begin
		struct=bt->read_collated(target_type, target_run, fpobjc=fpobjc)
	endif
	nstruct=n_elements(struct)

	fitsfile=self->tile_density_file(target_type, target_run, btflags, $
		fpobjc=fpobjc, noradcheck=noradcheck, $
		dperp0=dperp0, ihi=ihi, rlo=rlo, triangle_frac=triangle_frac, $
		outdir=outdir)

	; In this case we will re-measure the ts flags
	if n_elements(dperp0) ne 0 $
			or n_elements(ihi) ne 0 or n_elements(rlo) ne 0 then begin
		bl=obj_new('bosstarget_lrg',/commiss)
		lim=bl->limits()
		if n_elements(dperp0) ne 0 then begin
			lim.dperp0=dperp0
		endif
		if n_elements(ihi) ne 0 then begin
			lim.ihi=ihi
		endif
		if n_elements(rlo) ne 0 then begin
			lim.rmaglim[0]=rlo
		endif

		res = bl->select(struct, /struct, pars=lim, /commissioning)

		boss_target1 = res.boss_target1
	endif else begin
		boss_target1 = struct.boss_target1
	endelse

	logic = self->orflags(boss_target1, btflags)

	; this is already applied by ::sgc_gather now, so it's redundant
	if target_type eq 'qso' then begin
		lohiz = sdss_flagval('boss_target1','qso_known_lohiz')
		logic = logic and (boss_target1 and lohiz) eq 0
	endif

	; random sample the triangle at the specified rate.  Will just set 
	; logic to zero for those objects
	if n_elements(triangle_frac) ne 0 then begin
		print,'Random sampling bermuda triangle with rate: ',triangle_frac,$
			f='(a,g0)'
		wtri = where(self->orflags(boss_target1,'gal_triangle'), ntri)
		if ntri eq 0 then message,'No triangle galaxies found'
		if triangle_frac eq 0.0 then begin
			print,'Removing ',ntri,'/',ntri,' from triangle',$
				format='(a,i0,a,i0,a)'
			logic[wtri] = 0
		endif else if triangle_frac eq 1.0 then begin
			print,'Removing ',0,'/',ntri,' from triangle',$
				format='(a,i0,a,i0,a)'
		endif else begin
			; random sample at triangle_frac
			s=sort(randomu(seed,ntri))
			; note using 1-triangle_frac here
			nremove = ntri*(1.0-triangle_frac)
			print,'Removing ',nremove,'/',ntri,' from triangle',$
				format='(a,i0,a,i0,a)'
			not_keep_ind = wtri[s[0:nremove-1]]
			logic[not_keep_ind] = 0
		endelse
	endif

	keep=where(logic)

	; order important since the tiles go over the pole
	tstruct = struct[keep]

	density_struct=self->sloane_tile_density($
		tstruct, sloane=sloane_struct)

	if not file_test(outdir) then file_mkdir, outdir

	print,'Writing to file: ',fitsfile
	mwrfits, density_struct, fitsfile, /create
end

function esboss::tile_density_file, target_type, target_run, btflags, $
		fpobjc=fpobjc, noradcheck=noradcheck, $
		dperp0=dperp0, ihi=ihi, rlo=rlo, $
		triangle_frac=triangle_frac, $
		outdir=outdir

	bt=obj_new('bosstarget')
	tmp_name=bt->target_file(target_type, target_run, /collate, $
		fpobjc=fpobjc)

	tmp_name=file_basename(tmp_name)

	btflags_string = strjoin(btflags, '-')
	extra=repstr(btflags_string,'_','-')

	if keyword_set(noradcheck) then begin
		extra=extra+'-noradcheck'
	endif

	if target_type eq 'lrg' then begin
		if n_elements(dperp0) ne 0 then begin
			extra += '-dperp0-'+string(dperp0,f='(f0.3)')
		endif
		if n_elements(ihi) ne 0 then begin
			extra+='-ihi-'+string(ihi,f='(f0.3)')
		endif
		if n_elements(rlo) ne 0 then begin
			extra+='-rlo-'+string(rlo,f='(f0.3)')
		endif
		if n_elements(triangle_frac) ne 0 then begin
			extra+='-trifrac-'+string(triangle_frac,f='(f0.3)')
		endif
	endif else if target_type eq 'qso' then begin
	endif

	extra += "-density"
	fitsfile = repstr(tmp_name, '-collate.fits', '-'+extra+'.fits')

	if n_elements(outdir) eq 0 then begin
		outdir = bt->target_dir(target_run)
		outdir = path_join(outdir, 'density')
	endif


	;fitsfile = target_run+'-'+target_type'-sgc-tile-density-'+extra+'.fits')
	fitsfile = path_join(outdir, fitsfile)

	return, fitsfile

end


function esboss::sloane_tile_density, struct,  $
		noradcheck=noradcheck, sloane_struct=sloane_struct
	; By default it will only use tiles that do not hit boundary
	; send /noradcheck to override this

	sloane_struct=self->sloane_tiles()

	self->sloane_match, struct, mstruct, mtile

	ntile=sloane_struct.ntile
	density_struct = {$
		noradcheck: keyword_set(noradcheck), $
		keep: lonarr(ntile), $
		count: lon64arr(ntile), $
		density: dblarr(ntile),$
		median:0d, $
		mean:0d, $
		sdev:0d}

	density_struct = create_struct(sloane_struct, density_struct)

	bs=binner(mtile, rev=rev, min=0, max=ntile-1)
	for i=0L, n_elements(bs.hist)-1 do begin

		if rev[i] ne rev[i+1] then begin
			w=rev[ rev[i]:rev[i+1]-1 ]
			nw = n_elements(w)
			density_struct.count[i] = nw
			density_struct.density[i] = nw/sloane_struct.tile_area

		endif
	endfor

	if not keyword_set(noradcheck) then begin
		wsgc = self->boss_bounds_select($
			clambda=density_struct.clambda, ceta=density_struct.ceta, nsgc, $
			radius=density_struct.tile_radius, /complete, area='SGC')
		wsgc_simple = self->boss_bounds_select($
			clambda=density_struct.clambda, ceta=density_struct.ceta, nsimp, $
			radius=density_struct.tile_radius, /complete, area='SGC_SIMPLE')

		indices=[wsgc, wsgc_simple]
		indices=indices[rem_dup(indices)]
	endif else begin
		indices = self->boss_bounds_select($
			clambda=density_struct.clambda, ceta=density_struct.ceta, nsgc, $
			/complete, area='SGC')
	endelse

	density_struct.keep[indices] = 1
	density_struct.median = median(density_struct.density[indices])	
	density_struct.mean = mean(density_struct.density[indices])	
	density_struct.sdev = sdev(density_struct.density[indices])	

	return, density_struct

end




pro esboss::make_run_html, tsstruct, target_type

	html_dir = '~/public_html/bosstarget/fchart/'

	nrun = n_elements(rem_dup(tsstruct.run))
	ncol = n_elements(rem_dup(tsstruct.camcol))
	if nrun gt 1 or ncol gt 1 then message,'Only run/camcol one at a time'

	runstr=run2string(tsstruct[0].run)
	colstr=string(tsstruct[0].camcol,f='(I0)')


	wtarget = self->btselect(tsstruct, target_type, count=ntarget)

	seed = 1134

	html_file = 'fchart-'+runstr+'-'+colstr+'-'+ntostr(seed)+'-'+name+'.html'
	html_file = path_join(html_dir, html_file)

	nrand = 100 < ntarget

	s=sort(randomu(seed, ntarget))

	ind = wtarget[s[0:nrand-1]]


	sdss_fchart_table, tsstruct[ind].ra, tsstruct[ind].dec, html_file, $
		/add_atlas, struct=tsstruct[ind]

end


pro esboss::make_run_html_run, reload=reload

	common plot_run_density_run_blk, t745, t1345

	dir='/mount/early1/bosstarget'

	if n_elements(t745) eq 0 or keyword_set(reload) then begin
		bt=obj_new('bosstarget')
		t745 = bt->read('lrg', 745, 1, target_dir=dir, /collate)
		t1345 = bt->read('lrg', 1345, 1, target_dir=dir, /collate)
	
		obj_destroy, bt
	endif

	self->make_run_html, t1345, ['gal_loz']
	self->make_run_html, t1345, ['gal_cmass']
	self->make_run_html, t1345, ['gal_cmass_notred']

	self->make_run_html, t745, ['gal_loz']
	self->make_run_html, t745, ['gal_cmass']
	self->make_run_html, t745, ['gal_cmass_notred']

end

function esboss::testbed_pars, runid
	case runid of
		'2009-06-08-testbed': begin
			if 0 then begin
				bt=obj_new('bosstarget')
				bt->runlist, runs, reruns
				rl=sdss_runlist(/science)
				match, rl.run, runs, mrl, mfl

				w82=where(rl[mrl].stripe eq 82)
				runs = rl[mrl[w82]].run
				print,runs
			endif
			runs = 1755
			return, {$
				runs: runs, $
				target_types: ['lrg'], $
				process_types: ['fpobjc'] $
			}
		end
		else: message,'Unknown runid: '+string(runid)
	endcase
end

pro esboss::run_sgc_testbed, runid, _extra=_extra

	target_dir=self->sgc_target_dir(runid)

	file_mkdir, target_dir

	pars=self->testbed_pars(runid)

	;extra_logic = self->sgc_extra_logic()

	bt=obj_new('bosstarget')

	for it=0L, n_elements(pars.target_types)-1 do begin
		target_type=pars.target_types[it]

		for pi=0L, n_elements(pars.process_types)-1 do begin
			ptype = pars.process_types[pi]
			if ptype eq 'sweep' then begin
				bt->process_datasweeps, $
					target_type, target_dir=target_dir, $
					runs=pars.runs, extra_logic=extra_logic, $
					_extra=_extra
			endif else if ptype eq 'fpobjc' then begin
				bt->process_fpobjc, $
					target_type, target_dir=target_dir, $
					runs=pars.runs, extra_logic=extra_logic, $
					_extra=_extra
			endif else begin
				message,'process_type must be sweep or fpobjc'
			endelse
		endfor
	endfor

	obj_destroy, bt

end

pro esboss::collate_sgc_testbed, subdir

	target_dir=self->sgc_target_dir(subdir)

	file_mkdir, target_dir

	pars=self->testbed_pars()

	bt=obj_new('bosstarget')

	for i=0L, n_elements(pars.runs)-1 do begin
		run=pars.runs[i]
		for col=1,6 do begin
			for it=0L, n_elements(pars.types)-1 do begin
				type=types[it]
				t=bt->read(type, run, col, /collate, $
					target_dir=target_dir, filename=fname)

				outf = file_basename(fname)
				outf = path_join(outdir, outf)
				outf_ascii = repstr(outf, 'fits', 'st')

				print,'Writing fits file: ',outf
				mwrfits, t, outf, /create
				print,'Writing idlstruct file: ',outf_ascii
				write_idlstruct, t, outf_ascii, /ascii
				print,'----------------------------------------------'
			endfor
		endfor
	endfor
	obj_destroy, bt

end

pro esboss::test_bounds_select, t, radius=radius

	if n_elements(radius) ne 0 then begin
		psfile='~/tmp/sgc-lrg-footprint-complete-radcheck.eps'
	endif else begin
		psfile='~/tmp/sgc-lrg-footprint-complete.eps'
	endelse
	keep=self->boss_bounds_select(struct=t,nkeep,$
		/doplot,yrange=[105,165],ystyle=3, $
		psfile=psfile, $
		/xgrid,/ygrid,xticklen=1,yticklen=1,/complete,$
		/reload,radius=radius, area='SGC')

end


pro esboss::boss_bounds_plot, areaname=areaname
	dir=getenv("BOSSTARGET_DIR")
	if dir eq '' then message,'BOSSTARGET not set up'
	dir = path_join(dir, sub=['geometry','boss_survey.par'])

	bounds = yanny_readone(file,/anon)


	if n_elements(areaname) eq 0 then begin
		areaname='SGC'
	endif

	w=where(bounds.areaname eq areaname, nb)

	padding=5
	yrange = [min(bounds[w].cetaMin)-padding, $
		      max(bounds[w].cetaMax)+padding]
	xrange = [min(bounds[w].clambdaMin)-padding, $
		      max(bounds[w].clambdaMax)+padding]

	plot, [0], /nodata, $
		yrange=yrange, xrange=xrange, $
		xstyle=1, ystyle=1, $
		xtitle=textoidl('\lambda_c'), ytitle=textoidl('\eta_c'), $
		/iso
	for i=0L, nb-1 do begin
		plot_box, $
			bounds[w[i]].clambdaMin, bounds[w[i]].clambdaMax, $
			bounds[w[i]].cetaMin, bounds[w[i]].cetaMax, $
			color=c2i('blue'), thick=1
	endfor


end

function esboss::boss_bounds_select, $
		struct=struct, clambda=clambda, ceta=ceta, $
		nkeep, radius=radius_input, $
		complete=complete, areaname=areaname, $
		reload=reload, $
		doplot=doplot, psfile=psfile, verbose=verbose, $
		_extra=_extra
	
	; if radius sent, then do two passes:
	;	make sure the central points are contained in one of the boxes
	;	for those contained, check the point is radius from the boundary	
	common boss_sgc_area_select, boss_bounds, boss_bounds_complete

	if n_elements(boss_bounds) eq 0 or keyword_set(reload) then begin
		dir=getenv("BOSSTARGET_DIR")
		if dir eq '' then message,'BOSSTARGET not set up'
		dir = path_join(dir, 'geometry')
		file = path_join(dir, 'boss_survey.par')

		boss_bounds = yanny_readone(file,/anon)
	endif
	if n_elements(boss_bounds_complete) eq 0 or keyword_set(reload) then begin
		dir = '~/idl.lib/data'
		file = path_join(dir, 'boss_survey_complete.par')

		boss_bounds_complete = yanny_readone(file, /anon)
	endif

	if keyword_set(complete) then begin
		bounds=boss_bounds_complete
	endif else begin
		bounds=boss_bounds
	endelse

	if n_elements(struct) ne 0 then begin
		npoints = n_elements(struct)

		if tag_exist(struct, 'clambda') then begin
			clambda=struct.clambda
			ceta=struct.ceta
		endif else begin
			eq2csurvey, struct.ra, struct.dec, clambda, ceta
		endelse

	endif else begin
		npoints=n_elements(clambda)
		neta=n_elements(ceta)
		if npoints eq 0 or neta eq 0 then begin
			message,'Enter either struct or clambda/eta'
		endif
	endelse
	nrad = n_elements(radius_input)


	
	if n_elements(areaname) ne 0 then begin
		if keyword_set(verbose) then print,'Using region: ',areaname
		w=where(bounds.areaname eq areaname, nw)
		if nw eq 0 then message,'areaname not found: '+string(areaname)
	endif else begin
		if keyword_set(complete) then begin
			if nrad eq 0 then begin
				if keyword_set(verbose) then print,'Using region: SGC'
				w=where(bounds.areaname eq 'SGC', nw)
			endif else begin
				if keyword_set(verbose) then print,'Using region: SGC_SIMPLE'
				w=where(bounds.areaname eq 'SGC_SIMPLE', nw)
			endelse
		endif else begin
			if keyword_set(verbose) then print,'Using all regions'
			nw=n_elements(bounds)
			w=lindgen(nw)
		endelse
	endelse




	format = $
		'(%"(clambda ge %f and clambda le %f and ceta ge %f and ceta le %f)")'

	region_index = replicate(-9999, npoints)
	; loop over the bounding boxes
	for i=0L, nw-1 do begin
		bound = bounds[w[i]]
		wbound=where( $
			(clambda ge bound.clambdaMin) $
			and (clambda le bound.clambdaMax) $
			and (ceta ge bound.cetaMin) $
			and (ceta le bound.cetaMax), nbound)
		if nbound ne 0 then begin
			region_index[wbound] = w[i]
		endif

		ws = string( $
			bound.clambdaMin, bound.clambdaMax, $
			bound.cetaMin, bound.cetaMax, $
			format=format)
		if keyword_set(verbose) then begin
			print,ws
			print,'In region '+strn(w[i])+'# found: ',nbound
		endif
	endfor


	if nrad gt 0 then begin
		if nrad eq npoints then begin
			radius = radius_input
		endif else if nrad eq 1 then begin
			radius = replicate(radius_input[0], npoints)
		endif else begin
			message,'Radius must be a scalar or same length as points'
		endelse
		; now loop through and apply a bounds radius-from-edge check
		h=histogram(region_index, min=0, max=n_elements(bounds), $
			rev=rev)	
		for i=0L, n_elements(h)-1 do begin
			if rev[i] ne rev[i+1] then begin
				wb=rev[ rev[i]:rev[i+1]-1 ]

				; perform radius bounds check
				bound = bounds[i]
				wgood = where( $
					(clambda[wb] ge (bound.clambdaMin+radius[wb]) ) $
					and (clambda[wb] le (bound.clambdaMax-radius[wb]) ) $
					and (ceta[wb] ge (bound.cetaMin+radius[wb]) ) $
					and (ceta[wb] le (bound.cetaMax-radius[wb]) ), $
					ngood, comp=wbad, ncomp=nbad)

				if nbad ne 0 then begin
					if keyword_set(verbose) then begin
						print,'Bound '+strn(i)+$
							' threw out '+ntostr(nbad)+' for radius check'
					endif
					;print,bound
					region_index[wb[wbad]] = -9999
				endif
			endif
		endfor
	endif

	keep = where(region_index ge 0, nkeep)

	if keyword_set(doplot) and nkeep gt 0 then begin
		if n_elements(psfile) ne 0 then begin
			begplot, psfile, /encap, /color, xsize=11, ysize=8.5
		endif

		plot, clambda, ceta, psym=3, /ynozero, $
			xtitle=textoidl('\lambda_c'), ytitle=textoidl('\eta_c'), $
			_extra=_extra
		oplot, clambda[keep], ceta[keep], psym=3, color=c2i('green')

		for i=0L, nw-1 do begin
			plot_box, $
			bounds[w[i]].clambdaMin, bounds[w[i]].clambdaMax, $
			bounds[w[i]].cetaMin, bounds[w[i]].cetaMax, $
			color=c2i('blue'), thick=1
		endfor

		if n_elements(psfile) ne 0 then begin
			endplot
		endif
	endif


	return, keep
end

pro esboss::plot_qso_sgc_footprints, struct

	sgc_bounds = self->sgc_bounds(/generous)

	if n_elements(struct) eq 0 then begin
		struct = self->sgc_collate_read('qso')
	endif

	rashift = struct.ra-180
	eq2csurvey, struct.ra, struct.dec, clambda, ceta
	glactc, struct.ra, struct.dec, 2000.0, l, b, 1, /degree

	bounds_logic = $
		(clambda gt sgc_bounds.clambda_range[0]) $
		and (clambda lt sgc_bounds.clambda_range[1]) $
		and (ceta gt sgc_bounds.ceta_range[0]) $
		and (ceta lt sgc_bounds.ceta_range[1]) 


	;bounds_logic = struct.ra ne 10000

	core = sdss_flagval('boss_target1', 'qso_core')
	bonus = sdss_flagval('boss_target1', 'qso_bonus')

	core_logic = (struct.boss_target1 and core) ne 0
	bonus_logic = (struct.boss_target1 and bonus) ne 0

	wall = where(bounds_logic $
		and ( core_logic or bonus_logic) )
	wcore = where(bounds_logic and core_logic )
	wbonus = where(bounds_logic and bonus_logic )

	; first plot up everything in lambda/eta
	psdir = '~/public_html/bosstarget/qso-footprint'	
	file_mkdir, psdir


	begplot, path_join(psdir, 'sgc-sweeps-qso-core-bonus-lameta.eps'), $
		/encap, xsize=11, ysize=8.5

	plot, clambda[wall], ceta[wall], psym=3, $
		/ynozero, $
		xtitle=textoidl('\lambda_c'), ytitle=textoidl('\eta_c')
	legend, ['qso core or bonus'], /right
	endplot

	begplot, path_join(psdir, 'sgc-sweeps-qso-core-bonus-lb.eps'), $
		/encap, xsize=11, ysize=8.5

	plot, l[wall], b[wall], psym=3, $
		/ynozero, $
		xtitle=textoidl('l'), ytitle=textoidl('b')
	legend, ['qso core or bonus'], /right
	endplot



	begplot, path_join(psdir, 'sgc-sweeps-qso-core-lameta.eps'), $
		/encap, xsize=11, ysize=8.5

	plot, clambda[wcore], ceta[wcore], psym=3, $
		/ynozero, $
		xtitle=textoidl('\lambda_c'), ytitle=textoidl('\eta_c')
	legend, ['qso core'], /right
	endplot

	begplot, path_join(psdir, 'sgc-sweeps-qso-core-lb.eps'), $
		/encap, xsize=11, ysize=8.5

	plot, l[wcore], b[wcore], psym=3, $
		/ynozero, $
		xtitle=textoidl('l'), ytitle=textoidl('b')
	legend, ['qso core'], /right
	endplot



	begplot, path_join(psdir, 'sgc-sweeps-qso-bonus-lameta.eps'), $
		/encap, xsize=11, ysize=8.5

	plot, clambda[wbonus], ceta[wbonus], psym=3, $
		/ynozero, $
		xtitle=textoidl('\lambda_c'), ytitle=textoidl('\eta_c')
	legend, ['qso bonus'], /right
	endplot

	begplot, path_join(psdir, 'sgc-sweeps-qso-bonus-lb.eps'), $
		/encap, xsize=11, ysize=8.5

	plot, l[wbonus], b[wbonus], psym=3, $
		/ynozero, $
		xtitle=textoidl('l'), ytitle=textoidl('b')
	legend, ['qso bonus'], /right
	endplot







	if 0 then begin	
		begplot, path_join(psdir, 'sgc-sweeps-qso-radec.eps'), $
			/encap, xsize=11, ysize=8.5

		plot, rashift[wall], struct[wall].dec, psym=3, $
			/ynozero, $
			xtitle=textoidl('ra-180'), ytitle=textoidl('dec')
		legend, ['qso core or bonus'], /right
		endplot
	endif


end

pro esboss::plot_footprint, target_type, target_run, areaname, $
		extra_name=extra_name, $
		dorand=dorand, $
		fracuse=fracuse, $
		target_flags=target_flags, $
		survey=survey, $
		str=str
	
	psdir=getenv('BOSS_TARGET')
	psdir=path_join(psdir, [target_run,'plots'])
	psfile=['bosstarget',target_type,target_run]
	if n_elements(extra_name) ne 0 then begin
		psfile=[psfile,extra_name]
	endif
	psfile=[psfile,areaname+'pos']

	case areaname of
		'all': begin
			doshift=1
			xsize=15
			ysize=8.5
		end
		'ngc': begin
			doshift=0
			xsize=11
			ysize=8.5
			xrange=[100,270]
			xstyle=3
			yrange=[-10,90]
			ystyle=3
		end
		'sgc': begin
			doshift=1
			xsize=11
			ysize=8.5
			xrange=[210,330]
			xstyle=3
			yrange=[-15,40]
			ystyle=3
		end
		else: message,'unknown area name: '+areaname
	endcase

	if doshift then begin
		psfile=[psfile,'rashift']
	endif

	psfile=strjoin(psfile,'-')+'.eps'
	psfile=path_join(psdir, psfile)
	;print,psfile
	;return




	if n_elements(str) eq 0 then begin
		bt=obj_new('bosstarget')
		str=bt->read_collated(target_type,target_run,extra_name=extra_name)
	endif


	if n_elements(target_flags) eq 0 then begin
		w=self->btselect(str,target_type)
	endif else begin
		w=self->btselect(str,target_flags)
	endelse

	begplot,psfile,xsize=xsize,ysize=ysize,/encap, /color


	if doshift then begin
		ra = shiftra(str.ra, 90)
		xtitle='RA-90'
	endif else begin
		ra=str.ra
		xtitle='RA'
	endelse
	dec=str.dec
	ytitle='DEC'

	psym=3
	plot, ra[w], dec[w], psym=psym, $
		xtitle=xtitle, ytitle=ytitle, $
		xrange=xrange, xstyle=xstyle, $
		yrange=yrange, ystyle=ystyle


	if target_type eq 'qso' then begin
		nonphot = where(str[w].photometric eq 0, nbad)
		if nbad ne 0 then begin
			pplot, ra[w[nonphot]], dec[w[nonphot]], psym=3, /over, color='red'
			legend,['targets','non-photometric'],/right,psym=[8,8], $
				color=[!p.color,c2i('red')]
		endif
	endif

	endplot
	print,'Converting to png'
	spawn,'converter -d 100 '+psfile

end

pro esboss::test2, all, sweeps, targets_notinsweeps 
	dir=path_join(getenv('BOSS_TARGET'),'esheldon')

	;sgc=path_join(dir, 'sgc_test/bosstarget-lrg-sgc-collate.fits')
	sgc=path_join(dir, 'sgc_test2/bosstarget-lrg-sgc-collate.fits')
	sgca=path_join(dir, 'sgc_test_all/bosstarget-lrg-sgc-collate.fits')


	if n_tags(sweeps) eq 0 then begin
		print,sgc
		sweeps=mrdfits(sgc,1)
	endif
	if n_tags(all) eq 0 then begin
		print,sgca
		all=mrdfits(sgca,1)
	endif


	primary = sdss_flagval('resolve_status','survey_primary')
	loz=sdss_flagval('boss_target1','gal_loz')
	cmass=sdss_flagval('boss_target1','gal_cmass')

	outdir='~/public_html/bosstarget/lrg-footprint'
	if 0 then begin
		wprigal=where((sweeps.resolve_status and primary) ne 0)
		waprigal=where((all.objc_type eq 3) $
			and ((all.resolve_status and primary) ne 0))

		self->plot_footprint, all[waprigal], $
			psfile=path_join(outdir,'sgc-fpObjc-prigal-lameta.eps'), $
			leg='fpObjc objc_type=3 and primary', /dorand, fracuse=0.3

		self->plot_footprint, sweeps[wprigal], $
			psfile=path_join(outdir, 'sgc-sweeps-prigal-lameta.eps'), $
			leg='sweeps objc_type=3 and primary', /dorand, fracuse=0.3
	endif


	w=where($
		((sweeps.resolve_status and primary) ne 0) $
		and ( (sweeps.boss_target1 and loz) ne 0 $
		      or (sweeps.boss_target1 and cmass) ne 0 ))
	wa=where($
		(all.objc_type eq 3) $
		and ((all.resolve_status and primary) ne 0) $
		and ( (all.boss_target1 and loz) ne 0 $
		      or (all.boss_target1 and cmass) ne 0 ))

	self->plot_footprint, all[wa], $
		psfile=path_join(outdir,'sgc-fpObjc-targets-lameta.eps'), $
		leg='fpObjc objc_type=3 and primary and (loz or cmass)'

	self->plot_footprint, sweeps[w], $
		psfile=path_join(outdir,'sgc-sweeps-targets-lameta.eps'), $
		leg='sweeps objc_type=3 and primary and (loz or cmass)'


	sphoto_match, all[wa], sweeps[w], mall, msweep

	if n_elements(mall) ne n_elements(wa) then begin
		twa = wa
		remove, mall, twa
		targets_notinsweeps = all[twa]

		print,'not in sweeps: ',n_elements(twa),'/',n_elements(wa),$
			format='(a,I0,a,I0)'

		self->plot_footprint, all[twa], $
			psfile=path_join(outdir,'sgc-fpObjc-targets-notinsweeps-lameta.eps'), $
			leg='fpObjc objc_type=3 and primary (loz or cmass) notinsweeps'

		bl=obj_new('bosstarget_lrg')
		begplot,path_join(outdir,'sgc-fpObjc-targets-notinsweeps-rpsf-maghist.eps'),xsize=11,ysize=8.5,/encap

		rpsf=reform(  22.5-2.5*alog10( all.psfflux[2] > 0.001)  )
		plothist, rpsf[twa], min=15, max=24, bin=0.5, $
			xtitle='r-band psf mag'
		legend, 'fpObjc objc_type=3 and primary and (loz or cmass) notinsweeps'
		endplot



		begplot,path_join(outdir,'sgc-fpObjc-targets-notinsweeps-cmodel-maghist.eps'),xsize=11,ysize=8.5,/encap

		bl=obj_new('bosstarget_lrg')
		cmodel=bl->make_cmodelmag(all)
		plothist, cmodel[3,twa], min=15, max=24, bin=0.5, $
			xtitle='i-band cmodel'
		legend, 'fpObjc objc_type=3 and primary and (loz or cmass) notinsweeps'
		endplot


		begplot,path_join(outdir,'sgc-fpObjc-targets-notinsweeps-modelmag-maghist-allband.eps'), xsize=11, ysize=8.5, /encap
		!p.multi=[0,3,2] 
		!p.charsize=!p.charsize*2
		!x.margin=[3,3]

		bl=obj_new('bosstarget_lrg')
		ms = bl->magstruct(targets_notinsweeps)
		for band=0,4 do begin
			plothist,ms.modelmag[band], $
				xtitle='modelmag['+!colors[band]+']', bin=0.1
			if band eq 0 then legend,'not in sweeps', charsize=!p.charsize/2
		endfor
		endplot
		!x.margin = [10,3]
	endif


end

function esboss::_struct_select, str, wstring, nw, verbose=verbose
	command='w=where('+wstring+', nw)'
	if keyword_set(verbose) then begin
		print,command
	endif
	if not execute(command) then begin
		message,'Could not execute command: '+wstring
	endif

	return,w
end

pro esboss::test, fpobjc, sweep, fpobjc_gal_primary_target_flux,sweep_gal_primary_target_flux , cs=cs

	;run=5636
	;run=994
	run=1010
	camcol=1
	rerun=137

	photo_sweep=getenv('PHOTO_SWEEP')
	photo_calib=getenv('PHOTO_CALIB')
	photo_resolve=getenv('PHOTO_RESOLVE')

	if n_tags(fpobjc) eq 0 then begin
		fpobjc = sdss_readobj(run, camcol, rerun=rerun)
	endif
	if n_tags(sweep) eq 0 then begin
		sweep = sweep_readobj(run, camcol, rerun=rerun,type='gal')
	endif

	band=2

	;fpobjc_gal = where(fpobjc.objc_type eq 3 and fpobjc.psfflux[band] gt 2)
	;sweep_gal = where(sweep.objc_type eq 3 and sweep.psfflux[band] gt 2)

	comm='str.objc_type eq 3 and str.psfflux['+strn(band)+'] gt 2'
	fpobjc_gal = self->_struct_select(fpobjc, comm, /verbose)
	sweep_gal = self->_struct_select(sweep, comm)


	primary=sdss_flagval('resolve_status', 'survey_primary')
	fpobjc_gal_primary = where($
		fpobjc.objc_type eq 3 $
		and fpobjc.psfflux[band] gt 2 $
		and (fpobjc.resolve_status and primary) gt 0)
	sweep_gal_primary = where($
		sweep.objc_type eq 3 $
		and sweep.psfflux[band] gt 2 $
		and (sweep.resolve_status and primary) gt 0)

	sphoto_match, fpobjc[fpobjc_gal_primary], sweep[sweep_gal_primary], $
		matched_fpobjc, matched_sweep


	matched_fpobjc = fpobjc_gal_primary[matched_fpobjc]
	matched_sweep = sweep_gal_primary[matched_sweep]

	cs = compare_struct(fpobjc[matched_fpobjc], sweep[matched_sweep])

	for i=0L, n_elements(cs)-1 do begin
		comm='maxdiff=max(abs( fpobjc[matched_fpobjc]'+cs[i].field+' - '+$
			'sweep[matched_sweep]'+cs[i].field+'))'
		if not execute(comm) then begin
			message,'Failed to execute command '+comm
		endif
		print,'max diff for "'+cs[i].field+'": '+strn(maxdiff)

		comm='maxrat=max(abs( fpobjc[matched_fpobjc]'+cs[i].field+' / '+$
			'sweep[matched_sweep]'+cs[i].field+'), min=minrat)'
		if not execute(comm) then begin
			message,'Failed to execute command '+comm
		endif
		print,'max ratio for "'+cs[i].field+'": '+strn(maxrat)
		print,'min ratio for "'+cs[i].field+'": '+strn(minrat)



	endfor

	tfpobjc=fpobjc_gal_primary

	
;	run=1010
;	camcol=1
;	field=58
;	id=354

	help,photo_calib
	help,photo_resolve
	help,photo_sweep
	help,run,rerun
	help,fpobjc, sweep
	help,fpobjc_gal,sweep_gal
	help,fpobjc_gal_primary,sweep_gal_primary

	if n_elements(matched_fpobjc) ne n_elements(tfpobjc) then begin
		remove, matched_fpobjc, tfpobjc
		;matched_fpobjc=fpobjc_gal_primary[matched_fpobjc]
		;matched_sweep=sweep_gal_primary[matched_sweep]
		help,matched_fpobjc, matched_sweep

		bl=obj_new('bosstarget_lrg')
		cmodel = bl->make_cmodelmag(fpobjc)
		begplot,'~/tmp/test-cmodel-fpobjc-notinsweeps-'+strn(band)+'.eps',xsize=11,ysize=8.5,/encap
		plothist, cmodel[3,tfpobjc], bin=0.5, min=15,max=24
		endplot
	endif else begin
		print,'all match'
	endelse

	bl=obj_new('bosstarget_lrg')

	fpres=bl->select(fpobjc)
	swres=bl->select(sweep)

	loz=sdss_flagval('boss_target1', 'gal_loz')
	cmass=sdss_flagval('boss_target1', 'gal_cmass')

	fpobjc_gal_primary_target_flux = where($
		fpobjc.objc_type eq 3 $
		and fpobjc.psfflux[band] gt 2 $
		and (fpobjc.resolve_status and primary) ne 0 $
		and (((fpres and loz) ne 0) or ((fpres and cmass) ne 0 )) )

	fpobjc_gal_primary_cmass_flux = where($
		fpobjc.objc_type eq 3 $
		and fpobjc.psfflux[band] gt 2 $
		and (fpobjc.resolve_status and primary) ne 0 $
		and ((fpres and cmass) ne 0 )  )

	fpobjc_gal_primary_loz_flux = where($
		fpobjc.objc_type eq 3 $
		and fpobjc.psfflux[band] gt 2 $
		and (fpobjc.resolve_status and primary) ne 0 $
		and ((fpres and loz) ne 0 )  )


	

	sweep_gal_primary_target_flux = where($
		sweep.objc_type eq 3 $
		and sweep.psfflux[band] gt 2 $
		and (sweep.resolve_status and primary) ne 0 $
		and (((swres and loz) ne 0) or ((swres and cmass) ne 0 )) )

	sweep_gal_primary_cmass_flux = where($
		sweep.objc_type eq 3 $
		and sweep.psfflux[band] gt 2 $
		and (sweep.resolve_status and primary) ne 0 $
		and ((swres and cmass) ne 0) )

	sweep_gal_primary_loz_flux = where($
		sweep.objc_type eq 3 $
		and sweep.psfflux[band] gt 2 $
		and (sweep.resolve_status and primary) ne 0 $
		and ((swres and loz) ne 0) )



	help,fpobjc_gal_primary_target_flux, $
		fpobjc_gal_primary_cmass_flux, $
		fpobjc_gal_primary_loz_flux, $
		sweep_gal_primary_target_flux ,$
		sweep_gal_primary_cmass_flux , $
		sweep_gal_primary_loz_flux 
end

pro esboss::setup_default
	setenv, 'PHOTO_RESOLVE='+getenv('PHOTO_RESOLVE_DEFAULT')
	setenv, 'PHOTO_SWEEP='+getenv('PHOTO_SWEEP_DEFAULT')
	setenv, 'PHOTO_CALIB='+getenv('PHOTO_CALIB_DEFAULT')
end
pro esboss::setup_sgc, reset=reset
	common sgc_common, $
		photo_sweep_original, photo_resolve_original, photo_calib_original

	if keyword_set(reset) then begin
		setenv, 'PHOTO_RESOLVE='+photo_resolve_original
		setenv, 'PHOTO_SWEEP='+photo_sweep_original
		setenv, 'PHOTO_CALIB='+photo_calib_original
	endif else begin
		photo_resolve_original = getenv('PHOTO_RESOLVE')
		photo_sweep_original = getenv('PHOTO_SWEEP')
		photo_calib_original = getenv('PHOTO_CALIB')
		setenv, 'PHOTO_RESOLVE='+getenv('PHOTO_RESOLVE_SGC')
		setenv, 'PHOTO_SWEEP='+getenv('PHOTO_SWEEP_SGC')
		setenv, 'PHOTO_CALIB='+getenv('PHOTO_CALIB_SGC')
	endelse
end

pro esboss::make_sgc_targets, target_type, subdir, $
		fpobjc=fpobjc, target_dir=target_dir
	if n_elements(target_type) eq 0 then begin
		on_error, 2
		message,'eb->make_sgc_targets, target_type, /fpobjc, target_dir='
	endif

	if n_elements(target_dir) eq 0 then begin
		target_dir = self->sgc_target_dir(subdir)
	endif
	if not file_test(target_dir) then begin
		file_mkdir, target_dir
	endif

	self->setup_sgc

	bt=obj_new('bosstarget')
	; make sure the cached runlist is the one for the above variables
	bt->cache_runlist, /force

	if keyword_set(fpobjc) then begin
		bt->process_fpobjc, target_type, target_dir=target_dir, /primary
	endif else begin
		bt->process_datasweeps, target_type, target_dir=target_dir
	endelse

	self->setup_sgc, /reset

end

function esboss::sgc_target_dir, subdir
	if n_elements(subdir) eq 0 then begin
		message,'usage: dir=eb->sgc_target_dir(subdir)'
	endif

	target_dir=getenv('BOSS_TARGET')

	elements=['esheldon']

	elements=['esheldon','sgc_test2',subdir]

	return, path_join(target_dir, subdir=elements)
end
function esboss::sgc_collate_file, type, subdir, job=job, $
		all=all, fpobjc=fpobjc, anyresolve=anyresolve
	target_dir = self->sgc_target_dir(subdir)

	addstr = type
	if keyword_set(all) then begin
		addstr = 'all-' + addstr
	endif else if keyword_set(fpobjc) then begin
		addstr = 'fpobjc-' + addstr
	endif

	outfile = path_join(target_dir, 'bosstarget-'+addstr+'-sgc-collate.fits')
	if keyword_set(anyresolve) then begin
		outfile=repstr(outfile,'collate','collate-anyresolve')
	endif

	if n_elements(job) ne 0 then begin
		outfile = repstr(outfile, '.fits', '-'+string(job,f='(i03)')+'.fits')
	endif
	return, outfile
end
function esboss::sgc_collate_read, target_type, subdir, all=all, fpobjc=fpobjc
	f=self->sgc_collate_file(target_type, subdir, all=all, fpobjc=fpobjc)
	print,'Reading file: ',f
	return, mrdfits(f, 1)
end


function esboss::qso_select_bypar, str, qso_types, pars, $
		count=count, usewhere=usewhere, res=res, $
		commissioning=commissioning
	if n_elements(str) eq 0 $
			or n_elements(qso_types) eq 0 $
			or n_elements(pars) eq 0 then begin
		message,'Usage: w=eb->qso_select_bypar(str, qso_types, pars, count=, /usewhere)', /inf
		message,'  This will re-select based on the new parameters'
	endif

	primary = sdss_flagval('resolve_status','survey_primary')
	if not keyword_set(usewhere) then begin
		bq=obj_new('bosstarget_qso')
		res = bq->select(str, pars=pars, /struct, /reselect, $
			commissioning=commissioning	)

		w=self->btselect(res, qso_types, /primary, count=count)

	endif else begin
		wstrings=self->qso_where_strings_bypar(pars)

		ws_primary=$
			string('(str.resolve_status and ',primary,') ne 0',f='(a,i0,a)')

		if in(qso_types, 'qso_core') then begin
			add_arrval, wstrings.ws_core, ws
		endif 
		if in(qso_types, 'qso_bonus') then begin
			add_arrval, wstrings.ws_bonus, ws
		endif 
		if in(qso_types, 'qso_known_midz') then begin
			add_arrval, wstrings.ws_known_midz, ws
		endif 
		
		if n_elements(ws) eq 0 then begin
			message,string('Unknown qso type(s): "'+strjoin(qso_types,'", "')+'"')
		endif

		ws = '('+strjoin(ws,' or ')+')'

		ws = ws_primary +' and '+wstrings.ws_known_lohiz+' and '+ws

		command = 'w = where(' + ws + ', count)'
		if not execute(command) then begin
			message,'could not execute command: '+command
		endif

	endelse

	return, w
end

function esboss::qso_where_strings_bypar, pars_input, verbose=verbose

	; get the default pars, combining kde and chi2
	bq=obj_new('bosstarget_qso')
	pars=bq->pars()

	; override defaults
	if n_elements(pars_input) ne 0 then begin
		struct_assign, pars_input, pars, /nozero
	endif

	; below we are anding the x2_star cut in order to
	; also include the chi^2.  In otherwords we are anding
	; kde and chi2 in each case to make composites that get
	; or'ed together

	; OK, now build the where strings
	gmag = '(22.5-2.5*alog10(str.psfflux[1] > 0.001) - str.extinction[1])'

	ws_core_faint=string($
		'str.kde_stardens_faint LT 10d^(',pars.logstardensmax_faint_restrictive,$
		') and str.kde_stardens_faint gt ',-999, $
		' and str.kde_qsodens_faint GE 10d^(',pars.logqsodensmin_faint_restrictive,$
		') and str.x2_star GT ',pars.x2star_restrictive,$
		' and '+gmag+' ge ',pars.gsplit, $
		f='(a,f0.3,a,i0,a,f0.3,a,f0.3,a,f0.3)')
	
	ws_core_bright=string($
		'str.kde_stardens_bright LT 10d^(',pars.logstardensmax_bright_restrictive,$
		') and str.kde_stardens_bright gt ',-999, $
		' and str.kde_qsodens_bright GE 10d^(',pars.logqsodensmin_bright_restrictive,$
		') and str.x2_star GT ',pars.x2star_restrictive,$
		' and '+gmag+' lt ',pars.gsplit, $
		f='(a,f0.3,a,i0,a,f0.3,a,f0.3,a,f0.3)')


	ws_core = '( ('+ws_core_bright+') or ('+ws_core_faint+') )'

	ws_bonus_faint=string($
		'str.kde_stardens_faint LT 10d^(',pars.logstardensmax_faint_permissive,$
		') and str.kde_stardens_faint gt ',-999, $
		' and str.kde_qsodens_faint GE 10d^(',pars.logqsodensmin_faint_permissive,$
		') and str.x2_star GT ',pars.x2star_permissive,$
		' and '+gmag+' ge ',pars.gsplit, $
		f='(a,f0.3,a,i0,a,f0.3,a,f0.3,a,f0.3)')

	ws_bonus_bright=string($
		'str.kde_stardens_bright LT 10d^(',pars.logstardensmax_bright_permissive,$
		') and str.kde_stardens_bright gt ',-999, $
		' and str.kde_qsodens_bright GE 10d^(',pars.logqsodensmin_bright_permissive,$
		') and str.x2_star GT ',pars.x2star_permissive,$
		' and '+gmag+' lt ',pars.gsplit, $
		f='(a,f0.3,a,i0,a,f0.3,a,f0.3,a,f0.3)')



	ws_bonus = '( ('+ws_bonus_bright+') or ('+ws_bonus_faint+') )'



	midz = string(sdss_flagval('boss_target1','qso_known_midz'), f='(I0)')
	ws_known_midz = '((str.boss_target1 and '+midz+') ne 0)'

	ws_any = '( ('+ws_bonus+') or ('+ws_core+') or ('+ws_known_midz+') )'

	; the above usually need to be combined with this
	lohiz = string(sdss_flagval('boss_target1','qso_known_lohiz'), f='(I0)')
	ws_known_lohiz = '(str.boss_target1 and '+lohiz+') eq 0'

	st={ws_any: ws_any, $
		ws_core: ws_core, $
		ws_bonus: ws_bonus, $
		ws_core_faint: ws_core_faint, $
		ws_core_bright: ws_core_bright, $
		ws_bonus_faint: ws_bonus_faint, $
		ws_bonus_bright: ws_bonus_bright, $
		ws_known_midz: ws_known_midz, $
		ws_known_lohiz: ws_known_lohiz}

	st=create_struct(st, pars)

	if keyword_set(verbose) then begin
		print,'core_bright: ',ws_core_bright
		print
		print,'core_faint: ',ws_core_faint
		print
		print,'bonus_bright: ',ws_bonus_bright
		print
		print,'bonus_faint: ',ws_bonus_faint
		print
		print,'known_midz: ',ws_known_midz
		print
		print,'any: ',ws_any
	endif

	return, st
end


pro esboss::sgc_gather, target_type, target_run, all=all, fpobjc=fpobjc, runs=runs, typecut=typecut, anyresolve=anyresolve, outfile=outfile, target_dir=target_dir, match_method=match_method

	; sgc:  means we are defaulting to commissioning style selection, 
	; including sampleing the triangle for lrgs

	if keyword_set(all) then begin
		collate=0
	endif else begin
		collate=1
	endelse

	self->setup_sgc

	bt=obj_new('bosstarget')
	; force using above environment variables
	bt->cache_runlist, /force

	primary = strn(sdss_flagval('resolve_status','survey_primary'))
	photometric = string(sdss_flagval('calib_status','photometric'),f='(i0)')


	if target_type eq 'lrg' then begin
		loz=strn(sdss_flagval('boss_target1', 'gal_loz'))
		cmass=strn(sdss_flagval('boss_target1', 'gal_cmass'))
		triangle=strn(sdss_flagval('boss_target1', 'gal_triangle'))
		lodperp=strn(sdss_flagval('boss_target1','gal_lodperp'))

		orflags=strjoin([loz,cmass,triangle,lodperp],'+')

		ws='(str.boss_target1 and ('+orflags+')) ne 0'

		if keyword_set(typecut) then begin
			ws = ws + ' and (str.objc_type eq 3)'
		endif

	endif else if target_type eq 'qso' then begin
			core = string(sdss_flagval('boss_target1','qso_core'), f='(I0)')
			bonus = string(sdss_flagval('boss_target1','qso_bonus'), f='(I0)')
			midz = string(sdss_flagval('boss_target1','qso_known_midz'), f='(I0)')
			lohiz = string(sdss_flagval('boss_target1','qso_known_lohiz'), f='(I0)')

			orflags=strjoin([core,bonus,midz],'+')

			ws='(str.boss_target1 and '+lohiz+') eq 0' + $
				' and (str.boss_target1 and '+orflags+') ne 0'

		if keyword_set(typecut) then begin
			ws = ws + ' and (str.objc_type eq 6)'
		endif
	endif

	if not keyword_set(anyresolve) then begin
		ws=ws+' and ((str.resolve_status and '+primary+') ne 0)'
	endif

	if n_elements(match_method) ne 0 then begin
		; just get everything that is primary for matches
		ws='(str.resolve_status and '+primary+') ne 0'
		collate=0
	endif else begin
		collate=1
	endelse
	bt->gather2file,target_type, target_run, where_string=ws,$
		collate=collate, columns=columns, /fast, runs=runs, $
		fpobjc=fpobjc, all=all, match_method=match_method, $
		everything=everything, $
		outfile=outfile
	
	self->setup_sgc, /reset
end



function esboss::make_bayes_input, struct
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
		flags2: larrval}	
		
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

pro esboss::run_bayes_prob, struct, probgal, probflags

	print,'Making input'
	inst=self->make_bayes_input(struct)
	compute_bayes_prob, inst, probgal, probflags

end



pro esboss::sgc_test_target_probgal_runlist, runs, reruns
	
	self->setup_sgc
	bt=obj_new('bosstarget')
	; force to make sure it picks up the above environment variables
	bt->runlist, runs, reruns, /force
	self->setup_sgc, /reset


	s=sort(runs)
	runs=runs[s]
	reruns=reruns[s]

	;info=sdss_runlist(runs)
	;w = where(info.stripe eq 81 or info.stripe eq 82 or info.stripe eq 83)
	;w = where(info.stripe eq 81 or info.stripe eq 82 or info.stripe eq 83)

	;runs=runs[w]
	;reruns=reruns[w]

	obj_destroy, bt

end
pro esboss::sgc_test_target_probgal, noclobber=noclobber, new=new
	; read *all* data, not just sweeps of gal or star, and fudge
	; up the objc_type so everything looks like a galaxy
	; then run the target code.
	; then we will look at cutting 
	;	objc_type == 3
	;   probgal > 0.8
	; and see how things look

	self->sgc_test_target_probgal_runlist, runs, junk

	; by default we use the old redux but get the runlist from the sgc
	; stuff because things are missing for the new

	extra='-oldredux'
	if keyword_set(new) then begin
		extra=''
		self->setup_sgc
	endif

	bt=obj_new('bosstarget')
	bt->cache_runlist,/force

	target_dir='/mount/early1/bosstarget-sgc-testbayes'+extra
	if not file_test(target_dir) then begin
		file_mkdir, target_dir
	endif
	target_type = 'lrg'

	bt->process_fpobjc, target_type, target_dir=target_dir, runs=runs, $
		noclobber=noclobber, /primary

	obj_destroy, bt

	if keyword_set(new) then begin
		self->setup_sgc,/reset
	endif

end

pro esboss::sgc_test_target_probgal_gather


	extra='-oldredux'
	self->sgc_test_target_probgal_runlist, runs, junk

	primary = strn(sdss_flagval('resolve_status','survey_primary'))
	loz=strn(sdss_flagval('boss_target1', 'gal_loz'))
	cmass=strn(sdss_flagval('boss_target1', 'gal_cmass'))
	ws='(str.resolve_status and '+primary+') ne 0'+$
		' and ( (str.boss_target1 and '+loz+') ne 0 or (str.boss_target1 and '+cmass+') ne 0)'


	target_dir='/mount/early1/bosstarget-sgc-testbayes'+extra

	self->setup_sgc
	bt=obj_new('bosstarget')
	res=bt->gather('lrg', target_dir=target_dir, $
		where_string=ws, runs=runs)
	self->setup_sgc,/reset

	outdir=target_dir
	outfile='bosstarget-sgc-testbayes'+extra+'-lrg-gather.fits'
	outfile = filepath(root=outdir, outfile)

	print,'Writing to file: ',outfile
	mwrfits, res, outfile, /create

end

pro esboss::sgc_test_plots, struct, psfile=psfile_input, leg=leg, $
		dorand=dorand, fracuse=fracuse

	if n_tags(struct) eq 0 then begin
		outfile = self->sgc_collate_file('lrg', all=all)
		struct = mrdfits(outfile, 1)
	endif

	eq2csurvey, struct.ra, struct.dec, lam, eta
	psdir='~/public_html/bosstarget/lrg-footprint'
	if not file_test(psdir) then begin
		file_mkdir, psdir
	endif
	psfile=path_join(psdir,file_basename(psfile_input))


	begplot,psfile,/color,/encap,xsize=11,ysize=8.5

	xtitle='clambda'
	ytitle='ceta'
	yrange=[80,170]
	ystyle=3
	psym=3

	if n_elements(dorand) ne 0 then begin
		plotrand, lam, eta, psym=psym, fracuse=fracuse,$ 
			xtitle=xtitle, ytitle=ytitle, yrange=yrange, ystyle=ystyle
	endif else begin
		plot, lam, eta, psym=psym, $
			xtitle=xtitle, ytitle=ytitle, yrange=yrange, ystyle=ystyle
	endelse

	if n_elements(leg) ne 0 then legend, leg

	endplot;,/trim;,/landfix


end
pro esboss::_bayes_plots, struct, cmodeli, type
	!p.multi=[0,0,2]
	wo = where(struct.objc_type eq 3, no)
	wp = where(struct.prob_gal_flags eq 0 $
		and struct.prob_gal gt 0.8 and struct.prob_gal le 1.0, np)

	nostr=strn(no)
	npstr=strn(np)
	plothist, cmodeli, bin=0.1, $
		/ylog, yrange=[0.9, 1.e6], ystyle=3, $
		xtitle='cmodel_counts[3]'
	plothist, cmodeli[wo], bin=0.1, /overplot, color='red'
	plothist, cmodeli[wp], bin=0.1, /overplot, color='blue', $
		line=2

	if n_elements(type) ne 0 then begin
		legend,type,/left
	endif


	if type eq 'gal_loz' then bin=0.05 else bin=0.01
	plothist, struct.prob_gal, bin=bin, min=0,max=1,$
		/ylog, yrange=[0.9, 6.e6], ystyle=3,$
		xtitle='prob_gal'
	plothist, struct[wo].prob_gal, bin=bin, min=0,max=1,$
		/overplot, color='red'
	plothist, struct[wp].prob_gal, bin=bin, min=0,max=1,$
		/overplot, color='blue', line=2

	legend,$
		['gal+stars', $
		'objc_type = 3 ('+nostr+')',$
		'prob_gal > 0.8 ('+npstr+')'],/center,/top,$
		color=[!p.color,c2i('red'),c2i('blue')],line=[0,0,2]

	!p.multi=0
end


pro esboss::sgc_test_target_probgal_plots, struct, cmodeli, new=new

	extra='-oldredux'
	if keyword_set(new) then begin
		self->setup_sgc
		extra=''
	endif


	target_dir=getenv('BOSS_TARGET')
	target_dir=filepath(root=target_dir, 'esheldon/sgc_test')

	outfile = filepath(root=target_dir, 'bosstarget-lrg-sgc-collate.fits')

	if n_tags(struct) eq 0 then begin
		struct=mrdfits(outfile, 1)
	endif

	bl=obj_new('bosstarget_lrg')
	if n_elements(cmodeli) ne n_elements(struct) then begin
		print,'Creating cmodel counts'
		cmodel = bl->make_cmodelmag(struct)
		cmodel = cmodel - struct.extinction
		cmodeli = reform(cmodel[3,*])
	endif
	wobjc_type3=where(struct.objc_type eq 3 and cmodeli lt 17,nw)
	if nw ne 0 then begin
		help,wobjc_type3
		struct[wobjc_type3].prob_gal = 1.0
		struct[wobjc_type3].prob_gal_flags = 0
	endif else begin
		message,'no objc_type == 3 and cmodel i < 17'
	endelse


	both_logic = $
		((struct.boss_target1 and sdss_flagval('boss_target1','gal_loz')) ne 0) $
		or $
		((struct.boss_target1 and sdss_flagval('boss_target1','gal_cmass')) ne 0)

	if 1 then begin
		eq2csurvey, struct.ra, struct.dec, lam, eta
		psfile=$
			  '~/public_html/bosstarget/lrg-sgsep/sgcap-sgsep-bayestest'+ $
			  extra+'-objc-type-lameta.eps'

		pngfile=$
			  '~/public_html/bosstarget/lrg-sgsep/sgcap-sgsep-bayestest'+ $
			  extra+'-objc-type-lameta.png'
		doz=0

		if doz then begin
			setupplot, 'Z'
			device, set_resolution=[1280,1024]
		endif else begin
			begplot,psfile,/color,/encap,xsize=11,ysize=8.5
		endelse

		wobjc=where(both_logic and struct.objc_type eq 3)
		plot, lam[wobjc], eta[wobjc], psym=3,$
			xtitle='clambda', ytitle='ceta',$
			yrange=[80,170],ystyle=3
		legend, 'objc_type = 3 and (gal_loz or gal_cmass)'

		if not doz then begin
			endplot;,/trim;,/landfix
		endif else begin
			print,'Writing ',pngfile
			write_image, pngfile, 'png', tvrd()
		endelse


		psfile=$
			  '~/public_html/bosstarget/lrg-sgsep/sgcap-sgsep-bayestest'+ $
			  extra+'-probgal-lameta.eps'
		pngfile=$
			  '~/public_html/bosstarget/lrg-sgsep/sgcap-sgsep-bayestest'+ $
			  extra+'-probgal-lameta.png'

		if not doz then begin
			begplot,psfile,/color,/encap,xsize=11,ysize=8.5
		endif
		wprob=where(both_logic $
			and struct.prob_gal gt 0.8 and struct.prob_gal le 1.0)
		plot, lam[wprob], eta[wprob], psym=3,$
			xtitle='clambda', ytitle='ceta',$
			yrange=[80,170],ystyle=3

		legend, 'prob_gal > 0.8 and (gal_loz or gal_cmass)'
		if not doz then begin
			endplot;,/trim;,/landfix
		endif else begin
			print,'Writing ',pngfile
			write_image, pngfile, 'png', tvrd()
		endelse

		if doz then begin
			setupplot,'X'
		endif

	endif




	types=['gal_loz-or-gal_cmass','gal_loz','gal_cmass']
	nt=n_elements(types)
	for i=0L, nt-1  do begin

		type = types[i]
		psfile=$
		  '~/public_html/bosstarget/lrg-sgsep/sgcap-sgsep-bayestest'+ $
		  extra+'-'+type+'.eps'
		begplot,psfile,/color,/encap,ysize=11,xsize=8.5


		
		if type eq types[0] then begin
			w=where(both_logic)
		endif else begin
			btflag = sdss_flagval('boss_target1',type)
			w = where((struct.boss_target1 and btflag) ne 0)
		endelse
		
		self->_bayes_plots,struct[w], cmodeli[w], type
		if i ne nt-1 then key=prompt_kbrd('hit enter')

		endplot
	endfor	


end

pro esboss::gather_std
	bt=obj_new('bosstarget')
	res=bt->gather('std', target_dir='/mount/early1/bosstarget', $
		/collate)
	file=path_join(target_dir, 'bosstarget-std-collate.fits')
	mwrfits, res, file, /create
end


function esboss::splitlist, list, nsplit=nsplit, nper=nper

	nlist = n_elements(list)
	if n_elements(nper) ne 0 then begin
		nsplit = (nlist/nper) + (nlist mod nper)
	endif else if n_elements(nsplit) ne 0 then begin
		nper = nlist/nsplit
		nleft = nlist mod nsplit
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
function esboss::riemann_nodeinfo
	return, {nnodes: 14, corespernode: 8, mempernode: 8}
end
pro esboss::sgc_runlist, runs, reruns, nsplit=nsplit, nper=nper
	self->setup_sgc

	bt=obj_new('bosstarget')
	; make sure the cached runlist is the one for the above variables
	bt->cache_runlist, /force

	bt->runlist, runs, reruns

	if n_elements(nsplit) ne 0 or n_elements(nper) ne 0 then begin
		runs = self->splitlist(runs, nsplit=nsplit, nper=nper)
		reruns = self->splitlist(reruns, nsplit=nsplit,  nper=nper)
	endif

	self->setup_sgc, /reset

	obj_destroy, bt

end

function esboss::sgc_bounds, padding=padding

	dir=getenv("BOSSTARGET_DIR")
	if dir eq '' then message,'BOSSTARGET not set up'
	dir = path_join(dir, 'geometry')
	file = path_join(dir, 'boss_survey.par')

	bounds = yanny_readone(file,/anon)
	areaname='SGC'
	w=where(bounds.areaname eq areaname, nb)

	if n_elements(padding) eq 0 then begin
		padding=0
	endif
	ceta_range = [min(bounds[w].cetaMin)-padding, $
		          max(bounds[w].cetaMax)+padding]
	clambda_range = [min(bounds[w].clambdaMin)-padding, $
		             max(bounds[w].clambdaMax)+padding]


	return, $
		{clambda_range: clambda_range, $ 
		ceta_range: ceta_range}
end

function esboss::riemann_default_njobs
	return, 50
end
function esboss::sgc_extra_logic
	bounds=self->sgc_bounds(padding=5)
	clmin = strn(bounds.clambda_range[0])
	clmax = strn(bounds.clambda_range[1])
	cemin = strn(bounds.ceta_range[0])
	cemax = strn(bounds.ceta_range[1])

	extra_logic =  $
		'(str.clambda gt '+clmin+') and (str.clambda lt '+clmax+') '+$
		'and (str.ceta gt '+cemin+') and (str.ceta lt '+cemax+')'

	return, extra_logic
end
pro esboss::make_partial_sgc_targets, target_type, target_run, jobnum, $
		nper=nper, $
		run_index=run_index, $
		all=all, $
		fpobjc=fpobjc, $
		_extra=_extra


	if n_elements(target_type) eq 0 or n_elements(jobnum) eq 0 then begin
		on_error, 2
		message,'eb->make_partial_sgc_targets, target_type, target_run, jobnum, nper=2, run_index=, /all, _extra='
	endif

	;if target_type eq 'qso' or target_type eq 'std' then begin
	;	extra_logic = self->sgc_extra_logic()
	;endif
	if n_elements(nper) eq 0 then begin
		nper=2
	endif
	self->sgc_runlist, runptrs, rerunptrs, nper=nper
	runs = *runptrs[jobnum]
	reruns = *rerunptrs[jobnum]

	if n_elements(run_index) ne 0 then begin
		runs=runs[run_index]
		reruns=reruns[run_index]
	endif

	splog,'Using runs: ',runs

	self->setup_sgc

	bt=obj_new('bosstarget')
	; make sure the cached runlist is the one for the above variables
	bt->cache_runlist, /force

	bt->process_runs, target_type, target_run,  $
		runs=runs, extra_logic=extra_logic, fpobjc=fpobjc, all=all, $
		_extra=_extra

	self->setup_sgc, /reset

	obj_destroy, bt

end

pro esboss::sgc_gather_partial, target_type, target_run, jobnum, nper=nper, $
		run_index=run_index, all=all, fpobjc=fpobjc, $
		combine=combine,ascii=ascii, $ ; these with regard to combine only
		match_method=match_method, $
		outfile=outfile, $
		_extra=_extra

	if n_elements(target_type) eq 0 or (n_elements(jobnum) eq 0 and not keyword_set(combine)) then begin
		on_error, 2
		message,'eb->sgc_gather_partial, target_type, target_run, jobnum, nper=2, run_index=, /all, /combine, /ascii, _extra=', /inf
		message,'Send /combine to combine the individual files.  ',/inf
		message,'  /ascii is with regard to the combined file only'
	endif

	;if target_type eq 'qso' or target_type eq 'std' then begin
	;	extra_logic = self->sgc_extra_logic()
	;endif
	if n_elements(nper) eq 0 then begin
		nper=2
	endif
	self->sgc_runlist, runptrs, rerunptrs, nper=nper
	njobs=n_elements(runptrs)

	bt=obj_new('bosstarget')
	if keyword_set(combine) then begin
		ptr_free, runptrs, rerunptrs
		outfile=bt->target_file(target_type, target_run, $
			all=all, fpobjc=fpobjc, match_method=match_method, /collate)
		flist=strarr(njobs)
		for job=0L, njobs-1 do begin
			flist[job] = $
				repstr(outfile,'.fits','-'+string(job,f='(i03)')+'.fits')
		endfor
		for ii=0L, n_elements(flist)-1 do begin
			if file_test(flist[ii]) then begin
				hdr=headfits(flist[ii],ext=1)
				help,hdr
				break
			endif
		endfor
		fxhclean,hdr
		sxdelpar, hdr, 'COMMENT'
		tot=mrdfits_multi(flist)
		if n_tags(tot) eq 0 then message,'Failed to read'
		ntot=n_elements(tot)
		if keyword_set(ascii) then begin
			outfile = repstr(outfile, '.fits', '-tab.st')
			splog,'Writing ',ntot,' to file: ',outfile,format='(a,i0,a)'
			write_idlstruct, tot, outfile, /ascii
		endif else begin
			splog,'Writing ',ntot,' to file: ',outfile,format='(a,i0,a)'
			mwrfits, tot, outfile, hdr, /create
		endelse

		return
	endif

	runs = *runptrs[jobnum]
	reruns = *rerunptrs[jobnum]
	ptr_free, runptrs, rerunptrs

	if n_elements(run_index) ne 0 then begin
		runs=runs[run_index]
		reruns=reruns[run_index]
	endif

	splog, 'Using runs: ',runs


	if n_elements(outfile) eq 0 then begin
		outfile=bt->target_file(target_type, target_run, $
			all=all, fpobjc=fpobjc, /collate, match_method=match_method)
		outfile = repstr(outfile,'.fits','-'+string(jobnum,f='(i03)')+'.fits')
	endif
	splog,'outfile: ',outfile,format='(a,a)'

	self->sgc_gather, $
		target_type, target_run, $
		fpobjc=fpobjc, runs=runs, $
		match_method=match_method, $
		outfile=outfile


	obj_destroy, bt

end


function esboss::sgc_slow_runs
	runs = $
		[4836, 4879, 7717, 7773, 7778, $   ; slowest
		 1659, 4207, 4822, 7712, 7777, 7824, $   
		 1739, 4192, 5640, 7713, 7781] ; faster

	 return,runs
end


; This version requires a version of bosstarget to be given
pro esboss::create_pbs_bycamcol, target_type, target_run, $
		bosstarget_v=bosstarget_v, $
		photoop_v=photoop_v, $
		idlutils_v=idlutils_v, $
		photo_calib=photo_calib, $
		photo_resolve=photo_resolve, $
		photo_sweep=photo_sweep, $
		pars=pars, $
		runs=runs2do, $
		oldcache=oldcache, $
		ignore_resolve=ignore_resolve, $
		commissioning=commissioning, $
		comm2=comm2 ; this is qso only currently


	if n_elements(bosstarget_v) eq 0 then begin
		bosstarget_v = "-r /home/esheldon/exports/bosstarget-work"
	endif
	if n_elements(photoop_v) eq 0 then photoop_v = "v1_9_4"
	if n_elements(idlutils_v) eq 0 then begin
		idlutils_v = "-r /home/esheldon/exports/idlutils-work"
	endif

	if n_elements(commissioning) eq 0 then commissioning=0
	if n_elements(comm2) eq 0 then comm2=0

	add_arrval,"setup tree", setups 
	add_arrval,"setup photoop "+photoop_v, setups
	add_arrval,"setup idlutils "+idlutils_v, setups
	add_arrval,"setup bosstarget "+bosstarget_v, setups

	
	if n_elements(photo_sweep) eq 0 then PHOTO_SWEEP='/clusterfs/riemann/raid007/bosswork/groups/boss/sweeps/2009-11-16.v2'
	if n_elements(photo_resolve) eq 0 then PHOTO_RESOLVE='/clusterfs/riemann/raid006/bosswork/groups/boss/resolve/2009-11-16'
	if n_elements(photo_calib) eq 0 then PHOTO_CALIB='/clusterfs/riemann/raid007/bosswork/groups/boss/calib/2009-06-14/calibs/fall09i'



	add_arrval, 'export PHOTO_SWEEP='+photo_sweep, setups
	add_arrval, 'export PHOTO_RESOLVE='+photo_resolve, setups
	add_arrval, 'export PHOTO_CALIB='+photo_calib, setups


	setups = strjoin(setups, ' && ')

	sweep_old=getenv('PHOTO_SWEEP')
	resolve_old=getenv('PHOTO_RESOLVE')
	setenv, 'PHOTO_SWEEP='+PHOTO_SWEEP
	setenv, 'PHOTO_RESOLVE='+PHOTO_RESOLVE
	bt=obj_new('bosstarget')
	bt->cache_runlist, /force
	if n_elements(runs2do) ne 0 then begin
		bt->match_runlist, runs2do, runs, reruns
	endif else begin
		bt->runlist, runs, reruns
	endelse
	setenv, 'PHOTO_SWEEP='+sweep_old
	setenv, 'PHOTO_RESOLVE='+resolve_old

	print,'Found: ',n_elements(runs),' runs',f='(a,i0,a)'


	pbs_dir='/home/esheldon/pbs/'+target_type+'/'+target_run
	file_mkdir, pbs_dir

	qsub_file = path_join(pbs_dir, 'submit-'+target_type+'-bycamcol')
	fbase = target_type

	if keyword_set(commissioning) then begin
		qsub_file+='-comm'
		fbase+='-comm'
	endif else if keyword_set(comm2) then begin
		qsub_file+='-comm2'
		fbase+='-comm2'
	endif

	qsub_file+='.sh'
	openw, qsub_lun, qsub_file, /get_lun


	nrun=n_elements(runs)
	ntot=nrun*6
	ii=0L
	for i=0L, nrun-1 do begin
		run = runs[i]
		rerun = reruns[i]

		rstr = run2string(run)

		for camcol=1,6 do begin
			cstr=string(camcol,f='(i0)')

			job_name = target_type+'-'+rstr+'-'+cstr

			pbs_file = repstr(job_name, target_type, fbase)+'.pbs'
			pbs_file=filepath(root=pbs_dir, pbs_file)

			target_command = $
				string(f='(%"%s, %s, %s, run=%d, camcol=%d")', $
				"    bt->process_runs", $
				"'"+target_type+"'", "'"+target_run+"'", run, camcol)


			;if keyword_set(commissioning) then begin
			;	target_command += ", /commissioning"
			;endif else if keyword_set(comm2) then begin
			;	target_command += ", /comm2"
			;endif
			idl_commands="bt=obj_new('bosstarget'"
			if keyword_set(commissioning) then begin
				idl_commands += ",/commissioning"
			endif else if keyword_set(comm2) then begin
				idl_commands += ",/comm2"
			endif
			if keyword_set(ignore_resolve) then begin
				idl_commands += ",/ignore_resolve"
			endif
			if keyword_set(oldcache) and target_type eq 'qso' then begin
				idl_commands += ",/oldcache"
			endif
			idl_commands += ")"

			if n_elements(pars) ne 0 then begin
				idl_commands = [idl_commands, 'pars='+tostring(pars)]
				target_command += ', pars=pars'
			endif 

			idl_commands=[idl_commands,target_command]

			pbs_riemann_idl, $
				pbs_file, idl_commands, setup=setups, job_name=job_name


			printf, qsub_lun, $
				'echo -n "',ii+1,'/',ntot,' ',pbs_file,' "',$
				format='(a,i0,a,i0,a,a,a)'
			printf, qsub_lun, 'qsub '+pbs_file
			ii=ii+1
		endfor
	endfor


	free_lun, qsub_lun

end






; This version requires a version of bosstarget to be given
pro esboss::create_pbs, target_type, target_run,  $
		bosstarget_v=bosstarget_v, $
		photoop_v=photoop_v, $
		idlutils_v=idlutils_v, $
		photo_calib=photo_calib, $
		photo_resolve=photo_resolve, $
		photo_sweep=photo_sweep, $
		$
		pars=pars, reselect=reselect, extra_name=extra_name, $
		$
		dotarget=dotarget,$
		dogather=dogather,$
		nper=nper, $
		fpobjc=fpobjc, $
		match_method=match_method, $
		commissioning=commissioning, $
		comm2=comm2, $
		noverify=noverify, $
		where_string=where_string, $
		$
		walltime=walltime

	if n_elements(dotarget) eq 0 and n_elements(dogather) eq 0 then begin
		dotarget=1
		dogather=1
	endif
	if n_elements(nper) eq 0 then begin
		nper=2
	endif


	if n_elements(fpobjc) eq 0 then fpobjc=0
	if n_elements(commissioning) eq 0 then commissioning=0
	if n_elements(comm2) eq 0 then comm2=0



	if n_elements(bosstarget_v) eq 0 then begin
		bosstarget_v = "-r /home/esheldon/exports/bosstarget-work"
	endif

	if n_elements(photoop_v) eq 0 then photoop_v = "v1_9_4"
	if n_elements(idlutils_v) eq 0 then begin
		idlutils_v = "-r /home/esheldon/exports/idlutils-work"
	endif


	;if n_elements(pars) ne 0 and n_elements(extra_name) eq 0 then begin
	;	message,'You must send an extra_name= with pars'
	;endif


	;add_arrval,"setup sas bosswork", setups 
	add_arrval,"setup tree", setups 
	add_arrval,"setup photoop "+photoop_v, setups
	add_arrval,"setup idlutils "+idlutils_v, setups
	add_arrval,"setup bosstarget "+bosstarget_v, setups



	if n_elements(photo_sweep) eq 0 then PHOTO_SWEEP='/clusterfs/riemann/raid007/bosswork/groups/boss/sweeps/2009-11-16.v2'
	if n_elements(photo_resolve) eq 0 then PHOTO_RESOLVE='/clusterfs/riemann/raid006/bosswork/groups/boss/resolve/2009-11-16'
	if n_elements(photo_calib) eq 0 then PHOTO_CALIB='/clusterfs/riemann/raid007/bosswork/groups/boss/calib/2009-06-14/calibs/fall09i'





	add_arrval, 'export PHOTO_SWEEP='+photo_sweep, setups
	add_arrval, 'export PHOTO_RESOLVE='+photo_resolve, setups
	add_arrval, 'export PHOTO_CALIB='+photo_calib, setups


	setups = strjoin(setups, ' && ')

	sweep_old=getenv('PHOTO_SWEEP')
	resolve_old=getenv('PHOTO_RESOLVE')
	setenv, 'PHOTO_SWEEP='+PHOTO_SWEEP
	setenv, 'PHOTO_RESOLVE='+PHOTO_RESOLVE

	bt=obj_new('bosstarget')
	bt->cache_runlist, /force
	bt->split_runlist, runptrs, rerunptrs, nper=nper
	njobs = n_elements(runptrs)
	ptr_free, runptrs, rerunptrs

	setenv, 'PHOTO_SWEEP='+sweep_old
	setenv, 'PHOTO_RESOLVE='+resolve_old


	pbs_dir='/home/esheldon/pbs/'+target_type+'/'+target_run
	file_mkdir, pbs_dir

	qsub_file = path_join(pbs_dir, 'submit-'+target_type)
	check_file = path_join(pbs_dir, 'check.sh')
	combine_file = path_join(pbs_dir, 'combine.sh')
	if n_elements(extra_name) ne 0 then begin
		combine_file = repstr(combine_file, '.sh', '-'+extra_name+'.sh')
	endif

	fbase = target_type

	if fpobjc then begin
		qsub_file+='-fpobjc'
		fbase+='-fpobjc'
	endif
	if keyword_set(commissioning) then begin
		qsub_file+='-comm'
		fbase+='-comm'
	endif else if keyword_set(comm2) then begin
		qsub_file+='-comm2'
		fbase+='-comm2'
	endif


	if n_elements(match_method) ne 0 then begin
		qsub_file += '-match-gather'
		fbase += '-match-gather'
	endif else begin
		if keyword_set(dotarget) then begin
			qsub_file += '-target'
			fbase += '-target'
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

	numstr = string(njobs-1, format='(I03)')
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
	printf, lun, '  photo_sweep="'+photo_sweep+'"'
	printf, lun, '  bt=obj_new("bosstarget")'
	printf, lun, $
		"  bt->gather_partial, '",target_type,"','",target_run,"'",+$
		", photo_sweep=photo_sweep, /combine",extra_gather, $
		f='(a,a,a,a,a,a,a)'
	printf, lun, 'EOF'
	free_lun, lun



	for job=0L, njobs-1 do begin

		jobstr = string(job,format='(I03)')

		jfbase = fbase + '-'+jobstr
		
		pbs_file=filepath(root=pbs_dir, jfbase+'.pbs')

		target_command = $
			string(f='(%"%s, %s, %s, %d, nper=%d")', $
			"    bt->process_partial_runs", $
			"'"+target_type+"'", "'"+target_run+"'", job, nper)

		gather_command = string(f='(%"%s, %s, %s, %d, nper=%d")', $
			"    bt->gather_partial", $
			"'"+target_type+"'", "'"+target_run+"'", job, nper)

		if keyword_set(fpobjc) then begin
			target_command += ", /fpobjc"
			gather_command += ", /fpobjc"
		endif
		if n_elements(match_method) ne 0 then begin
			target_command += ", match_method='"+match_method+"'"
			gather_command += ", match_method='"+match_method+"'"
		endif 

		;if keyword_set(commissioning) then begin
		;	target_command += ", /commissioning"
		;	gather_command += ", /commissioning"
		;endif else if keyword_set(comm2) then begin
		;	target_command += ", /comm2"
		;endif

		if keyword_set(commissioning) then begin
			idl_commands = "bt=obj_new('bosstarget',/commissioning)"
		endif else if keyword_set(comm2) then begin
			idl_commands = "bt=obj_new('bosstarget',/comm2)"
		endif else begin
			idl_commands = "bt=obj_new('bosstarget')"
		endelse




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
			target_command += ", pars=pars"
			; turn this off for std since was *adding* boss_target1 in
			; really dumb
			;gather_command += ", pars=pars"
		endif 
		if n_elements(extra_name) ne 0 then begin
			;gather_command += ", pars=pars, extra_name='"+extra_name+"'"
			gather_command += ", extra_name='"+extra_name+"'"
		endif
			

		if keyword_set(dotarget) then begin
			idl_commands=[idl_commands,target_command]
		endif
		if keyword_set(dogather) then begin
			idl_commands=[idl_commands,gather_command]
		endif

		job_name = target_type+"-"+jobstr

		pbs_riemann_idl, pbs_file, idl_commands, setup=setups, job_name=job_name, $
			walltime=walltime
	endfor
end



pro esboss::create_sgc_pbs, target_type, target_run, $
		dotarget=dotarget,$
		dogather=dogather,$
		sleep=sleep, $
		nper=nper, $
		fpobjc=fpobjc, $
		match_method=match_method, $
		commissioning=commissioning

	if n_elements(dotarget) eq 0 and n_elements(dogather) eq 0 then begin
		dotarget=1
		dogather=1
	endif
	if n_elements(nper) eq 0 then begin
		nper=2
	endif
	if n_elements(sleep) ne 0 then begin
		sleepstr=strn(sleep)
	endif


	if n_elements(fpobjc) eq 0 then fpobjc=0
	if n_elements(commissioning) eq 0 then commissioning=0

	self->sgc_runlist, runptrs, rerunptrs, nper=nper
	njobs = n_elements(runptrs)
	ptr_free, runptrs, rerunptrs

	pbs_dir='/home/esheldon/pbs/'+target_type+'/'+target_run
	file_mkdir, pbs_dir

	qsub_file = path_join(pbs_dir, 'submit-'+target_type)
	fbase = 'sgc-'+target_type

	if fpobjc then begin
		qsub_file+='-fpobjc'
		fbase+='-fpobjc'
	endif
	if keyword_set(commissioning) then begin
		qsub_file+='-comm'
		fbase+='-comm'
	endif

	if n_elements(match_method) ne 0 then begin
		qsub_file += '-match-gather'
		fbase += '-match-gather'
	endif else begin
		qsub_file += '-target-gather'
		fbase += '-target-gather'
	endelse

	qsub_file+='.sh'
	openw, lun, qsub_file, /get_lun

	numstr = string(njobs-1, format='(I03)')
    printf, lun
    printf, lun, 'for i in `seq -w 0 '+numstr+'`; do' 
    printf, lun, '    file='+fbase+'-${i}.pbs'
    printf, lun, '    echo "qsub $file"'
    printf, lun, '    qsub $file'
	if n_elements(sleep) ne 0 then begin
		printf, lun, '    sleep '+sleepstr
	endif
    printf, lun, 'done'
	free_lun, lun

	for job=0L, njobs-1 do begin

		jobstr = string(job,format='(I03)')

		jfbase = fbase + '-'+jobstr
		
		pbs_file=filepath(root=pbs_dir, jfbase+'.pbs')

		target_command = $
			string(f='(%"%s, %s, %s, %d")', $
			"    eb->make_partial_sgc_targets", $
			"'"+target_type+"'", "'"+target_run+"'", job)

		gather_command = string(f='(%"%s, %s, %s, %d")', $
			"    eb->sgc_gather_partial", $
			"'"+target_type+"'", "'"+target_run+"'", job)

		if keyword_set(fpobjc) then begin
			target_command += ", /fpobjc"
			gather_command += ", /fpobjc"
		endif
		if n_elements(match_method) ne 0 then begin
			target_command += ", match_method='"+match_method+"'"
			gather_command += ", match_method='"+match_method+"'"
		endif
		if keyword_set(commissioning) then begin
			target_command += ", /commissioning"
		endif


		idl_commands = "eb=obj_new('esboss')"
		 if keyword_set(dotarget) then begin
			 idl_commands=[idl_commands,target_command]
		 endif
		 if keyword_set(dogather) then begin
			 idl_commands=[idl_commands,gather_command]
		 endif

		job_name = target_type+"-"+jobstr

		pbs_riemann_idl, pbs_file, idl_commands, job_name=job_name
	endfor
end


function esboss::get_all_sweep_filenames, type
	bt=obj_new('bosstarget')
	bt->runlist, runs, reruns
	ntot = n_elements(runs)*6

	flist = strarr(ntot)
	ii=0
	for i=0L, n_elements(runs)-1 do begin
		run=runs[i]
		rerun=reruns[i]
		for camcol=1,6 do begin
			sweep_file = sdss_name('calibObj.'+type,run,camcol, rerun=rerun)
			foundfile = (file_search(sweep_file+'*'))[0]
			if foundfile eq '' then begin
				print,'File not found: ',sweep_file
			endif else begin
				flist[ii] = foundfile
			endelse

			ii = ii+1
		endfor
	endfor

	return, flist
end


pro esboss::test_new_qso,objs,oldres,newres
	bt=obj_new('bosstarget')
	bt->runlist, runs, reruns

	n=35
	if n_elements(objs) eq 0 then begin
		objs=sweep_readobj(runs[n], 3, rerun=reruns[n], type='star')
		oldres = bt->read('qso', runs[n], 3, target_dir='/clusterfs/riemann/raid006/bosswork/groups/boss/target/esheldon/sgc_test2', /collate)
	endif

	bq=obj_new('bosstarget_qso')
	if n_elements(newres) eq 0 then begin
		newres = bq->select(objs, /struct)
		mwrfits, newres, '~/tmp/testqso.fits', /create
	endif

	core=sdss_flagval('boss_target1','qso_core')
	bonus=sdss_flagval('boss_target1','qso_bonus')
	midz=sdss_flagval('boss_target1','qso_known_midz')
	lohiz=sdss_flagval('boss_target1','qso_known_lohiz')
	;midz=sdss_flagval('boss_target1','qso_known_mid')
	;lohiz=sdss_flagval('boss_target1','qso_known_lohi')
	photometric = sdss_flagval('calib_status','photometric')

	wold_corebonus=where( $
		(oldres.boss_target1 and (core+bonus)) ne 0 $
		and (oldres.calib_status and photometric) ne 0, nold )	

	wnew_corebonus=where( $
		(newres.boss_target1 and (core + bonus)) ne 0 $
	)


	wnew_midz_and_corebonus = where($
		(newres.boss_target1 and midz) ne 0 $
		and (newres.boss_target1 and (core+bonus)) ne 0)
	help,wnew_midz_and_corebonus
	print,wnew_midz_and_corebonus[sort(wnew_midz_and_corebonus)]

	wnew_lohiz_and_corebonus = where($
		(newres.boss_target1 and lohiz) ne 0 $
		and (newres.boss_target1 and (core+bonus)) ne 0)
	help,wnew_lohiz_and_corebonus
	print,wnew_lohiz_and_corebonus[sort(wnew_lohiz_and_corebonus)]

	wnew_lohiz_and_corebonusmidz = where($
		(newres.boss_target1 and lohiz) ne 0 $
		and (newres.boss_target1 and (core+bonus+midz)) ne 0)
	help,wnew_lohiz_and_corebonusmidz
	print,wnew_lohiz_and_corebonusmidz[sort(wnew_lohiz_and_corebonusmidz)]



	wnew=where( $
		(newres.boss_target1 and lohiz) eq 0 $
		and (newres.boss_target1 and (core + bonus + midz)) ne 0 $
	)

	help, wold_corebonus, wnew_corebonus
	match, wold_corebonus, wnew_corebonus, $
		match_corebonus_old, match_corebonus_new, /sort
	help,match_corebonus_old, match_corebonus_new

	help, wnew
	match, wnew_corebonus, wnew, mncb, mn, /sort
	help, mncb, mn
	remove, mncb, wnew_corebonus
	remove, mn, wnew
	print,wnew_corebonus
	print,wnew


end




pro esboss::compare_concentration, $
		sweep, sweep_iconc, fpobjc_iconc
	; compare outputs from 
	;  standard gal datasweep run
	;  adding i concentration > 0.3 to datasweep run
	;  running on fpObjc with i concentration > 0.3

	outdir='/clusterfs/riemann/raid006/bosswork/groups/boss/target/esheldon/sgc_test2/compare-concentration'

	if n_elements(sweep) eq 0 then begin
		sweep=self->sgc_collate_read('lrg','2009-06-02')
		sweep_iconc = self->sgc_collate_read('lrg','2009-06-08-testbed')
		fpobjc_iconc = self->sgc_collate_read('lrg','2009-06-08-testbed', $
			/fpobjc)
	endif


	target_flags=['gal_loz','gal_cmass']
	wsweep = self->btselect(sweep, 'lrg', target_flags, count=nsweep)
	wsweep_iconc = self->btselect(sweep_iconc, 'lrg', nsweep_iconc)
	wfpobjc_iconc = self->btselect(fpobjc_iconc, 'lrg', nfpobjc_iconc)

	help,sweep,wsweep
	help,sweep_iconc,wsweep_iconc
	help,fpobjc_iconc, wfpobjc_iconc
	print

	print,'matching sweep to sweep_iconc'
	sphoto_match, sweep[wsweep], sweep_iconc[wsweep_iconc], $
		m_s2si, m_si2s
	help,wsweep,wsweep_iconc,m_s2si, m_si2s

	print
	print,'matching sweep_iconc to fpobjc_iconc'
	sphoto_match, sweep_iconc[wsweep_iconc], fpobjc_iconc[wfpobjc_iconc], $
		m_si2fi, m_fi2si
	help,wsweep_iconc,wfpobjc_iconc,m_si2fi, m_fi2si

	print,'in fpObjc but not sweep for iconc: ',(nfpobjc_iconc-nsweep_iconc)/float(nfpobjc_iconc)*100,'%',f='(a,g0.3,a)'
		


	; get the objects thrown out from the sweep run when demanding
	; iconc > 0.3
	; remove the matches
	sweep_removed_by_iconc = wsweep
	remove, m_s2si, sweep_removed_by_iconc 
	help, sweep_removed_by_iconc 


	winterp_sweep = where((sweep[wsweep].objc_flags and $
		sdss_flagval('object1','interp')) ne 0)
	help,winterp_sweep
	winterp_removed = where((sweep[sweep_removed_by_iconc].objc_flags and $
		sdss_flagval('object1','interp')) ne 0)
	help,winterp_removed

end




pro esboss::compare_concentration1755, $
		sweep, sweep_iconc, fpobjc_iconc
	; compare outputs from 
	;  standard gal datasweep run
	;  adding i concentration > 0.3 to datasweep run
	;  running on fpObjc with i concentration > 0.3

	dir='/clusterfs/riemann/raid006/bosswork/groups/boss/target/esheldon/sgc_test2/compare-concentration'
	sweep_file= $
		path_join(dir,'bosstarget-lrg-sgc-collate-anyresolve-001755.fits')
	sweep_iconc_file= $
	  path_join(dir,'bosstarget-iconc-lrg-sgc-collate-anyresolve-001755.fits')
	fpObjc_iconc_file= $
		path_join(dir,$
			'bosstarget-iconc-fpobjc-lrg-sgc-collate-anyresolve-001755.fits')

	if n_elements(sweep) eq 0 then begin
		sweep=mrdfits(sweep_file, 1)
		sweep_iconc = mrdfits(sweep_iconc_file,1)
		fpObjc_iconc = mrdfits(fpObjc_iconc_file,1)
	endif


	wsweep = self->btselect(sweep, 'lrg')
	wsweep_iconc = self->btselect(sweep_iconc, 'lrg')
	wfpobjc_iconc = self->btselect(fpobjc_iconc, 'lrg')


	print,'matching sweep to sweep_iconc'
	sphoto_match, sweep[wsweep], sweep_iconc[wsweep_iconc], $
		m_s2si, m_si2s
	help,wsweep,wsweep_iconc,m_s2si, m_si2s

	print
	print,'matching sweep_iconc to fpobjc_iconc'
	sphoto_match, sweep_iconc[wsweep_iconc], fpobjc_iconc[wfpobjc_iconc], $
		m_si2fi, m_fi2si
	help,wsweep,wfpobjc_iconc,m_si2fi, m_fi2si


	; get the objects thrown out from the sweep run when demanding
	; iconc > 0.3
	; remove the matches
	unmatched_sweep = wsweep
	remove, m_s2si, unmatched_sweep
	help, unmatched_sweep

	html_file = path_join(dir, 'sweep-removed-by-iconc-001755.html')
	print,'Creating html: ',html_file
	sdss_fchart_table, $
		sweep[unmatched_sweep].ra, sweep[unmatched_sweep].dec, $
		html_file, scale=0.396/2

	missed_fits=path_join(dir, 'sweep-removed-by-iconc-001755.fits')
	print,'Writing file: ',missed_fits
	mwrfits, sweep[unmatched_sweep], missed_fits, /create
end


pro esboss::test_bermuda_triangle, frac, struct=struct, mst=mst, $
		commissioning=commissioning, dops=dops

	; test the bermuda triangle.  Sparse sample it and see how it looks

	;commissioning=0

	bt=obj_new('bosstarget')
	bl=obj_new('bosstarget_lrg', commissioning=commissioning)
	lim=bl->limits()

	bt->runlist, runs, reruns

	;run=1755
	;rerun=137

	runs=runs[0:10]
	reruns=reruns[0:10]

	target_run='2009-06-10-comm'

	psdir = bt->target_dir(target_run)
	psdir = path_join(psdir, 'plots')

	if keyword_set(dops) then begin
		psfile=target_run+'-test-missed-triangle'
		if keyword_set(commissioning) then psfile+='-comm'
		psfile+='.eps'
		psfile=path_join(psdir,psfile)
		begplot,psfile,/color,/land;,/encap
		!p.charsize=1
		;!p.thick=2
		;!x.thick=2
		;!y.thick=2
	endif

	nstruct=n_elements(struct)
	if nstruct eq 0 then begin
		struct=bt->gather('lrg',target_run,runs=runs,/everything,/collate,/fast)
		nstruct=n_elements(struct)

		mst=bl->magstruct(struct)

		; the commissioning here doesn't matter since we don't use those
		; params
		logic=bl->flag_logic(struct, mst, commissioning=commissioning)
		logic = logic and mst.fibermag[3,*] gt lim.ilow_fib $
			and mst.fibermag[3,*] lt lim.ihi_fib
		logic = logic $
			AND (mst.cmodelmag[3,*] GT lim.ilow) $
			AND (mst.cmodelmag[3,*] LT lim.ihi)

		; this is just to keep the list small
		logic = logic and (mst.dperp gt -0.1) and (mst.dperp lt 1.5)

		w=where(logic)
		struct=struct[w]
		mst=mst[w]

	endif

	if !d.name eq 'X' then begin
		!p.background = c2i('white')
		!p.color = c2i('black')
	endif else begin
	endelse

	cmass_color = 'red'
	;triangle_color = 'blue'
	sparse_color ='blue'

	;pplot, /iso, $
	plotrand, /iso, $
		mst.cmodelmag[3], mst.dperp, psym=3, $
		ystyle=3, xstyle=3, $
		xtitle='cmodelmag[3]', ytitle='dperp'

	hiz_logic = bl->target_logic(mst, 'gal_hiz', $
		commissioning=commissioning)
	cmass_logic = bl->target_logic(mst, 'gal_cmass', $
		commissioning=commissioning)

	wcmass = where(cmass_logic)
	pplot, /overplot, $
		mst[wcmass].cmodelmag[3], mst[wcmass].dperp, psym=3, color=cmass_color




	; now define the "missing triangle"
	wtriangle = where(hiz_logic and not cmass_logic, ntriangle)
	if !d.name eq 'X' then color='green' else color='blue'
	;pplot, /overplot, $
	;	mst[wtriangle].cmodelmag[3], mst[wtriangle].dperp, psym=3, $
	;	color=triangle_color


	; sparse sample it with fraction "frac"
	nfrac = long(frac*ntriangle)
	r = randomu(seed, ntriangle)
	s=sort(r)
	wsparse = wtriangle[s[0:nfrac-1]]

	pplot, /overplot, $
		mst[wsparse].cmodelmag[3], mst[wsparse].dperp,  $
		color=sparse_color, psym=3

	mess = $
		['gal_cmass', $
		'Triangle sampling: '+string(frac,f='(f0.2)')]
	legend,mess,psym=8, color=[cmass_color, sparse_color]
	if keyword_set(commissioning) then begin
		legend,'Commissioning', /right
	endif


	if keyword_set(dops) then endplot,/landfix

end

pro esboss::lrg_cmodelmag_vs_dperp, struct, fracuse=fracuse, $
		psfile=psfile

	if n_elements(psfile) ne 0 then begin
		begplot,psfile,/landscape,/color
	endif
	
	plotrand, struct.cmodelmag[3], struct.dperp, psym=3, $
		fracuse=fracuse, indices=indices, $
		xrange=[12,20.1], xstyle=1, $
		yrange=[0.2, 1.8], ystyle=1, $
		xtitle='cmodelmag[3]', ytitle='dperp'

	wloz=self->btselect(struct[indices], 'gal_loz', count=nloz)
	wcmass=self->btselect(struct[indices], 'gal_cmass', count=ncmass)
	wtriangle=self->btselect(struct[indices], 'gal_triangle', count=ntriangle)
	wlodperp=self->btselect(struct[indices], 'gal_lodperp', count=nlodperp)

	wloz=indices[wloz]
	wcmass=indices[wcmass]
	wtriangle=indices[wtriangle]
	wlodperp=indices[wlodperp]

	cloz=c2i('green')
	ccmass=c2i('blue')
	ctriangle=c2i('red')
	clodperp=c2i('magenta')

	pplot, struct[wloz].cmodelmag[3], struct[wloz].dperp, psym=3, $
		color=cloz, /over
	pplot, struct[wcmass].cmodelmag[3], struct[wcmass].dperp, psym=3, $
		color=ccmass, /over
	pplot, struct[wtriangle].cmodelmag[3], struct[wtriangle].dperp, psym=3, $
		color=ctriangle, /over
	pplot, struct[wlodperp].cmodelmag[3], struct[wlodperp].dperp, psym=3, $
		color=clodperp, /over

	plegend, $
		['gal_loz','gal_cmass','gal_triangle','gal_lodperp'], $
		psym=8, $
		color=[cloz,ccmass,ctriangle,clodperp]

	if n_elements(psfile) ne 0 then begin
		endplot,/landfix
	endif
end

pro esboss::compare_sweep_versions, sold, snew

	common compare_flist_common, flistold, flistnew

	base=getenv('PHOTO_SWEEP_BASE')
	oldbase=path_join(base,'full_02apr06')
	newbase=path_join(base,'sgc_test2_uber')

	rerun='137'
	oldrerunbase=path_join(oldbase,rerun)
	newrerunbase=path_join(newbase,rerun)
	if n_elements(flistold) eq 0 then begin
		flistold = file_search(oldrerunbase,'calibObj*star.fits.gz')
		flistnew = file_search(newrerunbase,'calibObj*star.fits.gz')
	endif

	if n_elements(sold) eq 0 then begin
		; read some random set of runs
		match, file_basename(flistold), file_basename(flistnew), mfo, mfn, /sort
		rind = sort(randomu(seed,n_elements(mfo)))
		toread_old = flistold[mfo[rind[0:10]]]
		toread_new = flistnew[mfn[rind[0:10]]]

		sold=mrdfits_multi(toread_old)
		snew=mrdfits_multi(toread_new)
	endif

	spherematch, snew.ra, snew.dec, sold.ra, sold.dec, 1d/3600d, mn, mo, maxmatch=1

	begplot,'~/tmp/compare_sweep_versions-'+rerun+'.eps', xsize=8.5, ysize=8, /encap

	plothist, sold[mo].psfflux[1]/snew[mn].psfflux[1], bin=0.001,  $
		/ylog,yrange=[0.9, 1.2e4], xstyle=3, $
		xtitle='old_psfflux[1]/new_psfflux[1]'
	endplot
end


pro esboss::plot_ukidss_colors, struct
	bq=obj_new('bosstarget_qso')
	lups = bq->get_lups(struct, /deredden)

	g=reform(lups[1,*]) & i=reform(lups[3,*]) & k = struct.ukidss_k
	imk = i-k 
	gmi = g-i

	begplot,'~/tmp/ukidss_colors.eps', /encap, /color


	w=where(struct.ukidss_id ge 0)
	wc=self->btselect(struct[w], 'qso_core')
	wb=self->btselect(struct[w], 'qso_bonus')

	wc=w[wc]
	wb=w[wb]

	starflag = sdss_flagval('boss_target1','qso_ukidss')
	wb_notstars = where( (struct[wb].boss_target1 and starflag) eq 0, $
		complement=wb_stars)
	wb_notstars = wb[wb_notstars]
	wb_stars = wb[wb_stars]

	psym=8
	symsize=0.5
	ccolor='blue'
	bcolor_notstar = 'darkgreen'
	bcolor_star = 'red'

	pplot, [0], /nodata, aspect=1, $
		xrange=[0,4.5], yrange=[-0.5,2], xstyle=3, ystyle=3, $
		xtitle='i-K', ytitle='g-i'

	pplot, imk[wc], gmi[wc], $
		psym=psym, symsize=symsize, /overplot, $
		color=ccolor
	pplot, imk[wb_notstars], gmi[wb_notstars], $
		psym=psym, symsize=symsize, /overplot, $
		color=bcolor_notstar
	pplot, imk[wb_stars], gmi[wb_stars], $
		psym=psym, symsize=symsize, /overplot, $
		color=bcolor_star

	plegend, ['qso_core','qso_bonus (not star)','qso_bonus (star)'], $
		/right, psym=8, color=[ccolor, bcolor_notstar, bcolor_star]

	endplot
end

pro esboss::write_all_sweep_radec
	bt=obj_new('bosstarget')
	bt->runlist, runs, reruns

	bq=obj_new('bosstarget_qso')
	pars=bq->pars()

	glim=22.0 + 0.1
	rlim=21.85 + 0.1

	outdir=getenv('BOSS_TARGET')
	outdir=path_join(outdir,['esheldon','sgc_allpos'])
	outfile=path_join(outdir,'sdss_sgc_allpos.st')
	splog,'Will write to file: ',outfile,form='(a,a)'

	if file_test(outfile) then begin
		file_delete, outfile
	endif

	stdef={run:0,rerun:0,camcol:0,field:0,id:0,ra:0d,dec:0d}
	primary=sdss_flagval('resolve_status','survey_primary')
	types = ['star','gal']
	for i=0L, n_elements(runs)-1 do begin
		for camcol=1,6 do begin
			
			for ti=0L, n_elements(types)-1 do begin
				splog,'run=',runs[i],' rerun=',reruns[i],' camcol=',camcol,$
					' type="',type,'"',form='(a,i0,a,i0,a,i0,a,a,a)'
				type=types[ti]
				t=sweep_readobj(runs[i],camcol,rerun=reruns[i],type=type)
				w=where((t.resolve_status and primary) ne 0, nw)
				;lups=bq->get_lups(t, /deredden)	
				;w=where((t.resolve_status and primary) ne 0 and $
				;	((lups[1,*] le glim) or (lups[2,*] le rlim)), nw)
				if nw ne 0 then begin
					splog,'  writing: ',nw,form='(a,i0)'
					outst = replicate(stdef, nw)
					struct_assign, t[w], outst, /nozero
					write_idlstruct, outst, outfile, /append
				endif
			endfor
		endfor
	endfor
end

pro esboss::trim_allpos_to_stripe82, limited=limited, t=t
	outdir=getenv('BOSS_TARGET')
	outdir=path_join(outdir,['esheldon','sgc_allpos'])
	if keyword_set(limited) then begin
		fend='_maglim'
	endif else begin
		fend=''
	endelse

	infile=path_join(outdir,'sdss_sgc_allpos'+fend+'.st')
	outfile=path_join(outdir,'sdss_sgc_allpos_stripe82'+fend+'.fits')

	if n_elements(t) eq 0 then begin
		print,'Reading: ',infile
		stop
		t=read_idlstruct(infile)
	endif
	eq2csurvey, t.ra, t.dec, clambda, ceta

	bq=obj_new('bosstarget_qso')
	wdec=where(t.dec le 1.25 and t.ra ge -1.25, ndec)

	st=replicate({clambda:0d, ceta:0d, dec:0d}, ndec)
	st.clambda=clambda[wdec]
	st.ceta=ceta[wdec]

	bmask=bq->bounds_bitmask(st)
	w82=where(bmask eq 0, n82)
	t82=t[w82]
	print,'Writing: ',outfile
	mwrfits, t82, outfile, /create
end




pro esboss::check_nn_densities, target_run, $
		nn4=nn4, nn5=nn5, qso=qso, kde_coadd=kde_coadd

	bq=obj_new('bosstarget_qso')
	bt=obj_new('bosstarget')

	if n_elements(qso) eq 0 then qso=bt->read_collated('qso',target_run)
	if n_elements(nn4) eq 0 then begin
		nn4=bq->nn_read(version='V4')
		nn5=bq->nn_read(version='V5')
	endif
	if n_elements(kde_coadd) eq 0 then kde_coadd=bq->kde_coadd_read()

	qc=self->btselect(qso, 'qso_core')
	qb=self->btselect(qso, 'qso_bonus')

	rad=2d/3600d

	spherematch, nn4.ra, nn4.dec, kde_coadd.ra, kde_coadd.dec, rad,$
		m_nn4_to_kde, m_kde_to_nn4, maxmatch=1
	spherematch, nn5.ra, nn5.dec, kde_coadd.ra, kde_coadd.dec, rad,$
		m_nn5_to_kde, m_kde_to_nn5, maxmatch=1

	spherematch, nn4.ra, nn4.dec, qso[qc].ra, qso[qc].dec, rad,$
		m_nn4_to_core, m_core_to_nn4, maxmatch=1
	spherematch, nn5.ra, nn5.dec, qso[qc].ra, qso[qc].dec, rad,$
		m_nn5_to_core, m_core_to_nn5, maxmatch=1

	spherematch, nn4.ra, nn4.dec, qso[qb].ra, qso[qb].dec, rad,$
		m_nn4_to_bonus, m_bonus_to_nn4, maxmatch=1
	spherematch, nn5.ra, nn5.dec, qso[qb].ra, qso[qb].dec, rad, $
		m_nn5_to_bonus, m_bonus_to_nn5, maxmatch=1

	combined_matches4 = $
		[m_nn4_to_kde, m_nn4_to_core, m_nn4_to_bonus]
	combined_matches4=combined_matches4[rem_dup(combined_matches4)]
	combined_matches5 = $
		[m_nn5_to_kde, m_nn5_to_core, m_nn5_to_bonus]
	combined_matches5=combined_matches5[rem_dup(combined_matches5)]

	nn4_density = self->bosstile_density(nn4, tiles=tiles, /median)
	nn4_kde_density = self->bosstile_density(nn4[m_nn4_to_kde], /median)
	nn4_core_density = self->bosstile_density(nn4[m_nn4_to_core], /median)
	nn4_bonus_density = self->bosstile_density(nn4[m_nn4_to_bonus], /median)
	nn4_combined_density = self->bosstile_density(nn4[combined_matches4], /median)

	nn5_density = self->bosstile_density(nn5, tiles=tiles, /median)
	nn5_kde_density = self->bosstile_density(nn5[m_nn5_to_kde], /median)
	nn5_core_density = self->bosstile_density(nn5[m_nn5_to_core], /median)
	nn5_bonus_density = self->bosstile_density(nn5[m_nn5_to_bonus], /median)
	nn5_combined_density = self->bosstile_density(nn5[combined_matches5], /median)


	f='(a-30,i-10,a-30,g)'

	print,'# nn4: ',n_elements(nn4),' density: ',nn4_density,f=f
	print,'# nn4 kde coadd matches: ',n_elements(m_nn4_to_kde),' remaining density: ',nn4_density-nn4_kde_density,f=f
	print,'# nn4 core matches: ',n_elements(m_nn4_to_core),' remaining density: ',nn4_density-nn4_core_density,f=f
	print,'# nn4 bonus matches: ',n_elements(m_nn4_to_bonus),' remaining density: ',nn4_density-nn4_bonus_density,f=f
	print,'# nn4 net matches: ',n_elements(combined_matches4),' remaining density: ',nn4_density-nn4_combined_density,f=f

	print,'# nn5: ',n_elements(nn5),' density: ',nn5_density,f=f
	print,'# nn5 kde coadd matches: ',n_elements(m_nn5_to_kde),' remaining density: ',nn5_density-nn5_kde_density,f=f
	print,'# nn5 core matches: ',n_elements(m_nn5_to_core),' remaining density: ',nn5_density-nn5_core_density,f=f
	print,'# nn5 bonus matches: ',n_elements(m_nn5_to_bonus),' remaining density: ',nn5_density-nn5_bonus_density,f=f
	print,'# nn5 net matches: ',n_elements(combined_matches5),' remaining density: ',nn5_density-nn5_combined_density,f=f


	print
end

;+
; NAME:
;   bosstile_sub_regions
; PURPOSE:
;   Break a mangle window up into regions for QSO targeting
; CALLING SEQUENCE:
;   qso_regions, geometry, polygons=, regions= [, /plot]
; INPUTS:
;   geometry - [Ngeom] geometry defining area to consider
;   polygons - [Npoly] polygons for describing regions
;   regions - [Nregion] list of regions, with tags:
;                 .REGION - index number (starting with zero)
;                 .AREA - area in sq deg of region
;                 .PSTART - first polygon in polygons associated with region
;                 .PEND - last polygon in polygons associated with region
; OPTIONAL KEYWORDS:
;   /doplot - bring up a plot window for user to inspect
; COMMENTS:
;   The code tries to break the geometry up into regions about 5.5 sq
;     deg in size, defined by the SDSSPix 4th level pixelization. Any
;     leftover at the edge that are < half that area are attached to a
;     larger neighbor.
;   The "regions" structure is the final list of regions to use for 
;     any averaging.  The "polygons" structure contains the polygons
;     defining the regions (a region may be described by more than one
;     polygon)
;   The REGION tag in polygons structure acts as an index into
;     "regions", indicating which region the polygon belongs to.
;   In addition, polygons is returned sorted by region number, and 
;     PSTART and PEND in the regions structure indicate the range in
;     the polygons array the region polygons can be found in.
;   Therefore, if you want to check whether a set of RAs and Decs are 
;     in region number "i" you can execute:
; 
;     isin= is_in_window(polygons[regions[i].pstart:regions[i].pend], $
;                        ra= ra, dec= dec)
;
; REVISION HISTORY:
;   1-Dec-2009 MRB, NYU 
;-

function esboss::bosstile_sub_regions_cache_file, chunks, maxres

	boss_target = getenv('BOSS_TARGET')
	dir=filepath(root=boss_target,sub='esheldon','region-cache')

	if not file_test(dir) then begin
		file_mkdir, dir
	endif

	fname=['region','cache']

	for i=0L, n_elements(chunks)-1 do begin
		chstr_tmp=strtrim(string(chunks[i]),2)

		fname = [fname, chstr_tmp]
	endfor

	fname = [fname, 'maxres'+strtrim(string(maxres),2)]
	fname = strjoin(fname,'-')+'.ply'

	return, filepath(root=dir, fname)

end

function esboss::bosstile_sub_regions_cache_read, chunks, maxres
	fname = self->bosstile_sub_regions_cache_file(chunks, maxres)
	print,'reading sub region polygon cache: ',fname
	read_mangle_polygons, fname, polygons
	return, polygons
end

pro esboss::calculate_bosstile_sub_regions, $
		chunks=chunks, $
		geometry=geometry, $
		maxres=maxres, $
		$
		polygons=polygons, $
		regions=regions, $
		doplot=doplot


	print,'reading bounds poly for chunks(s): ',chunks
	bt=obj_new('bosstarget')
	;geometry = self->bosstile_poly_read(chunks, /bounds,/verbose)
	geometry = bt->chunk_polygon_read(chunks,/bounds,/verbose)

	if n_elements(maxres) eq 0 then maxres=5
	maxres=long64(maxres)

	print,'Using maxres=',maxres



	; first check for cache
	cache_file = self->bosstile_sub_regions_cache_file(chunks, maxres)
	if file_test(cache_file) and not keyword_set(force) then begin
		polygons = self->bosstile_sub_regions_cache_read(chunks, maxres)
	endif else begin

		full= 5.5*(4.0/maxres)^2 ;; full area of each pixel is 5.5 sq deg

		;; do initial pixelize 
		write_mangle_polygons, '/tmp/geom.ply', geometry
		geometry.weight=1.
		spawn, 'snap /tmp/geom.ply /tmp/snap.ply'
		spawn, 'balkanize /tmp/snap.ply /tmp/balkan.ply'
		spawn, 'unify /tmp/balkan.ply /tmp/unif.ply'
		; pixelize 
		;  d=sdspix
		;  0 polygons per pixel?
		;  4 max resolution of 4
		command=string('pixelize -Pd0,',maxres,$
					   ' /tmp/unif.ply /tmp/pixel.ply', $
					   f='(a,i0,a)')
		spawn, command
		;spawn, 'pixelize -Pd0,4 /tmp/unif.ply /tmp/pixel.ply'
		spawn, 'unify /tmp/pixel.ply /tmp/init.ply'
		read_mangle_polygons, '/tmp/init.ply', init
		init.weight= findgen(n_elements(init))+1.

		;; now attach small polygons to nearest large one
		area= init.str*(180./!DPI)^2
		xx=vmid(init)
		x_to_angles, xx, th, phi
		ra= th
		dec= 90.-phi
		ilo= where(area lt 0.5*full, nlo)
		ihi= where(area ge 0.5*full, nhi)
		if(nlo gt 0 and nhi gt 0) then begin
			for i=0L, nlo-1L do begin
				spherematch, $
					ra[ilo[i]], dec[ilo[i]], ra[ihi], dec[ihi], 5., m1, m2
				if(m1[0] ge 0) then $
					init[ilo[i]].weight= init[ihi[m2]].weight
			endfor
		endif
		write_mangle_polygons, '/tmp/join.ply', init
		spawn, 'snap /tmp/join.ply /tmp/snap.ply'
		spawn, 'balkanize /tmp/snap.ply /tmp/balkan.ply'
		spawn, 'unify /tmp/balkan.ply /tmp/pixel.ply'
		read_mangle_polygons, '/tmp/pixel.ply', polygons

		print,'writing cache file: ',cache_file
		write_mangle_polygons, cache_file, polygons
	endelse

	weights= polygons.weight
	regstr= replicate({region:-1L}, n_elements(polygons))
	isort= sort(weights)
	iuniq= uniq(weights[isort])
	regtmp={region:-1L, area:0.D, pstart:-1L, pend:-1L}
	regions= replicate(regtmp, n_elements(iuniq))
	ist=0L
	for i=0L, n_elements(iuniq)-1L do begin
		ind= iuniq[i]
		wcurr= isort[ist:ind]
		regstr[wcurr].region= i
		regions[i].region=i
		regions[i].area=total(polygons[wcurr].str)*(180.D/!DPI)^2
		regions[i].pstart= ist
		regions[i].pend= ind
		ist= ind+1L
	endfor

	polygons= struct_addtags(polygons, 'region', '-1L')
	polygons.region = regstr.region
	polygons.weight=1.
	polygons= polygons[isort]

	seed=11
	colors= randomu(seed, long(max(regions.region))+1L)*200+40
	color= colors[long(polygons.region)]

	if keyword_set(doplot) then begin
		plot_poly, geometry, /fill, color=255, offset=100.
		plot_poly, polygons, color=color, offset=100., /over, /fill
		if n_elements(init) ne 0 then begin
			plot_poly, init, out=1, color=180, offset=100., /over
		endif
		plot_poly, polygons, out=3, color=80, offset=100., /over
	endif

end


pro esboss::cache_bosstile_sub_regions, chunks, maxres

	common bosstile_subreg_density_block, $
		chunks_cache, maxres_cache, $
		polygons, regions, geometry

	recache = 0

	ncache = n_elements(chunks_cache)
	if ncache eq 0 then begin
		recache = 1
	endif else begin

		if maxres ne maxres_cache then begin
			; maxres has changed: recache
			recache = 1
		endif else begin
			if n_elements(chunks) ne ncache then begin
				recache = 1
			endif else begin
				match, chunks, chunks_cache, minput, mcache, /sort
				if minput[0] eq -1 or n_elements(minput) ne ncache then begin
					recache = 1
				endif
			endelse
		endelse
	endelse

	if recache then begin
		if n_elements(polygons) ne 0 then destruct_polygon, polygons
		if n_elements(geometry) ne 0 then destruct_polygon, geometry
		delvarx, polygons, geometry

		print,'cacheing sub regions'
		self->calculate_bosstile_sub_regions, $
			chunks=chunks, polygons=polygons, regions=regions, $
			maxres=maxres, geometry=geometry
	endif

end

pro esboss::bosstile_sub_regions, chunks, maxres, $
		polygons_out, regions_out, geometry_out


	common bosstile_subreg_density_block, $
		chunks_cache, maxres_cache, $
		polygons, regions, geometry

	; memory manageded internally
	self->cache_bosstile_sub_regions, chunks, maxres

	polygons_out = polygons
	regions_out = regions
	geometry_out = geometry
end



; 
; shirley's "ored" cuts
;

function esboss::orcut_table, type=type
	; a table with various "target densities".  They
	; don't quite match up to reality so we will explore
	;file=expand_path('~/idl.lib/data/shirley-cuts.st')
	file=expand_path('~/idl.lib/data/orcuts-table.st')
	cuts = read_idlstruct(file)
	if  n_elements(type) ne 0 then begin
		case type of
			'high': w=where(cuts.target_density gt 30)
			'low':  w=where(cuts.target_density lt 30)
			else: message,'Send high or low'
		endcase
		cuts=cuts[w]
	endif
	return, cuts
end

function esboss::orcut_logic, str, type=type, cuts=cuts
	if n_elements(cuts) eq 0 then begin
		cuts=self->orcut_table(type=type)
		w=where(cuts.target_density eq 60)
		cuts = cuts[w]
	endif
	return, $
		(str.kde_prob gt cuts.kde_prob_cut) $
		or (str.like_ratio gt cuts.like_ratio_cut) $
		or (str.nn_xnn gt cuts.nn_xnn_cut)
end
function esboss::orcut_select, str, nw, type=type, cuts=cuts
	logic = self->orcut_logic(str, cuts=cuts, type=type)
	w=where(logic, nw)
	return, w
end

function esboss::orcut_select_all, str, count, cuts=cuts, type=type
	; also include special logic for known/first
	good = lonarr(n_elements(str))

	; always include known and first
	flags=sdss_flagval('boss_target1',['qso_known_midz','qso_first_boss'])
	w=where( (str.boss_target1 and flags) ne 0, nw)
	if nw ne 0 then begin
		good[w] = 1
	endif

	worcut = self->orcut_select(str, norcut, cuts=cuts, type=type)
	if norcut ne 0 then begin
		good[worcut] = 1
	endif

	return, where( good, count )
end

function esboss::fit_and_eval, x, y, degree, x2eval, yfit=yfit
	; fit a polynomial and evaluate it at the requested x value
	pfit=poly_fit(x, y, degree)
	yfit = poly(x, pfit)
	yval=poly(x2eval, pfit)

	return, yval
end

pro esboss::orcut_plot, type, psfile=psfile, interpval=interpval
	cuts = self->orcut_table(type=type)

	if n_elements(psfile) ne 0 then begin
		begplot, psfile,/color,/encap,xsize=8.5,ysize=11
		;!p.charsize=1
	endif

	erase & multiplot, [1,3], /square
	pplot, cuts.target_density, cuts.kde_prob_cut, $
		yticklen=0.04, $
		psym=-8, xstyle=3, ystyle=3, ytitle='kde prob cut'

	if n_elements(interpval) ne 0 then begin
		yval=self->fit_and_eval($
			cuts.target_density, cuts.kde_prob_cut, 3, interpval, yfit=yfit)

		pplot, cuts.target_density, yfit, color='red',/over
		pplot, [0,interpval],[yval,yval],line=2,color='darkgreen',/over
		pplot, [interpval,interpval],[0,yval],line=2,color='darkgreen',/over
		plegend, string(f='("cut_val: ",f0.5)',yval),/bottom,/left
	endif


	multiplot
	pplot, cuts.target_density, cuts.like_ratio_cut, $
		yticklen=0.04, $
		psym=-8, xstyle=3, ystyle=3, ytitle='like_ratio cut'

	if n_elements(interpval) ne 0 then begin
		;if type eq 'high' then begin
		if 1 then begin
			order=2
			yval=self->fit_and_eval($
				cuts.target_density, cuts.like_ratio_cut, order, interpval, $
				yfit=yfit)

			pplot, cuts.target_density, yfit, color='red',/over
		endif else begin
			yval = interpol(cuts.like_ratio_cut, cuts.target_density, interpval)
		endelse
		pplot, [0,interpval],[yval,yval],line=2,color='darkgreen',/over
		pplot, [interpval,interpval],[0,yval],line=2,color='darkgreen',/over
		plegend, string(f='("cut_val: ",f0.3)',yval),/bottom,/left
	endif



	multiplot
	pplot, cuts.target_density, cuts.nn_xnn_cut, $
		yticklen=0.04, $
		psym=-8, xstyle=3, ystyle=3, ytitle='nn_xnn cut', $
		xtitle='original density'

	if n_elements(interpval) ne 0 then begin
		if type eq 'high' then order=3 else order=2
		yval=self->fit_and_eval($
			cuts.target_density, cuts.nn_xnn_cut, order, interpval, yfit=yfit)

		pplot, cuts.target_density, yfit, color='red',/over
		pplot, [0,interpval],[yval,yval],line=2,color='darkgreen',/over
		pplot, [interpval,interpval],[0,yval],line=2,color='darkgreen',/over
		plegend, string(f='("cut_val: ",f0.4)',yval),/bottom,/left
	endif




	multiplot,/reset
	if n_elements(psfile) ne 0 then endplot


	return

end

pro esboss::explore_orcuts, type, full=full
	; explore the ored cuts shirley sent on a subset of the area
	; type is 'high' for the area around 60 and 'low' for around 20
	chunk=2
	target_run = '2009-12-10-newlike2-norank'

	; get her cuts
	cuts = self->orcut_table(type=type)
	if type eq 'high' then begin
		goal=60.0 
		range=[35,90]
	endif else begin
		goal=20.0
		range=[13,26]
	endelse

	; get chunk2 data
	bt=obj_new('bosstarget')
	str=bt->read_collated('qso',target_run)

	; plotting
	targdir=filepath( $
		root=getenv('BOSS_TARGET'), subdir=[target_run,'plots'], 'tune-orcuts')
	if not file_test(targdir) then file_mkdir, targdir


	psfile=string(f='(a,"-chunk",i0,"-tune-orcut-",a,".eps")',target_run,chunk,type)
	dfile=string(f='(a,"-chunk",i0,"-densmap-",a,".eps")',target_run,chunk,type)
	psfile=filepath(root=targdir,psfile)
	dfile=filepath(root=targdir,dfile)

	; select everything that we called a quasar in preselection (this 
	; probably won't make any cuts actually)
	flags=$
		['qso_nn',$
		 'qso_like',$
		 'qso_kde', $
		 'qso_known_midz',$
		 'qso_first_boss']
	wall=self->btselect(str, flags)
	str=str[wall]

	self->plot_subreg_density,str.ra,str.dec,chunk,dfile

	begplot, psfile, /encap, /color, xsize=8.5,ysize=8.0
	!p.charsize=2

	ncuts = n_elements(cuts)
	actual_density = fltarr(ncuts)
	for i=0L, ncuts-1 do begin
		f='(a,f0.3)'
		print,'original density "goal" was: ',cuts[i].target_density, f=f
		print,'kde_prob_cut: ',cuts[i].kde_prob_cut, f=f
		print,'like_ratio_cut: ', cuts[i].like_ratio_cut, f=f
		print,'nn_xnn_cut: ',cuts[i].nn_xnn_cut,f=f
		print

		; This will pass known/first right through the cut
		w=self->orcut_select_all(str, cuts=cuts[i])

		; get density in smaller chunks
		ds=self->bosstile_subreg_density(str[w].ra,str[w].dec,chunk)

		; just use the flatter area
		wclam = where( ds.dstruct.clambda gt -41 and ds.dstruct.count gt 0 $
			and ds.dstruct.area gt 0.2)

		count = total( ds.dstruct[wclam].count )
		area = total( ds.dstruct[wclam].area )

		;if keyword_set(full) then begin
		;	self->compare_densities, chunk, 'qso', run, $
		;		outdir='/home/esheldon/tmp', str=str[w]
		;endif

		actual_density[i] = count/area
		print,'Total Density: ',actual_density[i]

		print,'----------------------------------------------------------'
	endfor

	pplot, cuts.target_density, actual_density, $
		aratio=1, $
		/ynozero, /iso, $
		xtitle='target density', ytitle='sub-region density', $
		xrange=range, yrange=range, xstyle=3, ystyle=3

	pplot, cuts.target_density, actual_density,/over, psym=8	
	pplot,[0,100],[0,100],/over, line=1

	lfit=poly_fit(cuts.target_density, actual_density,1)
	yfit = poly(cuts.target_density, lfit)
	pplot, cuts.target_density, yfit, color='red',/over


	val=interpol(cuts.target_density, yfit, goal)
	pplot, [0,val], [goal,goal], line=2, color='darkgreen', /over
	pplot, [val,val], [0,goal], line=2, color='darkgreen',/over

	plegend, $
		[string(f='("interp: ",f0.2)',val)]
	
	colprint,cuts.target_density,actual_density

	endplot
	
end



function esboss::gett, struct, name
	if strlowcase(name) eq 'ra' then begin
		c=shiftra(struct.ra, /wrap)
	endif else begin
		w=where(tag_names(struct) eq strupcase(name), nw)
		c=struct.(w)
	endelse
	return, c
end




pro esboss::compare_densities, $
		chunks, target_type, target_run, $
		reselect=reselect, $
		fixup=fixup, $
		extra_name=extra_name, $
		coord=coord, $
		xrange=xrange, yrange=yrange, $
		goal=goal, $
		outdir=outdir, $
		$
		str=str

	if n_params() lt 3 then begin
		on_error,2
		print,'Usage: '
		print,'  eb->compare_densities, chunks, target_type, target_run, '
		print,'     /reselect, /fixup, '
		print,'     extra_name=, coord=, outdir=, '
		print,'     xrange=, yrange=, '
		print,'     str='
		print
		message,'Halting'
	endif

	if n_elements(coord) eq 0 then coord='ra'
	if coord eq 'clambda' then csort=1 else csort=0
	if coord eq 'ra' then rasort=1 else rasort=0

	if n_elements(outdir) eq 0 then begin
		outdir=getenv('BOSS_TARGET')
		outdir=path_join(outdir, [target_run,'plots'])
	endif
	file_mkdir, outdir

	bt=obj_new('bosstarget')
	btypeobj=obj_new('bosstarget_'+target_type)
	if n_elements(str) eq 0 or keyword_set(reselect) then begin
		str=bt->read_collated(target_type, target_run,extra_name=extra_name)
	endif

	if keyword_set(reselect) then begin
		if chunks eq 2 then comm2=1
		if chunks eq 1 then commissioning=1
		str.boss_target1 = btypeobj->select(str, $
			comm2=comm2, commissioning=commissioning, /reselect)
	endif

	if keyword_set(extra_name) then begin
		ex='-'+extra_name
	endif else begin
		extra_name=''
		ex=''
	endelse

	chstr=strjoin(string(chunks,f='(i0)'), '-')
	psfile=path_join($
		outdir, $
		target_run+'-'+target_type+ex+'-chunk'+chstr+'-compare-densities.eps')
	if keyword_set(fixup) then psfile=repstr(psfile, '.eps', '-fixup.eps')
	if keyword_set(reselect) then psfile=repstr(psfile, '.eps', '-reselect.eps')



	; just to get the matches
	dstruct=$
		self->bosstile_density(str.ra,str.dec,chunks, $
			matched_objects=matched_objects)






	begplot,psfile,/color,/land

	if n_elements(xrange) eq 0 then begin
		;xrange=[-180,180]
		xv = self->gett(str[matched_objects],coord)
		xrange=[min(xv),max(xv)]
	endif

	case target_type of
		'qso': begin
			;types = [$
			;	'qso_core','qso_bonus','qso_like','qso_nn','qso_known_midz',$
			;	'qso_first_boss']
			types=['qso_nn','qso_like','qso_kde',$
				   'qso_core_main',$
				   'qso_bonus_main','qso_known_midz','qso_first_boss']
			if n_elements(yrange) eq 0 then yrange=[0,140]
		end
		'lrg': begin
			types = ['gal_loz', 'gal_cmass','gal_cmass_sparse']
			if n_elements(yrange) eq 0 then yrange=[0,240]
		end
		'std': begin
			types = ['std_fstar']
			if n_elements(yrange) eq 0 then yrange=[0,10]
		end
		else: message,"target_type must be 'qso','lrg','std'"
	endcase


	print,'xrange=',xrange
	print,'yrange=',yrange
	plot,[0],/nodata,  $
		xrange=xrange,yrange=yrange, xstyle=3, $
		xtitle=strupcase(coord), ytitle='#/square degree'

	if n_elements(goal) eq 0 then begin
		if target_type eq 'qso' then begin
			goal=60 
		endif else begin
			goal=155
			if n_elements(extra_name) ne 0 then begin
				if extra_name eq 'noknown' then begin
					goal=130
				endif
			endif
		endelse
	endif

	pplot, /over, [-180, 180], [goal,goal], line=2, color='grey75'

	ntypes=n_elements(types)

	nntot=0
	for i=0L, ntypes-1 do begin
		val=sdss_flagval('boss_target1',types[i])
		;ww=where(str.boss_target1 eq val, nn)
		ww=where((str.boss_target1 and val) ne 0, nn)
		if nn ne 0 then begin
			nntot+=1
		endif
	endfor
	nntot += nntot
	;colors=make_rainbow(11)
	colors=make_rainbow(nntot)
	colors=reverse(colors)
	icolor=0

	tilepoly=self->bosstile_poly_read(chunks)
	total_area=total(tilepoly.str)*(180d/!dpi)^2
	destruct_polygon, tilepoly

	print
	print,'Densities:   Total Area: ', total_area,' sq deg',f='(a,f0.2,a)'
	print,'--------------------------------------'
	;print,'target type','     count','    count',$
	;	form='(a-20,a-10,a-10)'


	for i=0L, ntypes-1 do begin
		name=types[i]

		case target_type of
			'lrg': begin
				; lrg have some unused types that we don't want to consider
				; so restrict to these types
				anylogic=self->btlogic(str, name, fixup=fixup, /silent)

				val=sdss_flagval('boss_target1', name)
				allvals = sdss_flagval('boss_target1', types)
				onlylogic = $
					(str.boss_target1 and val) ne 0 $
					and (str.boss_target1 and allvals) ne allvals
			end
			else: begin
				anylogic = $
					self->btlogic(str, name, fixup=fixup, /silent)
				onlylogic = $
					self->btlogic(str, name, fixup=fixup, /only, /silent)

			end

		endcase

		wany=where(anylogic, nany)
		wonly=where(onlylogic, nonly)

		dens_any=0.
		dens_only=0.

		if nonly ne 0 then begin
				dens_only_struct=$
					self->bosstile_density($
						str[wonly].ra,str[wonly].dec, chunks, $
						rasort=rasort, csort=csort, $
						matched_objects=only_matched_objects, $
						nmatched_objects=nonly, $
						mean_density=dens_only)
		endif

		if nany ne 0 then begin
			dens_any_struct=$
				self->bosstile_density($
					str[wany].ra,str[wany].dec, chunks, $
					rasort=rasort, csort=csort, $
					matched_objects=any_matched_objects, $
					nmatched_objects=nany, $
					mean_density=dens_any)
			;print,dens_any_struct.density
			c=self->gett(dens_any_struct,coord)
			s=sort(c)
			pplot, c[s], $
				dens_any_struct[s].density, /over, $
				color=colors[icolor]
			c=self->gett(dens_any_struct,coord)
			s=sort(c)
			pplot, c[s], $
				dens_any_struct[s].density, psym=8,/over, $
				color=colors[icolor]

			thisname=strlowcase(name)+': '+string(dens_any,f='(f0.1)')
			add_arrval, thisname, names
			add_arrval, icolor, icolorkeep
			add_arrval, 0, lines
			icolor+=1
		endif


		print,$
			name,$
			'any: ',$
			nany,$
			dens_any,$
			'only: ',$
			nonly,$
			dens_only, $
			form='(a-20,a,i-7,f5.1,"/sqdeg    ",a,i-7,f5.1,"/sqdeg")'
		
	endfor


	print
	print,'Composite Densities:   Total Area: ', total_area,' sq deg',f='(a,f0.2,a)'
	print,'--------------------------------------'


	cline=1


	case target_type of
		'qso': begin
			comp_orflags=types[0]
			types2add = types[1:n_elements(types)-1]
		end
		'lrg':  begin
			comp_orflags=types[0]
			types2add = types[1:n_elements(types)-1]
		end
		'std': begin
			; do nothing
		end
		else: message,"target_type must be 'qso','lrg','std'"
	endcase


	; add qso_bonus to this if it starts contributing
	for i=0L, n_elements(types2add)-1 do begin
		wone =self->btselect(str, types2add[i], count=count,/silent)
		if count gt 0 then begin
			add_arrval, types2add[i], comp_orflags


			w = self->btselect(str, comp_orflags, /silent, fixup=fixup)

			ds= self->bosstile_density($
					str[w].ra,str[w].dec,chunks, $
					rasort=rasort, $
					csort=csort,$
					nmatched_objects=nany, $
					mean_density=dens_any)

			c=self->gett(ds,coord)
			s=sort(c)
			pplot, c[s], ds[s].density, /over, $
				color=colors[icolor], line=cline
			pplot, c[s], ds[s].density, psym=8,/over, $
				color=colors[icolor]

			comp_name = strjoin(comp_orflags, ' or ')

			thisname=comp_name+': '+string(dens_any, f='(f0.1)')
			add_arrval, thisname, names
			add_arrval, icolor, icolorkeep
			add_arrval, cline, lines
			icolor+=1
			print,'    ', strjoin(comp_orflags, ' or ')
			print,nany,dens_any
		endif
	endfor

	plegend, $
		reverse(names), /right, $
		colors=reverse(colors[icolorkeep]), $
		psym=8,$
		symsize=replicate(1.2, n_elements(names)), $
		$;lines=reverse(lines), $
		charsize=0.9

	endplot,/landfix



	return


end

pro esboss::plot_chunk_poly, chunks, bounds=bounds, offset=offset, $
		xsize=xsize, ysize=ysize, $
		xrange=xrange, yrange=yrange, $
		_extra=_extra, $
		filename=filename

	if n_elements(filename) ne 0 then begin
		if n_elements(xsize) eq 0 then xsize=8.5
		if n_elements(ysize) eq 0 then ysize=8.0

		begplot,filename,xsize=xsize,ysize=ysize,/color,/encap
	endif

	nchunk=n_elements(chunks)
	for i=0L, nchunk-1 do begin
		ply = self->bosstile_poly_read(chunks[i], bounds=bounds)
		add_arrval, ply, ply_all
	endfor

	!x.style=3
	plot_poly, ply_all, outline_thick=2, offset=offset, color=!p.color, $
		xrange=xrange, yrange=yrange, $
		_extra=_extra
	!x.style=0

	if nchunk gt 1 then begin
		colors = make_rainbow(nchunk)
		for i=0L, nchunk-1 do begin
			ply = self->bosstile_poly_read(chunks[i], bounds=bounds)

			plot_poly, ply, outline_thick=2, offset=offset, color=colors[i],/over

			destruct_polygon, ply
		endfor
		plegend, string(chunks, f='("chunk",i0)'), color=colors, line=0, /right
	endif else begin
		plegend, string(chunks, f='("chunk",i0)'), /right
	endelse


	destruct_polygon, ply_all

	if n_elements(filename) ne 0 then begin
		endplot
	endif
end


function esboss::read_shirley_spec, matched=matched
	fall='/home/shirleyho/research/BOSS_spectra/qso_selection/data_files/BOSS_Quasars_1PCplus.fits'
	fmatch = '/home/shirleyho/research/BOSS_spectra/qso_selection/data_files/BOSS_Quasars_1PCplus.reducflds.fits'

	if not keyword_set(matched) then begin
		return,mrdfits(fall,1)
	endif else begin
		return,mrdfits(fmatch,1)
	endelse
end

pro esboss::plot_in_chunks, x, y, chunksize, over_index=over_index, $
		xtitle=xtitle, ytitle=ytitle, title=title, $
		tilepoly=tilepoly, offset=offset, $
		minx=minx, maxx=maxx,$
		_extra=_extra

	; plot in chunks of x, size chunksize, using the /iso keyword
	; to keep correct aspect

	if n_elements(minx) eq 0 then minx=min(x)
	if n_elements(maxx) eq 0 then maxx=max(x)
	
	chunkrat = (maxx-minx)/float(chunksize)
	nchunks = long64(chunkrat)
	if chunkrat gt nchunks then begin
		nchunks += 1
	endif

	;!p.charsize=1

	!p.multi=[0,0,nchunks]
	print,'nchunks=',nchunks,f='(a,i0)'

	;erase & multiplot, [1, nchunks], $
	;	mxtitle=xtitle, mytitle=ytitle;, ygap=0.03

	didtitle=0
	for i=0L, nchunks-1 do begin
		begx = minx+i*chunksize
		endx = minx+(i+1)*chunksize
		print,begx,endx

		w=where(x ge begx and x le endx, nw)
		if nw ne 0 then begin
			delvarx, tit
			if not didtitle then tit=title
			didtitle=1
			plot, x[w], y[w], psym=3, xrange=[begx,endx], xstyle=1,iso=1, $
				xtitle=xtitle, ytitle=ytitle, title=tit, _extra=_extra, $
				xticklen=0.04
		endif

		if n_elements(over_index) ne 0 then begin
			w=where(x[over_index] ge begx and x[over_index] le endx, nw)
			if nw ne 0 then begin
				w=over_index[w]
				pplot, x[w], y[w], psym=3, color='blue', /overplot
			endif
		endif

		if n_elements(tilepoly) ne 0 then begin
		plot_poly, tilepoly, /over,outline_thick=2, color=!p.color, $
			offset=offset;, xrange=[begx,endx]
		endif

		;multiplot
	endfor

	;multiplot, /reset
	!p.multi=0
end

pro esboss::plot_targets_and_tiles, $
		chunks, target_type, target_run, target_flags, $
		outdir=outdir, $
		str=str, $
		inwindow=inwindow, $
		observed=observed, $
		overstruct=overstruct, $
		minra=minra, maxra=maxra, $
		reselect=reselect, $
		extra_name=extra_name, _extra=_extra

	if n_params() lt 4 then begin
		on_error,2
		message,'usage: plot_targets_and_tiles, chunks, target_type, '+$
			'target_run, target_flags, outdir=, str=, /reselect, '+$
			'extra_name=, _extra='
	endif

	if n_elements(outdir) eq 0 then begin
		outdir=getenv('BOSS_TARGET')
		outdir=path_join(outdir, [target_run,'plots'])
	endif
	file_mkdir, outdir

	bt=obj_new('bosstarget')
	if n_elements(str) eq 0 or keyword_set(reselect) then begin
		str=bt->read_collated(target_type, target_run,extra_name=extra_name)
	endif
	if keyword_set(reselect) then begin
		if chunks eq 2 then comm2=1
		if chunks eq 1 then commissioning=1
		btypeobj=obj_new('bosstarget_'+target_type)
		str.boss_target1 = btypeobj->select(str, $
			comm2=comm2, commissioning=commissioning, /reselect)
	endif


	if target_type eq 'std' then begin
		all_psym=8
		intile_psym=8
		all_symsize=0.25
		intile_symsize=0.75 
	endif else if target_type eq 'qso' then begin
		all_psym=8
		all_symsize=0.1
		intile_psym=8
		intile_symsize=0.20
	endif else begin
		all_psym=3
		intile_psym=8
		intile_symsize=0.15
	endelse


	flagnames=strjoin(target_flags, '-')
	flagnames=repstr(flagnames, '_','-')

	fname=target_run+'-target-and-tile-radec-'+flagnames+'.eps'
	if n_elements(extra_name) ne 0 then begin
		fname=repstr(fname, '.eps', '-'+extra_name+'.eps')
	endif
	psfile=path_join(outdir, fname)


	w=self->btselect(str, target_flags)

	splog,'Checking inwindow'
	tilepoly=self->bosstile_poly_read(chunks)
	boundpoly=self->bosstile_poly_read(chunks, /bounds)

	ra  = str.ra
	dec = str.dec

	inwindow_overall = $
		where(is_in_window(boundpoly, ra=ra[w], dec=dec[w]), nwin)
	w=w[inwindow_overall]

	if n_elements(inwindow) eq 0 then begin
		inwindow=where(is_in_window(tilepoly, ra=ra[w], dec=dec[w]), nwin)
		;inwindow = w[inwindow]
	endif


	if chunks eq 1 then begin
		offset=90.0
		ra = shiftra(ra, offset)

		begplot, psfile, /encap, /color, ysize=11, xsize=8.5
		chunksize=20
		title = 'chunk '+strn(chunks)
		if n_elements(extra_name) ne 0 then title += ' "'+extra_name+'"'
		self->plot_in_chunks, ra[w],dec[w], chunksize, $
			xtitle='RA-90', ytitle='DEC', $
			title=title, $
			over_index=inwindow, $
			minx=minra, maxx=maxra, $
			tilepoly=tilepoly, offset=offset, $
			_extra=_extra

		endplot
		return
	endif

	begplot, psfile, /encap, /color, xsize=11, ysize=8.5
	plot, ra[w], dec[w], psym=all_psym, symsize=all_symsize, $
		xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, $
		xtitle=textoidl('RA'), ytitle=textoidl('DEC'), $
		/iso

	target_color='blue'
	pplot, /over, $
		ra[w[inwindow]], dec[w[inwindow]], $
		psym=intile_psym, symsize=intile_symsize, color=target_color, $
		_extra=_extra


	plot_poly, tilepoly, /over,outline_thick=2, color=!p.color, offset=offset
	destruct_polygon, tilepoly
	destruct_polygon, boundpoly

	tiles=self->bosstile_read(chunks)
	center_psym=7
	center_color='red'
	pplot, tiles.ra, tiles.dec, psym=7, color=center_color, /overplot

	messages = ['targets','tile centers']
	psym=[intile_psym, center_psym]
	colors=[target_color, center_color]

	if keyword_set(observed) then begin
		tiles=self->read_observed_plates(chunk=chunks)
		pplot, tiles.ra, tiles.dec, psym=center_psym, color='orange', /overplot

		messages=[messages,'observed tiles']
		psym=[psym, center_psym]
		colors=[colors,'orange']
	endif

	if n_elements(overstruct) ne 0 then begin
		spectra_psym=4
		spectra_color='darkgreen'
		pplot, overstruct.ra, overstruct.dec, /overplot, $
			psym=spectra_psym, color=spectra_color, thick=0.5, $
			_extra=_extra

		messages=[messages,'spectra']
		psym=[psym, spectra_psym]
		colors=[colors,spectra_color]
	endif

	;plegend, ['tiles','observed'], psym=7, color=['red','orange'], /right
	plegend, messages, psym=psym, color=colors,/right,/bottom
	endplot

end





function esboss::fstar_dist_cuts, nstep
	nstep=100
	cuts = arrscl(findgen(nstep), 0.04, 0.10)
	return, cuts
end
pro esboss::tune_fstar_std, struct=struct, density=density

	bt=obj_new('bosstarget')
	if n_elements(struct) eq 0 then begin
		struct=bt->read_collated('std','2009-07-30c')
	endif

	std = obj_new('bosstarget_std')
	lim=std->standards_limits()

	mags = 22.5 - 2.5*alog10( struct.psfflux > 0.001 ) - struct.extinction
	umg = reform( mags[0,*] - mags[1,*] )
	gmr = reform( mags[1,*] - mags[2,*] )
	rmi = reform( mags[2,*] - mags[3,*] )
	imz = reform( mags[3,*] - mags[4,*] )
	dist = reform( sqrt( $
		(umg-lim.umg_val)^2 $
		+ (gmr-lim.gmr_val)^2 $
		+ (rmi-lim.rmi_val)^2 $
		+ (imz-lim.imz_val)^2 ) )

	cuts = self->fstar_dist_cuts(nstep)

	density=dblarr(nstep)
	
	fstar_flag=sdss_flagval('boss_target1','std_fstar')
	frac_mindens = dblarr(nstep)
	for i=0L, nstep-1 do begin

		dcut = cuts[i]
		w=where( (struct.boss_target1 and fstar_flag) ne 0 $
			and dist lt dcut, nw)

		tdens_struct=self->bosstile_density(struct[w])

		
		wmindens = where(tdens_struct.count ge 16, nmindens)
		frac_mindens[i] = float(nmindens)/n_elements(tdens_struct)
		
		dens = median(tdens_struct.density)
		print, 'cut: ',dcut,' density: ',dens,f='(a,f0.4,a,f0.3)'


		add_arrval, tdens_struct, density_struct
		density[i] = dens
	endfor


	mindens = 16d/(!dpi*self->bosstile_radius()^2)
	begplot,'~/tmp/test_std_density.eps', /color, xsize=8.5, ysize=8


	cutval=interpol(cuts, density, mindens)

	pplot, cuts, density, xtitle='color distance', ytitle='density #/sq deg', $
		aspect=1
	pplot, [0,cutval],[mindens,mindens], /overplot,color='red'
	pplot, [cutval,cutval],[0,mindens], /overplot,color='red'

	endplot

	begplot,'~/tmp/test_std_density_frac.eps', /color, xsize=8.5, ysize=8
	pplot, cuts, frac_mindens, $
		xtitle='color distance', ytitle='Fraction with nstar >= 16', aspect=1
	endplot


end



pro esboss::mmt_match_fpobjc_pbs
	mmtdir = getenv("BOSSTARGET_DIR")
	mmtdir=path_join(mmtdir, "data")

	mmtname = path_join(mmtdir, "mmt-plate-data.dat")

	stdef = {ra:0d, dec:0d, rad:0d}
	read_struct, mmtname, stdef, plateinfo, skiplines=1

	outdir='~/pbs/mmt-fpobjc-match'
	for fieldid=0L, n_elements(plateinfo)-1 do begin
		fstr1 = string(fieldid, f='(i02)')
		pbs_file=path_join(outdir, 'mmt-fpobjc-match-'+fstr1+'.pbs')

		print,pbs_file

		fstr2=string(fieldid,f='(i0)')
		job_name='mmt-fpobjc-'+fstr1
		idl_commands="boss_mmt_match, 'star', /fpobjc, fieldid="+fstr2
		pbs_riemann_idl, pbs_file, idl_commands, setup=setups, job_name=job_name
	endfor

end

pro esboss::compare_like_new, str=str
	bt=obj_new('bosstarget')
	target_run = '2009-12-08-newlike-loose'
	if n_elements(str) eq 0 then str=bt->read_collated('qso',target_run)

	w=where(str.like_ratio gt 1.e-30 and str.like_ratio_old gt 1.e-30)

	loglike=alog10( str[w].like_ratio )
	loglike_old=alog10( str[w].like_ratio_old )

	psfile='like-compare-new.eps'
	begplot,psfile,/color,/encap,xsize=8.5,ysize=8
	plothist, loglike,  bin=0.1,max=3,xtitle='log(like_ratio)', $
		xrange=[-2,3],xstyle=1
	plothist, loglike_old, bin=0.1,max=3,/over,color='blue'
	plegend,['old','new'],line=0,color=[!p.color,c2i('blue')],/right
	endplot

	psfile='like-compare-new-zoom.eps'
	begplot,psfile,/color,/encap,xsize=8.5,ysize=8
	plothist, loglike,  bin=0.1,min=-2,max=3,xtitle='log(like_ratio)', $
		xrange=[-2,3],xstyle=1
	plothist, loglike_old, bin=0.1,min=-2,max=3,/over,color='blue'
	plegend,['old','new'],line=0,color=[!p.color,c2i('blue')],/right
	endplot
end




pro esboss::boss1_test, targs, lrg_intermed, qso_intermed, lrg, qso, std
	; test Michael's target lists compared to the input

	target_run='comm'

	dir='~/boss1-test'
	targfile=path_join(dir,'final-boss1.fits')
	lrg1file=path_join(dir,'lrg-boss1.fits')
	qso1file=path_join(dir,'qso-boss1.fits')
	std1file=path_join(dir,'std-boss1.fits')

	if n_elements(targs) eq 0 then begin
		targs=mrdfits(targfile,1)
		lrg_intermed=mrdfits(lrg1file,1)
		qso_intermed=mrdfits(qso1file,1)

		bt=obj_new('bosstarget')
		qso=bt->read_collated('qso',target_run)
		std=bt->read_collated('std',target_run)
		lrg=bt->read_collated('lrg',target_run, extra='noknown')


		; so we will get matches
		targs.sourcetype = strtrim(targs.sourcetype,2)
		targs.rerun = strtrim(targs.rerun,2)

		lrg.rerun=strtrim(lrg.rerun,2)
		lrg_intermed.rerun=strtrim(lrg_intermed.rerun,2)
		qso.rerun=strtrim(qso.rerun,2)
		qso_intermed.rerun=strtrim(qso_intermed.rerun,2)
	endif


	; first compare what we think was input to the intermediate file
	begplot,path_join(dir, 'lrg-boss1-test.ps'),/land, /color
	xtit='RA' & ytit='DEC'

	; compare to the intermediate file
	symsize=0.15
	pplot, shiftra(lrg.ra,/wrap), lrg.dec, psym=8, symsize=symsize, $
		xtit=xtit, ytit=ytit,$
		yrange=[-1.5,2], ystyle=3
	pplot, shiftra(lrg_intermed.ra,/wrap), lrg_intermed.dec, /over, $
		psym=8, symsize=symsize, $
		color='red'
	legend,['lrg '+target_run, 'lrg-boss1.fits'], $
		psym=8, color=[!p.color,c2i('red')]
	endplot, /landfix


	begplot,path_join(dir, 'qso-boss1-test.ps'),/land, /color
	symsize=0.25
	pplot, shiftra(qso.ra,/wrap), qso.dec, psym=8, symsize=symsize, $
		xtit=xtit, ytit=ytit,$
		yrange=[-1.5,2], ystyle=3
	pplot, shiftra(qso_intermed.ra,/wrap), qso_intermed.dec, /over, $
		psym=8, symsize=symsize, $
		color='red'
	legend,['qso '+target_run, 'qso-boss1.fits'], color=[!p.color,c2i('red')], $
		psym=8
	endplot,/landfix


	; compare to the final file
	sourcetype=strtrim(targs.sourcetype, 2)

	begplot,path_join(dir, 'lrg-final-boss1-test.ps'),/land, /color
	wlrg=where(strmatch(sourcetype,'LRG'))
	symsize=0.15
	pplot, shiftra(lrg.ra,/wrap), lrg.dec, psym=8, symsize=symsize, $
		xtit=xtit, ytit=ytit,$
		yrange=[-1.5,2], ystyle=3
	pplot, shiftra(targs[wlrg].ra,/wrap), targs[wlrg].dec, /over, $
		psym=8, symsize=symsize, $
		color='red'
	legend,['lrg '+target_run, 'lrg final-boss1.fits'], $
		psym=8, color=[!p.color,c2i('red')]
	endplot,/landfix


	begplot,path_join(dir, 'qso-final-boss1-test.ps'),/land, /color
	wqso=where(strmatch(sourcetype,'QSO'))
	symsize=0.25
	pplot, shiftra(qso.ra,/wrap), qso.dec, psym=8, symsize=symsize, $
		xtit=xtit, ytit=ytit,$
		yrange=[-1.5,2], ystyle=3
	pplot, shiftra(targs[wqso].ra,/wrap), targs[wqso].dec, /over, $
		psym=8, symsize=symsize, $
		color='red'
	legend,['qso '+target_run, 'qso final-boss1.fits'], $
		color=[!p.color,c2i('red')], $
		psym=8

	endplot, /landfix





	; compare structures
	lrg_cs = compare_struct(lrg, lrg_intermed)
	if n_elements(lrg_cs) ne 1 or lrg_cs[0].ndiff ne 0 then begin
		splog,'Differences found between lrg-'+target_run+' and lrg-boss1.fits'
		help,lrg_cs,/str
	endif else begin
		splog,'No differences found between lrg-'+target_run+' and lrg-boss1.fits'
	endelse

	sphoto_match, lrg, targs[wlrg], mlrg, mtarglrg
	lrg_final_cs = compare_struct(lrg[mlrg], targs[wlrg[mtarglrg]])
	if n_elements(lrg_final_cs) ne 1 or lrg_final_cs[0].ndiff ne 0 then begin
		splog,'Differences found between lrg-'+target_run+' and lrg targets in final-boss1.fits'
		help,lrg_final_cs,/str
	endif else begin
		splog,'No differences found between lrg-'+target_run+' and lrg targets in final-boss1.fits'
	endelse




	qso_cs = compare_struct(qso, qso_intermed)
	if n_elements(qso_cs) ne 1 or qso_cs[0].ndiff ne 0 then begin
		splog,'Differences found between qso-'+target_run+' and qso-boss1.fits'
		help,qso_cs,/str
	endif else begin
		splog,'No differences found between qso-'+target_run+' and qso-boss1.fits'
	endelse

	sphoto_match, qso, targs[wqso], mqso, mtargqso
	qso_final_cs = compare_struct(qso[mqso], targs[wqso[mtargqso]])
	if n_elements(qso_final_cs) ne 1 or qso_final_cs[0].ndiff ne 0 then begin
		splog,'Differences found between qso-'+target_run+' and qso targets in final-boss1.fits'
		help,qso_final_cs,/str
	endif else begin
		splog,'No differences found between qso-'+target_run+' and qso targets in final-boss1.fits'
	endelse



end






pro esboss::create_loosecoadd_pbs
	pars={kde_coadd_dotrim:0}

	self->create_pbs, $
		'qso', '2009-08-02c', pars=pars, $
		extra_name='loosecoadd',/comm, /dogather

end
pro esboss::create_qso_gather_pbs_reselect

	;bosstarget_v="-r /home/esheldon/exports/bosstarget-v1_0_2"
	bosstarget_v="-r /home/esheldon/exports/bosstarget-work"
	idlutils_v="-r /home/esheldon/exports/idlutils-v5_4_7"
	photoop_v="v1_9_4"
	extra_name = "with-ukidss"

	orflagnames=[$
		'qso_core',$
		'qso_kde_coadd',$
		'qso_known_midz',$
		'qso_nn',$
		'qso_like',$
		'qso_ukidss']

	flagvals=sdss_flagval('boss_target1',orflagnames)
	flagstr=string(flagvals,f='(i0)')

	where_string="((str.boss_target1 and 8192) eq 0) and ((str.boss_target1 and "+flagstr+") ne 0) and ((str.resolve_status and 256) ne 0)"

	self->create_pbs, $
		'qso', '2009-08-05c',  /comm, /dogather, /reselect, $
		bosstarget_v=bosstarget_v, $
		idlutils_v=idlutils_v, $
		photoop_v=photoop_v, $
		where_string=where_string, $
		extra_name=extra_name
end

pro esboss::create_lrg_gather_pbs_reselect_noknown, target_run


	bosstarget_v="-r /home/esheldon/exports/bosstarget-v1_0_4"
	;bosstarget_v="-r /home/esheldon/exports/bosstarget-work"
	idlutils_v="-r /home/esheldon/exports/idlutils-v5_4_7"
	photoop_v="v1_9_4"

	bt=obj_new('bosstarget')

	ws=bt->default_where_string('lrg')
	known=string(sdss_flagval('boss_target1','sdss_gal_known'),f='(i0)')
	ws += ' and ((str.boss_target1 and '+known+') eq 0)'
	self->create_pbs, $
		'lrg', target_run,  /comm, /dogather, /reselect, $
		bosstarget_v=bosstarget_v, $
		idlutils_v=idlutils_v, $
		photoop_v=photoop_v, $
		where_string=ws, $
		extra_name='noknown'
end


pro esboss::create_fullrun_pbs, target_type, target_run, $
		dogather=dogather
	;bosstarget_v="-r /home/esheldon/exports/bosstarget-v1_0_4"
	bosstarget_v="-r /home/esheldon/exports/bosstarget-work"
	idlutils_v="-r /home/esheldon/exports/idlutils-v5_4_7"
	photoop_v="v1_9_4"


	if target_type eq 'qso' and not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, target_type, target_run,  $
			/commissioning,  $
			bosstarget_v=bosstarget_v, $
			idlutils_v=idlutils_v, $
			photoop_v=photoop_v
	endif else begin
		self->create_pbs, target_type, target_run,  $
			dogather=dogather, $
			/commissioning,  $
			bosstarget_v=bosstarget_v, $
			idlutils_v=idlutils_v, $
			photoop_v=photoop_v
	endelse
end

pro esboss::create_gather_pbs_bitmask, target_run
	; 
	; collate a sample only cut by bitmask
	;
	if n_elements(target_run) eq 0 then begin
		message,'Please enter a target_run'
	endif
	if target_run eq 'comm' then begin
		commissioning=1
		PHOTO_SWEEP=$
			"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-06-14"
	endif else if target_run eq 'comm2' then begin
		comm2=1
		PHOTO_SWEEP=$
			"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"
	endif else begin
		message,'support target_run = '+target_run
	endelse
	target_type='qso'
	extra_name='bitmask'
	where_string='str.bitmask eq 0'

	bosstarget_v="-r /home/esheldon/exports/bosstarget-work"
	idlutils_v="-r /home/esheldon/exports/idlutils-v5_4_7"
	photoop_v="v1_9_4"

	self->create_pbs, target_type, target_run,  $
		where_string=where_string, $
		extra_name=extra_name, $
		/dogather, $
		commissioning=commissioning,  $
		comm2=comm2,$
		photo_sweep=photo_sweep, $
		bosstarget_v=bosstarget_v, $
		idlutils_v=idlutils_v, $
		photoop_v=photoop_v

end

; Wrappers for a particular setup

; Another Nikhil hack
pro esboss::create_pbs_qso_sweep20091116_local, dogather=dogather
	; use local codes with the new photo_sweep for post-commissioning
	; testing
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid2006/bosswork/groups/boss/sweeps/2009-11-16"

	;pars={likelihood_thresh: 1d-30}
	trun='2009-12-07m'
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			pars=pars, $
			photo_sweep=PHOTO_SWEEP
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather
	endelse
end



; this is Nikhil's temporary hack
pro esboss::create_pbs_qso_sweep20091127_local, dogather=dogather
	; use local codes with the new photo_sweep for post-commissioning
	; testing
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid2006/bosswork/groups/boss/sweeps/2009-11-27.tmp"
	IDLUTILS_V="v5_4_9"

	trun='2009-11-29'
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V, $
			/dogather
	endelse
end


pro esboss::create_pbs_qso_sweep20090928_postcomm_local, dogather=dogather
	; use local codes with the new photo_sweep for post-commissioning
	; testing
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"
	IDLUTILS_V="v5_4_9"

	trun='2009-11-28'
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V, $
			/dogather
	endelse
end

pro esboss::create_pbs_qso_sweep20090928_100sqdeg, dogather=dogather
	; loose cuts so we can go for ~100/sq degree.  This will most likely
	; have to be redone when the retraining starts.
	;
	; use local codes with the new photo_sweep for post-commissioning
	; testing
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"
	IDLUTILS_V="v5_4_9"

	trun='2009-12-01-100'
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V, $
			/dogather
	endelse
end

pro esboss::create_pbs_qso_sweep20090928_loose_stardens, dogather=dogather
	; tune with x2star at 3 but more permissive in stardens. Then we will
	; figure out what gives more than 80/sq deg
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"
	IDLUTILS_V="v5_4_9"

	pars = {$
		likelihood_thresh: 0.02, $
		nn_xnn_thresh:     0.4, $
		x2star_permissive: 3.0, $
		logstardensmax_bright_permissive: 0.5, $
		logstardensmax_faint_permissive: 0.3 $
	}

	trun='2009-12-01-stardens'
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V, $
			pars=pars
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V, $
			/dogather
	endelse
end



pro esboss::create_pbs_qso_sweep20090928_loose, dogather=dogather
	; loose cuts so we can go for ~100/sq degree.  This will most likely
	; have to be redone when the retraining starts.
	;
	; use local codes with the new photo_sweep for post-commissioning
	; testing
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"
	IDLUTILS_V="v5_4_9"

	pars = {$
		x2star_permissive: 0.1, $
		likelihood_thresh: 0.02, $
		nn_xnn_thresh:     0.4  $
	}

	trun='2009-12-02-loosenn2'
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V, $
			pars=pars
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V, $
			/dogather
	endelse
end

pro esboss::create_pbs_qso_sweep20090614_chunk1_newcode, dogather=dogather

	; run these sweeps through the new code
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-06-14"
	IDLUTILS_V="v5_4_9"

	trun='2009-12-02-chunk1-newcode'
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V, $
			/dogather
	endelse

end

pro esboss::create_pbs_std_20091116
	; pick up the missing runs
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"

	; same run as the qso run on this area
	target_run='2009-12-10f'
	self->create_pbs, 'std', target_run, photo_sweep=PHOTO_SWEEP

end





pro esboss::create_pbs_qso_20091116_finish, dogather=dogather
	; pick up the missing runs
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"

	target_run='2009-12-10f'

	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', target_run, $
			photo_sweep=PHOTO_SWEEP
	endif else begin
		self->create_pbs, 'qso', target_run, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather
	endelse
end



pro esboss::create_pbs_qso_sweep20090928_newlike2_loose, dogather=dogather
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"

	trun='2009-12-10-newlike2-loose'
	pars={norank:1, likelihood_thresh: 0}

	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			pars=pars, $
			photo_sweep=PHOTO_SWEEP
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather
	endelse
end
pro esboss::create_pbs_qso_sweep20090928_newlike2_norank, dogather=dogather
	; this one has the tuned likelihoods but no rank flags
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"

	target_run='2009-12-10-newlike2-norank'

	pars={norank:1}
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', target_run, $
			pars=pars, $
			photo_sweep=PHOTO_SWEEP
	endif else begin
		self->create_pbs, 'qso', target_run, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather
	endelse
end
pro esboss::create_pbs_qso_newlike2, dogather=dogather
	; This one has the tuned rank cuts.
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"

	target_run='2009-12-10-newlike2'

	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', target_run, $
			pars=pars, $
			photo_sweep=PHOTO_SWEEP
	endif else begin
		self->create_pbs, 'qso', target_run, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather
	endelse
end

pro esboss::create_pbs_qso_newlike1, dogather=dogather
	; This one has the tuned rank cuts.
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-06-14"

	target_run='2009-12-10-newlike1'

	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', target_run, $
			photo_sweep=PHOTO_SWEEP
	endif else begin
		self->create_pbs, 'qso', target_run, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather
	endelse
end




pro esboss::create_pbs_qso_newlike34, dogather=dogather
	; chunk1 on new likelihood tuned to 0.97 cut
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"

	trun='2009-12-10-newlike34'

	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather
	endelse
end



pro esboss::create_pbs_qso_sweep20090614_newlike1, dogather=dogather
	; chunk1 on new likelihood tuned to 0.97 cut
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-06-14"

	trun='2009-12-09-newlike1'
	;pars={likelihood_thresh: 1d-30}
	;pars={oldlike: 1}

	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			pars=pars, $
			photo_sweep=PHOTO_SWEEP
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather
	endelse
end



pro esboss::create_pbs_qso_sweep20090614, dogather=dogather
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-06-14"

	trun='2009-12-09-rank1'
	;pars={likelihood_thresh: 1d-30}
	pars={oldlike: 1}

	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			pars=pars, $
			photo_sweep=PHOTO_SWEEP
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather
	endelse
end


pro esboss::create_pbs_qso_sweep20090614_stripe82_pairs
	; just some N/S pairs on stripe 82 for alternatives to primary. These
	; are good-to-bad
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-06-14"

	runs=[4207,3388, $ ; good seeing
		  4858,4203,$  ; median seeing
		  4198,3434]   ; bad seeing
	trun='2009-11-30-stripe82pairs'
	pars={likelihood_thresh: 1d-30}
	self->create_pbs_bycamcol, 'qso', trun, $
		runs=runs, $
		pars=pars, $
		photo_sweep=PHOTO_SWEEP, $
		/ignore_resolve
end



pro esboss::create_pbs_qso_sweep20090928_newcode_oldpars, dogather=dogather
	; use local codes with the new photo_sweep for post-commissioning
	; testing
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"
	IDLUTILS_V="v5_4_9"

	trun='2009-11-30-oldpars'
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V, $
			/comm2
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			idlutils_v=IDLUTILS_V, $
			/dogather, $
			/comm2
	endelse
end

pro esboss::create_pbs_qso_sweep20090928_newlike, dogather=dogather
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"

	;pars={likelihood_thresh: 1d-30}
	; turn off ranking selection for now, since we still need shirley's 
	; numbers
	pars={norank:1}
	trun='2009-12-08-newlike'
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			pars=pars
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather
	endelse
end




pro esboss::create_pbs_lrg_newlrg34, noknown=noknown
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"

	run='2009-12-10-newlrg34'
	if not keyword_set(noknown) then begin
		self->create_pbs, 'lrg', run, photo_sweep=PHOTO_SWEEP
	endif else begin

		bt=obj_new('bosstarget')
		ws=bt->default_where_string('lrg')
		known=string(sdss_flagval('boss_target1','sdss_gal_known'),f='(i0)')
		ws += ' and ((str.boss_target1 and '+known+') eq 0)'

		self->create_pbs, 'lrg', run, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather, /reselect, $
			where_string=ws, $
			extra_name='noknown'
	endelse
end


pro esboss::create_pbs_lrg_newlrg1, noknown=noknown
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-06-14"

	run='2009-12-10-newlrg1'
	if not keyword_set(noknown) then begin
		self->create_pbs, 'lrg', run, photo_sweep=PHOTO_SWEEP
	endif else begin

		bt=obj_new('bosstarget')
		ws=bt->default_where_string('lrg')
		known=string(sdss_flagval('boss_target1','sdss_gal_known'),f='(i0)')
		ws += ' and ((str.boss_target1 and '+known+') eq 0)'

		self->create_pbs, 'lrg', run, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather, /reselect, $
			where_string=ws, $
			extra_name='noknown'
	endelse
end

pro esboss::create_pbs_lrg_newlrg2, noknown=noknown
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"

	run='2009-12-10-newlrg2'
	if not keyword_set(noknown) then begin
		self->create_pbs, 'lrg', run, photo_sweep=PHOTO_SWEEP
	endif else begin

		bt=obj_new('bosstarget')
		ws=bt->default_where_string('lrg')
		known=string(sdss_flagval('boss_target1','sdss_gal_known'),f='(i0)')
		ws += ' and ((str.boss_target1 and '+known+') eq 0)'

		self->create_pbs, 'lrg', run, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather, /reselect, $
			where_string=ws, $
			extra_name='noknown'
	endelse
end
pro esboss::create_pbs_lrg_newlrg_loose_sparse, noknown=noknown
	; no cuts on sparse sampling
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"

	run='2009-12-10-newlrg2-loose-sparse'
	pars={ramod_thresh:0}
	if not keyword_set(noknown) then begin
		self->create_pbs, 'lrg', run, photo_sweep=PHOTO_SWEEP, $
			pars=pars
	endif else begin

		bt=obj_new('bosstarget')
		ws=bt->default_where_string('lrg')
		known=string(sdss_flagval('boss_target1','sdss_gal_known'),f='(i0)')
		ws += ' and ((str.boss_target1 and '+known+') eq 0)'

		self->create_pbs, 'lrg', run, $
			photo_sweep=PHOTO_SWEEP, $
			/dogather, /reselect, $
			where_string=ws, $
			extra_name='noknown'
	endelse
end




pro esboss::create_pbs_lrg_sweep20090928_comm2, noknown=noknown
	; use local codes with the new photo_sweep for ngc commissioning
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"
	IDLUTILS_V="v5_4_9"
	BOSSTARGET_V="v1_1_2"

	run='comm2'
	if not keyword_set(noknown) then begin
		self->create_pbs, 'lrg', run, $
			/commissioning, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V
	endif else begin

		bt=obj_new('bosstarget')
		ws=bt->default_where_string('lrg')
		known=string(sdss_flagval('boss_target1','sdss_gal_known'),f='(i0)')
		ws += ' and ((str.boss_target1 and '+known+') eq 0)'

		self->create_pbs, 'lrg', run, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V, $
			/commissioning, $
			/dogather, /reselect, $
			where_string=ws, $
			extra_name='noknown'
	endelse
end



pro esboss::create_pbs_qso_sweep20090928_comm2, dogather=dogather
	; use local codes with the new photo_sweep for ngc commissioning
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"
	IDLUTILS_V="v5_4_9"
	BOSSTARGET_V="v1_1_2"

	trun='comm2'
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V, $
			/comm2
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V, $
			/comm2, /dogather
	endelse
end


pro esboss::create_pbs_qso_sweep20090928_main003, dogather=dogather
	; use local codes with the new photo_sweep for ngc commissioning
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"
	IDLUTILS_V="v5_4_11"
	BOSSTARGET_V="v2_0_0"

	trun='comm2'
	pars={nocalib:1}
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V, $
			pars=pars
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V, $
			pars=pars, $
			/dogather
	endelse
end



pro esboss::create_pbs_std_sweep20090928_local
	; use local codes with the new photo_sweep for ngc commissioning
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"
	IDLUTILS_V="v5_4_9"
	BOSSTARGET_V="v1_1_2"

	run='comm2'
	;run='2009-10-07-test'
	;run='2009-10-10-test'
	self->create_pbs, 'std', run, $
		photo_sweep=PHOTO_SWEEP, $
		bosstarget_v=BOSSTARGET_V, $
		idlutils_v=IDLUTILS_V
end


pro esboss::create_pbs_lrg_main002_nocalib, noknown=noknown
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"
	IDLUTILS_V="v5_4_11"
	;BOSSTARGET_V="v2_0_0"

	run='main002'
	pars={nocalib:1}
	if not keyword_set(noknown) then begin
		self->create_pbs, 'lrg', run, $
			pars=pars, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V
	endif else begin

		bt=obj_new('bosstarget')
		ws=bt->default_where_string('lrg')
		known=string(sdss_flagval('boss_target1','sdss_gal_known'),f='(i0)')
		ws += ' and ((str.boss_target1 and '+known+') eq 0)'

		self->create_pbs, 'lrg', run, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V, $
			pars=pars, $
			/dogather, /reselect, $
			where_string=ws, $
			extra_name='noknown'
	endelse
end

pro esboss::create_pbs_qso_main002_nocalib, dogather=dogather
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"
	IDLUTILS_V="v5_4_11"
	;BOSSTARGET_V="v2_0_0"

	trun='main002'

	pars={nocalib:1}
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			pars=pars, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V
	endif else begin
		self->create_pbs, 'qso', trun, $
			pars=pars, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V, $
			/dogather
	endelse
end


pro esboss::create_pbs_lrg_main001, noknown=noknown
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"
	IDLUTILS_V="v5_4_11"
	BOSSTARGET_V="v2_0_0"

	run='main001'
	if not keyword_set(noknown) then begin
		self->create_pbs, 'lrg', run, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V
	endif else begin

		bt=obj_new('bosstarget')
		ws=bt->default_where_string('lrg')
		known=string(sdss_flagval('boss_target1','sdss_gal_known'),f='(i0)')
		ws += ' and ((str.boss_target1 and '+known+') eq 0)'

		self->create_pbs, 'lrg', run, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V, $
			/dogather, /reselect, $
			where_string=ws, $
			extra_name='noknown'
	endelse
end
pro esboss::create_pbs_qso_main001, dogather=dogather
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"
	IDLUTILS_V="v5_4_11"
	BOSSTARGET_V="v2_0_0"

	trun='main001'
	if not keyword_set(dogather) then begin
		self->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V
	endif else begin
		self->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V, $
			/dogather
	endelse
end
pro esboss::create_pbs_std_main001
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"
	IDLUTILS_V="v5_4_11"
	BOSSTARGET_V="v2_0_0"

	; same run as the qso run on this area
	target_run='main001'
	self->create_pbs, 'std', target_run, $
		photo_sweep=PHOTO_SWEEP, $
		bosstarget_v=BOSSTARGET_V, $
		idlutils_v=IDLUTILS_V
end




;+
; NAME:
;   plot_poly
; PURPOSE:
;   plots a mangle polygon (or passes back what you need to plot it)
; CALLING SEQUENCE:
;   plot_poly, poly [, offset=, xrange=, yrange=, filename=, /fill, $
;      /nooutline, xsize=, ysize=, /over, color=, minside=, dangle=, $
;      outline_thick=, splot=, /aitoff ]
; INPUTS:
;   poly - mangle polygon (e.g. one created by construct_polygon())
; OPTIONAL INPUTS:
;   offset - offset to apply to ra
;   xrange, yrange - ranges to pass to plot
;   filename, xsize, ysize - PS file to output to (and its sizes)
;   minsize, dangle - pass to gverts
; OPTIONAL KEYWORDS:
;   /fill - fill
;   /over - over plot on current device
;   /nooutine - do not outline
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   01-Oct-2002  Written by MRB (NYU)
;	Added ability to send xrange= *or* yrange=.  Added _extra= keyword.
;	2010-01-13 - Erin Sheldon, BNL
;-
;------------------------------------------------------------------------------
pro esboss::plot_poly,poly,offset=offset,xrange=xrange,yrange=yrange, $
              filename=filename,fill=fill,nooutline=nooutline, $
              xsize=xsize, ysize=ysize, over=over, color=color, $
              minside=minside, dangle=dangle, outline_thick=outline_thick, $
              splot=splot, soplot=soplot, aitoff=aitoff, noplot=noplot, $
			  _extra=_extra

if(not keyword_set(offset)) then offset=0.
if(not keyword_set(outline_thick)) then outline_thick=0.001
if(not keyword_set(minside)) then minside=7
if(not keyword_set(dangle)) then dangle=0.05
if(n_elements(color) eq 0) then $
  use_color=fix(floor(64.0+127.0*randomu(seed,n_elements(poly))))
if(n_elements(color) eq 1) then use_color=replicate(color,n_elements(poly))
if(n_elements(color) gt 1) then use_color=color

if(keyword_set(filename)) then begin
    set_print,filename=filename,xsize=xsize,ysize=ysize
endif

;if(keyword_set(fill)) then $
  ;loadct,0
if(keyword_set(offset)) then $
  xtitle='ra-'+strtrim(string(offset),2) $
else $
  xtitle='ra'

if n_elements(xrange) eq 0 then do_xrange=1 else do_xrange=0
if n_elements(yrange) eq 0 then do_yrange=1 else do_yrange=0

if (do_xrange or do_yrange) and not keyword_set(over) then begin

    if do_yrange then yrange=[1.e+10,-1.e+10]
    if do_xrange then xrange=[1.e+10,-1.e+10]
    for i=0L, n_elements(poly)-1L do begin
        if(garea(poly[i]) gt 0.) then begin
            verts=0
            gverts,poly[i],angle=angle, verts=verts, edges=edges, ends=ends, $
              dangle=dangle,minside=minside,/plot
            if(keyword_set(verts)) then begin
                x_to_angles,verts,ra,theta
                ramoffset=(ra-offset) 
                indx=where(ramoffset lt 0.,count)
                if(count gt 0) then ramoffset[indx]=ramoffset[indx]+360.
                dec=reform(90.-theta,n_elements(theta))
                ramoffset=reform(ramoffset,n_elements(ramoffset))
				if do_yrange then begin
					yrange[0]=min([yrange[0],dec])
					yrange[1]=max([yrange[1],dec])
				endif
				if do_xrange then begin
					xrange[0]=min([xrange[0],ramoffset])
					xrange[1]=max([xrange[1],ramoffset])
				endif
            endif
        endif
    endfor
endif

if(NOT keyword_set(over)) then begin
    if keyword_set(splot) then begin
		splot,[0],[0],/nodata,xtitle=xtitle,ytitle='dec', $
			xrange=xrange,yrange=yrange, _extra=_extra
	endif else begin
		plot,[0],[0],/nodata,xtitle=xtitle,ytitle='dec', $
			xrange=xrange,yrange=yrange, _extra=_extra
	endelse
endif

for i=0L, n_elements(poly)-1L do begin
    if(garea(poly[i]) gt 0.) then begin
        verts=0
        gverts,poly[i],angle=angle, verts=verts, edges=edges, ends=ends, $
          dangle=dangle,minside=minside,/plot
        if(keyword_set(verts)) then begin
            x_to_angles,verts,ra,theta
            ramoffset=(ra-offset) 
            indx=where(ramoffset lt 0.,count)
            if(count gt 0) then ramoffset[indx]=ramoffset[indx]+360.
            dec=90.-theta
            xx=ramoffset
            yy=dec
            if(keyword_set(aitoff)) then begin
                aitoff, ramoffset, dec, xx, yy
            endif
            if(keyword_set(fill)) then begin
                if(keyword_set(splot)) then $
                  spolyfill,xx,yy,color=use_color[i],noclip=0 $
                else $
                  polyfill,xx,yy,color=use_color[i],noclip=0
            endif
            if(NOT keyword_set(nooutline)) then begin
                if(keyword_set(splot)) then $
                  soplot,xx,yy,color=use_color[i], $
                  thick=outline_thick $
                else $
                  oplot,xx,yy,color=use_color[i], $
                  thick=outline_thick 
            endif
        endif
    endif
endfor

if(keyword_set(filename)) then begin
    end_print
endif

end


pro esboss::compare_jessica

	; Read in the targeting file
	erinfile = '/home/jessica/boss/boss-qso-erin.fits'
	print,'reading like stat file: ',erinfile
	collateinfo = mrdfits(erinfile, 2)


	newfile = '/clusterfs/riemann/raid006/bosswork/groups/boss/target/2010-03-02-lnnv/bosstarget-qso-2010-03-02-lnnv-collate.fits'
	print,'reading collated file: ',newfile
	newtable = mrdfits(newfile, 2)


	bq=obj_new('bosstarget_qso')

	print,'running sphoto_match'
	sphoto_match, newtable, collateinfo, mnewtable, mcollateinfo
	;spherematch, collateinfo.ra, collateinfo.dec, newtable.ra, $
	;	newtable.dec, 2./3600, i3, i4, d34, maxmatch=1

	;liketruthmatch = newtable[i4]
	;likechi = collateinfo[i3]

	liketruthmatch = newtable[mnewtable]
	likechi = collateinfo[mcollateinfo]
	
	help,likechi

	mclike_ratio = bq->calculate_mcvalue_rmag_like_ratio( likechi, liketruthmatch )


	outstruct = {run:0,rerun:'',camcol:0,field:0,id:0L,ra:0d,dec:0d,mclike_ratio:0.0}
	outstruct = replicate(outstruct, n_elements(likechi))
	copy_struct, liketruthmatch, outstruct
	outstruct.mclike_ratio = mclike_ratio

	outfile='~/tmp/compare-jessica.fits'
	print,'writing outfile: ',outfile
	mwrfits,outstruct,outfile,/create

end


pro esboss__define
	struct = {$
		esboss, $
		verbose:0 $
	}
end
