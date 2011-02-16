function firstlobe::init, matchtype=matchtype
	common firstlobe_block, cached_matchtype, firstradio
    if n_elements(cached_matchtype) eq 0 then begin
        cached_matchtype = 'none'
    endif

    self->set_pars, matchtype=matchtype
	return, 1
end

pro firstlobe::set_pars, matchtype=matchtype
    if n_elements(matchtype) ne 0 then begin
        self.matchtype = matchtype
    endif
    if self.matchtype ne 'midpoint' then begin
        self.matchtype = 'standard'
    endif
end


function firstlobe::match, objs

    ntot = n_elements(objs)
	if ntot eq 0 then begin
		on_error, 2
		splog,'Usage: mstruct=bl->match(objs, matchtype="standard")'
		message,'Halting'
	endif

    self->set_pars, matchtype=matchtype

	common firstlobe_block, cached_matchtype, firstradio

	self->do_cache

    ; the only check we will make is primary
	primaryflag = sdss_flagval('RESOLVE_STATUS','SURVEY_PRIMARY')
	w = where((objs.resolve_status AND primaryflag) ne 0, nuse)
    if nuse eq 0 then begin
        splog,'No primary objects'
        return, -1
    endif else begin
        splog,'Found ',nuse,'/',ntot,' primary',format='(a,i0,a,i0,a)'
    endelse

	nobj=n_elements(objs)
	boss_target1=lon64arr(nobj)
	firstradio_struct = replicate({firstradio_id:-9999L, firstradio_dist:-9999.9},  nobj)

    if self.matchtype eq 'midpoint' then begin
        matchrad = double(firstradio.srad) ;radians
    endif else begin
        matchrad = 2d/3600d ; degrees
    endelse

	;spherematch, $
	;	objs[w].ra, objs[w].dec, firstradio.ra, firstradio.dec, matchrad, $
	;	objs_match, first_match, distances,$
	;	maxmatch=1
	htm_match, $
		firstradio.ra, firstradio.dec, objs[w].ra, objs[w].dec,  matchrad, $
		first_match, objs_match, distances,$
		maxmatch=1



    if first_match[0] ne -1 then begin
        w=w[objs_match]

        tresult = objs[w]
        result = $
            struct_addtags(tresult, ['first_id','first_matchrad'],['0L','0d'])

        result.first_id = first_match
        result.first_matchrad = distances
        nmatch = n_elements(result)
    endif else begin
        nmatch=0L
        result=-1
    endelse


	print
	splog,'Found ',nmatch,' FIRST matches',format='(a,i0,a)'
	print

	return, result
end

function firstlobe::struct, n
    struct={ $
        run:0, $
        rerun: '', $
        camcol: 0, $
        field: 0, $
        id: 0, $
        first_id: 0L, $
        first_matchrad: 0d $
    }

    if n_elements(n) ne 0 then begin
        struct=replicate(struct, n)
    endif
    return, struct
end



function firstlobe::file

    if self.matchtype eq 'midpoint' then begin
        dir=getenv('BOSS_TARGET')
        if dir eq '' then message,'$BOSS_TARGET not set'
        file=filepath(root=dir, sub=['esheldon','firstlobe'],'midpoints_final.fits.gz')
    endif else begin
        dir=filepath(root=getenv("BOSSTARGET_DIR"), "data")
        if dir eq '' then message,'$BOSSTARGET_DIR not set'
        file = filepath(root=dir, "first_08jul16.fits")
    endelse
	return, file
end


pro firstlobe::do_cache

	common firstlobe_block, cached_matchtype, firstradio

	if n_elements(firstradio) eq 0 or cached_matchtype ne self.matchtype then begin

		file=self->file()
		splog,'Caching in memory: ',file,form='(a,a)'
		firstradio=mrdfits(file,1)
		if n_tags(firstradio) eq 0 then begin
			message,string('Could not read file: ',file,f='(a,a)')
		endif
        cached_matchtype = self.matchtype
	endif

end



pro firstlobe::make_pbs, matchtype, dogather=dogather
    photo_sweep = '/clusterfs/riemann/netapp/groups/boss/sweeps/2010-01-11'
    photo_resolve = '/clusterfs/riemann/raid007/bosswork/groups/boss/resolve/2010-01-11'
    extra_setups='setup bosstarget -r ~esheldon/exports/bosstarget-work'


    if matchtype eq 'midpoint' then begin
        command = 'fl=obj_new("firstlobe",matchtype="midpoint") & res=fl->match(objs)'
        star_proctype = 'flstar-midpoint'
        gal_proctype = 'flgal-midpoint'
    endif else begin
        command = 'fl=obj_new("firstlobe") & res=fl->match(objs)'
        star_proctype = 'flstar'
        gal_proctype = 'flgal'
    endelse

    procrun = '01'
    sw=obj_new('sweeps','star')
    sw->create_pbs, star_proctype, procrun, command, $
        photo_sweep=photo_sweep, $
        photo_resolve=photo_resolve, $
        extra_setups=extra_setups, $
        dogather=dogather, $
        doproc=doproc

    sw=obj_new('sweeps','gal')
    sw->create_pbs, gal_proctype, procrun, command, $
        photo_sweep=photo_sweep, $
        photo_resolve=photo_resolve, $
        extra_setups=extra_setups, $
        dogather=dogather, $
        doproc=doproc
end




pro firstlobe__define
	struct = {$
		firstlobe, $
        matchtype:'' $
	}
end


