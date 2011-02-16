function runboss::init
	return, 1
end


;
; Bug fix in std using nocalib
;

pro runboss::create_pbs_std_main004_nocalib
	eb=obj_new('esboss')
	photo_sweep=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"
	idlutils_v="v5_4_11"
	bosstarget_v="-r /home/esheldon/exports/bosstarget-v2_0_2"

	; same run as the qso run on this area
	target_run='main004'
	pars={nocalib:1}
	eb->create_pbs, 'std', target_run, $
		pars=pars, $
		photo_sweep=photo_sweep, $
		bosstarget_v=bosstarget_v, $
		idlutils_v=idlutils_v
end



;
; run of chunk2 with new code
;

pro runboss::create_pbs_lrg_main003, noknown=noknown
	; rerun chunk2 with new code
	eb=obj_new('esboss')
	photo_sweep=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"
	idlutils_v="v5_4_11"
	bosstarget_v="v2_0_0"

	run='main003'
	if not keyword_set(noknown) then begin
		eb->create_pbs, 'lrg', run, $
			pars=pars, $
			photo_sweep=photo_sweep, $
			bosstarget_v=bosstarget_v, $
			idlutils_v=idlutils_v
	endif else begin

		bt=obj_new('bosstarget')
		ws=bt->default_where_string('lrg')
		known=string(sdss_flagval('boss_target1','sdss_gal_known'),f='(i0)')
		ws += ' and ((str.boss_target1 and '+known+') eq 0)'

		eb->create_pbs, 'lrg', run, $
			photo_sweep=photo_sweep, $
			bosstarget_v=bosstarget_v, $
			idlutils_v=idlutils_v, $
			pars=pars, $
			/dogather, /reselect, $
			where_string=ws, $
			extra_name='noknown'
	endelse
end


pro runboss::create_pbs_qso_main003, dogather=dogather
	eb=obj_new('esboss')
	; rerun chunk2 with new code
	photo_sweep=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"
	idlutils_v="v5_4_11"
	bosstarget_v="v2_0_0"

	trun='main003'
	if not keyword_set(dogather) then begin
		eb->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=photo_sweep, $
			bosstarget_v=bosstarget_v, $
			idlutils_v=idlutils_v, $
			pars=pars
	endif else begin
		eb->create_pbs, 'qso', trun, $
			photo_sweep=photo_sweep, $
			bosstarget_v=bosstarget_v, $
			idlutils_v=idlutils_v, $
			pars=pars, $
			/dogather
	endelse
end
pro runboss::create_pbs_std_main003
	eb=obj_new('esboss')
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-09-28"
	IDLUTILS_V="v5_4_11"
	BOSSTARGET_V="v2_0_0"

	; same run as the qso run on this area
	target_run='main003'
	eb->create_pbs, 'std', target_run, $
		photo_sweep=PHOTO_SWEEP, $
		bosstarget_v=BOSSTARGET_V, $
		idlutils_v=IDLUTILS_V
end



;
; mai002, chunks 3-4 *without* non-photometric objects removed
;


pro runboss::create_pbs_lrg_main002_nocalib, noknown=noknown
	eb=obj_new('esboss')
	photo_sweep=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"
	idlutils_v="v5_4_11"
	bosstarget_v="v2_0_1"

	run='main002'
	pars={nocalib:1}
	if not keyword_set(noknown) then begin
		eb->create_pbs, 'lrg', run, $
			pars=pars, $
			photo_sweep=photo_sweep, $
			bosstarget_v=bosstarget_v, $
			idlutils_v=idlutils_v
	endif else begin

		bt=obj_new('bosstarget')
		ws=bt->default_where_string('lrg')
		known=string(sdss_flagval('boss_target1','sdss_gal_known'),f='(i0)')
		ws += ' and ((str.boss_target1 and '+known+') eq 0)'

		eb->create_pbs, 'lrg', run, $
			photo_sweep=photo_sweep, $
			bosstarget_v=bosstarget_v, $
			idlutils_v=idlutils_v, $
			pars=pars, $
			/dogather, /reselect, $
			where_string=ws, $
			extra_name='noknown'
	endelse
end

pro runboss::create_pbs_qso_main002_nocalib, dogather=dogather
	eb=obj_new('esboss')
	photo_sweep=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"
	idlutils_v="v5_4_11"
	bosstarget_v="v2_0_1"

	trun='main002'

	pars={nocalib:1}
	if not keyword_set(dogather) then begin
		eb->create_pbs_bycamcol, 'qso', trun, $
			pars=pars, $
			photo_sweep=photo_sweep, $
			bosstarget_v=bosstarget_v, $
			idlutils_v=idlutils_v
	endif else begin
		eb->create_pbs, 'qso', trun, $
			pars=pars, $
			photo_sweep=photo_sweep, $
			bosstarget_v=bosstarget_v, $
			idlutils_v=idlutils_v, $
			/dogather
	endelse
end


;
; mai001, chunks 3-4 with non-photometric objects removed
;


pro runboss::create_pbs_lrg_main001, noknown=noknown
	eb=obj_new('esboss')
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"
	IDLUTILS_V="v5_4_11"
	BOSSTARGET_V="v2_0_0"

	run='main001'
	if not keyword_set(noknown) then begin
		eb->create_pbs, 'lrg', run, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V
	endif else begin

		bt=obj_new('bosstarget')
		ws=bt->default_where_string('lrg')
		known=string(sdss_flagval('boss_target1','sdss_gal_known'),f='(i0)')
		ws += ' and ((str.boss_target1 and '+known+') eq 0)'

		eb->create_pbs, 'lrg', run, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V, $
			/dogather, /reselect, $
			where_string=ws, $
			extra_name='noknown'
	endelse
end
pro runboss::create_pbs_qso_main001, dogather=dogather
	eb=obj_new('esboss')
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"
	IDLUTILS_V="v5_4_11"
	BOSSTARGET_V="v2_0_0"

	trun='main001'
	if not keyword_set(dogather) then begin
		eb->create_pbs_bycamcol, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V
	endif else begin
		eb->create_pbs, 'qso', trun, $
			photo_sweep=PHOTO_SWEEP, $
			bosstarget_v=BOSSTARGET_V, $
			idlutils_v=IDLUTILS_V, $
			/dogather
	endelse
end
pro runboss::create_pbs_std_main001
	eb=obj_new('esboss')
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-11-16"
	IDLUTILS_V="v5_4_11"
	BOSSTARGET_V="v2_0_0"

	; same run as the qso run on this area
	target_run='main001'
	eb->create_pbs, 'std', target_run, $
		photo_sweep=PHOTO_SWEEP, $
		bosstarget_v=BOSSTARGET_V, $
		idlutils_v=IDLUTILS_V
end





pro runboss__define
	struct = {$
		runboss, $
		crap:0 $
	}
end
