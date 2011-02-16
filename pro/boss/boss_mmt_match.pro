

pro boss_mmt_match, type, all=all, match_firstrass=match_firstrass, $
		fpobjc=fpobjc, fieldid=fieldid

	if n_elements(type) eq 0 then begin
		print,'usage:  boss_mmt_match, type'
		print,'type should be "star" or "gal"'
		on_error, 2
		message,'halting'
	endif

	; for testing target selection
	bq=obj_new('bosstarget_qso')

	rass_matchdist = 60.0/3600.0
	first_matchdist = 2.0/3600.0

	mmtdir = getenv("BOSSTARGET_DIR")
	mmtdir=path_join(mmtdir, "data")

	mmtname = path_join(mmtdir, "mmt-plate-data.dat")

	boss_target = getenv("BOSS_TARGET")
	;outdir = path_join(boss_target, subdir=["esheldon","mmt"])
	outdir = path_join(boss_target, "mmt")

	if keyword_set(fpobjc) then begin
		mtype='fpObjc'
	endif else begin
		mtype='sweep'
	endelse
	if keyword_set(all) then begin
		outdir = path_join(outdir, mtype+'-match-'+type+'-all')
		addstr = 'all-'
	endif else begin
		outdir = path_join(outdir, mtype+'-match-'+type)
		addstr = ''
	endelse

	if keyword_set(match_firstrass) then begin
		outdir = path_join(outdir, '-firstrass')
		addstr=addstr+'-firstrass'
	endif

	file_mkdir, outdir

	stdef = {ra:0d, dec:0d, rad:0d}
	read_struct, mmtname, stdef, plateinfo, skiplines=1

	if n_elements(fieldid) ne 0 then begin
		startnum=fieldid
		endnum=fieldid
	endif else begin
		startnum=0
		endnum=n_elements(plateinfo)-1
	endelse
	for i=startnum, startnum do begin
		if i gt 0 then begin
			print,'-----------------------------------------------------'
		endif
		outname = 'mmt-'+mtype+'-match-'+addstr+type+ntostr(i,f='(i02)')+'.fits'
		epsname = 'mmt-'+mtype+'-match-'+addstr+type+ntostr(i,f='(i02)')+'.eps'

		outname = path_join(outdir, outname)
		epsname = path_join(outdir, epsname)
		print,'Will write to file: ',outname

		print,'Reading '+mtype+' round :'
		print,'    ra: ',plateinfo[i].ra
		print,'    dec:',plateinfo[i].dec
		print,'    rad:',plateinfo[i].rad

		data=0
		if not keyword_set(fpobjc) then begin
			data = sdss_sweep_circle($
				plateinfo[i].ra, $
				plateinfo[i].dec, $
				plateinfo[i].rad, all=all, type=type)
		endif else begin
			data = sdss_circle($
				plateinfo[i].ra, $
				plateinfo[i].dec, $
				plateinfo[i].rad, rerun=301, all=all, objtype=type)
		endelse
		if n_tags(data) eq 0 then begin
			print,'NO MATCHES!!'
			continue
		endif
		help,data

		if not keyword_set(match_firstrass) then begin
			print,'Writing to file: ',outname
			mwrfits, data, outname, /create
		endif else begin
			; now match to first and rass
			print,'matching rass to radius: ',rass_matchdist
			rassmatch = es_catmatch($
				data.ra, data.dec, matchdist=rass_matchdist, cat='RASS')
			wr=where(rassmatch.rass_matchdist ge 0, nwr)
			print,'Found '+ntostr(nwr)+' rass matches'

			print,'matching first to radius: ',first_matchdist
			firstmatch = es_catmatch($
				data.ra, data.dec, matchdist=first_matchdist, cat='FIRST')
			wf=where(firstmatch.first_matchdist ge 0, nwf)
			print,'Found '+ntostr(nwf)+' first matches'

			outdata = create_struct(data[0], rassmatch[0], firstmatch[0])

			outdata = create_struct(outdata, $
				{x2_phot:-9999.,x2_star:-9999.,z_phot:-9999.,$
					qso_badflags:0L,bosstarget1:0L} )

			outdata = replicate(outdata, n_elements(data))


			struct_assign, data, outdata, /nozero
			struct_assign, rassmatch, outdata, /nozero
			struct_assign, firstmatch, outdata, /nozero

			; add target selection for qso
			if type eq 'star' then begin
				allow_move=0
			endif else begin
				allow_move=1
			endelse
			res = bq->select(outdata, /struct, allow_move=allow_move)
			outdata.x2_phot = res.x2_phot
			outdata.x2_star = res.x2_star
			outdata.z_phot = res.z_phot
			outdata.qso_badflags = res.bitmask
			outdata.bosstarget1 = res.bosstarget1


			print,'Writing to file: ',outname
			mwrfits, outdata, outname, /create

			begplot,epsname,/encap,/color,xsize=8.5
			!p.charsize=1
			!p.thick=1
			!x.thick=1
			!y.thick=1
			pplot, outdata.ra, outdata.dec, psym=3, xstyle=3, ystyle=3, $
				xtitle='RA', ytitle='DEC', aspect=1

			for j=0L,nwr-1 do begin
				tvcircle, $
					rass_matchdist, $
					outdata[wr[j]].rass_ra, $
					outdata[wr[j]].rass_dec, $
					/data, color=c2i('red')
			endfor

			if nwf gt 0 then begin
				ora = outdata[wf].first_ra
				odec = outdata[wf].first_dec
				if nwf eq 1 then begin
					ora = [ora]
					odec = [odec]
				endif
				pplot, ora, odec, psym=8, color=c2i('blue'), /overplot
			endif

			legend, ['rass','first'], color=[c2i('red'),c2i('blue')], $
				line=0, /right


			endplot,/trim
		endelse


	endfor

end
