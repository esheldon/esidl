function get_shirley_plate_data
	plate_data = {platedata, plateid:0L, ra:0d, dec:0d}
	plate_data = replicate(plate_data, 6)

	plate_data[0] = {platedata, plateid:3673, ra:117.14, dec:46.63}
	plate_data[1] = {platedata, plateid:3656, ra:110.84 , dec:  41.97}
	plate_data[2] = {platedata, plateid:3657,  ra:112.23, dec:   39.36}
	plate_data[3] = {platedata, plateid:3659,  ra:113.33, dec:   42.67}
	plate_data[4] = {platedata, plateid:3660,  ra:113.91, dec:   44.32}
	plate_data[5] = {platedata, plateid:3661,   ra:114.08, dec:  38.38}

	return,plate_data
end

function recoverprobs_matchplates, t
	plate_data = get_shirley_plate_data()

	plate_radius = 1.49d

	print,'Matching plates within ',plate_radius,' degrees',f='(a,f0.2,a)'
	spherematch, plate_data.ra, plate_data.dec, t.ra, t.dec, $
		plate_radius, $
		pm, tm, maxmatch=0

	; remove duplicates
	if tm[0] ne -1 then begin
		tm = tm[ rem_dup(tm) ]
	endif
	return, tm
end
pro prepare_recoverprobs, struct=struct, remake=remake, combine=combine,$
		fix=fix
	target_run="comm2"
	bt=obj_new('bosstarget')

	file=bt->target_file('qso',target_run,/collate,extra_name='bitmask')
	dir=file_dirname(file)
	dir = path_join(dir, 'recoverprobs')


	if keyword_set(combine) then begin
		if keyword_set(fix) then begin
			orig_files=file_search(dir,'chunk[0-9][0-9][0-9].fits')
			files=file_search(dir, 'chunk*rankprob-patrankV4.fits-dups')

			orig=mrdfits_multi(orig_files)
			comb=mrdfits_multi(files)

			pid = photoid(orig)
			keep = rem_dup(pid)



			print,'Keeping ',n_elements(keep),'/',n_elements(comb),' unique',$
				f='(a,i0,a,i0,a)'
			comb = comb[keep]

		endif else begin
			files=file_search(dir, 'chunk*rankprob-patrankV4.fits')
			comb=mrdfits_multi(files)
		endelse

		combined_file=path_join(dir,target_run+'-rankprob-patrankV6.fits')
		print,'Will write to combined file: ',combined_file
		; from Adam's code
		; Sort on the various parameters and append the rank

		comb.like_rank=sort(sort(comb.like_ratio))
		comb.like_rank_pat=sort(sort(comb.like_ratio_pat))
		comb.kde_rank=sort(sort(comb.kde_prob))
		comb.kde_rank_pat=sort(sort(comb.kde_prob_pat))
		comb.nn_rank=sort(sort(comb.nn_xnn))

		print,'Writing to combined file: ',combined_file
		mwrfits, comb, combined_file,/create
		return
	endif


	if n_elements(struct) eq 0 then begin
		t=bt->read_collated('qso',target_run,extra_name='bitmask',file=file)
	endif

	fout=file_basename(file)
	fout=repstr(fout, '.fits', '-platematch.fits')
	fout=path_join(dir, fout)
	psfile=repstr(fout, '.fits', '.eps')

	if not file_test(fout) then begin
		print,'Will write overall match file: ',fout


		; match to the plates
		matches = recoverprobs_matchplates(struct)

		t=struct[matches]

		print,'Writing ',n_elements(t),' to overall file: ',fout,f='(a,i0,a,a)'
		mwrfits, t, fout, /create

		begplot,psfile,/encap
		plot, t.ra, t.dec, psym=3, /iso, /ynozero
		endplot
	endif else begin
		print,'Reading overall match file: ',fout
		t=mrdfits(fout,1)
	endelse



	; break into 50 chunks

	bosstarget_v = "-r /home/esheldon/exports/bosstarget-work"
	photoop_v = "v1_9_4"
	idlutils_v="v5_4_9"


	add_arrval,"setup tree", setups 
	add_arrval,"setup photoop "+photoop_v, setups
	add_arrval,"setup idlutils "+idlutils_v, setups
	add_arrval,"setup bosstarget "+bosstarget_v, setups

	setups = strjoin(setups, ' && ')

	nt = n_elements(t)
	nmain = 100
	chunksize = nt/nmain

	nleft = nt mod chunksize
	if nleft gt 0 then nchunks = nmain + 1

	for i=0L, nchunks-1 do begin
		ind1 = i*chunksize
		ind2 = (i+1)*chunksize-1
		if ind2 ge nt then begin
			ind2=nt-1
		endif
		print,ind1,ind2

		jobname=string('chunk',i,f='(a,i03)')
		chunkfile=jobname+'.fits'
		chunkfile=path_join(dir, chunkfile)

		print,'Writing chunk file: ',chunkfile
		mwrfits, t[ind1:ind2], chunkfile, /create

		pbsfile=repstr(chunkfile,'.fits','.pbs')
		idl_command="bosstarget_qso_recoverprobs,'"+chunkfile+"',ext=1,/fast"

		print,'Writing pbs file: ',pbsfile
		pbs_riemann_idl, pbsfile, idl_command, setup=setups, job_name=jobname

	endfor

end


