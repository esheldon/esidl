function type2runpair, type
	case type of
		'good':   return, [4207,3388]
		'median': return, [4858,4203]
		'bad':    return, [4198,3434]
		else: message,string('Bad type: ',type,f='(a,a)')
	endcase
end

pro run_gather, target_run, type
	bt=obj_new('bosstarget')
	sold=getenv('PHOTO_SWEEP')
	PHOTO_SWEEP=$
		"/clusterfs/riemann/raid006/bosswork/groups/boss/sweeps/2009-06-14"
	setenv,'PHOTO_SWEEP='+photo_sweep

	runs = type2runpair(type)

	bt->gather2file, 'qso', target_run, runs=runs, /nomatch, /run_primary, $
		extra_name=type
	setenv,'PHOTO_SWEEP='+sold
end

pro process_stripe82_runpairs, type, gather=gather

	target_run='2009-11-30-stripe82pairs'
	chunk=1

	bt=obj_new('bosstarget')
	eb=obj_new('esboss')

	if keyword_set(gather) then begin
		run_gather, target_run, type
	endif

	eb->make_fullcache_and_inwindow, 'qso', target_run, chunk, $
		extra_name=type, /restrict

end
