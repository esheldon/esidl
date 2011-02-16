pro tmp3, n, xnn, znn, xnnslow, znnslow

	common mytmp_block, t

	if n_elements(t) eq 0 then begin
		t=sweep_readobj(4264,3,type='star',rerun=301)
	endif
	nobj=n_elements(t)

	if n_elements(n) eq 0 then n=nobj

	bq=obj_new('bosstarget_qso')

	mag=bq->get_lups(t, /deredden)

	help,mag
	i=10
	c=2

	print,'IDL: for obj=',i,' cindex=',c,' psfmag=',mag[c,i]
	print,'IDL: for obj=',i,' cindex=',c,' psfmag=',mag[i*5 + c]

	print,'fast'
	tm1=systime(1)
	bq->nn_run, t[0:n-1], xnn, znn
	ptime,systime(1)-tm1


	print,'slow'
	tm1=systime(1)
	bq->nn_run_old, t[0:n-1], xnnslow, znnslow
	ptime,systime(1)-tm1


	;colprint, xnn, znn, xnnslow, znnslow

end

