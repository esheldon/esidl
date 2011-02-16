pro print_compare_struct_header
	print
	print,$
		'tag name',$
		'max(abs(diff))',$
		'med(diff)',$
		'sdev(diff)', $
		'max(abs(%diff))', $
		'med(%diff)', $
		'sdev(%diff)', $
		f='(A30, 6A20)'

end
pro print_compare_struct, structA, structB, full=full, dstruct=dstruct

	if n_elements(dstruct) eq 0 then begin
		dstruct = compare_struct(structA, structB)
	endif

	if n_elements(dstruct) eq 1 and dstruct[0].ndiff eq 0 then begin
		print,'The matching tags for these structs are the same.'
		return
	endif

	ndiff = n_elements(dstruct)

	if not keyword_set(full) then print_compare_struct_header

	for i=0L, ndiff-1 do begin
		tname = dstruct[i].field
		ta=dstruct[i].tag_num_a
		tb=dstruct[i].tag_num_b

		diff = structA.(ta) - structB.(tb)
		pdiff = diff / structA.(ta) * 100
		med_diff = median(diff)
		med_pdiff = median(pdiff)
		stddev_diff = stddev( diff, /nan)
		stddev_pdiff = stddev( pdiff, /nan)

		if keyword_set(full) then print_compare_struct_header
		print, $
			tname,$
			max(abs(diff),/nan),$
			med_diff, $
			stddev_diff, $
			max(abs(pdiff),/nan), $
			med_pdiff, $
			stddev_pdiff, $
			f='(A30, 6g20)'

		if keyword_set(full) then begin
			print
			print,tname+':  Element-by-Element.  hit a key (q to quit)'
			key=get_kbrd(1)
			if key eq 'q' then return

			print,'diff', '%diff', f='(2A20)'
			forprint, diff, pdiff, f='(2g20)'
		endif
	endfor

end
