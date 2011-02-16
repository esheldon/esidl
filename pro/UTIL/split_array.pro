;+
; NAME:
;   split_array
;
; PURPOSE:
;   split an array into chunks.  The sub-arrays are returned as
;   elements in a pointer array.
;
; CALLING SEQUENCE:
;   ptrlist = split_array(array, nsplit=, nper=)
;
; INPUTS:
;   array: Any array.
;
; KEYWORDS:
;   nsplit: The number of chunks to return.
;   nper:   The number of elements per chunk.
;
; OUTPUTS:
;   A pointer array with the sub-arrays.
;
; EXAMPLE:
;   array = [1,5,8,2,3,5,6,6,7]
;   plist = split_array(array, nper=3)
;   help, plist
;
; MODIFICATION HISTORY:
;   Creation: Erin Sheldon, BNL, 2011-02-01
;
;-
function split_array, list, nsplit=nsplit, nper=nper
    if n_elements(list) eq 0 then begin
        on_error, 2
        print,'usage: ptrlist = split_array(array, nsplit=, nper=)'
        print,'   you must send nsplit or nper'
        message,'halting'
    endif

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

