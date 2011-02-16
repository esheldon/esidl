pro match_pos, photoobj,photoindex, sexcat, match_struct, matches

if (n_params() eq 0) then begin
	print,'Syntax: match_pos, photoobj, photoindex, sexcat, match_struct,matches'
	return
endif


tolerance = 6.0

w = where(abs( photoobj.rowc[2] - sexcat.y_image ) le tolerance and $
		abs(photoobj.colc[2] - sexcat.x_image) le tolerance)

if (n_elements(w) eq 1) then begin
	if (w[0] eq -1) then  begin 
		return
	endif else if (n_elements(matches) eq 1 and matches[0].photoindex eq -1) then begin
		matches.photoindex = photoindex
		matches.sexindex = w[0]
	endif else begin
		tmp = replicate(match_struct,1)
		tmp.photoindex = photoindex
		tmp.sexindex = w[0]
		matches = [matches,tmp]
	endelse
endif else begin
	print,'More that one match for ',photoindex
endelse

return
end