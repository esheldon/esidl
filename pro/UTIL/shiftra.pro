;+
; NAME:
;	shiftra
;
; CALLING SEQUENCE:
;	1) ra_new = shiftra(ra, ra_shift)
;	or
;	2) ra_new = shiftra(ra,/wrap)
;
; PURPOSE:
;	mode 1)
;		Shift ra by a given amount, making sure to wrap around the 0-360 
;		boundary
;			E.g. ra = ra-shift   if ra < 0 ra+=360
;	mode 2)
;		If /wrap is sent make points with ra 180-360 instead run from -180-0
;
;
; INPUTS:
;	ra: in degrees.  Can be an array.
;	shift: in degrees
; KEYWORDS:
;	/wrap: If set make ra run from -180,180 with no shifting otherwise. Ignore
;		the shift argument
;
; OUTPUTS:
;	The shifted ra
; MODIFICATION HISTORY:
;	2009-06-12: From other code, Erin Sheldon, BNL
;	2009-07-30: Added /wrap keywrod.  Erin Sheldon, BNL
;
;-
function shiftra, ra, ra_shift, wrap=wrap
	if n_elements(ra) eq 0 or $
			(not keyword_set(wrap) and n_elements(ra_shift) eq 0) then begin
		on_error,2
		print,'Usage: new_ra = shiftra(ra, shift, /wrap)'
		message,'Halting'
	endif

	if keyword_set(wrap) then begin
		ranew = ra
		w=where(ra gt 180d, nw)
		if nw ne 0 then begin
			ranew[w] -= 360
		endif
	endif else begin
		ranew = ra - ra_shift
		w=where(ranew lt 0.0, nw)
		if nw ne 0 then ranew[w] += 360
	endelse
	return, ranew
end
