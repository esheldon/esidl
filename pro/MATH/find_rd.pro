pro find_rd, bt, rd

;finds the ratio of (Rd/Re) assuming that Id = Ie
; Where rd=r devauc and rd = r exponential

if (n_params() eq 0) then begin
	print, '-syntax: find_rd, bt, rd'
	return
endif

rd = sqrt( (1.0/0.28)*(1 - bt)/(bt) )

return
end
