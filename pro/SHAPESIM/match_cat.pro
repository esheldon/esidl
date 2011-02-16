pro match_cat, photocat, sexcat, goodphoto, matches,stars=stars

if (n_params() eq 0) then begin
	print,'Syntax: match_cat, photocat, sexcat, goodphoto,matches,stars=stars'
	return
endif

print, 'Number of photo objects: ',n_elements(photocat)
print, 'Number of sextractor objects: ',n_elements(sexcat)


;Make strong cuts, don't know how sextracter deblends

if keyword_set(stars) then begin
	extract_stars, photocat, 2, goodphoto
endif else begin
	cut_bad_photo, photocat, indexes
	goodphoto = photocat[indexes]
endelse



print,'Number of good photo objects found: ',n_elements(goodphoto)

match_struct = create_struct(name='mtch', 'photoindex',-1,'sexindex',-1)
matches = match_struct

for i=0,n_elements(goodphoto)-1 do begin
	match_pos, goodphoto[i],i, sexcat, match_struct, matches
endfor

print,n_elements(matches), '  matches found'

return
end


