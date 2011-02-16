FUNCTION sdss_fchart_url, ra, dec, scale=scale, width=width, height=height,$
                          opt=opt, query=query, navigator=navigator
                          
	defsize=400
	if n_elements(width) eq 0 then begin
		if n_elements(height) ne 0 then width=height else width=defsize
	endif
	if n_elements(height) eq 0 then begin
		height=width
	endif

	width_string = strtrim(string(width, f='(i)'), 2)
	height_string = strtrim(string(height, f='(i)'), 2)

	if n_elements(scale) eq 0 then begin
		scale_string = '0.396127'
	endif else begin
		scale_string = strtrim(string(scale[0]), 2)
	endelse


	ra_string = strtrim( string(ra[0], format='(f15.8)'), 2)
	dec_string = strtrim( string(dec[0], format='(f15.8)'), 2)

	if n_elements(opt) eq 0 then opt = 'G'
	if n_elements(query) eq 0 then query = ''
  
	if keyword_set(navigator) then begin 
		base_url = 'http://casjobs.sdss.org/dr7/en/tools/chart/navi.asp'
	endif else begin 
		base_url = 'http://casjobs.sdss.org/ImgCutoutDR7/getjpeg.aspx'
	endelse 

	url = base_url + '?'+$
		'ra='+ra_string         +$
		'&dec='+dec_string       +$
		'&scale='+scale_string   +$
		'&width='+width_string   +$
		'&height='+height_string +$
		'&opt='+opt              +$
		'&query='+query
		
  return,url

END 
