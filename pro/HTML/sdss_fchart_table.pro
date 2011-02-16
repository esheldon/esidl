pro sdss_fchart_table, ra, dec, html_file, add_atlas=add_atlas, struct=st,  scale=scale, width=width, height=height, linksonly=linksonly, rebinfac=rebinfac, overwrite=overwrite

	if n_params() lt 3 then begin
		print,'usage: '
		print,'  sdss_fchart_table, ralist, declist, html_file, /add_atlas, struct=, cutout=, /linksonly'
		on_error, 2
		message,'halting'
	endif

	if n_elements(rebinfac) eq 0 then rebinfac=1

	openw, lun, html_file, /get_lun

	htmlbase=file_basename(html_file)

	if keyword_set(add_atlas) then begin
		image_dir = html_file+'.images'
		if not fexist(image_dir) then file_mkdir, image_dir
	endif

	printf, lun, '<html>'
	printf, lun, '<body>'

	printf, lun, '<table border=1>'
	
	header = "  <tr><th>"
	if keyword_set(linksonly) then begin
		elements = ["RA","DEC","fchart","navigate"]
	endif else begin
		elements = ["RA","DEC","fchart (click navigate)"]
	endelse
	if keyword_set(add_atlas) then begin
		elements = [elements,"atlas image","field"]
	endif

	header = header + strjoin(elements, "</th><th>") + "</th></tr>"
	printf, lun, header

	for i=0L, n_elements(ra)-1 do begin

		furl = sdss_fchart_url(ra[i], dec[i], $
			scale=scale, width=width, height=height)
		nurl = sdss_fchart_url(ra[i], dec[i],/nav, $
			scale=scale, width=width, height=height)

		if n_elements(cutout) eq 0 then begin
			cutout = width
		endif

		line = "  <tr><td>"

		if keyword_set(linksonly) then begin
			elements = [ntostr(ra[i]), ntostr(dec[i]), $
				"<a href='"+furl+"'>fchart</a>", $
				"<a href='"+nurl+"'>navigate</a>"]
		endif else begin
			elements = [ntostr(ra[i]), ntostr(dec[i]), $
				"<a href='"+nurl+"'><img src='"+furl+"'></a>"]
		endelse

		if keyword_set(add_atlas) then begin
			; don't re-create images if they exist unless /overwrite is sent
			afilebase=sdss_objname(st[i], front='atlas-', suffix='-all.jpg')
			filebase=sdss_objname(st[i], front='atlas-', suffix='.jpg')

			afile=path_join(image_dir, afilebase)
			file=path_join(image_dir, filebase)

			if (not fexist(afile)) or keyword_set(overwrite) then begin

				atlas2jpg, st[i], outdir=image_dir,/all,/sheldon,$
					cutout=cutout, rebin=rebinfac

			endif

			relfile = path_join('./'+htmlbase+'.images', filebase)
			relafile = path_join('./'+htmlbase+'.images', afilebase)
			if keyword_set(linksonly) then begin
				elements = [elements, $
					"<a href='./"+relfile+"'>atlas_image</a>", $
					"<a href='./"+relafile+"'>field</a>"]
			endif else begin
				elements = [elements, $
					"<img src='./"+relfile+"'>", $
					"<img src='./"+relafile+"'>"]
			endelse

		endif

		line = line + strjoin(elements, "</td<td>")+"</td></tr>"

		printf, lun, line
	endfor

	setupplot,'x'

	printf, lun, '</table>'
	printf, lun, '</body>'
	printf, lun, '</html>'

	free_lun, lun
	

end
