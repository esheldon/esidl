function _get_dims, scale, width, height, i
    if n_elements(scale) eq 1 then sc=scale[0] else sc=scale[i]
    if n_elements(width) eq 1 then w=width[0] else w=width[i]
    if n_elements(height) eq 1 then h=height[0] else h=height[i]
    return,{scale:ntostr(sc), width:ntostr(long(w)), height:ntostr(long(h))}
end

pro casimagelinks, ra, dec, scale=scale, width=width, height=height, extra_columns=columns, embed=embed, headmod=headmod, file=file

    if n_elements(scale) eq 0 then scale=0.85
    if n_elements(width) eq 0 then width=1024
    if n_elements(height) eq 0 then height=1024
    if n_elements(headmod) eq 0 then headmod=10

    main = 'http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?'
    
    if n_elements(file) ne 0 then begin
        openw, lun, file, /get_lun
    endif else begin
        lun=-1
    endelse

    printf, lun, '<html>'
    printf, lun, '<body>'
    printf, lun, '<table border=1>'

    hline = '  <tr><th>RA</th><th>DEC</th>'
    if n_elements(columns) ne 0 then begin
        tags=strlowcase(tag_names(columns)) & ntags=n_elements(tags)
        hline = hline + '<th>'+strjoin(tags, '</th><th>')+'</th>'
    endif else begin
        hline = hline + '<th>RA</th><th>DEC</th><th>link</th>'
    endelse
    hline = hline + '<th>link</th>'
    hline = hline+'</tr>'

    printf, lun, hline

    nra=n_elements(ra)
    for i=0L, nra-1 do begin
        rastr = ntostr(ra[i])
        decstr = ntostr(dec[i])

        t=_get_dims(scale, width, height, i)

        q='ra='+rastr
        q=q+'&dec='+decstr
        q=q+'&scale='+t.scale
        q=q+'&width='+t.width
        q=q+'&height='+t.height
        q=q+'&opt=G'

        url = main+q

        if keyword_set(embed) then begin
            url = '<img src="'+url+'">'
        endif else begin
            url = '<a href="'+url+'">image</a>'
        endelse

        line = '  <tr><td>'+rastr+'</td><td>'+decstr+'</td>'
        if n_elements(columns) ne 0 then begin
            delvarx, larr
            for j=0L, ntags-1 do begin
                add_arrval, ntostr( columns[i].(j) ), larr
            endfor
            line=line+'<td>' +strjoin(larr, '</td><td>') + '</td>'
        endif 

        line = line+'<td>'+url+'</td>'
        line = line + '</tr>'

        printf, lun, line

        if ( ((i+1) mod headmod) eq 0 ) then printf, lun, hline
    endfor

    printf, lun, '</table>'
    printf, lun, '</body>'
    printf, lun, '</html>'

    if lun ne -1 then free_lun, lun

end
