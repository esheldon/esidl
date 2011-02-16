function struct_stats, struct, tags, $
          weights=weights, index=index, status=status

    if n_params() lt 2 then begin
        print,'-Syntax: ss=struct_stats(struct, tags, weights=, index=, status=)'
        print
        message,'Halting'
    endif

    status=1
    ns=n_elements(struct)
    if n_elements(weights) eq 0 then weights=replicate(1.0, ns)
    if n_elements(index) eq 0 then index=lindgen(ns)

    nt=n_elements(tags)

    statst = {tagname: '', mean: 0d, sdev: 0d, err: 0d}
    statst = replicate(statst, nt)
   
    nkeep=0L 
    for i=0L, nt-1 do begin

        ;; strip trailing [index] or (index)
        bracket_pos  = stregex(tags[i], '[\[\(]')
        bracket_pos2 = stregex(tags[i], '[]\)]')
        if bracket_pos ne -1 then begin 

            tag = strmid(tags[i], 0, bracket_pos)

            ;; Now extract the element for our name
            element_str = '_'+strmid(tags[i], $
                                    bracket_pos+1, bracket_pos2-bracket_pos-1)
        endif else begin 
            tag = tags[i]
            element_str = ''
        endelse 

        if tag_exist(struct[0], tag) then begin 

            add_arrval, i, keep
            nkeep=nkeep+1

            wtot = total( weights[index] )
            mean_command = $
                'mean_tag = total( weights[index]*struct[index].'+tags[i]+' )/wtot'
            sdev_command = $
                'sdev_tag = total( weights[index]*(struct[index].'+tags[i]+'-mean_tag)^2 )/wtot'
            err_command = $
                'err_tag = sqrt( total( weights[index]^2*(struct[index].'+tags[i]+'-mean_tag)^2 )/wtot^2 > 0)'

            if not execute(mean_command) then message,'Unable to take mean of tag '+tags[i]
            if not execute(sdev_command) then message,'Unable to take sdev of tag '+tags[i]
            if not execute(err_command)  then message,'Unable to take err  of tag '+tags[i]

            
            statst[i].tagname = tag
            statst[i].mean = mean_tag
            statst[i].sdev = sdev_tag
            statst[i].err = err_tag

        endif else begin 
            message,'Tag '+ntostr(tag)+' does not exist in structure',/inf
        endelse 

    endfor

    if nkeep gt 0 then begin
        status=0
        statst = statst[keep]
        return, statst
    endif else begin
        return, -1
    endelse
    
end 
