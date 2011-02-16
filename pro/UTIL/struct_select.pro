function struct_select, str, where_string, nkeep, verbose=verbose

    nkeep = 0L
    nstr=n_elements(str) & nselect = n_elements(where_string)
    if nstr eq 0 or nselect eq 0 then begin 
        print,'-Syntax: keep = struct_select(str, where_string, nkeep, /verbose)'
        print
        print,'where_string is a string sent to the where() function.  It must'
        print,'  refer to "str"'
        print,'e.g. where_string = "str.z gt 0.05 and str.x le 23"'
        return,-1
    endif 
  
    select_command = 'keep = where('+where_string[0]+', nkeep)'
	if keyword_set(verbose) then begin
		print,'selecting: "'+select_command+"'"
	endif
    if not execute(select_command) then begin 
        message,'Error executing select statement'
    endif 
    return,keep

end 




