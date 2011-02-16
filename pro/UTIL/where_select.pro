function where_select, var, where_string, nkeep

    nkeep = 0L
    nvar=n_elements(var) & nselect = n_elements(where_string)
    if nvar eq 0 or nselect eq 0 then begin 
        print,'-Syntax: keep = where_select(var, where_string, nkeep)'
        print
        print,'where_string is a string sent to the where() function.  It must'
        print,'  refer to "var"'
        print,'e.g. where_string = "var.z gt 0.05 and var.x le 23"'
        return,-1
    endif 
  
    select_command = 'keep = where('+where_string[0]+', nkeep)'
    if not execute(select_command) then begin 
        message,'Error executing select statement'
    endif 
    return,keep

end 




