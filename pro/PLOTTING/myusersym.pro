pro myusersym, symname, symsize=symsize

    if n_params() lt 1 then begin 
        print,'myusersym, symname, symsize=symsize'
        print
        print,'valid symnames:'
        print,'  fill_triangle'
        print,'  open_circle'
        print,'  fill_circle'
        print,'  fill_diamond'
        print,'  vertical_line'
        return
    endif 

    if n_elements(symsize) eq 0 then symsize=1.0

    case strlowcase(symname) of 

        'fill_triangle': begin
            x=[-1,1,0,-1]*symsize
            y=[-sqrt(3)/2,-sqrt(3)/2,sqrt(3)/2,-sqrt(3)/2]*symsize
            usersym,x,y,/fill
        end 
        'open_circle': begin
            b=findgen(17) * (!pi*2/16.)
            usersym,cos(b)*symsize,sin(b)*symsize
        end 
        'fill_circle': begin
            b=findgen(17) * (!pi*2/16.)
            usersym,cos(b)*symsize,sin(b)*symsize,/fill
        end 
        'fill_diamond': begin
            x = [-1, 0, 1, 0, -1]*symsize
            y = [0, 1, 0, -1, 0]*symsize
            usersym,x,y,/fill
        end 
        'vertical_line': begin
            x = [0,0]
            y = [-1,1]*symsize
            usersym, x, y
        end
        else: print,"Don't know about symbol "+symname
    endcase 

end
