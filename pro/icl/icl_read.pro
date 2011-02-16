function icl_read, dtype, ftype, ra, dec, band, silent=silent, hdr=hdr, status=status

    if n_elements(dtype) eq 0 then begin
        print,'-Syntax: res = icl_read(dtype, ftype [, ra, dec, band, /silent, hdr=, status=])'
        print,' dtype = input|output'
        print,' For dtype=input the inputs fits file is read'
        print,' For dtype=output you must also enter an ftype and ra/dec/band to read'
        print,' ftype = bcg: all objects removed except BCG'
        print,'         all: all objects removed'
        print,'         none: no objects removed'
        print
        on_error, 2
        message,'Halting'
    endif

    icl = obj_new('icl')
    nd=n_elements(dtype) & nf=n_elements(ftype) 
    nra = n_elements(ra) & ndec = n_elements(dec)

    case strlowcase(dtype) of
        'input': begin
            file = icl->file('input',ftype)
            if not keyword_set(silent) then print,'Reading input structure: ',file
            result = mrdfits(file, 1, hdr, status=status)
        end
        'output': begin
            file = icl->output_file(ftype, ra, dec, band)
            if not keyword_set(silent) then print,'Reading output image: ',file
            result = mrdfits(file[0], 0, hdr, status=status)
        end
        else: message,'Unknown dtype: '+ntostr(dtype)+' try input|output'
    endcase

    obj_destroy, icl
    return, result
end
