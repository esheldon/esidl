function icl::init, sample
    if n_elements(sample) eq 0 then sample=1
    self.sample = sample
    return, 1
end

function icl::sample
    return, self.sample
end
function icl::sample_string, sample=sample
    if n_elements(sample) eq 0 then sample=self->sample()
    return, 'sample'+string(sample, f='(i02)')
end

function icl::rootdir
    return, esheldon_config('icl_root')
end
function icl::basedir, btype, sample=sample
    if n_elements(btype) eq 0 then begin
        print,'-Syntax: dir = icl->dir(base_type)'
        print,' base_type = input|output|sim|gonzeles'
        on_error, 2
        message,'Halting'
    endif

    root = self->rootdir()
    case strlowcase(btype) of
        'input':    return, ddir(root=root, ['icl', 'input'])
        'output':   return, ddir(root=root, ['icl', 'output'])
        'sim':      return, ddir(root=root, ['icl', 'sim'])
        'gonzales': return, ddir(root=root, ['icl', 'gonzales'])
        else: message,'Unsupported type: '+ntostr(btype)+'. Try input|output|sim|gonzales'
    endcase
end


function icl::dir, btype, ftype, sample=sample, subtype=subtype
    tf=self->file(btype, ftype, sample=sample, subtype=subtype, dir=dir)
    return, dir
end
function icl::file, btype, ftype, sample=sample, subtype=subtype, bin=bin, all=all, nodir=nodir, dir=dir

    if n_elements(btype) eq 0 or n_elements(ftype) eq 0 then begin
        print,'-Syntax: f=icl->file(btype, ftype, sample=, subtype=, bin=, /all, /nodir, dir=)'
        on_error, 2
        message,'Halting'
    endif
    root = self->rootdir()

    add_arrval, 'icl', dfinputs
    add_arrval, strlowcase(btype), dfinputs

    ; The sample
    if n_elements(sample) eq 0 then begin
        sample = self->sample()
    endif
    sstr = self->sample_string(sample=sample)
    add_arrval, sstr, dfinputs

    ; sub-samples start a new heirarchy
    if n_elements(subtype) ne 0 then begin
        add_arrval, 'sub', dfinputs
        add_arrval, subtype, dfinputs
    endif

    add_arrval, ftype, dfinputs

    ;dfinputs = '"'+strjoin(dfinputs, '", "')+'"'

    nbin = n_elements(bin)
    if nbin ne 0 then begin
        postfix = string(bin, format='(i02)')
    endif else begin
        nbin = 1
        postfix=''
    endelse

    for i=0L, nbin-1 do begin
        pf = postfix[i]
        command = 'tfile = dfile(root=root, dfinputs, postfix=pf, ext=".fits", nodir=nodir, dir=dir)'
        ;print,command
        if not execute(command) then message,'Could not create file name'
        add_arrval, tfile, files
    endfor

    return, files
end

; Converted tables from gonzales et al.
function icl::gonzales_file, gtype
    if n_elements(gtype) eq 0 then begin
        print,'-Syntax: dir = icl->gonzales_file(gtype)'
        print,' gtype = sersic|2dev'
        on_error, 2
        message,'Halting'
    endif
    dir = self->basedir('gonzales')
    case strlowcase(gtype) of
        'sersic': return,concat_dir(dir, 'sersic.dat')
        '2dev': return,concat_dir(dir, '2dev.dat')
        else: message,'Unknown file type: '+ntostr(gtype)+' try sersic|2dev'
    endcase
end


; input file
function icl::input_file
    dir = self->basedir('input')
    return, concat_dir(dir, 'maxbcg_icl_sample.fits')
end





; Output files and simulated icl files
function icl::band2cstr, band
    cstr = ['u','g','r','i','z']
    errmess = 'Band must be 0,1,2,3,4 or '+strjoin(cstr,',')
    if size(band,/tname) eq 'STRING' then begin
        w=where(cstr eq strlowcase(band[0]), nw)
        if nw eq 0 then message,errmess
        return, cstr[w]
    endif else begin
        b = fix(band[0])
        if b lt 0 or b gt 4 then message,errmess
        return, cstr[b]
    endelse
end
function icl::ftype2blanton, type
    case strlowcase(type[0]) of
        'bcg': return, 'cen-'
        'all': return, 'all-'
        'none': return, ''
        else: message,'Unknown type: '+ntostr(type)+' try bcg|all|none|sim'
    endcase
end

function icl::output_prefix, ra, dec
    return,'mBCG-'+hogg_iau_name(ra,dec,'')
end
function icl::output_dir, ra, dec, sample=sample, prefix=prefix

    nra=n_elements(ra) & ndec=n_elements(dec)
    if nra eq 0 or ndec eq 0 or (nra ne ndec) then begin
        print,'-Syntax: dir = icl->dir(ra, dec, prefix=)'
        print,'   ftype = bcg|all|none  a scalar'
        print,'   ra/dec same length'
        on_error, 2
        message,'Halting'
    endif

    ihr=long(ra/15.)
    idec=long(abs(dec)/2.)*2.
    dsign=replicate('p', ndec)

    w=where(dec lt 0.0, nw)
    if nw ne 0 then dsign[w] = 'm'

    basedir = self->basedir('output')
    sstr=self->sample_string(sample=sample)
    dir=path_join(basedir, sstr)
    dir=concat_dir(dir, string(ihr,f='(i2.2)')+'h' )
    dir=concat_dir(dir, dsign+strtrim(string(idec, f='(i2.2)'),2) )
    prefix=self->output_prefix(ra, dec)
    dir=concat_dir(dir,prefix)

    return, dir
   
end
function icl::output_file, ftype, ra, dec, band
    nt=n_elements(ftype) 
    nra = n_elements(ra) & ndec = n_elements(dec) 
    nband=n_elements(band)
    if nt ne 1 or nra eq 0 or ndec eq 0 or nband ne 1 or (nra ne ndec) then begin
        print,'-Syntax: of = icl->output_file(ftype, ra, dec, band)'
        print,' ftype = bcg: all objects removed except BCG'
        print,'         all: all objects removed'
        print,'         none: no objects removed'
        print,' ra,dec must be same length'
        print," band must be a scalar 0,1,2,3,4 or 'u','g','r','i','z'"
        on_error, 2
        message,'Halting'
    endif

    dir = self->output_dir(ra, dec, prefix=prefix)
    fpref = self->ftype2blanton(ftype)
    file = fpref+prefix+'-'+self->band2cstr(band)+'.fits.gz'

    file = concat_dir(dir, file)
    return, file
 
end



function icl::stackdir, sample, createdir=createdir
    dir = self->basedir('stack')
    if keyword_set(createdir) then file_mkdir, dir
    return, dir
end
function icl::stackfile, sample, bin=bin

end

function icl::goodind, ftype, ra, dec, band

    nra=n_elements(ra)
    
    goodind = intarr(nra)

    files = self->output_file(ftype, ra, dec, band)
    for i=0L, nra-1 do begin
        if fexist(files[i]) then goodind[i] = 1
    endfor

    goodind = where(goodind,count)
    print,'Found: '+ntostr(count)+'/'+ntostr(nra)

    return, goodind
end


; stack some clusters.  Note this only makes sense in a very narrow redshift
; slice
function icl::stackbypos, ftype, ra, dec, band, wsum=wsum, winverse=winverse

    nra=n_elements(ra)

    files = self->output_file(ftype, ra, dec, band)
    ; first go through and get all image sizes
    print,'Reading headers'
    for i=0L, nra-1 do begin
        hdr = headfits(files[i], err=err)
        if err ne '' then begin
            ;message,'Failed to read header from file: '+files[i], /inf
        endif else begin
            add_arrval, sxpar(hdr, 'naxis1'), sx 
            add_arrval, sxpar(hdr, 'naxis2'), sy 
        endelse 
    endfor

    ngood = n_elements(sx)
    sxmax = max(sx, min=sxmin)
    symax = max(sy, min=symin)

    print,'    Found: '+ntostr(ngood)+'/'+ntostr(nra)
    print,'    Max size: ',sxmax,symax
    print,'    Min size: ',sxmin,symin

    imsum = dblarr(sxmax, symax)
    wsum = dblarr(sxmax, symax)
    for i=0L, nra-1 do begin
        im = icl_read('output',ftype, ra[i], dec[i], band, /silent, status=imstatus)
        if imstatus eq 0 then begin
            sz = size(im, /dim)
            sxi = sz[0]
            syi = sz[1]

            ; Michael put the centers at n/2 instead of (n-1)/2
            xmin = sxmax/2 - sxi/2
            xmax = sxmax/2 + sxi/2 > 0 < sxmax-1
            ymin = symax/2 - syi/2
            ymax = symax/2 + syi/2 > 0 < symax-1

            ; create the weight image from the zeros in the image
            wim = im
            w=where(im ne 0, nw, comp=comp, ncomp=ncomp)
            if nw ne 0 then wim[w] = 1.0
            if ncomp ne 0 then wim[comp] = 0.0

            ; copy in data and weight
            imsum[xmin:xmax, ymin:ymax] = $
                imsum[xmin:xmax, ymin:ymax] + im[*,*]
            wsum[xmin:xmax, ymin:ymax] = $
                wsum[xmin:xmax, ymin:ymax] + wim[*,*]
        endif else begin
            ;message,'Could not read file: '+files[i], /inf
        endelse
    endfor   

    winverse = wsum
    winverse[*] = 0.0
    w=where(wsum ne 0, nw)
    if nw ne 0 then begin
        winverse[w] = 1.0/wsum[w]
    endif
    imout = imsum*winverse

    return, imout
end

pro icl__define
    struct = {              $
            icl,            $
            sample: 0       $
        }
end
