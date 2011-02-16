function zweight::init, sample
    if n_elements(sample) eq 0 then begin
        message,'You must enter a sample number'
        return, 0
    endif

    ps = { $
            sample: 0, $
            ssample: 'dr5', $
            psample: 'p6', $
            nneigh: 0 $
        }

    case sample of
        0: begin
            ps.nneigh = 20
        end
        1: begin
            ps.nneigh = 5
        end
        2: begin
            ps.psample = 'p6z30'
            ps.nneigh = 5
        end
        else: message,'Unknown sample number: '+ntostr(sample)
    endcase

    ps.sample = sample

    print_struct, ps
    self.pstruct = ptr_new(ps, /no_copy)
    self.project = 'weight'

    return, 1

end

function zweight::par_struct
    return, *self.pstruct
end






function zweight::root, plot=plot
    if keyword_set(plot) then begin
        root = esheldon_config('plot_dir')
        root = path_join(root, 'photoz/uchicago')
        return, root
    endif else begin
        return,sdssidl_config('photoz_dir')
    endelse
end

function zweight::ext, ftype
    if n_elements(ftype) eq 0 then return, '.dat'
    case ftype of
        'bindefs': return, '.st'
        'script': return, '.sh'
        'pbs': return, '.pbs'
        'zbinhist': return, '.st'
        'meanscinv': return, '.bin'
        else: return, '.dat'
    endcase
end

function zweight::samplename, dtype
    ps = self->par_struct()
    case strlowcase(dtype) of
        'pinput': el = ps.psample
        'sinput': el = ps.ssample
        'output': el=[ps.ssample, ps.psample, 'nn'+string(ps.nneigh,f='(i02)')]
        else: message,'Unknown dtype: '+ntostr(dtype)
    endcase 

    sname = strjoin(el, '-')
    return, sname
        
end

function zweight::dir, dtype, ftype, project=project, createdir=createdir
    tf=self->file(dtype, ftype, project=project, createdir=createdir, dir=dir)
    return, dir
end
function zweight::file, dtype, ftype, number=number, project=project, dir=dir, createdir=createdir

    on_error, 2
    if n_elements(dtype) eq 0 then begin
        print,'-Syntax: files = zw->file(dtype [,ftype, number=, project=, dir=, /createdir])'
        print,' dtype = pinput|sinput|output'
        message,'Halting'
    endif

    if n_elements(project) eq 0 then project=self.project

    ext=self->ext(ftype)
    root=self->root()

    sname = self->samplename(dtype)
    pathel = [dtype, sname]

    if n_elements(ftype) ne 0 then begin
        pathel = [pathel, ftype]
    endif

    file = datafile(pathel, root=root, project=project, $
                    number=number, ext=ext, dir=dir, createdir=createdir)
    return, file
end
 

function zweight::plotfile, ftype, number=number, dir=dir, project=project, createdir=createdir, encapsulated=encapsulated

;    on_error, 2
;    if n_elements(dtype) eq 0 then begin
;        print,'-Syntax: files = zw->file(dtype [,ftype, number=, project=, dir=, /createdir])'
;        print,' dtype = pinput|sinput|output'
;        message,'Halting'
;    endif

    if keyword_set(encapsulated) then begin
        ext = '.eps'
    endif else begin
        ext = '.ps'
    endelse
        
    root=self->root(/plot)

    if n_elements(project) eq 0 then project=self.project

    sname = self->samplename('output')
    pathel = sname

    if n_elements(ftype) ne 0 then begin
        pathel = [pathel, ftype]
    endif

    file = datafile(pathel, root=root, project=project, $
                    number=number, ext=ext, dir=dir, createdir=createdir)
    return, file


end



function zweight::input_structdef, n
    st = { $
            zs: 0.0, $
            zp: 0.0, $
            cm: fltarr(5) $
        }
    if n_elements(n) ne 0 then begin
        st = replicate(st, n)
    endif
    return, st
end
function zweight::output_structdef, n
    st = { $
            zs: 0.0, $
            zp: 0.0, $
            w: 0.0, $
            cm: fltarr(5) $
        }
    if n_elements(n) ne 0 then begin
        st = replicate(st, n)
    endif
    return, st
end


function zweight::datread, type, file
    if type eq 'input' then sdef = self->input_structdef() else sdef=self->output_structdef()
    n = numlines(file)
    read_struct, file, sdef, struct
    return, struct
end
function zweight::read, dtype, ftype, number=number, silent=silent
    f = self->file(dtype, ftype, number=number)
    if not keyword_set(silent) then print,'Reading file: ',f
    if n_elements(ftype) eq 0 then begin
        st = self->datread('input',f)
    endif else begin
        case ftype of
            'bindefs': st = read_idlstruct(f,silent=silent)
            'zbinhist': st = read_idlstruct(f,silent=silent)
            'raw': st = self->datread('output',f)
            else: message,'Dont support ftype = '+ntostr(ftype)+' yet'
        endcase
    endelse
    return, st
end







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Code for sub-samples
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; marcos says a million is about the right number
function zweight::ninput
    return, 1000000
end
function zweight::zrange
    return, [0.1, 0.9]
end

pro zweight::calculate_zbins, cat=cat
    if n_elements(cat) eq 0 then cat=self->read_cat()
    zrange = self->zrange()
    nperbin = self->ninput()

    print,'Binning into nperbin='+ntostr(nperbin)
    bs = binner(cat.zp, nperbin=nperbin, min=zrange[0], max=zrange[1])

    file = self->file('pinput', 'bindefs', /createdir)
    print,'Writing bindefs file: ',file
    write_idlstruct, bs, file, /ascii
end
pro zweight::zbins, zmin, zmax
    st = self->read('pinput', 'bindefs')
    zmin = st.xlow
    zmax = st.xhigh
end
function zweight::zbin_where_string, nbin

    self->zbins, zmin, zmax
    nbin = n_elements(zmin)
    for i=0L, nbin-1 do begin
        zminstr = ntostr(zmin[i])
        zmaxstr = ntostr(zmax[i])

        ts = '(var.zp ge '+zminstr+ $
              ' and var.zp lt '+zmaxstr+')'
        add_arrval, ts, ws
    endfor    
    return, ws
end

function zweight::where_string, sample
    if n_elements(sample) eq 0 then begin
        ps = self->par_struct()
        sample=ps.psample
    endif

    case sample of
        'p6z30': begin
            ws = self->zbin_where_string()
        end
        else: message,'Unknown psample: '+ntostr(sample)
    endcase

    return, ws
end

function zweight::nbin, sample
    ws=self->where_string(sample)
    return, n_elements(ws)
end





function zweight::read_cat
    ps=self->par_struct()
    case ps.psample of
        'p6': table = 'scat_princeton6'
        'p6z30': table = 'scat_princeton6'
        else: message,'Dont support psample '+ps.psample+' yet'
    endcase

    query = 'select photoz_z as zp, counts_model as cm from '+table
    print,query
    struct = pgsql_query(query)

    return, struct
end

pro zweight::setup, cat=cat

    ; read in the source catalog
    if n_elements(cat) eq 0 then begin
        cat = self->read_cat() 
    endif

    ; how many do we want in the input file?
    ninput = self->ninput()

    ws = self->where_string()
    nn=n_elements(ws)
    for i=0L, nn-1 do begin

        print,ws[i]
        keep = where_select(cat, ws[i], nkeep)

        if nkeep eq 0 then message,'no objects passed cut'

        ; Want at least 90% of requesed ninput
        nlowest = long(0.9*ninput)
        if nkeep lt ninput then begin
            if nkeep lt nlowest then begin
                message,'Found '+ntostr(nkeep)+$
                    ' but need at least 90% of '+ntostr(ninput)
            endif
        endif else begin
            nkeep = ninput
        endelse

        print,'Extracting first '+ntostr(nkeep)
        keep = keep[0:nkeep-1]
        tcat = cat[keep]

        file = self->file('pinput', number=i, /createdir)

        st = self->input_structdef(nkeep)
        copy_struct, tcat, st

        print,'Writing to file: ',file
        ascii_write, st, file, status=status
        if status ne 0 then message,'Could not write file'

    endfor 
end


function zweight::convert2mafalda, file
    root = self->root()
    replace = '/home/esheldon'
    nfile = repstr(file, root, replace)
    return, nfile
end
pro zweight::write_scripts, number=number, mafalda=mafalda

    ps = self->par_struct()
    files = self->file('output', 'script', number=number, /createdir)
    nf=n_elements(files)
    for i=0L,  nf-1 do begin

        file = files[i]
        comm = '/usr/bin/time -p weightgal'
        sf = self->file('sinput')
        pf = self->file('pinput',number=number[i])
        ; will have to run at least one ahead of time with dotrain=1
        ; for each spectro file. I only have one of these now 2007-04-25
        dotrain = '0'
        docolor = '0'
        fixvol = '0'
        nneigh = ntostr(ps.nneigh)
        output = self->file('output', 'raw', number=number[i], /createdir)

        if keyword_set(mafalda) then begin
            file = repstr(file, '.sh', '-mafalda.sh') 

            sf = self->convert2mafalda(sf) 
            pf = self->convert2mafalda(pf) 
            output = self->convert2mafalda(output) 
        endif


        print,'Writing to file: ',file
        openw, lun, file, /get_lun

        line = [comm,sf,pf,dotrain,docolor,fixvol,nneigh,output]
        line = strjoin(line, ' ')
        printf, lun, line


        address = 'erin.sheldon@gmail.com'
        printf,lun,'dt=`date`'
        mess = 'finished weightgal: '+file_basename(output)
        line = 'echo "'+mess+' ${dt}" | mail '+address+' -s "'+mess+'"'
        printf, lun, line

        free_lun, lun
    endfor
end

pro zweight::write_pbs, number=number

; just to be safe
    hstr = '80'

    pbsfiles = self->file('output', 'pbs', number=number, /createdir)
    scriptfiles = self->file('output', 'script', number=number, /createdir)
    nf = n_elements(pbsfiles)

    for i=0L, nf-1 do begin
        ; pbs file to write
        pbsfile = pbsfiles[i]
        print,pbsfile

        ; scriptfile name -on mafalda-
        scriptfile = scriptfiles[i]
        scriptfile = repstr(scriptfile, '.sh', '-mafalda.sh') 
        scriptfile = self->convert2mafalda(scriptfile)

        jobname = 'w'+ntostr(number[i])
        openw, lun, pbsfile, /get_lun

        printf, lun, "#!/bin/sh"
        printf, lun, "#PBS -V"               ; Keep environment
        printf, lun, "#PBS -l cput="+hstr+":0:0"   ; at most this many hours
        printf, lun, "#PBS -l mem=200mb"    
        printf, lun, "#PBS -N "+jobname      ; job name
        printf, lun, "#PBS -m ae"            ; send email on abort and end
        printf, lun
        printf, lun, "# The output file on scratch disk"
        printf, lun, "script="+scriptfile
        printf, lun, "output=${script}.log"
        printf, lun, 'bash $script &> $output'
        printf, lun

        free_lun, lun
    endfor



end




pro zweight::compare_legend, colors
    legend, ['p zphot', 'w spec zphot', 'w spec zspec'], $
        line=0, color=colors, $
        /right, box=0, charsize=1
end

; for comparing input/output.  When number is sent it currently assumes the
; binning is in zphot
pro zweight::compare, number, ws, ps, dops=dops, $
        nboot=nboot, zbin=zbin, zmin=zmin, zmax=zmax

    if n_elements(number) eq 0 then begin
        print,'-Syntax: zw->compare, number, [ws, ps, /dops, nboot=, zbin=, zmin=, zmax=]'
        on_error, 2
        message,'Halting'
    endif

    psfile = self->plotfile('compare', number=number, /createdir)
    if keyword_set(dops) then begplot,psfile,/color

    if n_elements(ws) eq 0 then begin
        ws = self->read('output','raw', number=number)
    endif
    if n_elements(ps) eq 0 then begin
        ps = self->read('pinput', number=number)
    endif


    ; Create index array for bootstrap samples
    ndata = n_elements(ws)

    ; For z binning
    if n_elements(zbin) eq 0 then zbin = 0.01
    if n_elements(zmin) eq 0 then zmin = 0.0
    if n_elements(zmax) eq 0 then zmax = 1.0

    if !d.name eq 'PS' then zpcolor=!blue else zpcolor=!green
    zpcolor=!darkgreen
    zscolor=!red


    ; This is just the zspec histogram with range in pz overplotted
    ; when number is sent
    yrange=[0,1.1]
    ylog=0
    plothist, ws.zs, weights=ws.w, peak=1, $
        aspect=!gratio, $
        ylog=ylog, yrange=yrange, ystyle=1, $
        xrange=[zmin,zmax], xstyle=1, $
        nboot=nboot, $
        bin=zbin, min=zmin, max=zmax, bhist=bhist, xtitle='zspec'
    
     if n_elements(number) ne 0 then begin
        self->zbins, zbmin, zbmax
        zbstr = string([zbmin[number],zbmax[number]], f='(F4.2)')
        mess = zbstr[0]+' < zphot < '+zbstr[1]
        legend, mess, /left, charsize=1

        oplot, [zbmin[number], zbmin[number]], [0, 10], line=2
        oplot, [zbmax[number], zbmax[number]], [0, 10], line=2
        legend, ['w zspec', 'zphot range'], line=[0,2], /right
    endif
 
    key=prompt_kbrd('hit a key')

    plothist, ps.zp, peak=1, yrange=[0,1.1], ystyle=1, $
        bin=zbin, min=zmin, max=zmax, $
        xtitle='z', ytitle='P(z)', aspect=!gratio, xrange=[0,zmax], /fill
    plothist, ws.zp, weights=ws.w, peak=1, nboot=nboot, /over,$
        bin=zbin, min=zmin, max=zmax, bhist=bhist, color=zpcolor, $
        line=0
    plothist, ws.zs, weights=ws.w, /overplot, peak=1, nboot=nboot, $
        bin=zbin, min=zmin, max=zmax, bhist=bhist, color=zscolor

    self->compare_legend, [!p.color, zpcolor, zscolor]

    if n_elements(number) ne 0 then begin
        bd = self->read('pinput', 'bindefs')
        self->zbins, zbmin, zbmax
        zbstr = string([zbmin[number],zbmax[number]], f='(F4.2)')
        mess = zbstr[0]+' < zphot < '+zbstr[1]
        legend, mess, /left, charsize=1
    endif

 
    key=prompt_kbrd('hit a key')

    !p.charsize=1
    erase & multiplot, [3,2], /square 

    mbin = 0.2
    mmin = 16. & mmax=28.
    xrange=[mmin,mmax]
    yrange=[0,1.2]

    for i=0,4 do begin
        if i gt 2 then xt='mag' else xt=''
        if i eq 0 or i eq 3 then yt='dN/dmag' else yt=''
        plothist, ps.cm[i], bin=mbin, min=mmin, max=mmax, peak=1, $
            xrange=xrange, xsty=1, yrange=yrange, ystyle=1, xtitle=xt, ytitle=yt
        plothist, ws.cm[i], bin=mbin, min=mmin, max=mmax, peak=1, $
            weights=ws.w, /overplot, color=zscolor, nboot=nboot
              
        legend, !colors[i], box=0, charsize=1
        if i ne 4 then multiplot
    endfor
    multiplot, /default

    key=prompt_kbrd('hit a key')

    cbin = [0.1,0.1,0.05,0.05]
    erase & multiplot, [2,2], /square
    for i=0,3 do begin
        if i gt 1 then xt='color' else xt=''
        if i eq 0 or i eq 2 then yt='dN/color' else yt=''

        color = ps.cm[i] - ps.cm[i+1]
        cmin=-2. & cmax=4.
        xrange=[cmin,cmax]
        yrange=[0,1.2]
        plothist, color, bin=cbin[i], min=cmin, max=cmax, peak=1, $
            xrange=xrange, xsty=1, yrange=yrange, ystyle=1, xtitle=xt, ytitle=yt
        color = ws.cm[i] - ws.cm[i+1]
        plothist, color, bin=cbin[i], min=cmin, max=cmax, peak=1, $
            weights=ws.w, /overplot, color=zscolor, nboot=nboot

        cstr = !colors[i]+' - '+!colors[i+1]
        legend, cstr, box=0, charsize=1
        if i ne 3 then multiplot
    endfor
    multiplot, /default

    if keyword_set(dops) then endplot
    return
 

end


; Create a file with all the weighted spec redshift histograms
pro zweight::create_zbinhist, zbin=zbin, zmin=zmin, zmax=zmax, nboot=nboot

    file = self->file('output', 'zbinhist', /createdir)
    print,'Will write to file: ',file
    if n_elements(zbin) eq 0 then zbin = 0.01
    if n_elements(zmin) eq 0 then zmin = 0.0
    if n_elements(zmax) eq 0 then zmax = 1.0
    if n_elements(nboot) eq 0 then nboot = 100

    self->zbins, zbmin, zbmax
    nbin = n_elements(zbmin)
    for i=0L, nbin-1 do begin

        print,'-----------------------------------------------------------'
        ps = self->read('pinput', number=i)
        ws = self->read('output','raw', number=i)
        bs = boot_hist(nboot, ws.zs, weight=ws.w, bin=zbin, min=zmin, max=zmax)

        if i eq 0 then begin
            st = create_struct('zpmean', 0.0, 'zpmin', 0.0, 'zpmax', 0.0, bs)
            binstruct = replicate(st, nbin)
        endif

        tbs = binstruct[i]
        copy_struct, bs, tbs

        tbs.zpmean = mean(ps.zp)
        tbs.zpmin = zbmin[i]
        tbs.zpmax = zbmax[i]
        binstruct[i] = tbs
    endfor

    print,'Writing to file: ',file
    write_idlstruct, binstruct, file

end

; for comparing input/output.  When number is sent it currently assumes the
; binning is in zphot
pro zweight::plothist_zbins, numbers, all=all,  $
        nboot=nboot, zbin=zbin, zmin=zmin, zmax=zmax, $
        norange=norange, nolegend=nolegend, color=color

    nn = n_elements(numbers)
    if nn eq 0 then begin
        print,'-Syntax: zw->plothist_zbins, numbers, [ nboot=, zbin=, zmin=, zmax=, /color]'
        on_error, 2
        message,'Halting'
    endif

    self->zbins, zbmin, zbmax

    if nn eq n_elements(zbmin) then begin
        nolegend=1
        norange=1
        nstr = 'all'
    endif else begin
        nstr = string(numbers, f='(i02)')
        nstr = strjoin(nstr, '-')
    endelse

    psfile = self->plotfile('hist-zbins', /createdir, /encap)
    psfile = repstr(psfile, '.eps', '-'+nstr+'.eps')
    if keyword_set(color) then psfile=repstr(psfile,'.eps','-color.eps')
    begplot,psfile,color=color,/encap

    ; For z binning
    if n_elements(zbin) eq 0 then zbin = 0.01
    if n_elements(zmin) eq 0 then zmin = 0.0
    if n_elements(zmax) eq 0 then zmax = 1.0


    charsize=1.5
    !p.charsize=1.5
    !x.thick=3
    !y.thick=3
    !p.thick=3

    yrange = [0, 25]

    if keyword_set(nolegend) then begin
        colors = make_rainbow(nn)
    endif else begin
        if keyword_set(color) then begin
            simpctable, colorlist=colors
        endif else begin
            loadct, 0
            colors = arrscl(findgen(nn), 0, 200)
        endelse
    endelse
    nclr = n_elements(colors)
    for i=0L, nn-1 do begin
        num = numbers[i]
        ws = self->read('output','raw', number=num)
        ps = self->read('pinput', number=num)

        if i ne 0 then overplot=1
        plothist, ws.zs, weights=ws.w, /norm, overplot=overplot, $
            aspect = !gratio, $
            bin=zbin, min=zmin, max=zmax, $
            ylog=ylog, yrange=yrange, ystyle=1, $
            xrange=[zmin,zmax], xstyle=1, $
            nboot=nboot, $
            xhist=xhist, yhist=yhist, histerr=histerr, $
            bhist=bhist, xtitle='zspec'

        print,max(yhist)
        pplot, xhist, yhist, yerr=histerr, psym=10, hat=0, color=colors[i], /overplot

        zbstr = string([zbmin[num],zbmax[num]], f='(F4.2)')
        mess = zbstr[0]+' < zphot < '+zbstr[1]
        add_arrval, mess, legstr

        if not keyword_set(norange) and i eq 0 then begin
            maxymin = interpol(yhist, xhist, zbmin[num])
            maxymax = interpol(yhist, xhist, zbmax[num])
            oplot, [zbmin[num], zbmin[num]], [0, maxymin], line=2
            oplot, [zbmax[num], zbmax[num]], [0, maxymax], line=2
            legend, ['w zspec', 'zphot range'], line=[0,2], /right, $
                charsize=charsize
        endif

    endfor

    if not keyword_set(nolegend) then begin
        legend, legstr, line=0, color=colors, /right, /center, $
            charsize=charsize*0.8
    endif
    endplot, /trim_bbox
end



function zweight::cleanup
    ptr_free, self.pstruct
    return, 1
end

pro zweight__define
    struct = {                 $
            zweight,           $
            sample: 0,         $
            project: '',       $
            pstruct: ptr_new() $
        }
            
end
