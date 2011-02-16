pro maxbcg_lensing_compare_samples, sample1, sample2, type, subtype=subtype, add_labels=add_labels

    if n_params() lt 3 then begin
        print,'-Syntax: maxbcg_lensing_compare_samples, sample1, sample2, type, subtype='
        on_error, 2
        message,'Halting'
    endif

    m1 = obj_new('maxbcg_lensing', sample1)
    m2 = obj_new('maxbcg_lensing', sample2)

    t1 = m1->lensread(type, subtype=subtype) 
    t2 = m2->lensread(type, subtype=subtype) 

    outdir = esheldon_config('plot_dir')
    outdir = concat_dir(outdir, 'maxbcg/compare_samples')
    fmt='(I02)'
    tstr = strlowcase(type)
    s1 = ntostr(sample1,form=fmt)
    s2 = ntostr(sample2,form=fmt)

    if n_elements(subtype) eq 0 then begin
        sstr='' 
    endif else begin
        ws=m1->where_string(subtype, labels=labels)
        sstr='-'+subtype
    endelse
    file = 'maxbcg-compare-'+tstr+'-'+s1+'-'+s2+sstr+'.eps'
    file = concat_dir(outdir, file)
    begplot, file, /color, /encapsulated
    !p.multi=0

    n=n_elements(t1)
    colors=make_rainbow(n)

    if keyword_set(add_labels) and n_elements(labels) ne 0 then begin
        xrange = [0.01, 500.0]
    endif else begin
        xrange = [0.01, 30.0]
    endelse
    yrange = [0.8,1.2]
    allrat = t2.sigma/t1.sigma
    sigma_clip, allrat, rmean, rsig
    yrange = rmean + rsig*2*[-1,1]
    ytit='sample'+s2+'/sample'+s1
    xtit=!mpcxtitle2
 
    rst = {ratio: t1[0].sigma, ratio_err: t1[0].sigma}
    rst = replicate(rst, n)
    nrad = n_elements(t1[0].meanr)

    for i=0L, n-1 do begin

        ratio = t2[i].sigma/t1[i].sigma
        rerr = ratio*sqrt( (t1[i].sigmaerr/t1[i].sigma)^2 + (t2[i].sigmaerr/t2[i].sigma)^2 )

        if i eq 0 then begin
            pplot, t1[i].meanr/1000, ratio, aspect=!gratio, /xlog, $
                xtit=xtit, ytit=ytit, $
                xrange=xrange, xstyle=3, $
                yrange=yrange, ystyle=3;, yerr=rerr
        endif else begin
            pplot, t1[i].meanr/1000, ratio, color=colors[i], /overplot;, yerr=rerr
        endelse

        rst[i].ratio = ratio
        rst[i].ratio_err = rerr 

        wmom, ratio, rerr, wm, ws, we
        print,wm,we

        add_arrval, wm, ratiomean
        add_arrval, we, ratioerr
        
    endfor
    obj_destroy, m1, m2

    if keyword_set(add_labels) then begin
        legend, labels, /right, color=colors, line=0, charsize=1
    endif

    for j=0L, nrad-1 do begin
        wmom, rst.ratio[j], rst.ratio_err[j], wm, ws, we
        add_arrval, wm, mratio
        add_arrval, we, mratio_err
    endfor

;    print,wwm, wwerr
    pplot, t1[0].meanr/1000.0, mratio, yerr=mratio_err, psym=-8, /over, thick=2

    print,'----------------------------------------------'
    wmom, ratiomean, ratioerr, wm, ws, we
    print,wm,we

    minr = min(t1[0].meanr/1000.0, max=maxr)
    miny = wm-we
    maxy = wm+we
    plot_box, minr, maxr, miny, maxy, /polyfill, /line_fill, orient=45

    fo = '(2(f4.2))'
    mess = [type, 'ratio = '+ntostr(wm,form=fo)+!csym.plusminus+ntostr(we,form=fo)]
    legend, mess, /bottom, /right, box=0

    endplot, /trim_bbox
end
