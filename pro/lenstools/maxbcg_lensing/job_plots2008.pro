pro job_plots2008, color=color
    ; This is a composite graph, probably move to separate file

    M2L = obj_new('maxbcg_m2l', [21,22], 4)
    m21 = obj_new('maxbcg_lensing', 21)

    cfac = 1.15
    
    if keyword_set(color) then begin
        plotfile = M2L->plotfile('ngals200_12', 'job-colloqium-m-r-color.eps')
    endif else begin
        plotfile = M2L->plotfile('ngals200_12', 'job-colloqium-m-r.eps')
    endelse
    begplot, plotfile, /encap, /color
    !p.charsize=!p.charsize*cfac

    t = mrdfits(M2L->m2lfile('ngals200_12'),1)

    col1_xtitle = textoidl('r [h^{-1} Mpc]')
    col2_xtitle = textoidl('N_{200}')
    col1_xrange = [0.015, 30.0]
    col2_xrange = [2,200]

    n=n_elements(t)

    if keyword_set(color) then begin
        colors = make_rainbow(n)
        ;colors = c2i(/sample)
    endif else begin
        ; greyscale
        loadct, 0
        colors = reverse(long(arrscl(findgen(n), 0, 225)))
    endelse

    labels = m21->labels('ngals200_12')

    ;
    ; M(<r)
    ;
    yrange = [5.e10, 1.e16]
    for i=0L, n-1 do begin

        ;if (i eq 0) or (i eq (n-1)) then begin
        ;if i eq 7 then begin
        ;    yerr=t[i].mass_err 
        ;    psym=-8
        ;    delvarx, line
        ;endif else begin
        ;    delvarx, yerr, psym
            line=0
        ;endelse
        if i eq 0 then begin
            pplot, t[i].r, t[i].mass, yerr=yerr, $
                aspect=1, $
                psym=psym, line=line, $
                xrange=col1_xrange, xstyle=3, /xlog, $
                yrange=yrange, ystyle=3, /ylog, $
                xticklen=0.04, yticklen=0.04, $
                ytitle = textoidl('M(<r) [h^{-1} M_{\odot}]'), $
                xtitle=col1_xtitle
        endif
        if i gt (n-5) then thick = 5 else thick=!p.thick
        pplot, t[i].r, t[i].mass, yerr=yerr, $
            psym=psym, line=line, thick=thick, $
            color=colors[i], /overplot, symsize=0.5, hat=0


    endfor
    pplot, t.r200, t.m200, psym=2, /overplot, symsize=1.15

    endplot, /trim

    ;
    ; M200 vs N200
    ;


    plotfile = M2L->plotfile('ngals200_12', 'job-colloqium-m-n.eps')
    begplot, plotfile, /encap, /color
    !p.charsize=!p.charsize*cfac

    yrange = [2.e12, 1.e15]
    ytitle = textoidl('M_{200} [h^{-1} M_{\odot}]')
    pplot, t.n200_red, t.m200, yerr=t.m200_err*2, psym=8, $
        aspect=1, $
        xrange=col2_xrange, xstyle=3, /xlog, $
        xticklen=0.04, yticklen=0.04, $
        yrange=yrange, ystyle=3+8, /ylog, symsize=1, hat=0, $
        ytitle=ytitle, xtitle=col2_xtitle
    ;axis, yaxis=1, yrange=yrange, ystyle=1, ytitle=ytitle, yticklen=0.04, $
    ;    ytickf='loglabels'

    ; add top axis
    ;rlin = 10.0^!y.crange
    ;toprange = fitst.norm*rlin^fitst.index
    ;axis, yaxis=1, xrange=toprange, xstyle=1, /xlog, $
    ;    xtitle=xtitle_top+' !c', charsize=2

    endplot, /trim

    ;key=prompt_kbrd('hit a key')

    ;
    ; M/L(<r)
    ;

    if keyword_set(color) then begin
        plotfile = M2L->plotfile('ngals200_12', 'job-colloqium-m2l-r-color.eps')
    endif else begin
        plotfile = M2L->plotfile('ngals200_12', 'job-colloqium-m2l-r.eps')
    endelse
    begplot, plotfile, /encap, /color
    !p.charsize=!p.charsize*cfac

    if keyword_set(color) then begin
        colors = make_rainbow(n)
        ;colors = c2i(/sample)
    endif else begin
        ; greyscale
        loadct, 0
        colors = reverse(long(arrscl(findgen(n), 0, 225)))
    endelse



    yrange = [5,1000]
    for i=0L, n-1 do begin

        ;if (i eq 0) or (i eq (n-1)) then begin
        ;if i eq 7 then begin
        ;    yerr=t[i].m2ltot_err 
        ;    psym=-8
        ;    delvarx, line
        ;endif else begin
        ;    delvarx, yerr, psym
            line=0
        ;endelse
        if i eq 0 then begin
            pplot, t[i].r, t[i].m2ltot, yerr=yerr, $
                aspect=1, $
                psym=psym, line=line, $
                xrange=col1_xrange, xstyle=3, /xlog, $
                yrange=yrange, ystyle=3, /ylog, $
                xticklen=0.04, yticklen=0.04, $
                xtickf='loglabels', ytickf='loglabels', $
                ytitle = textoidl('(M/L)(<r) [h M_{\odot}/L_{\odot}]'), $
                xtitle=col1_xtitle
        endif
        if i gt (n-5) then thick = 5 else thick=!p.thick
        pplot, t[i].r, t[i].m2ltot, yerr=yerr, $
            psym=psym, line=line, thick=thick, $
            color=colors[i], /overplot

    endfor

    pplot, t.r200, t.m2l200, psym=2, /overplot, symsize=1.15

    ;legend, labels, line=0, color=colors, /right, /bottom, charsize=0.5

    endplot, /trim



    ;key=prompt_kbrd('hit a key')

    ;
    ; M200/L200 vs N200
    ;

    plotfile = M2L->plotfile('ngals200_12', 'job-colloqium-m2l-n.eps')
    begplot, plotfile, /encap, /color
    !p.charsize=!p.charsize*cfac


    yrange = [0,600]
    pplot, t.n200_red, t.m2l200, yerr=t.m2l200_err, psym=8, $
        aspect=1, $
        xtitle=col2_xtitle, $
        xrange=col2_xrange, xstyle=3, /xlog, $
        xticklen=0.04, yticklen=0.04, $
        yrange=yrange, ystyle=3+8, symsize=1, hat=0, $
        ytitle = textoidl('(M/L)_{200} [h M_{\odot}/L_{\odot}]')

    ;pplot, t.n200_red*1.05, t.m2lmax, yerr=t.m2lmax_err, psym=4, $
    ;    /overplot, color='darkgreen'
 
    ;ytitle = textoidl('(M/L)_{200} [h M_{\odot}/L_{\odot}]')
    ;axis, yaxis=1, yrange=yrange, ystyle=1, ytitle=ytitle, yticklen=0.04, $
    ;    ytickf='(i)'


    endplot, /trim

    obj_destroy, M2L
end

