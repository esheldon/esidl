pro _get_mpvalue, n, mp, xwhen, ywhen

    case n of
        12: begin
            mp = [4,3]
        end
        16: begin
            mp = [4,4]
        end
        else: message,'bad n'
    endcase
    xwhen = (lindgen(n)/mp[0]) ge (mp[1]-1)
    ywhen = (lindgen(n) mod mp[0]) eq 0
end


;m200 = 2.923e14 h^{-1} Msun (r200/( h^{-1} Mpc)^3
pro maxbcg_lensing_plot_m2l, subtype

    if n_elements(subtype) eq 0 then begin
        print,'maxbcg_lensing_plot_m2l, subtype, /allfour'
        return
    endif

    ; hard wiring the samples for simplicity
    mass_sample = [21,22]
    lsample = 4

    ; build the basic dir and file names
    dir = '~/plots/m2l/maxbcg/'+subtype
    if not fexist(dir) then file_mkdir, dir

    mstr = strjoin( ntostr(mass_sample), '-')
    lstr = ntostr(lsample)

    front = 'm2l-'+subtype+'-m'+mstr+'-l'+lstr


    output_name = path_join(dir, front+'-m2l.fits')
    print,'Reading ',output_name
    m2lstruct = mrdfits(output_name,1)

    ;
    ; M(r) plots and fits overplotted
    ;

    n=n_elements(m2lstruct)
    _get_mpvalue,n, mp, xwhen, ywhen

    name = path_join(dir, front+'-massfits.eps')
    begplot, name, /color, /encap
    !p.charsize=1
    !p.thick=2
    !x.thick=2
    !y.thick=2
    erase & multiplot, mp, /square, $
        mxTitle=!mpcxtitle2, $
        myTitle='M(<r) [ h!U-1!N M'+sunsymbol()+' ]'

    for i=0L, n-1 do begin

        if xwhen[i] then xtickf='loglabels' else xtickf=''
        if ywhen[i] then ytickf='loglabels' else ytickf=''
        pplot, m2lstruct[i].r, m2lstruct[i].mass, yerr=m2lstruct[i].mass_err, $
            psym=8, symsize=0.5, hat=0, /xlog, /ylog, $
            xrange=[0.01, 60.0], yrange=[0.1,40000]*1.e12, $
            xstyle=3, ystyle=3, $
            xticklen=0.04, yticklen=0.04, $
            xtickf=xtickf, ytickf=ytickf

        pplot, m2lstruct[i].r, m2lstruct[i].nfw_yfit, color=!orange, /over
        pplot, m2lstruct[i].r, m2lstruct[i].linear_yfit, color=!blue, /over
        pplot, m2lstruct[i].r, m2lstruct[i].mass_yfit, color=!DarkGreen, /over

        pplot, [m2lstruct[i].r200], [m2lstruct[i].m200], $
            psym=8,  color=!red, /over
        multiplot
    endfor

    multiplot,/default
    endplot,/trim_bbox


return
    ;
    ; M/L at r200 vs. M200
    ;
    name = path_join(dir, front+'-m2l200-vs-m200.eps')
    begplot, name, /encap
    !p.multi=0

    xtitle='M!D200!N [h!U-1!N M'+sunsymbol()+']'
    ytitle = 'M/L!D200!N [h M'+sunsymbol()+'/L'+sunsymbol()+']'
    pplot, m2lstruct.m200, m2lstruct.m2l200, $
        xerr=m2lstruct.m200_err, yerr=m2lstruct.m2l200_err, psym=8, $
        aspect=1, /xlog, $
        xtitle=xtitle, ytitle=ytitle

    key = prompt_kbrd('hit a key')
    if key eq 'q' then return

    endplot,/trim



    ;
    ; M/L at the last point as a function of M200
    ;
    name = path_join(dir, front+'-m2llast-vs-m200.eps')
    begplot, name,/encap, /color


    xtitle='M!D200!N [h!U-1!N M'+sunsymbol()+']'
    ytitle = 'M/L!DLast!N [h M'+sunsymbol()+'/L'+sunsymbol()+']'

    yrange = [0,900]
    pplot, m2lstruct.m200, m2lstruct.m2lmax, $
        xerr=m2lstruct.m200_err, yerr=m2lstruct.m2lmax_err, psym=8, $
        aspect=1, /xlog, $
        xtitle=xtitle, ytitle=ytitle, yrange=yrange, ystyle=1

    wmom, m2lstruct.m2lmax, m2lstruct.m2lmax_err, wlm, ws, wlerr

    mlstr = ntostr(wlm,f='(F5.1)')+!csym.plusminus+ntostr(wlerr, f='(F4.1)')

    oplot, [1,1.e20], [wlm,wlm]+wlerr, line=2
    oplot, [1,1.e20], [wlm,wlm]
    oplot, [1,1.e20], [wlm,wlm]-wlerr, line=2

    legend, 'M/L Last Point: '+mlstr

    endplot,/trim

    print,'Mean M/L mean = '+ntostr(wlm)+!plusminus+ntostr(wlerr)



    ;
    ; Asymptotic M/L as a function of M200
    ;
    name = path_join(dir, front+'-m2lasym-vs-m200.eps')
    begplot, name,/encap, /color


    xtitle='M!D200!N [h!U-1!N M'+sunsymbol()+']'
    ytitle = 'M/L!DAsym!N [h M'+sunsymbol()+'/L'+sunsymbol()+']'

    wbad=where(m2lstruct.m2lasym gt 700,nbad,comp=comp)

    yerr = m2lstruct.m2lasym_err
    m2lasym = m2lstruct.m2lasym

    if nbad ne 0 then begin
        ;; get average conversion of regular error to mcmc error
        errfac = mean( yerr[comp]/m2lstruct[comp].m2lasym_err_0 )
        m2lasym[wbad] = m2lstruct[wbad].m2lasym_0
        yerr[wbad] = m2lstruct[wbad].m2lasym_err_0*errfac
    endif

    ;yerr = m2lmax_err
    ayrange = [100, 700]
    pplot, m2lstruct.m200, m2lasym, xerr=m2lstruct.m200_err, yerr=yerr, psym=8, $
        aspect=1, /xlog, $
        xtitle=xtitle, ytitle=ytitle, yrange=ayrange, ystyle=3
    wmom, m2lasym, yerr, wam, ws, waerr

    mastr = ntostr(wam,f='(F5.1)')+!csym.plusminus+ntostr(waerr, f='(F3.1)')

    oplot, [1,1.e20], [wam,wam]+waerr, line=2
    oplot, [1,1.e20], [wam,wam]
    oplot, [1,1.e20], [wam,wam]-waerr, line=2

    legend, 'M/L Asymptotic: '+mastr

    endplot,/trim

    print,'Mean M/L asymptotic = '+ntostr(wam)+!plusminus+ntostr(waerr)



    ;
    ; Both as a function of M200
    ;
    name = path_join(dir, front+'-m2lasym-last-vs-m200.eps')
    begplot, name,/encap, /color


    xtitle='M!D200!N [h!U-1!N M'+sunsymbol()+']'
    ytitle = 'M/L [h M'+sunsymbol()+'/L'+sunsymbol()+']'

    ;yerr = m2lmax_err
    pplot, m2lstruct.m200, m2lstruct.m2lmax, $
        xerr=m2lstruct.m200_err, yerr=m2lstruct.m2lmax_err, psym=8, $
        aspect=1, /xlog, $
        xtitle=xtitle, ytitle=ytitle, yrange=yrange, ystyle=1

    pplot, m2lstruct.m200*1.1, m2lasym, xerr=m2lstruct.m200_err, yerr=yerr, $
        psym=7,/overplot, color=!blue

    oplot, [1,1.e20], [wam,wam]+waerr, line=2, color=!blue
    oplot, [1,1.e20], [wam,wam], color=!blue
    oplot, [1,1.e20], [wam,wam]-waerr, line=2, color=!blue

    oplot, [1,1.e20], [wlm,wlm]+wlerr, line=2
    oplot, [1,1.e20], [wlm,wlm]
    oplot, [1,1.e20], [wlm,wlm]-wlerr, line=2

    amess = 'Asymptotic: '+ mastr
    lmess = 'Last Point: '+ mlstr
    legend, [amess,lmess],psym=[8,7],color=[!p.color,!blue]

    endplot,/trim



end 
