;;Plotting M/L for 16-bin lum split.


function _m2lstruct, nrad, nbin
    arrval=dblarr(nrad)
    arrval2=dblarr(nrad,nrad)
    st = {$
            r:arrval, $
            $
            mass: arrval, $
            mass_err: arrval, $
            mass_yfit: arrval, $
            nfw_yfit: arrval,$
            linear_yfit: arrval, $
            $
            lum: arrval, $
            lum_err: arrval, $
            lumtot: arrval, $
            lumtot_err: arrval, $
            $
            m2l: arrval, $
            m2l_err: arrval, $
            m2l_cov: arrval2, $
            m2ltot: arrval, $
            m2ltot_err: arrval, $
            m2ltot_cov: arrval2, $
            $
            r200: 0d, $
            n200_red: 0L, $ 
            l200_red: 0d, $
            l200: 0d, $
            l200_err: 0d, $
            m200: 0d, $
            m200_err: 0d, $
            m2l200: 0d, $
            m2l200_err: 0d, $
            m2lmax: 0d, $
            m2lmax_err: 0d, $
            $
            m2lasym: 0d, $
            m2lasym_err: 0d, $
            rs: 0d, $
            rs_err: 0d, $
            index: 0d, $
            index_err: 0d, $
            yfit: arrval, $
            $
            m2lasym_0: 0d, $
            m2lasym_err_0: 0d, $
            rs_0: 0d, $
            rs_err_0: 0d, $
            index_0: 0d, $
            index_err_0: 0d, $
            yfit_0: arrval $
        }
        
    st = replicate(st,  nbin)
    return, st

end

pro _compare200, m
    ; for calculations
    m21=obj_new('maxbcg_lensing', 21)
    nm = n_elements(m)
    print,' r200int r200nfw m200int m200nfw'
    for i=0L, nm-1 do begin
        mv = m21->massvals(m[i])
        st = m21->nfwfit(mv.r, mv.mass, mv.mass_err, 0.25)
        sh = m21->m200(m[i])
        print,sh.r200,st.r200,sh.m200,st.m200
    endfor

end

pro _calcm2l, m, m_cov, lin, ltotin, lin_cov, m2l, m2l_cov, m2ltot, m2ltot_cov
    
    ; First radius bin of luminosity is not in the mass. Get the ones we want
    num=n_elements(m)

    l = dblarr(num)
    ltot = l
    l_cov = dblarr(num,num)
    for i=1,num do begin
        l[i-1] = lin[i]
        ltot[i-1] = ltotin[i]
        for j=1,num do begin
            l_cov[i-1,j-1] = lin_cov[i,j]
        endfor
    endfor
        
    m2l = m/l
    m2ltot = m/ltot

    m2l_cov = dblarr(num, num)
    m2ltot_cov = m2l_cov
    for i=0L, num-1 do begin
        for j=0L, num-1 do begin
            m2l_cov[i,j] = m2l[i]*m2l[j]*( m_cov[i,j]/m[i]/m[j] + $
                                           l_cov[i,j]/l[i]/l[j] )
            m2ltot_cov[i,j] = m2ltot[i]*m2ltot[j]*( m_cov[i,j]/m[i]/m[j] + $
                                                    l_cov[i,j]/ltot[i]/ltot[j] )
        endfor
    endfor

end

function _fit_m200_and_m2l, mass_sample, lsample, subtype

    ; correlate and lensing objects
    cobj = obj_new('maxbcg_correlate', lsample)
    mobj = obj_new('maxbcg_lensing', mass_sample[0]) ; vers doesn't matter here

    ; read outputs for correlate and mass
    c = cobj->corr_read('dd','jackknife_invert', sub=subtype)
    cm = cobj->corr_read('dd','means', sub=subtype)

    m = mobj->lensread('jackknife_invert', sub=subtype, sample=mass_sample)

    n=n_elements(c)
    nr = n_elements(c[0].ir)

    w = ( lindgen(nr) )[1:nr-1]

    mark_color = !red
    ocolor = !darkgreen

    m2lytitle = 'M/L [h M'+sunsymbol()+'/L'+sunsymbol()+']'
    xtitle = !mpcxtitle3d
    if keyword_set(allfour) then !p.charsize=1
    for i=0L, n-1 do begin 

        print
        print,'bin = '+ntostr(i+1)+'/'+ntostr(n)


        radius = c[i].ir[1:nr-1]
        mass = m[i].mass_out*1.e12
        mcov = m[i].mass_out_cov*1.e12*1.e12
        mass_err = m[i].mass_out_err*1.e12
        lum = c[i].iradilum*1.e10
        lcov = c[i].iradilum_cov*1.e10*1.e10
        lum_err = c[i].iradilum_err*1.e10
        lumtot = lum + c[i].bcg_ilum_mean*1.e10
        lumtot_err = sqrt( lum_err^2 + (c[i].bcg_ilum_err*1.e10)^2 )

        mv = mobj->massvals(m[i])
        st = mobj->nfwfit(mv.r, mv.mass, mv.mass_err, 0.25)
        tr200 = st.r200
        tm200 = st.m200
        tm200_err = st.m200_err


        tl200 = interpol(lumtot, c[i].ir, tr200)
        tl200_err = interpol(lumtot_err, c[i].ir, tr200)

        tit = c[i].label+' '+!csym.times+' 10!U10!N L'+sunsymbol()


        _calcm2l, mass, mcov, lum, lumtot, lcov, m2l, m2l_cov, m2ltot, m2ltot_cov
        m2ltot_err = cov2diag(m2ltot_cov)
        m2l_err = cov2diag(m2l_cov)
        ;endplot
        ntrial = 1000000
        ;ntrial = 100000
        burnin = long( ntrial/10.0 )

        ;m2ltot_cov = diagonal_array(m2ltot_err^2)
        fitst = fitm2l_mcmc(radius, m2ltot, m2ltot_cov, ntrial, burnin, /doplot)
        ;stop
 
        nrat=n_elements(m2ltot)

        xrange= [0.015,30]

        if tr200 gt 0 then begin 

            if keyword_set(allfour) then begin

                !p.multi = [0,2,2]
                pplot, m[i].ir, mass, yerr=mass_err, $
                    psym=8, /xlog, /ylog, $
                    xtitle=xtitle, ytitle = 'M(<r) [h!U-1!NM'+sunsymbol()+']', $
                    title=tit, $
                    aspect=1, symsize=0.7, $
                    xrange=xrange, xstyle=3, yrange=[1.e11,5.e15], ystyle=3

          
                oplot, [0.001, 100], [tm200, tm200], color=mark_color
                oplot, [tr200, tr200], [0.01, 1.e19], color=mark_color

                ;; Cumulative lum, including BCG
                pplot, c[i].ir, lumtot, yerr=lum_err, $
                    psym=8, /xlog, /ylog, /nohat, $
                    xtitle=xtitle, ytitle = 'L(<r) [h!U-2!NL'+sunsymbol()+']', $
                    aspect=1, symsize=0.7, $
                    xrange=xrange, xstyle=3, yrange=[1.e8,1.e13], ystyle=3            

                ;; Cumulative lum, not including BCG
                pplot, c[i].ir, lum, yerr=lum_err, /nohat, color=ocolor, /overplot, $
                    xrange=xrange, xstyle=3, yrange=[5,2000], ystyle=3
                legend,['L(<r)','L(<r) + BCG'], psym=[1,8], color=[ocolor,!p.color],$
                    /right,/bottom,box=0,charsize=1


                oplot, [0.001, 100], [tl200, tl200], color=mark_color
                oplot, [tr200, tr200], [0.01, 1.e19], color=mark_color



                ;; Both together
                m2ltot_err = cov2diag(m2ltot_cov)
                pplot, radius, m2ltot, yerr=m2ltot_err, psym=8, /xlog, $
                    xtitle=xtitle, ytitle=m2lytitle, /ylog, $
                    aspect=1, symsize=0.7, $
                    xrange=xrange, xstyle=3, yrange=[5,2000], ystyle=3
                pplot, c[i].ir[w], ratio, yerr=ratio_err, color=ocolor, /overplot
                legend,['L(<r)','L(<r) + BCG'], psym=[1,8], color=[ocolor,!p.color],$
                    /right,/bottom,box=0,charsize=1

            endif

            add_arrval, fitst, fitstruct

            yrange = prange(m2ltot,m2ltot_err)
            yrange=[5,2000]
            yrange[0] = yrange[0] > 1
            pplot, radius, m2ltot, yerr=m2ltot_err, psym=8, /xlog, $
                yrange=yrange, ystyle=3,/ylog, xtitle=xtitle, ytitle=m2lytitle, $
                title=tit, $
                aspect=1, symsize=0.7
            pplot, radius, m2l, /overplot, color=!blue
            pplot, radius, fitst.yfit,/over
            pplot, radius, fitst.yfit_0,color=!DarkGreen, line=2, /over

            ww = where(m2ltot gt 0)
            tm2l200 = interpol(m2ltot[ww], radius[ww], tr200)
            tm2l200_err = interpol(m2ltot_err[ww], radius[ww], tr200)

            ;print,tm2l200,!plusminus,tm2l200_err

            if tm2l200 lt 0 then stop

            oplot, [tr200,tr200],[0.001,10000], color=mark_color
            oplot, [0.001,10000], [tm2l200,tm2l200], color=mark_color

            legend,['M/L!Dtot!N','M/L','Fit'], line=0, psym=[8,4,0], $
                color=[!p.color, !blue, !p.color], $
                /bottom,/right,box=0,/clear

            if i eq 0 then begin
                m2lstruct = _m2lstruct(n_elements(radius), n_elements(c))
            endif

            m2lstruct[i].r = radius

            m2lstruct[i].mass = mv.mass*1.e12
            m2lstruct[i].mass_err = mv.mass_err*1.e12
            m2lstruct[i].mass_yfit = st.mass_yfit
            m2lstruct[i].nfw_yfit = st.nfw_yfit
            m2lstruct[i].linear_yfit = st.linear_yfit

            m2lstruct[i].lum = lum[1:nr-1]
            m2lstruct[i].lum_err = lum_err[1:nr-1]
            m2lstruct[i].lumtot = lumtot[1:nr-1]
            m2lstruct[i].lumtot_err = lumtot_err[1:nr-1]

            m2lstruct[i].m2l = m2l
            m2lstruct[i].m2l_err = m2l_err
            m2lstruct[i].m2ltot = m2ltot 
            m2lstruct[i].m2ltot_err = m2ltot_err

            m2lstruct[i].r200 = tr200
            m2lstruct[i].n200_red = round(m[i].mean_ngals200)
            m2lstruct[i].l200_red = m[i].mean_ilum200
            m2lstruct[i].l200 = tl200
            m2lstruct[i].l200_err = tl200_err
            m2lstruct[i].m200 = tm200
            m2lstruct[i].m200_err = tm200_err
            m2lstruct[i].m2l200 = tm2l200
            m2lstruct[i].m2l200_err = tm2l200_err
            m2lstruct[i].m2lmax = m2l[nr-2]
            m2lstruct[i].m2lmax_err = m2l_err[nr-2]

            m2lstruct[i].m2lasym = fitst.m2l
            m2lstruct[i].m2lasym_err = fitst.m2l_err
            m2lstruct[i].rs = fitst.rs
            m2lstruct[i].rs_err = fitst.rs_err
            m2lstruct[i].index = fitst.index
            m2lstruct[i].index_err = fitst.index_err
            m2lstruct[i].yfit = fitst.yfit

            m2lstruct[i].m2lasym_0 = fitst.m2l_0
            m2lstruct[i].m2lasym_err_0 = fitst.m2l_err_0
            m2lstruct[i].rs_0 = fitst.rs_0
            m2lstruct[i].rs_err_0 = fitst.rs_err_0
            m2lstruct[i].index_0 = fitst.index_0
            m2lstruct[i].index_err_0 = fitst.index_err_0
            m2lstruct[i].yfit_0 = fitst.yfit_0

        endif 
      
    endfor 


    obj_destroy, cobj, mobj

    return, m2lstruct


end

pro maxbcg_lensing_calc_m2l, subtype, allfour=allfour

    if n_elements(subtype) eq 0 then begin
        print,'maxbcg_lensing_calc_m2l, subtype, /allfour'
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

    ; For each bin, fit for the NFW m200 (no mis-centering taken into
    ; account), calculate the M/L, and fit the basic M/L function
    ; using fitm2l_mcmc

    name = path_join(dir, front+'-each.ps')
    begplot, name=name, /color
    m2lstruct = _fit_m200_and_m2l(mass_sample, lsample, subtype)
    endplot

    output_name = path_join(dir, front+'-m2l.fits')
    print,'Writing fits file: ',output_name
    mwrfits, m2lstruct, output_name, /create

end 
