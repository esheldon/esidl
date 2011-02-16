function maxbcg_m2l::init, msample, lsample

    self.mass_sample = ptr_new(msample)
    self.lum_sample = ptr_new(lsample)
    return, 1
end

;
; samples and files
;

function maxbcg_m2l::mass_sample
    return, *self.mass_sample
end
function maxbcg_m2l::lum_sample
    return, *self.lum_sample
end

function maxbcg_m2l::mstr
    msample = self->mass_sample()
    return,strjoin( ntostr(msample), '-')
end
function maxbcg_m2l::lstr
    lsample = self->lum_sample()
    return,strjoin( ntostr(lsample), '-')
end

function maxbcg_m2l::dir, subtype, createdir=createdir
    dir = '~/plots/m2l/maxbcg/'+subtype
    if keyword_set(createdir) then begin
        if not fexist(dir) then file_mkdir, dir
    endif
    return, dir
end

function maxbcg_m2l::file_front, subtype
    mstr = self->mstr()
    lstr = self->lstr()
    return, 'm2l-'+subtype+'-m'+mstr+'-l'+lstr
end
function maxbcg_m2l::m2lfile, subtype
    dir = self->dir(subtype)
    front = self->file_front(subtype)
    return, path_join(dir, front+'-m2l.fits')
end

function maxbcg_m2l::plotfile, subtype, post, createdir=createdir
    dir = self->dir(subtype, createdir=createdir)
    front = self->file_front(subtype)
    return, path_join(dir, front+'-'+post)
end

function maxbcg_m2l::m2lstring
    ;return, !csym.delta_cap+'M/'+!csym.delta_cap+'L'
    return, textoidl('\DeltaM/\DeltaL')
end

function maxbcg_m2l::subtype_legend, subtype
    case subtype of
        'ngals200_12': return,'N200' 
        'ilum200_16': return,'L200' 
        else: message,'no subtype match for '+subtype
    endcase
end

;
; Methods to calculate the M/L for given sample and subtype
;

;
; The M/L data structure.  Contains much of what we need for
; the plots and tables.
;
function maxbcg_m2l::m2lstruct, nrad, nbin
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
            m2l_corr: arrval2, $
            m2ltot: arrval, $
            m2ltot_err: arrval, $
            m2ltot_cov: arrval2, $
            m2ltot_corr: arrval2, $
			$
            m2n: arrval, $
            m2n_err: arrval, $
            m2n_cov: arrval2, $
            m2n_corr: arrval2, $
            m2ntot: arrval, $
            m2ntot_err: arrval, $
            m2ntot_cov: arrval2, $
            m2ntot_corr: arrval2, $
            $
            meanlum: arrval, $
            meanlum_err: arrval, $
            $
            r200: 0d, $
            lbcg_red: 0d, $
            n200_red: 0d, $ 
            l200_red: 0d, $
            l200: 0d, $
            l200_err: 0d, $
            m200: 0d, $
            m200_err: 0d, $
            m2l200: 0d, $
            m2l200_err: 0d, $
            m2lmax: 0d, $
            m2lmax_err: 0d, $
            m2n200: 0d, $
            m2n200_err: 0d, $
            m2nmax: 0d, $
            m2nmax_err: 0d, $
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

; Here we assume the bin centers are the same and all arrays and matrices
; are the same shape
; This can also be used to calculate mass to number ratios
pro maxbcg_m2l::_calcm2l, $
        m, m_cov, l, ltot, l_cov, m2l, m2l_cov, m2ltot, m2ltot_cov
    
    m2l = m/l
    m2ltot = m/ltot

    num=n_elements(m)
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


function maxbcg_m2l::sub_matrix, mat, ind

	n=n_elements(ind)
	newmat = replicate(mat[0], n, n)
	for i=0L, n-1 do begin
		for j=0L, n-1 do begin
			newmat[i,j] = mat[ind[i], ind[j]]
		endfor
	endfor
	return, newmat
end

pro maxbcg_m2l::getboth, subtype, msample, lsample, m, l

	mobj = obj_new('maxbcg_lensing', msample[0])
	lobj = obj_new('maxbcg_correlate', lsample)

    l = lobj->corr_read('dd','jackknife_invert', sub=subtype)
    m = mobj->lensread('jackknife_invert', sub=subtype, sample=msample)
	obj_destroy, mobj, lobj

end



; You must enter the radial indices for mass and lum so that the radii
; essentially match.  No offset between the two is valid
function maxbcg_m2l::fit_m200_and_m2l, subtype, mradind, lradind

	nmradind = n_elements(mradind)
	nlradind = n_elements(lradind)
	if n_elements(subtype) eq 0 or nmradind eq 0 or nlradind eq 0 then begin
		print,'usage: st=m2l->fit_m200_and_m2l(subtype, mradind, lradind)'
		message,'Halting'
	endif
	if nmradind ne nlradind then begin
		message,'mradind and lradind must be equal length'
	endif

    mass_sample = self->mass_sample()
    lsample = self->lum_sample()
    ; correlate and lensing objects
    cobj = obj_new('maxbcg_correlate', lsample)
    mobj = obj_new('maxbcg_lensing', mass_sample[0]) ; vers doesn't matter here

    ; read outputs for correlate and mass
    c = cobj->corr_read('dd','jackknife_invert', sub=subtype)
    ;cm = cobj->corr_read('dd','means', sub=subtype)

    m = mobj->lensread('jackknife_invert', sub=subtype, sample=mass_sample)

    n=n_elements(c)

    mark_color = c2i('red')
    ocolor = c2i('darkgreen')

    m2lytitle = self->m2lstring()+'(< r) [h M'+sunsymbol()+'/L'+sunsymbol()+']'
    xtitle = textoidl('r [h^{-1}Mpc]')
    if keyword_set(allfour) then !p.charsize=1
    for i=0L, n-1 do begin 

        print
        print,'bin = '+ntostr(i+1)+'/'+ntostr(n)

        radius = c[i].ir[lradind]

		; mass
        mass = m[i].mass_out[mradind]*1.e12
        mcov = self->sub_matrix(m[i].mass_out_cov, mradind)*1.e12*1.e12
        mass_err = m[i].mass_out_err[mradind]*1.e12

		; luminosity
        lum = c[i].iradilum[lradind]*1.e10
        lcov = self->sub_matrix(c[i].iradilum_cov, lradind)*1.e10*1.e10
        lum_err = c[i].iradilum_err[lradind]*1.e10
        lumtot = lum + c[i].bcg_ilum_mean*1.e10
        lumtot_err = sqrt( lum_err^2 + (c[i].bcg_ilum_err*1.e10)^2 )


		; number
		; see the ::invert function which is adding 1 for the BCG; this 
		; is not consistent and should be changed
		num = c[i].iradcounts[lradind]
		ncov = self->sub_matrix(c[i].iradcounts_cov, lradind)
		num_err = c[i].iradcounts_err[lradind]
		numtot = 1+num
		numtot_err = num_err



        mv = mobj->massvals(m[i])
        st = mobj->nfwfit($
			mv.r[mradind], mv.mass[mradind], mv.mass_err[mradind], 0.25)
        tr200 = st.r200
        tm200 = st.m200
        tm200_err = st.m200_err


        tl200 = interpol(lumtot, c[i].ir, tr200)
        tl200_err = interpol(lumtot_err, c[i].ir, tr200)

        tit = c[i].label+' '+!csym.times+' 10!U10!N L'+sunsymbol()


        self->_calcm2l, mass, mcov, lum, lumtot, lcov, m2l, m2l_cov, m2ltot, m2ltot_cov
        m2ltot_err = cov2diag(m2ltot_cov)
        m2l_err = cov2diag(m2l_cov)

        self->_calcm2l, mass, mcov, num, numtot, ncov, m2n, m2n_cov, m2ntot, m2ntot_cov
        m2ntot_err = cov2diag(m2ntot_cov)
        m2n_err = cov2diag(m2n_cov)


        ntrial = 1000000
        burnin = long( ntrial/10.0 )

        ;m2ltot_cov = diagonal_array(m2ltot_err^2)
        fitst = fitm2l_mcmc(radius, m2ltot, m2ltot_cov, ntrial, burnin, /doplot)
        ;stop
 
        nrat=n_elements(m2ltot)

        xrange= [0.015,30]

        if tr200 gt 0 then begin 

            if keyword_set(allfour) then begin

                !p.multi = [0,2,2]
                pplot, radius, mass, yerr=mass_err, $
                    psym=8, /xlog, /ylog, $
                    xtitle=xtitle, ytitle = !csym.delta_cap+'M(<r) [h!U-1!NM'+sunsymbol()+']', $
                    title=tit, $
                    aspect=1, symsize=0.7, $
                    xrange=xrange, xstyle=3, yrange=[1.e11,5.e15], ystyle=3

          
                oplot, [0.001, 100], [tm200, tm200], color=mark_color
                oplot, [tr200, tr200], [0.01, 1.e19], color=mark_color

                ;; Cumulative lum, including BCG
                pplot, radius, lumtot, yerr=lum_err, $
                    psym=8, /xlog, /ylog, /nohat, $
                    xtitle=xtitle, ytitle = !csym.delta_cap+'L(<r) [h!U-2!NL'+sunsymbol()+']', $
                    aspect=1, symsize=0.7, $
                    xrange=xrange, xstyle=3, yrange=[1.e8,1.e13], ystyle=3            

                ;; Cumulative lum, not including BCG
                pplot, radius, lum, yerr=lum_err, /nohat, color=ocolor, /overplot, $
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
                pplot, radius, ratio, yerr=ratio_err, color=ocolor, /overplot
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
            pplot, radius, m2l, /overplot, color=c2i('blue')
            pplot, radius, fitst.yfit,/over
            pplot, radius, fitst.yfit_0,color=c2i('darkgreen'), line=2, /over

            ww = where(m2ltot gt 0)
            tm2l200 = interpol(m2ltot[ww], radius[ww], tr200)
            tm2l200_err = interpol(m2ltot_err[ww], radius[ww], tr200)
            if tm2l200 lt 0 then stop

            tm2n200 = interpol(m2ntot[ww], radius[ww], tr200)
            tm2n200_err = interpol(m2ntot_err[ww], radius[ww], tr200)
            if tm2n200 lt 0 then stop

            oplot, [tr200,tr200],[0.001,10000], color=mark_color
            oplot, [0.001,10000], [tm2l200,tm2l200], color=mark_color

            legend,['M/L!Dtot!N','M/L','Fit'], line=0, psym=[8,4,0], $
                color=[!p.color, c2i('blue'), !p.color], $
                /bottom,/right,box=0,/clear

            if i eq 0 then begin
                m2lstruct = self->m2lstruct(n_elements(radius), n_elements(c))
            endif

            cumstruct = cobj->cumulative_densities(c[i])

            mlum=cumstruct.radilumdens[lradind]/cumstruct.radnumdens[lradind]
            mlum_err=mlum*sqrt( $
                (cumstruct.radilumdens_err[lradind]/cumstruct.radilumdens[lradind])^2 + $
                (cumstruct.radnumdens_err[lradind]/cumstruct.radnumdens[lradind])^2 )
            m2lstruct[i].meanlum=mlum
            m2lstruct[i].meanlum_err=mlum_err

            m2lstruct[i].r = radius
			nr = n_elements(radius)

            m2lstruct[i].mass = mass
            m2lstruct[i].mass_err = mass_err
            m2lstruct[i].mass_yfit = st.mass_yfit
            m2lstruct[i].nfw_yfit = st.nfw_yfit
            m2lstruct[i].linear_yfit = st.linear_yfit

            m2lstruct[i].lum = lum
            m2lstruct[i].lum_err = lum_err
            m2lstruct[i].lumtot = lumtot
            m2lstruct[i].lumtot_err = lumtot_err

            m2lstruct[i].m2l = m2l
            m2lstruct[i].m2l_err = m2l_err
            m2lstruct[i].m2l_cov = m2l_cov
            m2lstruct[i].m2l_corr = cov2corr(m2l_cov)
            m2lstruct[i].m2ltot = m2ltot 
            m2lstruct[i].m2ltot_err = m2ltot_err
            m2lstruct[i].m2ltot_cov = m2ltot_cov
            m2lstruct[i].m2ltot_corr = cov2corr(m2ltot_cov)

            m2lstruct[i].m2n = m2n
            m2lstruct[i].m2n_err = m2n_err
            m2lstruct[i].m2n_cov = m2n_cov
            m2lstruct[i].m2n_corr = cov2corr(m2n_cov)
            m2lstruct[i].m2ntot = m2ntot 
            m2lstruct[i].m2ntot_err = m2ntot_err
            m2lstruct[i].m2ntot_cov = m2ntot_cov
            m2lstruct[i].m2ntot_corr = cov2corr(m2ntot_cov)


            m2lstruct[i].r200 = tr200
            m2lstruct[i].lbcg_red = m[i].mean_bcg_ilum
            m2lstruct[i].n200_red = m[i].mean_ngals200
            m2lstruct[i].l200_red = m[i].mean_ilum200
            m2lstruct[i].l200 = tl200
            m2lstruct[i].l200_err = tl200_err
            m2lstruct[i].m200 = tm200
            m2lstruct[i].m200_err = tm200_err

            m2lstruct[i].m2l200 = tm2l200
            m2lstruct[i].m2l200_err = tm2l200_err
            m2lstruct[i].m2lmax = m2l[nr-1]
            m2lstruct[i].m2lmax_err = m2l_err[nr-1]

            m2lstruct[i].m2n200 = tm2n200
            m2lstruct[i].m2n200_err = tm2n200_err
            m2lstruct[i].m2nmax = m2n[nr-1]
            m2lstruct[i].m2nmax_err = m2n_err[nr-1]

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

pro maxbcg_m2l::calcm2l, subtype, mradind, lradind

	nmradind = n_elements(mradind)
	nlradind = n_elements(lradind)
	if n_elements(subtype) eq 0 or nmradind eq 0 or nlradind eq 0 then begin

		print,'usage: m2l->calcm2l, subtype, mradind, lradind'
		message,'Halting'
	endif
	if nmradind ne nlradind then begin
		message,'mradind and lradind must be equal length'
	endif

    ; For each bin, fit for the NFW m200 (no mis-centering taken into
    ; account), calculate the M/L, and fit the basic M/L function
    ; using fitm2l_mcmc

    name = self->plotfile(subtype, 'each.ps')
    begplot, name=name, /color
    m2lstruct = self->fit_m200_and_m2l(subtype, mradind, lradind)
    endplot

    file = self->m2lfile(subtype)
    print,'Writing fits file: ',file
    mwrfits, m2lstruct, file, /create

end 


;
; Plotting stuff
;
function maxbcg_m2l::mpvalue, n, axis_struct=axis_struct

    case n of
        12: begin
            mp = [4,3]
        end
        16: begin
            mp = [4,4]
        end
        else: message,'bad n'
    endcase
    xtickshow = (lindgen(n)/mp[0]) ge (mp[1]-1)
    ytickshow = (lindgen(n) mod mp[0]) eq 0

    axis_struct = {                $
            xtickshow: xtickshow,  $
            ytickshow: ytickshow   $
        }
    return, mp
end

pro maxbcg_m2l::plot_masses, subtype, color=color

    m2lfile = self->m2lfile(subtype)
    m2lstruct = mrdfits(m2lfile, 1)

    m=obj_new('maxbcg','dr406')
    ws = m->where_string(subtype, labels=labels)
    obj_destroy, m

    post = 'massfits'
    if keyword_set(color) then post = post+'-color.eps' else post=post+'.eps'
    pfile = self->plotfile(subtype,post)
    begplot, pfile, /color, /encap

    n=n_elements(m2lstruct)
    mp = self->mpvalue(n, axis_struct=as)

    !p.charsize=1
    !p.thick=2
    !x.thick=2
    !y.thick=2
    !p.charthick=2
    mytitle=TeXToIDL('\DeltaM(<r) [h^{-1} M_{\odot}]')
    titsize=1.5
    erase & multiplot, mp, /square, $
        mxTitle=!mpcxtitle3d, mXtitSize=titsize, mXtitOffset=0.5, $
        myTitle=mytitle, myTitSize=titsize

    line_both = 0
    line_nfw = 3
    line_twohalo = 1
    if keyword_set(color) then begin
        clr_twohalo = c2i('blue')
        clr_both = c2i('DarkGreen')
        clr_nfw = c2i('orange')
        clr_r200 = c2i('red')
    endif else begin
        clr_twohalo = c2i('grey30')
        clr_both = !p.color
        clr_nfw = c2i('grey50')
        clr_r200 = c2i('black')
    endelse

    for i=0L, n-1 do begin

        if as.xtickshow[i] then xtickf='loglabels' else xtickf=''
        if as.ytickshow[i] then ytickf='loglabels' else ytickf=''
        pplot, m2lstruct[i].r, m2lstruct[i].mass, yerr=m2lstruct[i].mass_err, $
            psym=8, symsize=0.5, hat=0, /xlog, /ylog, $
            xrange=[0.011, 60.0], yrange=[0.1,40000]*1.e12, $
            xstyle=3, ystyle=3, $
            xticklen=0.04, yticklen=0.04, $
            xtickf=xtickf, ytickf=ytickf

        pplot, m2lstruct[i].r, m2lstruct[i].nfw_yfit, $
            color=clr_nfw, line=line_nfw, /over
        pplot, m2lstruct[i].r, m2lstruct[i].linear_yfit, $
            color=clr_twohalo, line=line_twohalo, /over
        pplot, m2lstruct[i].r, m2lstruct[i].mass_yfit, $
            color=clr_both, line=line_both, /over

        pplot, [m2lstruct[i].r200], [m2lstruct[i].m200], $
            psym=2,  color=clr_r200, /over

        legend, labels[i], charsize=0.7, margin=0
        multiplot
    endfor

    multiplot,/default
    endplot,/trim_bbox

end

; Array of plots
pro maxbcg_m2l::plot_luminosities, subtype

    m2lfile = self->m2lfile(subtype)
    m2lstruct = mrdfits(m2lfile, 1)

    pfile = self->plotfile(subtype,'lum.eps')
    begplot, pfile, /color, /encap

    n=n_elements(m2lstruct)
    mp = self->mpvalue(n, axis_struct=as)

    !p.charsize=1
    !p.thick=2
    !x.thick=2
    !y.thick=2
    erase & multiplot, mp, /square, $
        mxTitle=!mpcxtitle3d, $
        myTitle=!csym.delta_cap+'L!Di!N(< r) [ h!U-2!N M'+sunsymbol()+' ]'

    for i=0L, n-1 do begin

        if as.xtickshow[i] then xtickf='loglabels' else xtickf=''
        if as.ytickshow[i] then ytickf='loglabels' else ytickf=''
        pplot, m2lstruct[i].r, m2lstruct[i].lum, $
            yerr=m2lstruct[i].lum_err, $
            psym=8, symsize=0.5, hat=0, /xlog, /ylog, $
            xrange=[0.01, 60.0], yrange=[0.005,10000]*1.e10, $
            xstyle=3, ystyle=3, $
            xticklen=0.04, yticklen=0.04, $
            xtickf=xtickf, ytickf=ytickf

        pplot, m2lstruct[i].r, m2lstruct[i].lumtot, /overplot
        pplot, [m2lstruct[i].r200], [m2lstruct[i].l200], $
            psym=2,  color=!red, /over
        multiplot
    endfor

    multiplot,/default
    endplot,/trim_bbox

end


; All overplottede
pro maxbcg_m2l::plot_luminosities_over, subtype, tot=tot, color=color, both=both

    m2lfile = self->m2lfile(subtype)
    m2lstruct = mrdfits(m2lfile, 1)

    if keyword_set(color) then post = '-color.eps' else post='.eps'
    if keyword_set(tot) then begin
        pfile = self->plotfile(subtype,'lumtotover'+post)
    endif else if keyword_set(both) then begin
        pfile = self->plotfile(subtype,'lumbothover'+post)
    endif else begin
        pfile = self->plotfile(subtype,'lumover'+post)
    endelse

    begplot, pfile, color=color, /encap

    n=n_elements(m2lstruct)
    mp = self->mpvalue(n, axis_struct=as)

    !p.charsize=2
    !p.thick=3
    !x.thick=3
    !y.thick=3
    !p.charthick=3

    if keyword_set(color) then begin
        colors = make_rainbow(n)
    endif else begin
        loadct, 0
        colors = arrscl( findgen(n), 0, 200)
    endelse

    for i=0L, n-1 do begin

        if keyword_set(tot) then begin
            lum = m2lstruct[i].lumtot
            lumerr = m2lstruct[i].lumtot_err
        endif else begin
            lum = m2lstruct[i].lum
            lumerr = m2lstruct[i].lum_err
        endelse
        if i eq 0 then begin
            ytitle=textoidl('\DeltaL(<r) [h^{-2} M_{\odot}]')
            pplot, m2lstruct[i].r, lum, yerr=lumerr, $
                aspect=1, $
                hat=0, /xlog, /ylog, $
                xrange=[0.01, 60.0], yrange=[0.005,50000]*1.e10, $
                xstyle=3, ystyle=3, $
                xticklen=0.02, yticklen=0.02, $
                xtickf=xtickf, ytickf=ytickf, $
                xTitle=!mpcxtitle3d, $
                yTitle=yTitle
        endif else begin
            pplot, m2lstruct[i].r, lum, yerr=lumerr, color=colors[i], $
                hat=0, /overplot
        endelse

        if keyword_set(both) and not keyword_set(tot) then begin
            pplot, m2lstruct[i].r, m2lstruct[i].lumtot, $
                /overplot, color=colors[i], line=2
        endif
        if keyword_set(tot) or keyword_set(both) then begin
            pplot, [m2lstruct[i].r200], [m2lstruct[i].l200], $
                psym=2,  color=colors[i], /over
        endif
    endfor

    pplot, m2lstruct.r200, m2lstruct.l200, psym=-2, /overplot

    m=obj_new('maxbcg','dr406')
    ws = m->where_string(subtype, labels=labels)
    obj_destroy, m

    labels=reverse(labels)
    colors=reverse(colors)

    if n_elements(lcharsize) eq 0 then lcharsize=1
    legend, labels, charsize=lcharsize, lin=0, color=colors


    endplot,/trim_bbox

end

pro maxbcg_m2l::_tagnums, struct, tagname, tagnum, errtagnum
	if not tag_exist(struct, tagname, index=tagnum) then begin
		message,'tag "'+tagname+'" does not exist'
	endif
	if not tag_exist(struct, tagname+'_err', index=errtagnum) then begin
		message,'tag "'+tagname+'_err" does not exist'
	endif
end

; all overplotted
pro maxbcg_m2l::plot_m2l_over, subtype, color=color, number=number

	m2lfile = self->m2lfile(subtype)
	m2lstruct = mrdfits(m2lfile, 1)

	m=obj_new('maxbcg','dr406')
	ws = m->where_string(subtype, labels=labels)
	obj_destroy, m

	xrange = [0.5*min(m2lstruct.r), 1.5*max(m2lstruct.r)]
	if keyword_set(number) then begin
		ytittex = '\DeltaM/\DeltaN [h^{-1} M_{\odot}]'
		post='m2nover'
		self->_tagnums, m2lstruct, 'm2ntot', tag, etag
		yrange = prange( m2lstruct.(tag), m2lstruct.(etag) )
		yrange[0] = 5e10
		o200 = m2lstruct.m2n200
	endif else begin
		ytittex = '\DeltaM/\DeltaL(<r) [h M_{\odot}/L_{\odot}]'
		post = 'm2lover'
		self->_tagnums, m2lstruct, 'm2ltot', tag, etag
		yrange = prange( m2lstruct.(tag), m2lstruct.(etag) )
		yrange[0] = 5
		o200 = m2lstruct.m2l200
	endelse

	if keyword_set(color) then begin
		post = post+'-color.eps' 
		clr_nobcg = 'blue'
		clr_r200 = 'red'
	endif else begin
		post=post+'.eps'
		clr_nobcg = 'grey30'
		if !d.name eq 'PS' then clr_r200 = 'black' else clr_r200 = 'white'
	endelse
	name = self->plotfile(subtype, post)
	begplot, name, /encap, color=color

	n=n_elements(m2lstruct)
	mp = self->mpvalue(n, axis_struct=as)
	colors = make_rainbow(n)

	line_nobcg = 1

	xtitle=textoidl('r [h^{-1} Mpc]')
	ytitle=textoidl(ytittex)
	!p.charsize=!p.charsize*1.2

	xtickf='loglabels'
	ytickf='loglabels'
	;err = m2lstruct[0].(etag)
	pplot, m2lstruct[0].r, m2lstruct[0].(tag), $
		yerr=err, $
		hat=0, /xlog, /ylog, $
		xrange=xrange, yrange=yrange, $
		xstyle=3, ystyle=3, $
		xticklen=0.04, yticklen=0.04, $
		xtickf=xtickf, ytickf=ytickf, $
		xtitle=xtitle, ytitle=ytitle, $
		aspect=1
	for i=0L, n-1 do begin
		;err = m2lstruct[i].(etag)
		pplot, m2lstruct[i].r, m2lstruct[i].(tag), yerr=err,$
			color=colors[i], /overplot
		pplot, [m2lstruct[i].r200], [o200[i]], psym=7, color=clr_r200, /over
	endfor

	legend, reverse(labels), charsize=1.3, /right, /bottom, line=0, color=reverse(colors)

    endplot,/trim
end



; array of plots with M/L fits 
pro maxbcg_m2l::plot_m2lfits, subtype, color=color

    m2lfile = self->m2lfile(subtype)
    m2lstruct = mrdfits(m2lfile, 1)

    m=obj_new('maxbcg','dr406')
    ws = m->where_string(subtype, labels=labels)
    obj_destroy, m

    post = 'm2lfits'
    if keyword_set(color) then post = post+'-color.eps' else post=post+'.eps'
    name = self->plotfile(subtype, post)
    begplot, name, /encap, color=color
 
    n=n_elements(m2lstruct)
    mp = self->mpvalue(n, axis_struct=as)

    line_nobcg = 1
    if keyword_set(color) then begin
        clr_nobcg = c2i('blue')
        clr_r200 = c2i('red')
    endif else begin
        clr_nobcg = c2i('grey30')
        clr_r200 = !p.color
    endelse
    !p.charsize=1
    !p.thick=2
    !x.thick=2
    !y.thick=2
    !p.charthick=2
    mytitle=textoidl(self->m2lstring()+'(<r) [h M_{\odot}/L_{\odot}]')
    titsize=1.5
    erase & multiplot, mp, /square, $
        mxTitle=!mpcxtitle3d, mXtitSize=titsize, mXtitOffset=0.5, $
        myTitle=mytitle, mYtitSize=titsize

    for i=0L, n-1 do begin

        if as.xtickshow[i] then xtickf='loglabels' else xtickf=''
        if as.ytickshow[i] then ytickf='loglabels' else ytickf=''
        pplot, m2lstruct[i].r, m2lstruct[i].m2ltot, $
            yerr=m2lstruct[i].m2ltot_err, $
            psym=8, symsize=0.5, hat=0, /xlog, /ylog, $
            xrange=[0.011, 60.0], yrange=[7,2000], $
            xstyle=3, ystyle=3, $
            xticklen=0.04, yticklen=0.04, $
            xtickf=xtickf, ytickf=ytickf

        pplot, m2lstruct[i].r, m2lstruct[i].m2l, color=clr_nobcg, $
            line=line_nobcg, /over

        if (m2lstruct[i].m2lasym gt 700) or (subtype eq 'ngals200_12' and i eq (n-1)) then begin
            yfit = m2lstruct[i].yfit_0
        endif else begin
            yfit = m2lstruct[i].yfit
        endelse

        pplot, m2lstruct[i].r, yfit, /overplot
        pplot, [m2lstruct[i].r200], [m2lstruct[i].m2l200], $
            psym=2,  color=clr_r200, /over

        legend, labels[i], charsize=0.7, margin=0, /right, /bottom

        multiplot
    endfor

    multiplot,/default
 
    endplot,/trim
end

pro maxbcg_m2l::plot_m2l200_vs_m200, subtype, number=number

    ;
    ; M/L at r200 vs. M200
    ;

    m2lfile = self->m2lfile(subtype)
    m2lstruct = mrdfits(m2lfile, 1)

	if keyword_set(number) then begin
		post='m2n200-vs-m200.eps'
		datpost='m2n200-vs-m200.dat'
		self->_tagnums, m2lstruct, 'm2n200', tag, etag
		ytitle = '(\DeltaM/\DeltaN)_{200} [h^{-1} M_{\odot}]'
	endif else begin
		post='m2l200-vs-m200.eps'
		datpost='m2l200-vs-m200.dat'
		self->_tagnums, m2lstruct, 'm2l200', tag, etag
		ytitle = '(\DeltaM/\DeltaL)_{200} [h M_{\odot}/L_{\odot}]'
	endelse
    name = self->plotfile(subtype, post)
    datname = self->plotfile(subtype, datpost)
    begplot, name, /encap
    !p.multi=0
	!p.charsize=!p.charsize*1.2

    xtitle=textoidl('M_{200} [h^{-1} M_{\odot}]')
	ytitle = textoidl(ytitle)
    pplot, m2lstruct.m200, m2lstruct.(tag), $
        xerr=m2lstruct.m200_err, yerr=m2lstruct.(etag), psym=8, $
        aspect=1, /xlog, $
        xtitle=xtitle, ytitle=ytitle, $
		ystyle=2

    endplot,/trim

	print,'Writing data file: ',datname
	openw, lun, datname, /get
	printf, lun, '# Column1: M_{200} [ h^{-1} M_{\odot} ]'
	printf, lun, '# Column2: \sigma(M_{200}) [ h^{-1} M_{\odot} ]'
	printf, lun, '# Column3: (M/N)_{200} [ h^{-1} M_{\odot} ]'
	printf, lun, '# Column3: \sigma((M/N)_{200}) [ h^{-1} M_{\odot} ]'
	colprint, lun=lun, $
		m2lstruct.m200, m2lstruct.m200_err, m2lstruct.(tag), m2lstruct.(etag)

	free_lun, lun

end

function maxbcg_m2l::get_mvsn_mapping, m2lstruct
    ; get mass observable relation for axes
    x=m2lstruct.n200_red
    y=m2lstruct.m200/1.e12
    yerr=m2lstruct.m200_err/1.e12
    guess = [1.0, 1.0]
    fitpower, x, y, yerr, guess, yfit, a, aerr
    a[0]=a[0]*1.e12
    yfit=yfit*1.e12
    return, { $
        n200:m2lstruct.n200_red, $
        m200:m2lstruct.m200, $
        m200_err:m2lstruct.m200_err, $
        yfit:yfit, $
        norm:a[0], $
        index:a[1]}
end
pro maxbcg_m2l::plot_m2l200_vs_ngals200, subtype

    ;
    ; M/L at r200 vs. N200
    ;

    m2lfile = self->m2lfile(subtype)
    m2lstruct = mrdfits(m2lfile, 1)

    ; get mass observable relation for axes
    fitst = self->get_mvsn_mapping(m2lstruct)

    ;pplot, fitst.n200, fitst.m200, yerr=fitst.m200_err, psym=8, /xlog, /ylog
    ;pplot, fitst.n200, fitst.yfit, /overplot


    name = self->plotfile(subtype, 'm2l200-vs-n200.eps')
    textfile=repstr(name, '.eps', '.dat')
    begplot, name, /encap
    !p.multi=0

    xtitle_top = textoidl('\DeltaM_{200} [h^{-1} M_{\odot}]')
    xtitle_bottom = textoidl('N_{200}')
    ytitle = self->m2lstring()+textoidl('_{200} [h M_{\odot}/L_{\odot}]')

    pplot, m2lstruct.n200_red, m2lstruct.m2l200, $
        yerr=m2lstruct.m2l200_err, psym=8, $
        aspect=1, /xlog, xrange=[1,300], xstyle=3+8, $
        xtitle=xtitle_bottom, ytitle=ytitle, charsize=2

    ; add top axis
    rlin = 10.0^!x.crange
    toprange = fitst.norm*rlin^fitst.index
    axis, xaxis=1, xrange=toprange, xstyle=1, /xlog, $
        xtitle=xtitle_top+' !c', charsize=2

    if 0 then begin
        tk=self->read_tinker()
        xtk = interpol(fitst.n200, fitst.yfit, tk.m180)
        w=where(xtk gt 2.5)
        pplot, xtk[w]*1.25, tk[w].m180/tk[w].ltot_r/3.0, /overplot
    endif else begin
        tk=self->read_tinker_new()
        xtk = interpol(fitst.n200, fitst.yfit, tk.m200)
        xtk = xtk*1.5
        w=where(xtk gt min(fitst.n200) and xtk lt max(fitst.n200))
        pplot, xtk[w], tk[w].m2l_i, /overplot
    endelse

    print,'************************************'
    print,'Fitting M/L vs N200'
    fitpower, m2lstruct.n200_red, m2lstruct.m2l200, m2lstruct.m2l200_err, $
        [10.0, 1.0], yfit, a, asig
    ;pplot, m2lstruct.n200_red, yfit, line=2, /overplot

    print,'************************************'
    print,'Fitting M/L vs M200'
    fitpower, m2lstruct.m200/1.e12, m2lstruct.m2l200, m2lstruct.m2l200_err, $
        [10.0, 1.0], yfit, a, asig
 

    endplot,/trim

    print,'Writing text file for jeremy: ',textfile
    head='    n200_red            m200        m200_err         m2l_200     m2l_200_err'
    colprint, $
        header=head, $
        m2lstruct.n200_red, $
        m2lstruct.m200, m2lstruct.m200_err, $
        m2lstruct.m2l200, m2lstruct.m2l200_err, file=textfile

end

pro maxbcg_m2l::plot_m2l200_two, subtype1, subtype2

    ;
    ; M/L at r200 vs. M200
    ;

    m2lfile1 = self->m2lfile(subtype1)
    m2lstruct1 = mrdfits(m2lfile1, 1)
    m2lfile2 = self->m2lfile(subtype2)
    m2lstruct2 = mrdfits(m2lfile2, 1)

    name = self->plotfile(subtype1+'-'+subtype2, 'm2l200-vs-m200.eps',$
                            /createdir)
    begplot, name, /color, /encap
    !p.multi=0

    psym1 = 8
    psym2 = 4

    color1 = !p.color
    color2 = c2i('blue')


    xtitle='M!D200!N [h!U-1!N M'+sunsymbol()+']'
    ytitle = self->m2lstring()+'!D200!N [h M'+sunsymbol()+'/L'+sunsymbol()+']'
    pplot, m2lstruct1.m200, m2lstruct1.m2l200, $
        xerr=m2lstruct1.m200_err, yerr=m2lstruct1.m2l200_err, psym=psym1, $
        hat=0, $
        aspect=1, /xlog, $
        xtitle=xtitle, ytitle=ytitle
    pplot, m2lstruct2.m200, m2lstruct2.m2l200, $
        xerr=m2lstruct2.m200_err, yerr=m2lstruct2.m2l200_err, $
        hat=0, $
        psym=psym2, color=color2, $
        /overplot
    label1 = self->subtype_legend(subtype1)
    label2 = self->subtype_legend(subtype2)
    legend, [label1, label2]+' binning', psym=[psym1,psym2],$
        colors=[color1,color2]

    endplot,/trim

end

; Plot M/L last point vs n200
pro maxbcg_m2l::plot_m2llast_vs_ngals200, subtype

    m2lfile = self->m2lfile(subtype)
    m2lstruct = mrdfits(m2lfile, 1)

    ; get mass observable relation for axes
    fitst = self->get_mvsn_mapping(m2lstruct)

    ;pplot, fitst.n200, fitst.m200, yerr=fitst.m200_err, psym=8, /xlog, /ylog
    ;pplot, fitst.n200, fitst.yfit, /overplot


    name = self->plotfile(subtype, 'm2llast-vs-n200.eps')
    begplot, name, /encap

    xtitle_bottom = textoidl('N_{200}')
    xtitle_top = textoidl('\DeltaM_{200} [h^{-1} M_{\odot}]')
    ytitle = self->m2lstring()+textoidl('_{rmax} [h M_{\odot}/L_{\odot}]')


    yrange = [0,900]
    pplot, m2lstruct.n200_red, m2lstruct.m2lmax, $
        yerr=m2lstruct.m2lmax_err, psym=8, $
        aspect=1, /xlog, xrange=[1,300], xstyle=3+8, $
        xtitle=xtitle_bottom, ytitle=ytitle, yrange=yrange, ystyle=1, $
        charsize=2
    wmom, m2lstruct.m2lmax, m2lstruct.m2lmax_err, wlm, ws, wlerr

    mlstr = ntostr(wlm,f='(F5.1)')+!csym.plusminus+ntostr(wlerr, f='(F4.1)')

    help, wlm, wlerr
    oplot, [1,1000], [wlm,wlm]+wlerr, line=2
    oplot, [1,1000], [wlm,wlm]
    oplot, [1,1000], [wlm,wlm]-wlerr, line=2

    ; add top axis
    rlin = 10.0^!x.crange
    toprange = fitst.norm*rlin^fitst.index
    axis, xaxis=1, xrange=toprange, xstyle=1, /xlog, $
        xtitle=xtitle_top+' !c', charsize=2


    endplot,/trim

   
end

;
; Plot M/L last point, M/L asymptotic, and both overplotted
pro maxbcg_m2l::plot_m2lasym, subtype

    ;
    ; Asymptotic M/L as a function of M200
    ;

    label = (strsplit(subtype, '_',/extract))[0]
 
    m2lfile = self->m2lfile(subtype)
    m2lstruct = mrdfits(m2lfile, 1)

    name = self->plotfile(subtype, 'm2lasym-vs-m200.eps')
    begplot, name,/encap, /color


    xtitle='M!D200!N [h!U-1!N M'+sunsymbol()+']'
    ytitle = self->m2lstring()+'!DAsym!N [h M'+sunsymbol()+'/L'+sunsymbol()+']'

    wbad=where(m2lstruct.m2lasym gt 700,nbad,comp=comp)
    m2lasym_err = m2lstruct.m2lasym_err
    m2lasym = m2lstruct.m2lasym

    if nbad ne 0 then begin
        ;; get average conversion of regular error to mcmc error
        errfac = mean( m2lasym_err[comp]/m2lstruct[comp].m2lasym_err_0 )
        m2lasym[wbad] = m2lstruct[wbad].m2lasym_0
        m2lasym_err[wbad] = m2lstruct[wbad].m2lasym_err_0*errfac
    endif

    ayrange = [100, 700]
    pplot, m2lstruct.m200, m2lasym, xerr=m2lstruct.m200_err, yerr=m2lasym_err, psym=8, $
        aspect=1, /xlog, $
        xtitle=xtitle, ytitle=ytitle, yrange=ayrange, ystyle=3
    wmom, m2lasym, m2lasym_err, wam, ws, waerr

    mastr = ntostr(wam,f='(F5.1)')+!csym.plusminus+ntostr(waerr, f='(F3.1)')

    oplot, [1,1.e20], [wam,wam]+waerr, line=2
    oplot, [1,1.e20], [wam,wam]
    oplot, [1,1.e20], [wam,wam]-waerr, line=2

    legend, self->m2lstring()+' Asymptotic: '+mastr

    legend, label, /right, /clear
    endplot,/trim


    ;
    ; Last point
    ;

    name = self->plotfile(subtype, 'm2llast-vs-m200.eps')
    begplot, name,/encap, /color


    xtitle='M!D200!N [h!U-1!N M'+sunsymbol()+']'
    ytitle = self->m2lstring()+'!DLast!N [h M'+sunsymbol()+'/L'+sunsymbol()+']'

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

    legend, self->m2lstring()+' Last Point: '+mlstr

    legend, label, /right, /clear
    endplot,/trim


    name = self->plotfile(subtype, 'm2lasym-last-vs-m200.eps')
    begplot, name,/encap, /color
    blue = c2i('blue')

    ytitle = self->m2lstring()+' [h M'+sunsymbol()+'/L'+sunsymbol()+']'
    pplot, m2lstruct.m200, m2lstruct.m2lmax, $
        xerr=m2lstruct.m200_err, yerr=m2lstruct.m2lmax_err, psym=8, $
        aspect=1, /xlog, $
        xtitle=xtitle, ytitle=ytitle, yrange=yrange, ystyle=1

    pplot, m2lstruct.m200*1.1, m2lasym, xerr=m2lstruct.m200_err, $
        yerr=m2lasym_err, $
        psym=7,/overplot, color=blue

    wmom, m2lasym, m2lasym_err, wam, ws, waerr
    wmom, m2lstruct.m2lmax, m2lstruct.m2lmax_err, wlm, ws, wlerr
    mastr = ntostr(wam,f='(F5.1)')+!csym.plusminus+ntostr(waerr, f='(F3.1)')
    mlstr = ntostr(wlm,f='(F5.1)')+!csym.plusminus+ntostr(wlerr, f='(F4.1)')

    oplot, [1,1.e20], [wam,wam]+waerr, line=2, color=blue
    oplot, [1,1.e20], [wam,wam], color=blue
    oplot, [1,1.e20], [wam,wam]-waerr, line=2, color=blue

    oplot, [1,1.e20], [wlm,wlm]+wlerr, line=2
    oplot, [1,1.e20], [wlm,wlm]
    oplot, [1,1.e20], [wlm,wlm]-wlerr, line=2

    amess = 'Asymptotic: '+ mastr
    lmess = 'Last Point: '+ mlstr
    legend, [lmess,amess],psym=[8,7],color=[!p.color,blue]

    legend, label, /right, /clear
    endplot,/trim




end

pro maxbcg_m2l::plotall, subtype
    self->plot_masses, subtype

    self->plot_luminosities, subtype
    self->plot_luminosities_over, subtype
    self->plot_luminosities_over, subtype, /tot
    self->plot_luminosities_over, subtype, /both

    self->plot_m2lfits, subtype

    self->plot_m2l200, subtype
    self->plot_m2lasym, subtype

end


;
; The basic M/L tables
;

; convert the value of ngals200 or ilum200 to a string
function maxbcg_m2l::getsubstr, subtype, struct
    case subtype of
        'ngals200_12': begin
            tagname = 'mean_ngals200'
            format='(I0)'
        end
        'ilum200_16': begin
            tagname = 'mean_ilum200'
            format='(F5.1)'
        end
        else: message,'unknown type '+subtype
    endcase
    if not tag_exist(struct, tagname, ind=ind) then message,'tag not found'

    str = string(struct.(ind), format=format)
    return,str
end

; Create the name and label for the split varialbe, \nvir, \lvir, etc.
pro maxbcg_m2l::name_label, subtype, name, label, units
    case subtype of
        'ngals200_12': begin
            name = '\nvir'
            label = 'm2lngals'
            units = ''
        end
        'ilum200_16': begin
            name = '\lvir'
            label = 'm2llum'
            units = '($h^{-2} 10^{10} L_{\odot}$)'
        end
        else: message,'unknown type '+subtype
    endcase
end

function maxbcg_m2l::read_tinker_new
    ; columns now m200crit, ngal, M/L(r-band), M/L(i-band)
    dir = '/home/users/esheldon/oh/tinker/m2lpredict'
    file=path_join(dir, 'predict-2007-09-08.dat')
    sdef = {m200: 0.0, $
            junk: 0, $
            m2l_r: 0.0, $
            m2l_i: 0.0}
    read_struct, file, sdef, st, skip=1
    return, st

end

pro maxbcg_m2l::plot_tinker_new
    st=self->read_tinker_new()
    w=where(st.m200 gt 5.e11 and st.m200 lt 1.e15)
    pplot, st[w].m200, st[w].m2l_i, yrange=[0,500], /xlog, aspect=1, $
        xtitle=textoidl('M_{200}'), ytitle=textoidl('(M/L)_{200}')
end

function maxbcg_m2l::read_tinker
    dir = '/home/users/esheldon/oh/tinker/m2lpredict'
    file=path_join(dir, 'Q5.dat')
    sdef = {m180: 0.0, $
            ltot_r_lstar: 0.0, $
            fcen: 0.0, $
            ltot_r_lstar_err: 0.0, $
            fcen_err: 0.0}
    read_struct, file, sdef, st, skip=1

    mstar_r = -20.44
    ;lstar_r = calc_sdss_lumsolar(mstar_r, clr=2)
    lstar_r = sdss_am2lumsolar(mstar_r, clr=2)
    ;lstar_r = lstar_r[0]*1.e10

    sdef = {m180: 0.0, $
            ltot_r_lstar: 0.0, $
            ltot_r_lstar_err: 0.0, $
            ltot_r: 0.0, $
            ltot_r_err: 0.0, $
            fcen: 0.0, $
            fcen_err: 0.0}
    struct = replicate(sdef, n_elements(st))

    copy_struct, st, struct
    struct.ltot_r = struct.ltot_r_lstar*lstar_r
    struct.ltot_r_err = struct.ltot_r_lstar_err*lstar_r
    return, struct
end
pro maxbcg_m2l::plot_tinker
    st=self->read_tinker()
    plot, st.m180, st.m180/st.ltot_r, yrange=[0,2000], /xlog
end

; create the latex table
pro maxbcg_m2l::create_table, subtype, mlumrad=mlumrad, m2lstruct=m2lstruct

    if n_elements(subtype) eq 0 then begin
        print,'m2l->create_table, subtype, mlumrad=, m2lstruct='
        return
    endif

    m2lfile = self->m2lfile(subtype)
    m2lstruct = mrdfits(m2lfile, 1)

    if n_elements(mlumrad) eq 0 then mlumrad=10

    ; Find nearest to the mlumrad
    raddiff = abs(m2lstruct[0].r - mlumrad)
    md=min(raddiff, wrad)

    mass_sample = self->mass_sample()
    mobj = obj_new('maxbcg_lensing', mass_sample[0]) ; vers doesn't matter here
    
    m=mobj->lensread('jackknife_invert', sub=subtype, samp=mass_sample)
    n=n_elements(m2lstruct)

    ; asymptotic had some instability
    wbad=where(m2lstruct.m2lasym gt 700,nbad,comp=comp)

    m2la = m2lstruct.m2lasym
    m2la_err = m2lstruct.m2lasym_err

    rs = m2lstruct.rs
    rs_err = m2lstruct.rs_err
    index=m2lstruct.index
    index_err = m2lstruct.index_err

    if nbad ne 0 then begin
        ;; get average conversion of regular error to mcmc error
        errfac = mean( m2la_err[comp]/m2lstruct[comp].m2lasym_err_0 )
        m2la[wbad] = m2lstruct[wbad].m2lasym_0
        m2la_err[wbad] = m2lstruct[wbad].m2lasym_err_0*errfac
        errfac = mean( rs_err[comp]/m2lstruct[comp].rs_err_0 )
        rs[wbad] = m2lstruct[wbad].rs_0
        rs_err[wbad] = m2lstruct[wbad].rs_err_0*errfac
        errfac = mean( index_err[comp]/m2lstruct[comp].index_err_0 )
        index[wbad] = m2lstruct[wbad].index_0
        index_err[wbad] = m2lstruct[wbad].index_err_0*errfac
    endif

    self->name_label, subtype, name, label, units
    print
    print
    print,'\begin{deluxetable*}{ccccccccc}'
    print,'\tabletypesize{\scriptsize}'
    print,'\tablecaption{\deltamtol\ Statistics for '+name+'\ Bins \label{tab:'+label+'}}'
    print,'\tablewidth{0pt}'
    print,'\tablehead{'
    print,'    \colhead{$\langle$'+name+' $\rangle$}       &'
    print,'    \colhead{$\Delta L_{200}$}         &'
    print,'    \colhead{$\Delta M_{200}$}  &'
    print,'    \colhead{$\left(\frac{\Delta M}{\Delta L}\right)_{200}$}  &'
    print,'    \colhead{$\left(\frac{\Delta M}{\Delta L}\right)_{rmax}$}  &'
    print,'    \colhead{$\langle L^{gal} \rangle $} &'
    print,'    \colhead{\rhalf} &'
    print,'    \colhead{$\alpha$} &'
    print,'    \colhead{$\left(\frac{\Delta M}{\Delta L}\right)_{asym}$}\\  '
    print,'    '+units+' &'
    print,'    ($h^{-2} 10^{10} L_{\odot}$) &'
    print,'    ($h^{-1} 10^{12} M_\odot$) &'
    print,'    ($h M_{\odot}/L_{\odot}$) &'
    print,'    ($h M_{\odot}/L_{\odot}$) &'
    print,'    ($h^{-2} 10^{10} L_{\odot}$) &'
    print,'    ($h^{-1} $Mpc) &'
    print,'    &'
    print,'    ($h M_{\odot}/L_{\odot}$)'
    print,'}'
    print,'\'
    print,'\startdata'

    ; substr & l200+/-err & m200+/-err & M/L200+/-err & M/Llast+/-err & M/Lasymp+/-err
    for i=0L, n-1 do begin
        substr = self->getsubstr(subtype, m[i])

        l200=string(m2lstruct[i].l200/1.e10, format='(F6.2)')
        l200err=string(m2lstruct[i].l200_err/1.e10, format='(F6.2)')

        m200=string(m2lstruct[i].m200/1.e12, format='(F6.2)')
        m200err=string(m2lstruct[i].m200_err/1.e12, format='(F6.2)')

        m2l200=string(m2lstruct[i].m2l200, format='(I0)')
        m2l200err=string(m2lstruct[i].m2l200_err, format='(I0)')

        m2lmax=string(m2lstruct[i].m2lmax, format='(I0)')
        m2lmaxerr=string(m2lstruct[i].m2lmax_err, format='(I0)')


        rss=string(rs[i], format='(F5.2)')
        rserrs=string(rs_err[i], format='(F5.2)')

        indexs=string(index[i], format='(F5.2)')
        indexerrs=string(index_err[i], format='(F5.2)')

        m2lasym=string(m2la[i], format='(I0)')
        m2lasymerr=string(m2la_err[i], format='(I0)')

        mlum=string(m2lstruct[i].meanlum[wrad], format='(F0.2)')
        mlum_err=string(m2lstruct[i].meanlum_err[wrad], format='(F0.2)')
        lumline='$'+mlum+' \pm '+mlum_err+'$ & '

        print,$
            substr+' & '+$
            '$'+l200+' \pm '+l200err+'$ & '+$
            '$'+m200+' \pm '+m200err+'$ & '+$
            '$'+m2l200+' \pm '+m2l200err+'$ & '+$
            '$'+m2lmax+' \pm '+m2lmaxerr+'$ & '+$
            lumline+$
            '$'+rss+' \pm '+rserrs+'$ & '+$
            '$'+indexs+' \pm '+indexerrs+'$ & '+$
            '$'+m2lasym+' \pm '+m2lasymerr+'$ ', format='($,a)'

        if i eq (n-1) then print,'' else print,' \\'

    endfor
    print,'\enddata'
    if subtype eq 'ngals200_12' then begin


        print,'\tablecomments{Mass-to-light ratio statistics for clusters binned by richness \nvir.  \deltal\ is the excess light over the mean luminosity density of the universe, and \deltam\ is the excess mass over the mean mass density. Subscripts $200$ and $rmax$ indicate quantities within $r_{200}$ and the maximum radius \maxinvertrad, respectively.  The subscript $asym$ refers to an asymptotic value from fitting the model described in the text.  Other parameters of this fit are $r_{1/2}$, the radius at which the M/L reaches half its asymptotic value, and $\alpha$, a power law index.  The model is not a physical model, so these uncertainties should be considered lower limits.  $L^{gal}$ is the mean luminosity of the neighboring galaxies used to calculate \deltal. No attempt was made to accountfor mis-centering of the BCGs relative to the halo mass peak; this affects the \deltamvir\ but should have little effect on \mtolvir, and no effect on\mtolmax.  The mean richness is shown but the ranges can be found in \citet{SheldonLensing07}}'

;        print,'\tablecomments{\mtolmax\ is the mass-to-light ratioi within \maxinvertrad. The $\Delta M/\Delta L_{asym}$ is the result of fitting a model.  The model is not a physical model, so these uncertainties should be considered lower limits. No attempt was made to account for mis-centering of the BCGs relative to the halo mass peak; this affects the \deltamvir\ but should have little effect on \mtolvir, and no effect on \mtolmax.  Masses are in units of $h^{-1} M_\odot$, luminosities are measured in units of $h^{-2} 10^{10} L_{\odot}$ and mass-to-light ratios in units of $h M_\odot/L_\odot$.  \rhalf\ is in units of $h^{-1} $Mpc.}'

    endif else begin
        print,'\tablecomments{Same as table \ref{tab:m2lngals} for clusters binned by  \lvir.}'
    endelse
    print,'\end{deluxetable*}'


    wmom, m2lstruct.meanlum[wrad], m2lstruct.meanlum_err[wrad], $
        wm, ws, we
    print,'Mean Mean lum = '+ntostr(wm)+' '+!plusminus+' '+ntostr(we)
end

; make a csv table and tex table
pro maxbcg_m2l::create_paper_csv

    Nfile = self->m2lfile('ngals200_12')
    N = mrdfits(Nfile, 1)
    Lfile = self->m2lfile('ilum200_16')
    L = mrdfits(Lfile, 1)

    mass_sample = self->mass_sample()
    lsample = self->lum_sample()
    cobj = obj_new('maxbcg_correlate', lsample)
    mobj = obj_new('maxbcg_lensing', mass_sample[0]) ; vers doesn't matter here
 
    Nc = cobj->corr_read('dd','jackknife_invert', sub='ngals200_12')
    Lc = cobj->corr_read('dd','jackknife_invert', sub='ilum200_16')
    Nm = mobj->lensread('jackknife_invert', sub='ngals200_12', sample=mass_sample)
    Lm = mobj->lensread('jackknife_invert', sub='ilum200_16', sample=mass_sample)

    obj_destroy, cobj, mobj

    ; These for some of the statistics
    ;Nd = self->lensread('jackknife', sub='ngals200_12', sample=[21,22])
    ;Ld = self->lensread('jackknife', sub='ilum200_16', sample=[21,22])



    dir=self->lensdir('jackknife',sample=[21,22])
    dir = repstr(dir, 'jackknife', 'csv')

    if not fexist(dir) then file_mkdir, dir
    file = 'maxbcg_sample21-22_ngals200_12_ilum200_12_m2l.csv'
    texfile = 'maxbcg_sample21-22_ngals200_12_ilum200_12_m2l.tex'
    file = path_join(dir, file)
    texfile = path_join(dir, texfile)

    print,'Writing to file: ',file
    openw, lun, file, /get


    nbin = n_elements(N)
    nrad = n_elements(N[0].r)

    corrs = 'corr'+ntostr(lindgen(nrad))
    corrs = strjoin(corrs, ',')
    printf, lun, $
        'type,binval,rmpc,bcgfrac,m2l,m2l_err,'+corrs

    for ibin=0L, nbin-1 do begin
        ;mean_bcg_lum = Nc[ibin].bcg_ilum_mean*1.e10
        mean_bcg_lum = Nc[ibin].bcg_ilum_mean*1.e10
        for irad=0L, nrad-1 do begin
            binval = Nm[ibin].mean_ngals200
            if (binval lt 9) then binval=rnd(binval)
            bcgfrac = mean_bcg_lum/N[ibin].lumtot[irad]
            ;if (bcgfrac gt 1) then print,'frac = ',bcgfrac,' Nlum = ', N[ibin].lum[irad],' lumtot = ',N[ibin].lumtot[irad],' bcglum = ',mean_bcg_lum
            printf, lun, $
                'N',binval, N[ibin].r[irad], $
                bcgfrac, $
                N[ibin].m2ltot[irad], N[ibin].m2ltot_err[irad], $
                f='(a,",",e0,",",e0,",",e0,",",e0,",",e0,$)'
            for irad2=0L,nrad-1 do begin
                printf, lun, N[ibin].m2ltot_corr[irad,irad2],$
                    f='(",",e0,$)'
            endfor
            printf, lun
        endfor
    endfor

    nbin = n_elements(L)
    nrad = n_elements(L[0].r)
    for ibin=0L, nbin-1 do begin
        mean_bcg_lum = Lc[ibin].bcg_ilum_mean*1.e10
        ;mean_bcg_lum = Lm[ibin].mean_bcg_ilum*1.e10
        for irad=0L, nrad-1 do begin
            bcgfrac = mean_bcg_lum/L[ibin].lumtot[irad]
            ;if (bcgfrac gt 1) then print,'frac = ',bcgfrac,' Llum = ', L[ibin].lum[irad],' lumtot = ',L[ibin].lumtot[irad],' bcglum = ',mean_bcg_lum
            printf, lun, $
                'L',L[ibin].l200_red, L[ibin].r[irad], $
                bcgfrac, $
                L[ibin].m2ltot[irad], L[ibin].m2ltot_err[irad], $
                f='(a,",",e0,",",e0,",",e0,",",e0,",",e0,$)'
            for irad2=0L,nrad-1 do begin
                printf, lun, L[ibin].m2ltot_corr[irad,irad2],$
                    f='(",",e0,$)'
            endfor
            printf, lun
        endfor
    endfor


    free_lun, lun


    print,'Writing to file: ',texfile
    openw, lun, texfile, /get

    printf, lun,'\begin{deluxetable*}{cccccc}'
    printf, lun,'\tablecaption{Mass-to-light ratio data for \maxbcg\ Clusters \label{tab:m2l}}'
    printf, lun,'\tablewidth{0pt}'
    printf, lun,'\tablehead{'
    printf, lun,'    \colhead{$X$}       &'
    printf, lun,'    \colhead{$\langle X \rangle$}         &'
    printf, lun,'    \colhead{$r$}  &'
    printf, lun,'    \colhead{$L_{BCG}/\Delta L$} &'
    printf, lun,'    \colhead{\deltamtol}  &'
    printf, lun,'    \colhead{$\sigma$(\deltamtol)}\\ '
    printf, lun,'    (N or L) &'
    printf, lun,'     &'
    printf, lun,'    ($h^{-1}$ Mpc) &'
    printf, lun,'    &'
    printf, lun,'    ($h M_{\odot}/L_{\odot}$) &'
    printf, lun,'    ($h M_{\odot}/L_{\odot}$) '
    printf, lun,'}'
    printf, lun,'\'
    printf, lun,'\startdata'

    ;format = '(a," & ",e0," & ",e0," & ",e0," & ",e0,$)'
    format = '(a," & ",e10.3," & ",e10.3," & ",e10.3," & ",e10.3," & ",e10.3,$)'
    nprint = 10
    count=0L
    for ibin=0L, nbin-1 do begin
        mean_bcg_lum = N[ibin].lumtot[nrad-1] - N[ibin].lum[nrad-1]
        for irad=0L, nrad-1 do begin
            binval = N[ibin].n200_red
            if (binval lt 9) then binval=rnd(binval)
            bcgfrac = mean_bcg_lum/N[ibin].lumtot[irad]
            printf, lun, $
                'N',binval, N[ibin].r[irad], $
                bcgfrac,$
                N[ibin].m2ltot[irad], N[ibin].m2ltot_err[irad], $
                f=format
            count = count+1
            if count gt nprint then begin
                printf, lun
                break
            endif
            printf, lun, ' \\'
        endfor
        if count gt nprint then break
    endfor
    printf, lun, '\enddata'
    printf, lun, "\tablecomments{Mass-to-light ratio data corresponding to figure "
    printf, lun, "\ref{fig:ngals200_m2lfits}."
    printf, lun, "The first column indicates binning on either \nvir, labeled ``N'', "
    printf, lun, "or \lvir, labeled ``L''. The second column is the mean value of "
    printf, lun, "\nvir\ or \lvir\ for the bin, with units of either number or  "
    printf, lun, "$10^{10} h^{-2} L_{\odot}$.  Column 3 is the fraction of the total"
    printf, lun, "excess light $\Delta L$ contained in the BCG.  Note the fraction of "
    printf, lun, "light in the BCG can be slightly greater than unity due to noise.  This table "
    printf, lun, "is available in its entirety in the electronic edition of "
    printf, lun, "the {\it Astrophysical Journal}, including the full dimensionless "
    printf, lun, "correlation matrix for the mass-to-light ratio profile.  A portion is shown here for "
    printf, lun, "guidance regarding its form and content.}"
    printf, lun, '\end{deluxetable*}'


    free_lun, lun



end



pro maxbcg_m2l::calculate_omega
    z=0.25

    m2l_last = 349.0
    m2l_last_err = 51.0

    m2l_asym = 335.0
    m2l_asym_err = 9.0

    ; ldens in *physical* coords at z=0.25
    ldens = 3.14e8
    ldens_err = 0.1e8

    omega_last = m2l_last*ldens/(1.0+z)^3/rhocrit()
    omega_last_err = omega_last*sqrt( (m2l_last_err/m2l_last)^2 + $
                                      (ldens_err/ldens)^2 )
    omega_asym = m2l_asym*ldens/(1.0+z)^3/rhocrit()
    omega_asym_err = omega_asym*sqrt( (m2l_asym_err/m2l_asym)^2 + $
                                      (ldens_err/ldens)^2 )


    print,'Omega(last) = '+$
        ntostr(omega_last)+' '+!plusminus+' '+ntostr(omega_last_err)
    print,'Omega(asym) = '+$
        ntostr(omega_asym)+' '+!plusminus+' '+ntostr(omega_asym_err)
end


pro maxbcg_m2l::for_eduardo
	sub='ngals200_12'

	lsamp=4
	msamp=[21,22]
	mstr = strjoin( string(msamp,f='(i0)'), '_')
	lstr = strjoin( string(lsamp,f='(i0)'), '_')
	m=obj_new('maxbcg_lensing', msamp[0])
	l=obj_new('maxbcg_correlate', lsamp)


	outdir='~/tmp/for_eduardo'
	file_mkdir, outdir

	mout=path_join(outdir, 'lensing_sample'+mstr+'_'+sub+'.fits')
	lout=path_join(outdir, 'corr_sample'+lstr+'_'+sub+'.fits')


	mst = m->lensread('jackknife_invert', sub=sub, sample=msamp)
	mtags=$
		['mean_ngals200','sdev_ngals200','err_ngals200', $
		'meanr', 'sigma', 'sigmaerr', 'covariance']
	newtags = $
		['ngals200_mean','ngals200_sdev','ngals200_err', $
		 'rmpc', 'delta_sigma', 'delta_sigma_err', 'delta_sigma_cov']
		   
	mst = extract_tags(mst, mtags)
	mst = rename_tags(mst, mtags, newtags)
	mst.rmpc = mst.rmpc/1000.

	print,'writing: ',mout
	mwrfits, mst, mout, /create


	lst = l->corr_read('dd','jackknife_invert', sub=sub)
	ltags = ['r','radnumdens','radnumdens_err','radnumdens_cov']
	newtags = ['rmpc','numdens','numdens_err','numdens_cov']
	lst = extract_tags(lst, ltags)
	lst = rename_tags(lst, ltags, newtags)

	print,'writing: ',lout
	mwrfits, lst, lout, /create

end

pro maxbcg_m2l__define
    struct = { $
            maxbcg_m2l, $
            mass_sample: ptr_new(), $
            lum_sample: ptr_new(), $
            inherits maxbcg_lensing $
        }
end
