pro shapenoise_vagc, struct, corr, weighted=weighted

    if n_elements(struct) eq 0 then begin

        dir='/global/data/vagc-dr6/vagc0/sdss/parameters'
        pattern = path_join(dir, 'calibObj*.fits')
        files = file_search(pattern)
        ; read a few of them
        files = files[1000:2400]


        columns = ['m_e1', 'm_e1e1err', 'm_rr_cc', 'm_cr4', 'm_e1_psf', $
            'm_rr_cc_psf', 'm_cr4_psf', 'modelflux']
        struct = mrdfits_multi(files, columns=columns)

        v=obj_new('vagc')
        ;struct = v->read('vagc','calibobj',run=756,camcol=3)
        corr = v->getsmear(struct)
        obj_destroy,v
    endif


    cmin=0.05
    cmax=0.8
    w=where(corr[2,*] gt cmin and corr[2,*] lt cmax and struct.m_e1[2] ge -1 and struct.m_e1[2] le 1)

    !p.multi=[0,0,2]

    e1corr = struct[w].m_e1[2] - struct[w].m_e1_psf[2]*corr[2,w]
    e1err = struct[w].m_e1e1err[2]
    e1corr2 = e1corr/(1.0-corr[2,w])
    e1err2 = e1err/(1.0 - corr[2,w])

    rmag = 22.5-2.5*alog10(struct.modelflux[2])

    plothist,corr[2,w],bin=0.01, min=0, max=1.1
    plothist,struct[w].m_e1[2],bin=0.01
    plothist,e1corr,bin=0.01,/over,color=c2i('darkgreen')
    plothist,e1corr2,bin=0.01,/over,color=c2i('blue')

    w2 = where(corr[2,w] gt cmin and corr[2,w] lt cmin+0.05 and $
                rmag[w] gt 16.5 and rmag[w] lt 17.5)
    print,sdev(struct[w[w2]].m_e1[2])
    print,sdev(e1corr[w2])
    print,sdev(e1corr2[w2])

    !p.multi=0



    key=prompt_kbrd('hit a key')

    nperbin=10000

    if keyword_set(weighted) then wts = 1.0/(e1err^2 + 0.32^2)
    bs=binner(rmag[w], e1corr, nperbin=nperbin, weights=wts, rev=rev)
    if keyword_set(weighted) then wts = 1.0/(e1err2^2 + 0.32^2)
    bs2=binner(rmag[w], e1corr2, nperbin=nperbin, weights=wts, rev=rev)

    begplot,'~/tmp/shapenoise.eps',/color,/encapsulated
    pplot, bs2.xmean, bs2.ysdev, yerr=bs2.ysdev/sqrt(2*nperbin), $
        aspect=1, xrange=[14,19], yrange=[0.25,0.5], xsty=3, ysty=3, $
        xtitle='rmag', ytitle='e1 sdev'
    pplot, [0, 100], [0.32, 0.32], /over
    ;pplot, bs2.xmean, bs2.ysdev, yerr=bs2.ysdev/sqrt(2*nperbin), $
    ;    /over, color=c2i('darkgreen')
    endplot
    return


    key=prompt_kbrd('hit a key')

    rmag = 22.5-2.5*alog10(struct.modelflux[2])

    if keyword_set(weighted) then wts = 1.0/(e1err^2 + 0.32^2)
    bs=binner(corr[2,w], e1corr, nperbin=nperbin, weights=wts, rev=rev)
    if keyword_set(weighted) then wts = 1.0/(e1err2^2 + 0.32^2)
    bs2=binner(corr[2,w], e1corr2, nperbin=nperbin, weights=wts, rev=rev)

    pplot, bs.xmean, bs.ysdev, yerr=bs.ysdev/sqrt(2*nperbin), $
        aspect=1, yrange=[0.1,0.4], $
        xtitle='rsmear', ytitle='e1 sdev'
    pplot, bs2.xmean, bs2.ysdev, yerr=bs2.ysdev/sqrt(2*nperbin), $
        /over, color=c2i('darkgreen')
    colprint, bs.ysdev, bs2.ysdev


    return

    key=prompt_kbrd('hit a key')

    rmag = 22.5-2.5*alog10(struct.modelflux[2])

    wts = 1.0/(e1err^2 + 0.32^2)
    bs=binner(1.0/wts, e1corr, nperbin=nperbin, weights=wts, rev=rev)
    wts = 1.0/(e1err2^2 + 0.32^2)
    bs2=binner(1.0/wts, e1corr2, nperbin=nperbin, weights=wts, rev=rev)

    pplot, bs.xmean, bs.ysdev, yerr=bs.ysdev/sqrt(2*nperbin), $
        aspect=1, yrange=[0.1,0.4], $
        xtitle='err^2 + 0.32^2', ytitle='e1 sdev'
    pplot, bs2.xmean, bs2.ysdev, yerr=bs2.ysdev/sqrt(2*nperbin), $
        /over, color=c2i('darkgreen')
    colprint, bs.ysdev, bs2.ysdev






end
