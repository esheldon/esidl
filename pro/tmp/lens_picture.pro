
pro lens_picture,imbg,imbr,imbi

    sx = 2000L
    sy = 2000L

    make_gaussian,g,size=[sx,sy],$
	counts=1000.,fwhm=1200.



    g = g*(200.0/max(g))
    ;g = g+100.0
    add_noise,g,go,gain=10
    ;g=go-100.0

;    run = 752
;    camcol = 6
;    rerun=1
;    field=489
;    id = 134

    run=752
    camcol = 1
    rerun = 1
    field = 244
    id = 39

    fetch_dir,run,camcol,rerun,dir,atldir
    read_tsobj,dir,l,start=field,nframes=1
;    photo_match,l.run,l.rerun,l.camcol,l.field,l.id,$
;	752,1,5,62,23,m1,m2
    w=where(l.id EQ id)
    print,w
    get_atlas,l,w,dir=atldir,img=img,imr=imr,imi=imi,/nodisplay,maxsize=[1000,1000]
    xs=(size(img))(1)
    ys=(size(img))(2)

;    g = replicate(0., sx, sy)
;    index = lindgen(sx*sy)
;    x = index MOD sx
;    y = index/sx

    cen=[(sx-1.)/2.,(sy-1.)/2.]

    xl=(sx-xs)/2.
    xh=(sx+xs)/2.-1
    yl=(sy-ys)/2.
    yh=(sy+ys)/2.-1

;    radius = xs/2.
;    halor = 2.5*radius
;    print,halor
;return
;    rr = sqrt( (x - cen[0])^2 + (y - cen[1])^2 )

;    w = where(rr LE halor, nw)
    
    imbg=g
    imbr=g
    imbi=g
;    sigma_clip, imr, mr, sr,nsig=3,niter=4,/silent

;    imbg[index] = 50.*sr/sqrt(1. + (.2*rr/xs)^2)
;    imbr[index] = 50.*sr/sqrt(1. + (.2*rr/xs)^2)
;    imbi[index] = 50.*sr/sqrt(1. + (.2*rr/xs)^2)

;    imbg[index[w]] = 2.*sr
;    imbr[index[w]] = 2.*sr
;    imbi[index[w]] = 2.*sr

;    imbg(xl:xh,yl:yh)=imbg(xl:xh,yl:yh)+img
;    imbr(xl:xh,yl:yh)=imbr(xl:xh,yl:yh)+imr
;    imbi(xl:xh,yl:yh)=imbi(xl:xh,yl:yh)+imi

    imbg(xl:xh,yl:yh)=imbg(xl:xh,yl:yh)+(img-1000.)
    imbr(xl:xh,yl:yh)=imbr(xl:xh,yl:yh)+(imr-1000.)
    imbi(xl:xh,yl:yh)=imbi(xl:xh,yl:yh)+(imi-1000.)
    

    return
    end
