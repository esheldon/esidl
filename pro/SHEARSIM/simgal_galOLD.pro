PRO simgal_gal, sigma, zlens, zsource, binsize, $
                noise=noise, niter=niter,write=write, $
                h=h,omega=omega, yrange=yrange

IF n_params() EQ 0 THEN BEGIN
    print,'-Syntax: simgal_gal, sigma, zlens, zsource, binsize, noise=noise, niter=niter,write=write,h=h,omega=omega'
    return
ENDIF 

time = systime(1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

IF NOT keyword_set(write) THEN write = 0
IF NOT keyword_set(h) THEN h=.7
IF NOT keyword_set(omega) THEN omega=1.
IF n_elements(yrange) EQ 0 THEN yrange=[-.01, .01]

IF NOT keyword_set(noise) THEN BEGIN
    noise=0
    niter = 1
ENDIF ELSE BEGIN
    IF n_elements(niter) EQ 0 THEN niter = 1000
ENDELSE 


sdec = 1200./3600.              
sra  = 1200./3600.

density = 7500.

rmin = double(10.)              ;arcsec
rmin = rmin/3600.               ;degrees
rmax = double(600.)             ;arcsec
rmax = rmax/3600.               ;degrees
binsize=binsize/3600.           ;degrees

nbin = long( (rmax - rmin)/binsize )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output postscript file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

IF write THEN BEGIN 
    fend = '.ps'
    addnum = 1
    add = '_'+ntostr(addnum)
    dir = '/sdss4/data1/esheldon/SHEARSIM/'
    file = dir+'gal_gal'+add+fend
    WHILE exist(file) DO BEGIN
        addnum=addnum+1
        add = '_'+ntostr(addnum)
        file = dir+'gal_gal'+add+fend
    ENDWHILE 
    print,'Writing to file: ',file
ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; declare some arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

default = double(1.e10)
betan = replicate(default, niter, nbin)
berad = betan
betanerr = betan
beraderr = betan
bmeanr = betan

etan = dblarr(nbin)
erad = etan
etanerr = etan
eraderr = etan
meanr = etan

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Estimate shear
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

IF noise THEN BEGIN
    tdens=50000
    sis_shear, zlens, zsource, sigma, sdec, sra, cat, cen, $
      h=h, omegamat=omega, density=tdens

    cenx=cen[0]
    ceny=cen[1]
    
    tan_shear, cat.e1, cat.e2, cat.uncert, $
      cat.dec, cat.ra, cenx, ceny, $
      rmin, rmax, binsize, $
      oetan, oerad, oetanerr, oeraderr, omeanr

ENDIF 

FOR i=0L, niter-1 DO BEGIN 

;    print, format='($,A)', '.'
    IF (i MOD 100 EQ 0 ) THEN print,ntostr(i+1)+'/'+ntostr(niter)
    sis_shear, zlens, zsource, sigma, sdec, sra, cat, cen, $
      h=h, omegamat=omega,density=density, noise=noise

    cenx = cen[0]
    ceny = cen[1]

    ;; IMPORTANT to get x and y right!
    tan_shear, cat.e1, cat.e2, cat.uncert, $
      cat.dec, cat.ra, cenx, ceny, $
      rmin, rmax, binsize, $
      tmpetan, tmperad, tmpetanerr, tmperaderr, tmpmeanr

    betan[i,*] = tmpetan
    berad[i,*] = tmperad
    betanerr[i,*] = tmpetanerr
    beraderr[i,*] = tmperaderr
    bmeanr[i,*] = tmpmeanr
    IF NOT noise THEN print,tmpetan[0], tmpmeanr[0]*3600.
ENDFOR 

IF NOT noise THEN BEGIN
    w=where(betan[0,*] NE 1.e10)
    etan = betan[0,w]
    erad = berad[0,w]
;    etanerr = betanerr[0,w]
;    eraderr = beraderr[0,w]
    etanerr[*]=0.
    eraderr[*]=0.
    meanr = bmeanr[0,w]
    
ENDIF ELSE BEGIN
    vint = .32^2

    FOR nb=0, nbin-1 DO BEGIN
        w=where( betan[*,nb] NE 1.e10 )

;        help,w
        wmom, betan[w,nb], sqrt(vint+betanerr[w,nb]^2), wmean, wsig, werr
        etan[nb] = wmean
        etanerr[nb] = werr

        wmom, berad[w,nb], sqrt(vint+beraderr[w,nb]^2), wmean, wsig, werr
        erad[nb] = wmean
        eraderr[nb] = werr

        meanr[nb] = median( bmeanr[w,nb] )
    ENDFOR 
ENDELSE 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make some plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

yrange=prange(etan/2., erad/2., etanerr/2., eraderr/2.)

IF write THEN makeps, file, /noland
pold=!p.multi
!p.multi = [0,1,2]

tit = 'sigma = '+ntostr(sigma)+' Zlens = '+ntostr(zlens)+' Zsource = ' + $
      ntostr(zsource)
xt='Radius'
yt='"Tangential Shear"'
ploterr, meanr*3600., etan/2.0, etanerr/2.0, $
         psym=1,yrange=yrange,xtitle=xt,ytitle=yt, title=tit
oplot,[0,1200],[0,0]
IF noise THEN oplot, omeanr*3600., oetan/2.0

yt='"Radial Shear"'
ploterr, meanr*3600., erad/2.0, eraderr/2.0, $
         psym=1,yrange=yrange,xtitle=xt,ytitle=yt
oplot,[0,1200],[0,0]

!p.multi=pold
IF write THEN ep

ptime, systime(1)-time

return
END 
