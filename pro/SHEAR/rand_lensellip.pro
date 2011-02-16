PRO rand_lensellip, nrand

IF n_params() LT 1 THEN BEGIN
    print,'-Syntax: rand_lensellip, nrand, lcat, scat=scat'
    return
ENDIF 

outdir='/sdss4/data1/esheldon/GAL_GAL/MASSMAPS/'

lcat = mrdfits('/sdss4/data1/esheldon/CORRECTED/run752_756_lensgal_r_overlap.fit',1)
scat = mrdfits('/sdss4/data1/esheldon/CORRECTED/run752_756_srcgal_r_overlap.fit',1)

wl=where(lcat.uncert LT .64,nwl)
ws=where(scat.uncert LT .64,nws)
print,'Using ',ntostr(nwl),' lens galaxies'
print,'Using ',ntostr(nws),' source galaxies'


r1=752
r2=756
clr=2

newlcat = lcat[wl]
e1 = lcat[wl].e1
e2 = lcat[wl].e2

print
print,' Doing ',ntostr(nrand),' randomizations'
print
FOR i=0, nrand-1 DO BEGIN 

    theta = arrscl( randomu(seed, nwl), -!pi/2.,!pi/2., arrmin=0., arrmax=1. )
    cos2 = cos(2.*theta)
    sin2 = sin(2.*theta)

    newlcat.e1 =  cos2*e1 + sin2*e2
    newlcat.e2 = -sin2*e1 + cos2*e2

    kappa_map, r1, r2, clr, newlcat, scat=scat[ws], $
            /allign, $
            gridsize=100., $
            slength=10.0, $
            rfac = 10., $
            stepfac=2., $
            /write, $
            /surface, $
            outdir=outdir, $
            verbose = 0
ENDFOR 
return
END 
