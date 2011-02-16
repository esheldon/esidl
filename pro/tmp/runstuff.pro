PRO runstuff,allign=allign

run1 = 752
run2 = 756
clr  = 2
slength = 10.                   ;arcsec
rfac    = 10.
gridsize = 100.
stepfac = 2.

outdir='/sdss4/data1/esheldon/GAL_GAL/MASSMAPS/'
lcat = mrdfits('/sdss4/data1/esheldon/CORRECTED/run752_756_lensgal_r_overlap.fit',1)
scat = mrdfits('/sdss4/data1/esheldon/CORRECTED/run752_756_srcgal_r_overlap.fit',1)

wl=where(lcat.uncert LT .64,nwl)
ws=where(scat.uncert LT .64,nws)
print,'Using ',ntostr(nwl),' lens galaxies'
print,'Using ',ntostr(nws),' source galaxies'

kappa_map, run1, run2, clr, lcat[wl], $
             slength=slength, $
             rfac = rfac, $
             gridsize=gridsize, $
             stepfac=stepfac, $
             scat=scat[ws], $
             /surface, $
             outdir=outdir, $
             /write, allign=allign, verbose=1

return
END 
