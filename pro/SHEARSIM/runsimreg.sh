#!/bin/sh

idl<<EOF
indir = '/sdss4/data1/esheldon/CORRECTED/'
lcat = mrdfits(indir+'run752_756_lensgal_r_overlap.fit',1)
scat = mrdfits(indir+'run752_756_srcgal_r_overlap.fit',1)
minra = 180.
maxra = 200.

wl=where(lcat.ra le maxra and lcat.ra ge minra, nwl)
ws=where(scat.ra le maxra and scat.ra ge minra, nws)

print
print,'Ra range ',minra,maxra
print,'Using ',nwl,' of the lenses'
print,'Using ',nws,' of the sources'
print,'For the no-noise sim'
print

plot,lcat.ra,lcat.dec,psym=3

sigma = 170.
cutoff = 150.
eend = '1'

runsimreg, eend, lcat=lcat[wl], scat=scat[ws], sigma=sigma, /use, /nocut

eend = '2'
runsimreg, eend, lcat=lcat[wl], scat=scat[ws], sigma=sigma, cutoff=cutoff, /use

eend = '3'
runsimreg, eend, lcat=lcat, scat=scat, sigma=sigma, /use, /nocut, /snoise

eend = '4'
runsimreg, eend, lcat=lcat, scat=scat, sigma=sigma, cutoff=cutoff, /use, /snoise

EOF
