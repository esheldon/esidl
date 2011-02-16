#!/bin/sh

idl<<EOF
colors=['u','g','r','i','z']

clr=1
dir = '/sdss3/data1/corrected/corr752/1/combined/'
f=dir+'run752_756_lensgal_'+colors[clr]+'_overlap.fit'
lcat=mrdfits(f,1)

voronoi_shear, lcat, clr

EOF


idl<<EOF
colors=['u','g','r','i','z']

clr=2
dir = '/sdss3/data1/corrected/corr752/1/combined/'
f=dir+'run752_756_lensgal_'+colors[clr]+'_overlap.fit'
lcat=mrdfits(f,1)

voronoi_shear, lcat, clr

EOF


idl<<EOF
colors=['u','g','r','i','z']

clr=3
dir = '/sdss3/data1/corrected/corr752/1/combined/'
f=dir+'run752_756_lensgal_'+colors[clr]+'_overlap.fit'
lcat=mrdfits(f,1)

voronoi_shear, lcat, clr

EOF
