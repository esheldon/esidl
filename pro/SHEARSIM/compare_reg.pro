PRO compare_reg

;; model
sigma = 170.
cutoff = 150.
zs=.285
zl=.172

indir='/sdss4/data1/esheldon/MANYSIM/'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; no cutoff, small region of run
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

f1=indir+'sim11_22_r_N1.dat'
f2=indir+'simregshears_11_22_r_N1.txt'

rdobjshear,f1,str
readcol, f2, mr, sh, sherr, blah2, blah3, blah4, blah5

shmodel = shearsis(sigma, zs, zl, mr)

tit='No cutoff'
xtit='Projected Radius (arcsec)'
ytit='Shear'
line=[0,2,3]
message=['Input','Standard Method','Regression']

plot, str.meanr, str.shear, tit=tit,xtit=xtit,ytit=ytit, line=line[1]
oplot, mr, sh, line=line[2]
oplot, mr, shmodel, line=line[0]
oplot, [0,10000],[0,0]

legend, message, line=line, /right

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; cutoff, small region of run
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

key=get_kbrd(1)
f1=indir+'sim11_22_r_N2.dat'
f2=indir+'simregshears_11_22_r_N2.txt'

rdobjshear,f1,str
readcol, f2, mr, sh, sherr, blah2, blah3, blah4, blah5

shmodel = shearsis_trunc(sigma, cutoff, zs, zl, mr)

tit='Cutoff = 150 arcsec'
xtit='Projected Radius (arcsec)'
ytit='Shear'
line=[0,2,3]
message=['Input','Standard Method','Regression']

plot, str.meanr, str.shear, tit=tit,xtit=xtit,ytit=ytit, line=line[1]
oplot, mr, sh, line=line[2]
oplot, mr, shmodel, line=line[0]
oplot, [0,10000],[0,0]

legend, message, line=line, /right

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; No cutoff, entire region of run, noise
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

key=get_kbrd(1)
f1=indir+'sim11_22_r_N3.dat'
f2=indir+'simregshears_11_22_r_N3.txt'

rdobjshear,f1,str
readcol, f2, mr, sh, sherr, blah2, blah3, blah4, blah5

shmodel = shearsis(sigma, zs, zl, mr)

tit='No cutoff'
xtit='Projected Radius (arcsec)'
ytit='Shear'
line=[0,2,3]
message=['Input','Standard Method','Regression']

ploterr, str.meanr, str.shear, str.shearerr, $
  tit=tit,xtit=xtit,ytit=ytit, line=line[1]
oploterr, mr, sh, sherr, line=line[2]
oplot, mr, shmodel, line=line[0]
oplot, [0,10000],[0,0]

legend, message, line=line, /right

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; cutoff, small region of run, noise
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

key=get_kbrd(1)
f1=indir+'sim11_22_r_N4.dat'
f2=indir+'simregshears_11_22_r_N4.txt'

rdobjshear,f1,str
readcol, f2, mr, sh, sherr, blah2, blah3, blah4, blah5

shmodel = shearsis_trunc(sigma, cutoff, zs, zl, mr)

tit='Cutoff = 150 arcsec'
xtit='Projected Radius (arcsec)'
ytit='Shear'
line=[0,2,3]
message=['Input','Standard Method','Regression']

ploterr, str.meanr, str.shear, str.shearerr, $
  tit=tit,xtit=xtit,ytit=ytit, line=line[1]
oploterr, mr, sh, sherr, line=line[2]
oplot, mr, shmodel, line=line[0]
oplot, [0,10000],[0,0]

legend, message, line=line, /right

return
END 
