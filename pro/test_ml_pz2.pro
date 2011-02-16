PRO test_ml_pz2,num,delin,del_ml,sig_ml,zs_sig=zs_sig,sig_e=sig_e

IF n_params() EQ 0 THEN BEGIN
    print,'-syntax test_ml_pz,num'
    return
endif

IF n_elements(del) EQ 0 THEN Del=.0111
IF n_elements(zl_mean) EQ 0 THEN zl_mean=.1
IF n_elements(zs_mean) EQ 0 then zs_mean=.15
IF n_elements(zs_sig) EQ 0 THEN zs_sig=.05
IF n_elements(sig_e) EQ 0 THEN sig_e=.35

zl=(randomn(seed,num)+5)^2
zl=zl*zl_mean/mean(zl)

;zl = replicate(0.1, num)

zs=(randomn(seed,num)+5)^2
zs=zs*zs_mean/mean(zs)

zsig=randomu(seed,num)*.07+.03-.05+zs_sig

zbar=zs+randomn(seed,num)*zsig

siginv = sigmacritinv(zl, zs)

!p.multi=[0,0,3]

plothist,zl,bin=0.01,/norm
plothist,zs,bin=0.01,/overplot,color=!blue,/norm
;oplot,[zl[0], zl[0]], [0, 10000], color=!red
plothist,siginv,bin=.01
shear=siginv*Del

print,'mean shear ',mean(shear)
print,'max shear ',max(shear)
print,'min shear ',min(shear)

sige=replicate(1.0,num)*sig_e
;sige=randomn(seed,num)*.1+sig_e > 0.1
vare=sige^2
e=2*shear+randomn(seed,num)*sig_e
misc=fltarr(num)
visc=misc

perlast=-1
FOR i=0L, num-1 DO BEGIN
    mean_inv_sig_crit,zbar(i),zsig(i),zl(i),mi,vi
    misc(i)=mi
    visc(i)=vi

    IF num GT 100 then begin
        per=fix(100.0*i/float(num))
        if per gt perlast then begin
            perlast=per
            PRINT, FORMAT = '(3x,1(i3.2,a3,3x))', per, '% '+STRING(BYTE(141))
        endif
    ENDIF
    
ENDFOR 

IF num GT 100 THEN print,FORMAT = '(3x,i3.2)',100

t=.5*total(e*misc/vare)
b1=total(misc^2/vare)
b2=total(visc*(1.0-(e^2/vare))/vare)
b=b1+b2

del_ml=t/b
var_ml=.25/b
sig_ml=sqrt(var_ml)
delin=del

print,'del ',del
print,'del_ml',del_ml,'  +/-  ',sqrt(var_ml)

print,'diff ', (del_ml-del)/sqrt(var_ml), 'sigma'

print,b1,b2,b

sig_del_ml = sqrt(var_ml)
mind = del_ml - 3.5*sig_del_ml
maxd = del_ml + 3.5*sig_del_ml
dd = arrscl( findgen(1000), mind, maxd )
pd = gaussprob(dd, del_ml, sig_del_ml)
plot, dd, pd
oplot, [del, del], [0, 10000], color=!red

!p.multi=0

return
END











