PRO galpluscluster

simpctable,/add

rmin = 10./1000.
rmax = 3000./1000.
rr = arrscl(findgen(1000), rmin, rmax)

gnorm = 1.0
gcore=3./1000.
gr = .5
wg=where( abs(rr-gr) EQ min(abs(rr-gr)) )
;galaxy = gnorm*1./sqrt(1. + ( (rr-gr)/gcore )^2 )/gcore
galaxy = gnorm*1./(1. + ((rr-gr)/gcore)^2 )/gcore^2

cnorm=100.*gnorm
ccore=1./1000.
;cluster = cnorm*1./sqrt(1. + (rr/ccore)^2 )/ccore
cluster = cnorm*1./(1. + (rr/ccore)^2)/ccore^2

rho_0 = 200.
tot = cluster + galaxy + rho_0

rhobar = '!S'+!tsym.rho+'!R!A'+!tsym.minus+'!N'
rho0 = !tsym.rho+'!D0!N'
rho0c = !tsym.rho+'!S!D0!N!R!Ucluster!N'
rho0g = !tsym.rho+'!S!D0!N!R!Ugal!N'

title=!tsym.rho+' = '+rho0+$
     'r!U'+!tsym.minus+'2!N      '+rho0c+' = 100'+!tsym.times+rho0g
aplot,!gratio,rr,tot, xrange=[.1,2], xstyle=1,$
      xtitle='r [Mpc]', ytitle=!tsym.rho+'(r)',$
      title=title,/xlog,/ylog,xticklen=0.04,yticklen=0.04
add_labels,xtickv=['.2','.3','.5']
oplot,rr,cluster+rho_0, line=3
oplot, [.1,10],[rho_0, rho_0], line=2

diff = tot[wg[0]] - galaxy[wg[0]]

;oplot, rr, galaxy+diff, color=!orangered

max=max(galaxy+diff)
rad = 100./1000.
oplot, [gr-rad, gr-rad], [rho_0, max],color=!green
oplot, [gr+rad, gr+rad], [rho_0, max],color=!green

;rad = 25./1000.
;oplot, [gr-rad, gr-rad], [rho_0, max],color=!red
;oplot, [gr+rad, gr+rad], [rho_0, max],color=!red

rad = 300./1000.
oplot, [gr-rad, gr-rad], [rho_0, max],color=!blue
oplot, [gr+rad, gr+rad], [rho_0, max],color=!blue

legend,[rhobar,$
        'Cluster + '+rhobar,$
        'Cluster + '+rhobar +' + Galaxy'],$
       line=[2,3,0],/right,box=0,/clear,charsize=1,$
       thick=replicate(!p.thick,3)

legend,['r-r!Dgal!N = 100 kpc','r-r!Dgal!N = 300 kpc'],$
       color=[!green,!blue],line=[0,0],charsize=1,$
       thick=replicate(!p.thick,2),box=0

END 
