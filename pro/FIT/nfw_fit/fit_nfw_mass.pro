pro fit_nfw_mass,r,mass,mass_err,p,perror,yr=yr,ynoedge=ynoedge,xnoedge=xnoedge,tit=tit,xtit=xtit,ytit=ytit, $
                 mass_yfit=mass_yfit, nfw_yfit=nfw_yfit, bias_yfit=bias_yfit, $
                 positive_bias=positive_bias

if n_params() eq 0 then begin
    print,'-syntax fit_nfw_mass,r,mass,mass_err,p,perror,yr=,ynoedge=,xnoedge=,/positive_bias'
    return
endif

recache=0
defsysv,'!linear_xi_cache',Ex=ex
if ex eq 0 or recache then begin
    print,'calculating linear Xi'
    rc=10.0^(dindgen(1000)/200-3)
    linear_xi,rc,xic
    xi2j3,rc,xic,j3c
    print,'caching linear Xi as "!linear_xi_cache"'
    defsysv,'!linear_xi_cache',ptr_new(xic)
    defsysv,'!linear_r_cache',ptr_new(rc)
    defsysv,'!linear_j3_cache',ptr_new(j3c)
endif

if n_elements(yr) eq 0 then begin
    yr=[min(mass)/3.>.01,max(mass)*3.< 5000]*1.e12
endif


if n_elements(xtit) eq 0 then begin
    xtit='r (h!E-1!N Mpc)'
endif


if n_elements(ytit) eq 0 then begin
    sun=sunsymbol()
    munit=' ( h!E -1!N M'+sun+'!N )'
    ytit='M '+munit
endif

ytickf='loglabels'
xtickf='loglabels'

xr=[.01,40]
;yr=[.05,7e3]

xst=1
yst=1
if keyword_set(xnoedge) then begin
    xtit=''
    xtickf='nolabels'
    print,'xnoedge'
endif

if keyword_set(ynoedge) then begin
    ytit=''
    ytickf='nolabels'
    print,'xnoedge'
endif


pplot,r,mass*1.e12,yerr=mass_err*1.e12,/xlog,/ylog,yr=yr,psym=3,yst=yst,xtit=xtit,ytit=ytit,$
  ytickf=ytickf,xtickf=xtickf,xr=xr,xst=xst,tit=tit,aspect=1

pi=replicate({fixed:0L,limits:[0d,0d],limited:[0,0]},3)
pi[0].limits=[0.0,500.0]
pi[0].limited=[1,0]
pi[1].limits=[0.1,100.0]
pi[1].limited=[1,1]
pi[2].limits=[-100,500.0]
pi[2].limited=[1,1]

BL=1.0d
cguess=3.0d
r200g=0.3d
guess=[r200g,cguess,BL]

; Prior bias > 0
if keyword_set(positive_bias) then begin 
    pi[2].limits = [0,500.0]
endif 


func='nfw_mass_plus_lin'
y=mass
yerr=mass_err
p=mpfitfun(func,r,y,yerr,guess,perror=perror,covar=covar,parinfo=pi,/quiet)
mass_yfit=nfw_mass_plus_lin(r,p)*1.e12

nfw_yfit=nfw_mass_plus_lin(r,[p[0:1],0.0])*1.e12
oplot,r,nfw_yfit,color=c2i('darkgreen')
bias_yfit=nfw_mass_plus_lin(r,[1e-5,p[1:2]])*1.e12
oplot,r,bias_yfit,color=c2i('blue')
oplot,r,mass_yfit,color=c2i('magenta')


return
end
