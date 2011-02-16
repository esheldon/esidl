FUNCTION deproject::init
  return, 1
END 

;docstart::deproject::ds2drho_mass
; NAME: deproject::ds2drho_mass
;  
; PURPOSE:
; Convert DeltaSigma measurements to the 3D over-density
; \Delta\rho. Propagate its covariance matrix and compute and apply the
; endpoint correction. Uses an interpolation methods based on fitting
; powerlaws to an N-point neighborhood of each bin. Results in some 
; smoothing but removed bias rsulting from linear or polynomial 
; interpolation for data that resembles local powerlaws like most
; real data. Furthermore calculate aperture masses with both formulas.
; Propogate covariances matrices for all. Also calculate circular
; velocities V_c=sqrt(G M(r)/r) with error.
;
; CALLING SEQUENCE:
; ds2drho_mass,r,ds,dscov,drho,drhocov,mi,micov,mo,mocov,vel,velerr,velcov,m,merr,$
; /corr2d, numd=numd,numx=numx,nocov=nocov,sm=sm      
;
; INPUTS:
; r  - array of binned radii, supposed to be the radius such that ds=DeltaSigma(r)  
;     may need to make corrections for effects of binning. Should be
;     in Mpc/h for numerical reasons.
; ds - DeltaSigma
; dscov - covariance matrix of ds 
; OPTIONAL KEYWORD INPUTS:
; nocov - if set will not compute the covariance matrices
; nsigo - the number of sigma above 0 that the data must be to be
;         interpolated without an offset. Default 1.
; endslope - The negative endpoint DeltaSigma logarithmic
;            slope. Default 1.
; noplot - do not plotting
; corr2d - treat this as a 2D correlation function like Sigma rather
;          than DeltaSigma. The mass is then really Mass - Mass[0]
;          since the inner mass point is undetermined.
; OUTPUTS:
; drho - 3D over-density computed from ds
; drhocov - the covariance matrix of drho computed from dscov using the
;         Jacobian, numerically differentiated.
; mi -  aperture mass by the inside formula
; mo - aperture mass by the outside formula
; micov - covariance matrix of mi
; mocov -covariance matrix of mo
; m - a simple average of mi and mo (can try weighted average if you want)
; merr - merr error on m
; vel -  circular velocity sqrt(G M(r)/r) in km/s
; velerr - error on vel
; velcov - covariance matrix on velocity
; COMMENTS:
;
;  Assumes that rho_m is constant which is true in comoving
;  coordinates only so one should work with comoving.
;  Uses rho_crit but keeps Omega_m=1.0. Thus it really outputs 
;  drho_gm*Omega_m. One has to divide the output by Omega_m if one
;  really want drho_gm. Don't forget to divide drhoerr by Omega_m and also
;  drhocov by Omega_m^2. Doesn't effect mass formulas.
;  the output drho has one less element then r and ds since the last
;  point cannot be evulated from the data.
;  
; EXAMPLES:
; IDL> file="Erins-absolutely-rediculously-long-filename.fit"
; IDL> str=mrdfits(file,1)
; IDL> r=str.meanr/1000.0 ; put in Mpc/h 
; IDL> ds=str.sigma
; IDL> dscov=str.covariance
; IDL> ds2drho_mass,r,ds,dscov,drho,drhocov,mi,micov,mo,mocov,vel,velerr,velcov,m,merr
; METHOD:
;
; PROCEDURES CALLED:
;  ds2drho
;  ds_drho2mass
;  tvim2
;  cov2corr()
; REVISION HISTORY:
;   20-May-2004  Written by David Johnston, Princeton.
;   09-26-2005    Added the corr2d keyword and functionality- Dave Caltech
;  
;docend::deproject::ds2drho_mass

function deproject::invert, r, ds, dscov, nocov=nocov, endslope=endslope, nsigo=nsigo, noplot=noplot, endc=endc, linmap=linmap, corr2d=corr2d

  self->ds2drho_mass, r, ds, dscov,drho,drhocov,mi,micov,mo,mocov,vel,velcov,m,mcov,$
    nocov=nocov,endslope=endslope,nsigo=nsigo,noplot=noplot,endc=endc,linmap=linmap,corr2d=corr2d

  nr = n_elements(r)
  nr2 = nr-1
  di = lindgen(nr2)*(nr2+1)

  st = $
    {r: r, $
     dsigma: ds, $
     dsigma_cov: dscov, $
     ir: r[0:nr-2], $
     drho: drho, $
     drho_err: sqrt(drhocov[di]), $
     drho_cov: drhocov, $
     massin: mi, $
     massin_err: sqrt(micov[di]), $
     massin_cov: micov, $
     massout: mo, $
     massout_err: sqrt(mocov[di]), $
     massout_cov: mocov, $
     masscomb: m, $
     masscomb_err: sqrt(mcov[di]), $
     masscomb_cov: mcov, $
     vel: vel, $
     vel_err: sqrt(velcov[di]), $
     vel_cov: velcov}
  return, st

end 
pro deproject::ds2drho_mass,r,ds,dscov,drho,drhocov,mi,micov,mo,mocov,vel,velcov,m,mcov,$
               nocov=nocov,endslope=endslope,nsigo=nsigo,noplot=noplot,endc=endc,linmap=linmap,corr2d=corr2d         

;a wrapper for ds2drho and ds_drho2mass

if n_params() eq 0 then begin
print,'-syntax ds2drho_mass,r,ds,dscov,drho,drhocov,mi,micov,mo,mocov,vel,velcov,m,mcov'
print,'nocov=nocov,endslope=endslope,nsigo=nsigo,endc=endc'
return
endif

n=n_elements(r)
rmin=min(r)
rmax=max(r)
xr=[(.5*rmin)>.001,rmax]

if not keyword_set(noplot) and !d.name eq 'X' then window,0
self->ds2drho,r,ds,drho,drhocov,drhoerr,endc,dscov=dscov,$
  yr=yr,nocov=nocov,noplot=noplot,corrmat=corrmat,endslope=endslope,nsigo=nsigo,corr2d=corr2d
drhoerr=sqrt(drhocov(lindgen(n-1)*n))

ni = n_elements(drho)
pr = r[0:ni-1]

dserr=sqrt(dscov(lindgen(n)*(n+1)))
if not keyword_set(noplot) and !d.name eq 'X' then window,1

if keyword_set(corr2d) eq 0 then begin
    self->ds_drho2mass,r,ds,dserr,drho,drhoerr,mass_in,mass_out,yr=yr,noplot=noplot,endslope=endslope,nsigo=nsigo
endif else begin
    self->drho2mass,r,drho,drhoerr,mass,yr=yr,noplot=noplot,nsigo=nsigo
    mass_out=mass
    mass_in=mass
endelse

mati=fltarr(n,n-1)
mato=fltarr(n,n-1)
dsp=ds
dsm=ds

FOR i=0, n-1 DO BEGIN
    ;a small stop size
    h=abs(0.005*abs(ds(i))+dserr(i))
    dsp(i)=ds(i)+h
    dsm(i)=ds(i)-h
    self->ds2drho,r,dsp,drhop,dscov=dscov,/nocov,/noplot,endslope=endslope,nsigo=nsigo,corr2d=corr2d
    self->ds2drho,r,dsm,drhom,dscov=dscov,/nocov,/noplot,endslope=endslope,nsigo=nsigo,corr2d=corr2d
    if keyword_set(corr2d) eq 0 then begin
        self->ds_drho2mass,r,dsp,dserr,drhop,drhoerr,mass_in_p,mass_out_p,/noplot,endslope=endslope,nsigo=nsigo
        self->ds_drho2mass,r,dsm,dserr,drhom,drhoerr,mass_in_m,mass_out_m,/noplot,endslope=endslope,nsigo=nsigo
    endif else begin
        self->drho2mass,r,drhop,drhoerr,massp,/noplot,nsigo=nsigo
        mass_in_p=massp
        mass_out_p=massp
        self->drho2mass,r,drhom,drhoerr,massm,/noplot,nsigo=nsigo
        mass_in_m=massm
        mass_out_m=massm
    endelse

    deri=(mass_in_p-mass_in_m)/(2*h)
    dero=(mass_out_p-mass_out_m)/(2*h)
    mati(i,*)=deri
    mato(i,*)=dero
    dsp(i)=ds(i)
    dsm(i)=ds(i)
ENDFOR

linmap=mato                     ;a linear mapping of delta sigma to mass
mass_in_cov=mati##(dscov##transpose(mati))
mass_in_err=mass_in_cov(lindgen(n-1)*n)
wbad=where(mass_in_err lt 0,nbad)
if nbad gt 0 then begin
    print,"Warning mass_in_cov has negitive diagonals"
    mass_in_err(*)=-999.0
endif else mass_in_err=float(sqrt(mass_in_err))

mass_out_cov=mato##(dscov##transpose(mato))
mass_out_err=mass_out_cov(lindgen(n-1)*n)
wbad=where(mass_out_err lt 0,nbad)
if nbad gt 0 then begin
    print,"Warning mass_out_cov has negitive diagonals"
    mass_out_err(*)=-999.0
endif else mass_out_err=float(sqrt(mass_out_err))

dscor = cov2corr(dscov, /fix)
drhocor = cov2corr(drhocov, /fix)
mass_in_cor = cov2corr(mass_in_cov, /fix)
mass_out_cor = cov2corr(mass_out_cov, /fix)



pmold=!p.multi
pcharold = !p.charsize
if not keyword_set(noplot) then begin
    if !d.name eq 'X' then window,2
    !p.charsize=1
    loadct,0,/silent
    !p.multi = [0,2,2]
    tvim2,dscor,tit="DS",r=[-1,1],/scale
    tvim2,drhocor,tit="DRHO",r=[-1,1],/scale
    tvim2,mass_in_cor,tit="MASS IN",r=[-1,1],/scale
    tvim2,mass_out_cor,tit="MASS_OUT",r=[-1,1],/scale
    simpctable
    if !d.name eq 'X' then window,3
    !p.multi=[0,2,3]
endif

diag=lindgen(n-1)*n
err1=sqrt(mass_in_cov(diag))
err2=sqrt(mass_out_cov(diag))

if not keyword_set(noplot) then begin
    !p.charsize=2
    yr=[0.8*min(mass_in>mass_out)>.01,1.2*max(mass_in<mass_out)]
    pplot,pr,mass_in,/xlog,/ylog,psym=3,tit="MASS IN",xr=xr,/xst,yr=yr,/yst, $
      aspect=1, xtickf='loglabels', ytickf='loglabels'
    errplot,r,mass_in+err1,mass_in-err1
    
    pplot,pr,mass_out,/xlog,/ylog,psym=3,tit="MASS OUT",xr=xr,/xst,yr=yr,/yst, $
      aspect=1, xtickf='loglabels', ytickf='loglabels'
    errplot,r,mass_out+err2,mass_out-err2
endif

;just average them for simplicity
;AND averge their covariances matrices 
;they are very correlated so perhaps the errors don't get
;much smaller
;This probably isn't the best thing to do.

mass=(mass_in+mass_out)/2.0
;masserr=(err1+err2)/2.0
mcov=(mass_in_cov+mass_out_cov)/2.0
masserr=sqrt(mcov(diag))
mcor = cov2corr(mcov, /fix)

if not keyword_set(noplot) then begin
    pplot,pr,mass,psym=3,xr=xr,/xst,/xlog,/ylog,tit="MASS",yr=yr,/yst,$
      aspect=1, xtickf='loglabels', ytickf='loglabels'
    errplot,r,mass+masserr,mass-masserr
    
    pplot,pr,mass/masserr,/xlog,tit="S/N",xr=xr,/xst, aspect=1.62, $
      xtickf='loglabels'
    loadct,0,/silent
    tvim2,mcor,r=[-1,1],/scale
    simpctable
endif

mo=mass_out
mi=mass_in
mocov=mass_out_cov
micov=mass_in_cov
m=mass
merr=masserr

;rotational circular velocity
;assumes M is in units of 1x10^12 M_sun
;and r is in Mpc/h
;units of vel is km/s
;vel=sqrt(G*M(r)/r)
;The number I get is 65.583
vconst=65.583

massmin=0.0001
massu=mass > massmin
;don't allow negative mass for velocity calculation

vel=vconst*sqrt(massu/r)
velerr=0.5*vel*masserr/massu
vmvec=reform(vel/massu)
velcov=0.25*mcov*(vmvec##vmvec) ;outer product of vectors
if not keyword_set(noplot) then begin
    pplot,pr,vel,/xlog,psym=1,xtit="r Mpc/h",ytit="Rot. Vel.",$
      aspect=1.62, xtickf='loglabels'
    errplot,r,vel-velerr,vel+velerr
endif

!p.multi=pmold
!p.charsize=pcharold
return
end



;docstart::deproject::ds2drho
; NAME: deproject::ds2drho
;  
; PURPOSE:
; Convert DeltaSigma measurements to the 3D over-density deltaRho
; Propagate its covariance matrix and compute and apply the
; endpoint correction. Uses a hybrid interpolation methods based 
; which is a powerlaw except for low S/N data when it becomes
; a powerlaw plus offset.
;
; CALLING SEQUENCE:
;   ds2drho,r,ds,drho,dscov,drhocov,drhoerr,endc,$
;    num=num,yr=yr,nocov=nocov,noplot=noplot,corrmat=corrmat,sm=sm
;
; INPUTS:
; r  - array of binned radii, supposed to be the radius such that ds=DeltaSigma(r)  
;     may need to make corrections for effects of binning. Should be
;     in Mpc/h for numerical reasons.
; ds - DeltaSigma
;  
; OPTIONAL INPUTS:
; dscov - delta sigma covariance matrix, needed if drhocov is desired 
;  
; OPTIONAL KEYWORD INPUTS:
; yr - yrange for some plots that are made by default
; nocov - if set will not compute the covariance matrix of drho
; noplot - if set will not make plots
; nsigo - the number of sigma above 0 that the data must be to be
;         interpolated without an offset. Default 1.
; endslope - The negative endpoint DeltaSigma logarithmic
;            slope. Default 1.
; corr2d - treat this as a 2D correlation function like Sigma rather than DeltaSigma 
;
; OUTPUTS:
; drho - the 3D over-density equal to rho_bar * xi
;
; drhocov - the covariance matrix of drho computed from dscov using the
;         Jacobian, numerically differentiated.
; OPTIONAL KEYWORD OUTPUTS
; corrmat - the correlation matrix from drhocov
; COMMENTS:
;
;  the output rho has one less element then r and ds since the last
;  point cannot be determined from the data.
;  
; EXAMPLES:
; IDL> file="Erins-absolutely-rediculously-long-filename.fit"
; IDL> str=mrdfits(file,1)
; IDL> r=str.meanr/1000.0 ; put in Mpc/h 
; IDL> ds=str.sigma
; IDL> dscov=str.covariance
; IDL> ds2drho,r,ds,drho,drhocov,drhoerr,endc,dscov=dscov
; ;or if one doesn't have coavriance matrix, fake it
; IDL> n=n_elements(r)
; IDL> dscov=fltarr(n,n)
; IDL> diag=indgen(n)*(n+1)
; IDL> dserr=str.sigmaerr
; IDL> dscov(daig)=dserr^2
; IDL> ds2drho,r,ds,drho,drhocov,drhoerr,endc,dscov=dscov
;
; METHOD:
; Uses a method based on powerlaw interpolation between points when
; ds > nsigo*dserr. For low S/N data or negative data an offset is
; added to make this condition true before determining the 2 powerlaw
; paramters. Then the offset is subtracted off. So the resulting interpolation
; is powerlaw plus constant for low S/N data. For data that looks
; locally like a powerlaw this is unbiased for high S/N data and
; mostly unbiased for low S/N data. For low S/N data the interpolation
; is closer to linear interpolation.
; The covariance matrix is evaulated by a numerical differentiation
; computation of the Jacobian and then this can be used to propagate
; the dscov.  
; PROCEDURES CALLED:
;  xi_pl_int - a numerical integration wrapper
;  interp_hyb()
;
; REVISION HISTORY:
;   20-Sept-2004  Written by David Johnston, Princeton.
;   06-Jul-2005   Converted from ds2xi, never should have used xi, David Johnston
;   09-26-2005    Added the corr2d keyword and functionality- Dave 
;docend::deproject::ds2drho


pro deproject::ds2drho, r, ds, drho, drhocov, drhoerr, endc, $
             dscov=dscov, yr=yr, nocov=nocov, noplot=noplot, corrmat=corrmat,$
             nsigo=nsigo, endslope=endslope, corr2d=corr2d

if n_params() LT 3 then begin
    on_error, 2
    print,'dep->ds2drho, r, deltasig, drho, drhocov, drhoerr, endc, '
    print,'              dscov=, /nocov, corrmat=, '
    print,'              /corr2d, '
    print,'              endslope=, nsigo=, '
    print,'              yr=, /noplot'
    message,'Halting'
endif

if n_elements(nsigo) eq 0 then nsigo=1.0
;the number of sigma to offset before fitting the powerlaw
if n_elements(endslope) eq 0 then endslope=1.0
;the endpoint slope of delta sigma for the endpoint correction
if n_elements(endslope_npoints) eq 0 then endslope_npoints=3
;the number of points at the end to use to fit the endpoint amplitude

if max(r) gt 100 then begin
    print,"convert r to Mpc/h"
    print,"I have trouble with big numbers"
    ;fix this, not sure why, numerical overflow?
    return
endif

pmold=!p.multi
pcharold = !p.charsize

n=n_elements(r)

nrmax=30.0         ;has to do with where to do final integration completely analytically
                   ;start this number times rmax, go to
                   ;infinity. A tiny correction.     
rmax=max(r)
rmin=min(r)
xr=[(.5*rmin) > .001 ,rmax*2]

if n_elements(dscov) gt 0 then begin
    dserr=dscov(lindgen(n)*(n+1))
    wbad=where(dserr le 0,nbad)
    if nbad gt 0 then begin
        print,"dscov must be positive definite"
        !p.multi=pmold
        return
    endif
    dserr=float(sqrt(dserr))
endif else begin
    ;just fake the errors, constant S/N
    dserr=abs(0.01*ds)
endelse

if not keyword_set(noplot) then doplot=1 else doplot=0


if n_elements(yr) eq 0 then yr=[(0.9*min(ds-dserr))>0.01,max(ds+dserr)*1.1]

if doplot then begin
    !p.multi=[0,2,2]
    !p.charsize=1
    pplot,r,ds,/xlog,/ylog,psym=1,yr=yr,/yst,$
      xtickf='loglabels', ytickf='loglabels', $
      aspect=1, $
      xtit="R",ytit="DS",tit="Input Delta Sigma",xr=xr,/xst
    errplot,r,ds-dserr,ds+dserr
    nip=100
    ipoints=fltarr(nip,n-1)
                                ;log interpolated points for plotting
    for i=0, n-2 do begin
        imin=alog(r[i])
        irange=alog(r[i+1])-alog(r[i])
        ipoints[*,i]=exp(findgen(nip)*irange/(nip-1)+imin)
        oplot,$
            ipoints[*,i], $
            self->interp_hyb(r[i],r[i+1],ds[i],ds[i+1],dserr[i],dserr[i+1],ipoints[*,i],nsig=nsigo),$
          color=c2i('green'),psym=3
    endfor
endif


slope=fltarr(n)
amp=slope
off=slope

for i=0, n-2 do begin
    ;get the slope and amplitude and offset for this interval
    junk=self->interp_hyb(r[i],r[i+1],ds[i],ds[i+1],dserr[i],dserr[i+1],xjunk,sl=sli,amp=ampi,off=offi,nsig=nsigo)
    ;print,sli,ampi,offi,ds(i)
    slope(i)=-sli
    amp(i)=ampi
    off(i)=offi
endfor
slope(n-1)=endslope

wnfin=where(finite(slope) eq 0,wnf)
if wnf gt 0 then begin
    print,'some slopes not finite'
    stop
endif

drho=fltarr(n-1)
endc=drho

;get the last few points to fix the endpoint powerlaw amplitude
wuse=indgen(endslope_npoints)+n-endslope_npoints
xuse=r(wuse)^(-endslope)
yuse=ds(wuse)
vuse=dserr(wuse)^2
;inverse variance weighted and get the Linear Least Squares best fit
amp_end=total(yuse*xuse/vuse)/total(xuse*xuse/vuse)
amp(n-1)=amp_end
if doplot then begin
    rend=[r(wuse(0)),1e5]
    oplot,rend,amp_end*rend^(-endslope),color=c2i('magenta')
endif

;now we know the local PL representation of DS including endpoint 
;extrapolation model. Perfrom the numerical
;integrations for each bin and add them up to get drho.
;Use Johnston et al. (2004) formulas.

for j=0, n-2 do begin
    for i=j, n-2 do begin
                                ;numerical integrationn is in here
        int=xi_pl_int(r(i),r(i+1),r(j),slope(i)+1)
        if not keyword_set(corr2d) then begin
                                ;treat this as DeltaSigma, the regular case
            int2=xi_pl_int(r(i),r(i+1),r(j),0.0+1)
            drho_1=amp(i)*(2.0-slope(i))*int
            drho_2=-off(i)*(2.0-0.0)*int2
            drho(j)=drho(j)+drho_1+drho_2
        endif else begin
                                ;treat this as a 2D correlation function
            drho(j)=drho(j)+amp(i)*slope(i)*int
        endelse
        
                                ;print,drho_1,drho_2
    endfor
    
    if not (keyword_set(noendc)) then begin
        rrmax=nrmax*rmax
                                ;do numerical in for most of it
        int=xi_pl_int(rmax,rrmax,r(j),slope(n-1)+1)
        if not keyword_set(corr2d) then begin
            endc(j)=amp_end*(2.0-slope(n-1))*int
        endif else begin
            endc(j)=amp_end*slope(n-1)*int
        endelse
        
                                ;but last bit to infinity (tiny
                                ;correction), do analytically
                                ;(approximate)
        if not keyword_set(corr2d) then begin
            endc(j)=endc(j)+amp_end*(2.0-slope(n-1))/!pi * rrmax^(-(slope(n-1)+1))/(slope(n-1)+1)
        endif else begin
            endc(j)=endc(j)+amp_end*slope(n-1)/!pi * rrmax^(-(slope(n-1)+1))/(slope(n-1)+1)
        endelse
    endif
endfor

;actually apply correction
drho=drho+endc                      

if doplot then begin
     yr=[(0.5*min(drho))>0.01,max(drho)*2]
     pplot,r(0:n-2),(drho-endc)>1e-6,/xlog,/ylog,yr=yr,/yst,psym=1,$
       xtickf='loglabels', ytickf='loglabels', $
       aspect=1, $
       xtit="r",ytit="Drho",xr=xr,/xst,tit='Output Drho'
    oplot,r(0:n-2),drho,psym=1,color=c2i('magenta')
endif

if keyword_set(nocov) or n_elements(dscov) eq 0 then begin
    !p.multi=pmold
    return
endif
;or else continue and do the covariance matrix be recursive calls
;to evaluate numerical derivative matrix the Jacobian.

mat=fltarr(n,n-1)   ;the jacobian
dsp=ds
dsm=ds

;fill in the jacobian by symmetric numerical
;differentiation. Calls this program recursively. Simplest
;method to code but not fastest. 

FOR i=0, n-1 DO BEGIN
    ;a small number
    h=abs(0.005*(abs(ds(i))+abs(dserr(i))))
    ;offset these
    dsp(i)=ds(i)+h
    dsm(i)=ds(i)-h
    self->ds2drho,r,dsp,drhop,/nocov,/noplot,dscov=dscov,nsigo=nsigo,endslope=endslope,corr2d=corr2d
    self->ds2drho,r,dsm,drhom,/nocov,/noplot,dscov=dscov,nsigo=nsigo,endslope=endslope,corr2d=corr2d
    ;the derivative as a difference.
    der=(drhop-drhom)/(2*h)
    mat(i,*)=der
    ;puts back the offset values to originals
    dsp(i)=ds(i)
    dsm(i)=ds(i)
ENDFOR

;propogate covariance matrices

drhocov=mat##(dscov##transpose(mat))
drhoerr=drhocov(lindgen(n-1)*n)
wbad=where(drhoerr le 0,nbad)
if nbad gt 0 then begin
    print,"Warning drhocov has negitive diagonals"
    drhoerr(*)=-999.0
    !p.multi=pmold
    return
endif
drhoerr=float(sqrt(drhoerr))

if doplot then begin
    errplot,r,drho+drhoerr,drho-drhoerr
    corrmat=drhocov*0.0
    corrmatds=dscov*0.0
    for i=0, n-2 do for j=0, n-2 do corrmat(i,j)=drhocov(i,j)/sqrt(drhocov(i,i)*drhocov(j,j))
    for i=0, n-1 do for j=0, n-1 do corrmatds(i,j)=dscov(i,j)/sqrt(dscov(i,i)*dscov(j,j))
    loadct,0,/silent
    tvim2,corrmatds,tit='DS Correlation Matrix',r=[-1,1],/scale
    tvim2,corrmat,tit='Drho Correlation Matrix',r=[-1,1],/scale
    simpctable
endif

!p.multi=pmold
!p.charsize=pcharold

return
end



function deproject::interp_hyb,x0,x1,y0,y1,e0,e1,x,nsig=nsig,lin=lin,pl=pl,amp=amp,slope=slope,off=off
;this functions takes an interval x0,x1 and fuction vales at those
;points y0 y1 and errors at those points e1 e2 and figures out an
;appropriate interpolation for points x. The interpolation model can
;be specified with keywords as either, linear, powerlaw or hybrid
;(default)
;The hybrid interpolation is given by y=Amp* x^(-slope) -off. This has
;better stability than trying to fit powerlaw to points close to
;zero or even negative points (impossible). The offset is determined as the amount
;to add to y to make the point nsig*sigma above zero.

if n_elements(x) eq 0 then x=1.0
;just a dummy for when you just want the slopes and offsets and amp

if keyword_set(lin) then begin
    dx=x1-x0
    y=(y1*(x-x0) + y0*(x1-x))/dx
    slope=1.0
    amp=(y1-y0)/dx
    off=(amp*x0-y0)
    return,y
endif

if keyword_set(pl) then begin
    dx=alog(x1)-alog(x0)
    y=(alog(y1)*(alog(x)-alog(x0)) + alog(y0)*(alog(x1)-alog(x)))/dx
    slope=(alog(y1)-alog(y0))/dx
    amp=exp((alog(x1)*alog(y0)-alog(x0)*alog(y1))/dx)
    off=0.0
    return,exp(y)
endif

if n_elements(nsig) eq 0 then nsig=1
off=(-min([y0-nsig*e0,y1-nsig*e1])) > 0.0

dx=alog(x1)-alog(x0)
y=(alog(y1+off)*(alog(x)-alog(x0)) + alog(y0+off)*(alog(x1)-alog(x)))/dx
slope=(alog(y1+off)-alog(y0+off))/dx
amp=exp((alog(x1)*alog(y0+off)-alog(x0)*alog(y1+off))/dx)
return,exp(y)-off
end






;docstart::deproject::ds_drho2mass
; NAME: deproject::ds_drho2mass 
;  
; PURPOSE:
;  Convert DeltaSigma and its derived drho to aperture mass , using both
;  mass formulas.
;
; CALLING SEQUENCE:
; ds_drho2mass,r,ds,dserr,drho,drhoerr,mass_in,mass_out,yr=yr,noplot=noplot
;
; INPUTS:
; r  - array of binned radii, supposed to be the radius such that ds=DeltaSigma(r)  
;     may need to make corrections for effects of binning. Should be
;     in Mpc/h for numerical reasons.
; ds - DeltaSigma
; dserr - error on delta sigma
; drho - the 3D over-density from program ds2drho  
; drhoerr -  error on drho used for weighting 
; 
;
; OPTIONAL INPUTS:
; yr - yrange used for plotting 
; noplot -  if set won't plot
; nsigo - the number of sigma above 0 that the data must be to be
;         interpolated without an offset. Default 1.
; endslope - The negative endpoint DeltaSigma logarithmic
;            slope. Default 1. Used here as the drho endpoint slope
;            which is -(endslope+1).
;
; OUTPUTS:
; mass_in  - The inside mass formula 
; mass_out - The outside mass formula
; OPTIONAL KEYWORD OUTPUTS
;
; COMMENTS:
; For error estimates use the program ds2drho_mass_hyb, which calls this one.
;  
; EXAMPLES:
; IDL> file="Erins-absolutely-rediculously-long-filename.fit"
; IDL> str=mrdfits(file,1)
; IDL> r=str.meanr/1000.0 ; put in Mpc/h 
; IDL> ds=str.sigma
; IDL> dscov=str.covariance
; IDL> ds2drho,r,ds,drho,drhocov,drhoerr,endc,dscov=dscov
; IDL> ds_drho2mass,r,ds,dserr,drho,drhoerr,mass_in,mass_out
; METHOD:
; Uses a method based on powerlaw interpolation between points when
; drho > nsigo*drhoerr. For low S/N data or negative data an offset is
; added to make this condition true before determining the 2 powerlaw
; paramters. Then the offset is subtracted off. So the resulting interpolation
; is powerlaw plus constant for low S/N data. For data that looks
; locally like a powerlaw this is unbiased for high S/N data and
; mostly unbiased for low S/N data. For low S/N data the interpolation
; is closer to linear interpolation.
; PROCEDURES CALLED:
; mass_pl_int
;  
; REVISION HISTORY:
;   20-Sept-2004  Written by David Johnston, Princeton.
;   06-July-2005  Converted from ds_xi2mass.pro David Johnston
;docend::deproject::ds_drho2mass

pro deproject::ds_drho2mass,r,ds,dserr,drho,drhoerr,mass_in,mass_out,yr=yr,noplot=noplot,nsigo=nsigo,endslope=endslope

if n_params() eq 0 then begin
    print,'-syntax dep->ds_drho2mass,r,ds,dserr,drho,drhoerr,mass_in,mass_out,yr=yr,'
    print,'noplot=noplot,nsigo=nsigo,endslope=endslope'
    return
endif

n=n_elements(r)
ni= n_elements(drho)
pr = r[0:ni-1]
if n_elements(nsigo) eq 0 then nsigo=1.0
;the number of sigma to offset before fitting the powerlaw
if n_elements(endslope) eq 0 then endslope=1.0
;the endpoint slope of delta sigma for the endpoint correction
if n_elements(endslope_npoints) eq 0 then endslope_npoints=3
;the number of points at the end to use to fit the endpoint amplitude

rmax=max(r)
rmin=min(r)
xr=[(.5*rmin) > .001, rmax]

doplot=1
if keyword_set(noplot) then doplot=0

if n_elements(yr) eq 0 then yr=[0.8*min(drho-drhoerr)>.01,max(drho+drhoerr)*1.2]

pmold = !p.multi
pcharold = !p.charsize
if doplot then begin
    !p.multi=[0,2,2]
    !p.charsize = 1
    pplot,pr,drho,/xlog,/ylog,psym=1,yr=yr,/yst,$
      xtickf='loglabels', ytickf='loglabels', $
      aspect=1, $
      xtit="R",ytit="DRHO",tit="Input Drho",xr=xr,/xst
    errplot,r,drho-drhoerr,drho+drhoerr
    nip=100
    ipoints=fltarr(nip,n-1)
                                ;log interpolated points for plotting
    for i=0, n-3 do begin
        imin=alog(r[i])
        irange=alog(r[i+1])-alog(r[i])
        ipoints[*,i]=exp(findgen(nip)*irange/(nip-1)+imin)
        ;oplot,ipoints[*,i],self->interp_hyb(r[i],r[i+1],ds[i],ds[i+1],dserr[i],dserr[i+1],ipoints[*,i],/lin),$
        ;  color=!red,psym=3
        ;oplot,ipoints[*,i],self->interp_hyb(r[i],r[i+1],ds[i]>.1,ds[i+1]>.01,dserr[i],dserr[i+1],ipoints[*,i],/pl),$
        ;  color=!blue,psym=3
        oplot,ipoints[*,i],self->interp_hyb(r[i],r[i+1],drho[i],drho[i+1],drhoerr[i],drhoerr[i+1],ipoints[*,i],nsig=nsigo),$
          color=c2i('green'),psym=3
    endfor
endif


amp=fltarr(n)
slope=amp
off=amp

;get slopes and amps and offsets


for i=0, n-3 do begin
    junk=self->interp_hyb(r[i],r[i+1],drho[i],drho[i+1],drhoerr[i],drhoerr[i+1],xjunk,sl=sli,amp=ampi,off=offi,nsig=nsigo)
    ;print,sli,ampi,offi,drho(i)
    slope(i)=-sli
    amp(i)=ampi
    off(i)=offi
endfor

slope(n-1)=endslope+1          
slope(n-2)=endslope+1     

;get the last few points for the endpoint amp determination
;do Least Squares Fit for amplitude.
wuse=indgen(endslope_npoints)+n-endslope_npoints-1
xuse=r(wuse)^(-(endslope+1))
yuse=drho(wuse)
vuse=drhoerr(wuse)^2
amp_end=total(yuse*xuse/vuse)/total(xuse*xuse/vuse)
amp(n-2)=amp_end
amp(n-1)=amp_end
if doplot then begin
    rend=[r(wuse(0)),1e5]
    oplot,rend,amp_end*rend^(-(endslope+1)),color=c2i('magenta')
endif

mass_in=fltarr(n-1)
mass_out=mass_in
imin=0
Rs=r(imin)
;imax=n-1
mass_in_c=fltarr(n-1)
mass_out_c=fltarr(n-1)

;do numerical integrations for each interval and add them up
;using Johnston et al. (2004) Formulas.

for j=1, n-2 do begin
    for i=0, j-1 do begin
        x0=sqrt((r(i)/Rs)^2-1)
        x1=sqrt((r(i+1)/Rs)^2-1)
        int=mass_pl_int(X0,X1,slope(i),/in)
        ;int2_num=mass_pl_int(X0,X1,0.0,/in)
        ;second intergral can be done analytically 
        int2=(x1-x0)+ 2*(x1^3-x0^3)/3.0
        amplitude=amp(i)
        mass_in(j)=mass_in(j)+2*amplitude*int*Rs^(-slope(i)) - 2*off(i)*int2
    endfor
endfor

imax=n-3
for j=0, n-2 do begin
    for i=j, imax do begin
        x0=sqrt((r(i)/r(j))^2-1)
        x1=sqrt((r(i+1)/r(j))^2-1)
        int=mass_pl_int(X0,X1,slope(i),/out)
        int2=mass_pl_int(X0,X1,0,/out)
        x0_c=sqrt((r(i)/Rs)^2-1)
        x1_c=sqrt((r(i+1)/Rs)^2-1)
        int_c=mass_pl_int(X0_c,X1_c,slope(i),/out)
        ;int_c2_num=mass_pl_int(X0_c,X1_c,0,/out)
        ;do it analytically
        int_c2=(x1_c-x0_c)+2*(x1_c^3-x0_c^3)/3.0 - 2*((1+x1_c^2)^(3/2.0) - (1+x0_c^2)^(3/2.0))/3.0
        mass_out(j)=mass_out(j)+2*amp(i)*int*r(j)^(-slope(i)) - 2*off(i)*int2
        mass_in_c(j)=mass_in_c(j)+2*amp(i)*int_c*Rs^(-slope(i)) - 2*off(i)*int_c2
    endfor
    i=n-2
    X0=sqrt((r(i)/r(j))^2-1)
    ;do the xmax to infinity correction
    intoc=mass_pl_int(X0,(X0+1)*20d,slope(i),/out)
    mass_out_c(j)=2*amp(i)*intoc*r(j)^(-slope(i))
    ;also for mass_in_c
    X0_c=sqrt((r(i)/Rs)^2-1)
    int_c=mass_pl_int(X0_c,(X1_c+1)*20,slope(i),/out)
    mass_in_c(j)=mass_in_c(j)+2*amp(i)*int_c*Rs^(-slope(i))
endfor

mass_in=(mass_in+ds(imin)/Rs)*!pi*Rs^3
mass_out=(mass_out+mass_out_c+ds/r)*!pi*r^3
mass_in_c=mass_in_c*!pi*Rs^3
mass_in=mass_in+mass_in_c

if doplot then begin
    yr=[0.8*min(mass_in>mass_out)> .1,max(mass_in<mass_out)*1.2]
    pplot,pr,mass_in,/xlog,/ylog,$
      xtickf='loglabels', ytickf='loglabels', $
      aspect=1, $
      xtit="r",ytit="mass_in (r)",psym=1,xr=xr,/xst,yr=yr,/yst
    pplot,pr,mass_in,/xlog,/ylog,$
      xtickf='loglabels', ytickf='loglabels', $
      aspect=1, $
      xtit="r",ytit="mass_out (r)",psym=1,xr=xr,/xst,yr=yr,/yst
    oplot,pr,mass_out,psym=7,color=c2i('green')
    pplot,pr,(mass_in+mass_out)/2.0,/xlog,/ylog,xtit="r",$
      xtickf='loglabels', ytickf='loglabels', $
      aspect=1, $
      ytit="mass (r)",psym=1,xr=xr,/xst,yr=yr,/yst,tit='Average'
    ;;print,"mean mass_in",mean(mass_in)
    ;;print,"mean mass_out",mean(mass_out)
endif


!p.multi = pmold
!p.charsize = pcharold

return
end












;docstart::deproject::drho2mass
; NAME: drho2mass 
;  
; PURPOSE:
;  Convert drho to aperture mass.  This is called by ::ds2drho_mass, but
;  does not propogate the covariance matrix, so you should use that program
;  with the /corr2d keyword.
;
; CALLING SEQUENCE:
; drho2mass,r,drho,drhoerr,mass,yr=yr,noplot=noplot
;
; INPUTS:
; r  - array of binned radii, supposed to be the radius such that ds=DeltaSigma(r)  
;     may need to make corrections for effects of binning. Should be
;     in Mpc/h for numerical reasons.
; drho - the 3D over-density from program ds2drho  
; drhoerr -  error on drho used for weighting 
; 
; OPTIONAL INPUTS:
; yr - yrange used for plotting 
; noplot -  if set won't plot
; nsigo - the number of sigma above 0 that the data must be to be
;         interpolated without an offset. Default 1.
; OUTPUTS:
; mass - the aperture mass
;
; OPTIONAL KEYWORD OUTPUTS
;
; COMMENTS:
;  This program doesn't use DeltaSigma and so doesn't include the
;  inner mass, inside the first data point. The resultant mass[0] is 
;  always zero.
;
; METHOD:
; Uses a method based on powerlaw interpolation between points when
; drho > nsigo*drhoerr. For low S/N data or negative data an offset is
; added to make this condition true before determining the 2 powerlaw
; paramters. Then the offset is subtracted off. So the resulting interpolation
; is powerlaw plus constant for low S/N data. For data that looks
; locally like a powerlaw this is unbiased for high S/N data and
; mostly unbiased for low S/N data. For low S/N data the interpolation
; is closer to linear interpolation.
; PROCEDURES CALLED:
; mass_pl_int
;  
; REVISION HISTORY:
;   20-Sept-2004  Written by David Johnston, Princeton.
;   11-Sept-2006  Converted from ds_drho2mass.pro David Johnston Caltech
;docend::deproject::drho2mass

pro deproject::drho2mass,r,drho,drhoerr,mass,yr=yr,noplot=noplot,nsigo=nsigo,debug=debug

if n_params() eq 0 then begin
    print,'-syntax dep->drho2mass,r,drho,drhoerr,mass,yr=yr,'
    print,'noplot=noplot,nsigo=nsigo'
    return
endif

;n=n_elements(r)
n=n_elements(drho)
pr = r[0:n-1]

if n_elements(nsigo) eq 0 then nsigo=1.0
;the number of sigma to offset before fitting the powerlaw

rmax=max(r)
rmin=min(r)
xr=[(.5*rmin) > .001, rmax]

doplot=1
if keyword_set(noplot) then doplot=0

if n_elements(yr) eq 0 then yr=[0.8*min(drho-drhoerr)>.01,max(drho+drhoerr)*1.2]

if doplot then begin
    !p.multi=[0,2,2]
    pplot,pr,drho,/xlog,/ylog,psym=1,yr=yr,/yst,$
      xtickf='loglabels', ytickf='loglabels', $
      aspect=1, $
      xtit="R",ytit="DRHO",tit="Input Drho",xr=xr,/xst
    errplot,pr,drho-drhoerr,drho+drhoerr
    nip=100
    ipoints=fltarr(nip,n-1)
                                ;log interpolated points for plotting
    for i=0, n-2 do begin
        imin=alog(r[i])
        irange=alog(r[i+1])-alog(r[i])
        ipoints[*,i]=exp(findgen(nip)*irange/(nip-1)+imin)
        ;oplot,ipoints[*,i],self->interp_hyb(r[i],r[i+1],ds[i],ds[i+1],dserr[i],dserr[i+1],ipoints[*,i],/lin),$
        ;  color=!red,psym=3
        ;oplot,ipoints[*,i],self->interp_hyb(r[i],r[i+1],ds[i]>.1,ds[i+1]>.01,dserr[i],dserr[i+1],ipoints[*,i],/pl),$
        ;  color=!blue,psym=3
        oplot,ipoints[*,i],self->interp_hyb(r[i],r[i+1],drho[i],drho[i+1],drhoerr[i],drhoerr[i+1],ipoints[*,i],nsig=nsigo),color=!green,psym=3
    endfor
endif

amp=fltarr(n)
slope=amp
off=amp

;get slopes and amps and offsets

for i=0, n-2 do begin
    junk=self->interp_hyb(r[i],r[i+1],drho[i],drho[i+1],drhoerr[i],drhoerr[i+1],xjunk,sl=sli,amp=ampi,off=offi,nsig=nsigo)
    if finite(sli) eq 0 then begin
        print,'ERROR - slope not finite'
        print,sli,ampi,offi,drho(i)
        stop
    endif
    slope(i)=-sli
    amp(i)=ampi
    off(i)=offi
endfor
mass=fltarr(n-1)


;do numerical integrations for each interval and add them up

al=2.0-slope
epmin=0.001                     ;powerlaw slopes too close to -1 are logarithms

mi=fltarr(n)
for i=0, n-2 do begin
    ep=al[i]+1.0
    if abs(ep) lt epmin then begin
        ;first order approximation goes smoothly into ep=0
        mi[i]=Amp[i]*alog(r[i+1]/r[i])*(1+ep*alog(r[i]*r[i+1]))
    endif else begin
        mi[i]=Amp[i]/(1+al[i])*(r[i+1]^(1+al[i])-r[i]^(1+al[i]))
    endelse
    mi[i]=(mi[i]-(off[i]/3.0)*(r[i+1]^3-r[i]^3))*4*!pi
endfor

mass=total(shift(mi,1),/cum)

if doplot then begin
    yr=[0.8*min(mass)> .1,max(mass)*1.2]
    pplot,pr,mass,/xlog,/ylog,$
      xtickf='loglabels', ytickf='loglabels', $
      aspect=1, $
      xtit="r",ytit="mass (r)",psym=1,xr=xr,/xst,yr=yr,/yst
    ;;print,"mean mass",mean(mass)
endif

if keyword_set(debug) then stop

return
end






PRO deproject__define
  struct = { $
             deproject, $
             dummy: 0 $
           }
END 
