PRO test_max_like_pz2,num,dd,like_d,dmax,del

  ;; input deltasig
  Del=.221
  ;;Del = 10.0

  ;; lens redshift
  zl=.1
  
  ;; mean source redshift
  meanz=.4
  sigz=.05

  ;; mean inverse critical density
  siginv = sigmacritinv(zl, meanz)

  ;; mean shear
  mean_shear=Del*siginv
  
  ;; shear error
  sig_shear=.04
  ;;sig_shear = 0.05

  print,' zl ',zl,' meanz ',meanz,' sigz ',sigz
  print,' siginv ',siginv,' mean_shear ',mean_shear ,' sig_shear ',sig_shear,' Del ',Del

  ;; draw from source redshift distribution
  zobj=randomn(seed,num)*sigz+meanz

  ;; inverse critical density for each object
  siginv_obj=sigmacritinv(zl,zobj)

  ;; shear induced in each object, plus random error
  shear=Del*siginv_obj+randomn(seed,num)*sig_shear

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; don't know what this is for
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;lnorm=1.0/(sqrt(2.0*!pi)*sig_shear)
  ;;lexp=((shear-mean_shear)/sig_shear)^2 <40.0
  ;;l=lnorm*exp(-.5*lexp)

  ;; some values at which to do the redshift
  ;; "integral"

  ntest=200
;  ztest=findgen(ntest)/ntest
;  ztest=2.0*ztest

  minzstest = 0.0
  maxzstest = 1.0

  gauleg, minzstest, maxzstest, ntest, ztest, WWi

  ;;stop
  ;; now the values of critical density at those redshifts
  sigtest=sigmacritinv(zl, ztest)
  u=((ztest-meanz)/sigz)^2 < 40.0
  
  ;; plot up these points and test shears
  !p.multi=[0,1,4]
  !p.charsize = 2
  zgauss=exp(-.5*u)/sqrt(2.*!pi)/sigz
  plot,ztest,zgauss
  oplot,ztest,sigtest*Del

  ;; plotting some stuff

  lnorm=1.0/(sqrt(2.0*!pi)*sig_shear)
  lexp=((mean_shear-Del*sigtest)/sig_shear)^2 < 40.0
  l=lnorm*exp(-.5*lexp)
  plot,ztest,l
  lnorm=1.0/(sqrt(2.0*!pi)*sig_shear)
  lexp=((mean_shear+4*sig_shear-Del*sigtest)/sig_shear)^2 < 40.0
  l=lnorm*exp(-.5*lexp)
  oplot,ztest,l,color=130

  ;; now lets do the likelihood
  nbins=100
  like_d=fltarr(nbins)

  mindd = 0.2
  maxdd = 0.25
  ;;mindd = 8.0
  ;;maxdd = 12.0
  dd=arrscl( findgen(nbins), mindd, maxdd)

  FOR j=0, nbins-1 DO BEGIN
      
      ;; temporary holder for likelihood.  Will sum log like over all
      ;; objects
      like=0.0
      FOR i=0, num -1 DO BEGIN

          fac = 2.
          u=(ztest-zobj[i])^2/sigz^2/fac < 40.0
          zgauss=exp(-.5*u)/sqrt(2.*!pi*fac*sigz^2)

          ;; calculate likelihood at each test source redshift
          lnorm=1.0/(sqrt(2.0*!pi)*sig_shear)
          lexp=((shear(i)-dd[j]*sigtest)/sig_shear)^2 < 40.0
          l=lnorm*exp(-.5*lexp)
          
          ;; This is Dave's "integral" over z
          ;;likej=total(l*zgauss)

          likej = total(l*zgauss*WWi)
        
          ;; print out some stuff
          IF i eq 20 AND 1 THEN BEGIN
              IF j EQ 0 THEN plot,ztest,l,yr=[0,15] ELSE oplot,ztest,l
              junk=''
              
              print,'likej ',likej,' d ',dd[j],' del ',del
          ENDIF

          ;; add up the log likelihood
          like=like + alog10(likej)
          
      ENDFOR

      ;; log likelihood for this density contrast
      like_d(j)=like
  ENDFOR

  ;; plot it up
  like_d=like_d-max(like_d)
  plot,dd,10.0^(like_d),/yno
  oplot,[Del,Del],[-1,1]*1e15,color=!red,linest=1

  ;; plot up max likelihood point
  junk=max(like_d,dmax)
  dmax=dd(dmax)
  oplot,[dmax,dmax],[-1,1]*1e15,color=!green

  ;; don't know

  lnorm=1.0/(sqrt(2.0*!pi)*sig_shear)
  lexp=(((mean_shear-dmax*sigtest)/sig_shear)^2) < 90.0
  print,max(lexp)
  l=lnorm*exp(-.5*lexp)


  print,'dmax ',dmax,' Del ',del
  return
END



