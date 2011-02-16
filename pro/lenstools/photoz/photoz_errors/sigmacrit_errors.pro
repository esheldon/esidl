PRO sigmacrit_errors, zL, norm, pow, dops=dops

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: sigmacrit_errors, zL [, norm, pow, dops=dops]'
      return
  ENDIF 

  tt=systime(1)

  dir='/net/cheops1/data0/esheldon/pzerr2scriterr/'
  fitfile = dir + 'sigerr_zL'+ntostr(zL,5)+'.fit'
  IF keyword_set(dops) THEN BEGIN 
      psfile = dir + 'sigerr_zL'+ntostr(zL,5)+'.ps'
      begplot,name=psfile,/color
  ENDIF 

  ;; see how gaussian in zs looks in sigma_crit and 1./sigma_crit
  ;; this is done for a fixed zL. A power law is fit to the
  ;; dependence of the ratio zerr/sigma_crit_err as a function
  ;; of source redshift
  ;; (zerr/sigmacrit_err) = norm(zL)*zsource^( -pow(zL) )

  ;; fixed lens redshift

  nzmean = 15
  minzs = 0.25 > (zL+0.01)
  maxzs = 0.8
  zmeanarr = arrscl( findgen(nzmean), minzs, maxzs)

  nzerr = 10
  zerrarr = arrscl( findgen(nzerr), 0.01, 0.1 )
  sigfracerrarr = fltarr(nzmean, nzerr)
  jsigfracerrarr = sigfracerrarr
  jsigfracerrarr_err = sigfracerrarr

  meanrat = fltarr(nzmean)
  meanraterr = fltarr(nzmean)

  ;; # of zs for genrand
  nzs = 100000L

  fitamp = fltarr(nzmean, nzerr)
  fitcen = fltarr(nzmean, nzerr)
  fitalpha = fltarr(nzmean, nzerr)
  fitsig = fltarr(nzmean, nzerr)

  FOR j=0L, nzmean-1 DO BEGIN 

      ;; source redshift distribution
      zmean = zmeanarr[j]

      ;; Mean DLs, Ds, DL
      DL = angdist_lambda(zL)
      meanDs = angdist_lambda(zmean)
      meanDLs = angdist_lambda(zmean, zL)
            
      ninsig = 20.

      FOR i=0L, nzerr-1 DO BEGIN 
          
          zerr = zerrarr[i]
          zbin = 2.*zerr/ninsig
          zs = randomu(seed, nzs, /normal)*zerr + zmean

          zmin = zmean - 3.*zerr
          zmax = zmean + 3.*zerr

          sigcrit = fltarr(nzs)

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; plot up the source redshift distribution
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          !p.multi=[0,0,2]
          xrange = [zmin,zmax]
          xrange[0] = xrange[0] < 0
          xrange[1] = xrange[1] > 1
          plothist, zs, bin=zbin, /norm, xrange=xrange, $
                    ytitle='P(z!Ds!N)', xtitle='Z!Ds!N'
          ;;oplot,z,zgauss,psym=10,color=!red
          legend,['z!Ds!N = '+ntostr(zmean,4,/round),$
                  !csym.sigma+'(z!Ds!N) = '+ntostr(zerr,5,/round)]

          ;; now calculate sigma_crit, for these
          fac1 = .3474          ;g/cm^2
          fac2 = 1.663e3        ;Msolar/parsec^2
          
          mean_sigcrit_inv = sigmacritinv(zl, zmean)*1.e4
          sigcrit_inv = sigmacritinv(zs, zL)*1.e4
          sigcrit_inv = sigcrit_inv[where(finite(sigcrit_inv))] 

          wzsh = where(sigcrit_inv GT 1.e-3)
          
          ;; there are many outliers
          sigma_clip, sigcrit_inv, meansigclip, sigerrclip, nsig=4, niter=10, /silent
          ;;w=where(sigcrit_inv LT meansigclip + 10.*sigerrclip)
          meansiginv = mean(sigcrit_inv)
          siginverr = sdev(sigcrit)

          ;; plot up sigma_crit distribution
          minusesig = mean_sigcrit_inv - 4.*sigerrclip
          maxusesig = mean_sigcrit_inv + 6.*sigerrclip
          xtitle = !csym.sigma_cap+'!Dcrit!N!U-1!N [10!U-4!N pc!U2_N / M'+sunsymbol()+$
            ' pc!U-2!N]'
          ytitle = 'P('+!csym.sigma_cap+'!Dcrit!N!U-1!N)'
          sbin = 2.*sigerrclip/ninsig
          plothist, sigcrit_inv, xhist, yhist, bin=sbin, /norm, $
                    min=minusesig, max=maxusesig,$
                    title = 'z!DL!N = '+ntostr(zL,5), $
                    xtitle=xtitle, $
                    ytitle=ytitle
          
          oplot,[mean_sigcrit_inv,mean_sigcrit_inv],[0,1000],$
                color=!green,thick=5
          oplot,[meansiginv,meansiginv],[0,1000],color=!red
          medsc=median(sigcrit_inv)

          mean_sigmacritinv, zL, zmean, zerr, 200, m2, /silent
          m2 = m2*1.e4

          ;;m2 = total(zs[wzsh]*sigcrit[wzsh])/total(zs[wzsh])
          oplot, [m2, m2], [0, 1000], color=!orange

          oplot,[medsc,medsc],[0,1000],color=!blue
          legend,['From mean z!Ds!N',$
                  'Mean of distribution', $
                  'Median'],colors=[!green,!red, !blue],line=[0,0,0],/right,$
                 charsize=0.7

          ;; oplot 1-sigma range
          oplot,[mean_sigcrit_inv + siginverr, mean_sigcrit_inv + siginverr],[0,1000],color=!magenta
          oplot,[mean_sigcrit_inv - siginverr, mean_sigcrit_inv - siginverr],[0,1000],color=!magenta
          
;          print,meansig,!plusminus,sigerr
;          print,'Zmean = ',zmean,' Zerr = ',zerr,' Sigerr/1.e4 = ',sigerr
;          print,'(sigerr/1.e4)/zerr = ',sigerr/zerr
          
          sigfracerrarr[j,i] = siginverr/mean_sigcrit_inv

          w=where(sigcrit_inv LT maxusesig AND sigcrit_inv GT minusesig, nw)
          cjackknife, sigcrit_inv[w], nw, 1, datamean, dataerr,statistic=2
          jsigfracerrarr[j,i] = datamean[0]/mean_sigcrit_inv
          jsigfracerrarr_err[j,i] = dataerr[0]/mean_sigcrit_inv
;          print,datamean,dataerr
;          print
          
;          IF j EQ 0 AND i EQ 0 THEN BEGIN 
;              tamp = 4.7*((zmean/0.25)^0.25)
;              talpha = 3.8*(0.25/zmean)*( (zerr/0.01)^0.25 )
;              tsig = 0.54*( (0.25/zmean)^0.5 )*( (zerr/0.01)^0.25  )
;          ENDIF ELSE BEGIN 
              
;              IF tamp = fitamp[j-1, i]
              

;              IF i EQ 0 THEN BEGIN 
;                  tamp = fitamp[j-1,i]
;          ENDELSE 
;          aguess = [tamp, talpha, tsig]

;          wxhist = where(xhist GT 0.2)
;          delvarx,afit
;          afit = comfit2(xhist[wxhist],yhist[wxhist],aguess,$
;                         /lognormal,/noderivative,$
;                         yfit=yfit,iter=iter)

;          fitamp[j,i] = afit[0]
;          fitalpha[j,i] = abs(afit[1])
;          fitsig[j,i] = abs(afit[2])

;          fitmean = fitalpha[j,i] + exp(0.5*fitsig[j,i]^2)
;          fiterr  = sqrt( (exp(2.*fitsig[j,i]^2) - exp(fitsig[j,i]^2)) > 0.0 )

;          print,'amp =   ',fitamp[j,i]
;          print,'alpha = ',fitalpha[j,i]
;          print,'sig =   ',fitsig[j,i]

;          oplot, xhist[wxhist], yfit, color=!green
;          oplot, [fitmean,fitmean], [0, 10000.], color=!blue
;          oplot, [fitmean+fiterr,fitmean+fiterr], [0, 10000.], color=!blue,line=2
;          oplot, [fitmean-fiterr,fitmean-fiterr], [0, 10000.], color=!blue,line=2

          !p.multi=0
key=get_kbrd(1)
          
      ENDFOR 

      print
      print,'ZL = ',zL
      print,'Zs = ',zmean
      print,'       zserr         amp        alpha        sig'
      forprint,zerrarr,fitamp[j,*],fitalpha[j,*],fitsig[j,*]
;key=get_kbrd(1)

  ENDFOR 

stop

  !p.multi=[0,0,2]
  xrange=[0.9*min(zerrarr), 1.1*max(zerrarr)]
  FOR j=0L, nzmean-1 DO BEGIN 

      IF j EQ 0 THEN BEGIN 
          xtitle = !csym.sigma+'(z!Ds!N)'
          ytitle = !csym.sigma+'('+!csym.sigma_cap+'!Dcrit!N!U-1!N) / '+$
            !csym.sigma_cap+'!Dcrit!N!U-1!N'
          plot,zerrarr,sigfracerrarr[j,*],psym=4, $
               xrange=xrange, $
               title='z!DL!N = '+ntostr(zL,5),$
               xtitle=xtitle, ytitle=ytitle,charsize=1.1
          oploterror,zerrarr,jsigfracerrarr[j,*],jsigfracerrarr_err[j,*], psym=1, color=!red
      ENDIF ELSE BEGIN 
          oplot,zerrarr,sigfracerrarr[j,*],psym=4
          oploterror,zerrarr,jsigfracerrarr[j,*],jsigfracerrarr_err[j,*], psym=1, color=!red
      ENDELSE 
  ENDFOR 

  legend,'zs = '+ntostr(zmeanarr,5),charsize=0.7

  FOR j=0L, nzmean-1 DO BEGIN 

      IF j EQ 0 THEN BEGIN 
          xtitle = !csym.sigma+'(z!Ds!N)'
          ytitle = '[ '+!csym.sigma+'('+!csym.sigma_cap+'!Dcrit!N!U-1!N) / '+$
            !csym.sigma_cap+'!Dcrit!N!U-1!N ] / '+!csym.sigma+'(z!Ds!N)'
          plot,zerrarr,sigfracerrarr[j,*]/zerrarr,psym=4, xrange=xrange,$
               xtitle=xtitle,ytitle=ytitle,charsize=1.1
          rat = jsigfracerrarr[j,*]/zerrarr
          raterr = jsigfracerrarr_err[j,*]/zerrarr
          oploterror,zerrarr,rat,raterr, psym=1, color=!red
      ENDIF ELSE BEGIN 
          oplot,zerrarr,sigfracerrarr[j,*]/zerrarr,psym=4
          rat = jsigfracerrarr[j,*]/zerrarr
          raterr = jsigfracerrarr_err[j,*]/zerrarr
          oploterror,zerrarr,rat,raterr, psym=1, color=!red
      ENDELSE 
      meanrat[j] = mean(rat)
      meanraterr[j] = sdev(rat)
      print,'Mean rat = ',meanrat[j]

      oplot,zerrarr,replicate(meanrat[j], nzerr)

  ENDFOR 

  legend,'zs = '+ntostr(zmeanarr,5) + ' rat = '+ntostr(meanrat,5),charsize=0.7

  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  !p.multi=0


  str=create_struct('zerr', zerrarr, $
                    'sigfracerr', sigfracerrarr, $
                    'jsigfracerr', jsigfracerrarr, $
                    'jsigfracerr_err', jsigfracerrarr_err, $
                    'zsmean', zmeanarr, $
                    'meanrat', meanrat, $
                    'meanraterr', meanraterr, $
                    'norm', 0.0, $
                    'pow', 0.0)


  fitpower,str.zsmean,str.meanrat,str.meanraterr,[.1,-1],yfit,aa
  norm = aa[0]
  pow = abs(aa[1])
  str.norm = norm
  str.pow = pow

  yr=prange(str.meanrat,str.meanraterr, /slack)
  xr = [0.9*minzs, 1.1*maxzs]
  xtitle = 'mean z!Ds!N'
  ytitle = 'mean [ '+!csym.sigma+'('+!csym.sigma_cap+'!Dcrit!N!U-1!N) / '+$
            !csym.sigma_cap+'!Dcrit!N!U-1!N ] / '+!csym.sigma+'(z!Ds!N)'
  aploterror,!gratio,str.zsmean,str.meanrat,str.meanraterr,$
             /xlog,/ylog,psym=1,$
             yrange=yr,xrange=xr,xsty=1,ysty=1,$
             xtitle=xtitle,ytitle=ytitle,charsize=1.1

  add_labels,xtickv=[0.2,0.3,0.4,0.5,0.6,0.7,0.8],ytickv=[0.2,0.4,0.7,2,3,4],$
             charsize=1

  func = str.norm*str.zsmean^(-str.pow)
  oplot, str.zsmean, func, color=!red

  legend,['Norm = '+ntostr(str.norm),'Power = '+ntostr(str.pow)],/right,charsize=1

  IF keyword_set(dops) THEN endplot

  mwrfits, str, fitfile, /create

  ptime,systime(1)-tt

END 
