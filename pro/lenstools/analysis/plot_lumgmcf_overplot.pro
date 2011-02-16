PRO plot_lumgmcf_overplot, stripes, bybin=bybin, loglog=loglog

  ;; plot combined gmcf for each luminosity bin
  ;; old mckay et al. stuff

  errfac=1.32

  IF n_elements(stripes) EQ 0 THEN stripes=[10,36,37,42,43,82]
  nstripe = n_elements(stripes)
  stripestr=''
  FOR i=0L, nstripe-1 DO stripestr = stripestr+'stripe'+ntostr(stripes[i])+'_'
  print,'stripestr = ',stripestr

  types=['lum1_','lum2_','lum3_','lum4_']
  names=['lum1 (low)','lum2','lum3','lum4 (high)']
  colors=['u','g','r','i','z']
  nclr=n_elements(colors)

  basedir = '/sdss5/data0/lensout/'
  ;; assumes the combined of all stripes is under "first" stripe
  paramdir = basedir+'stripe'+ntostr(stripes[0])+'/sublum/' + colors +'/'

  ;; this is just to get mean luminosities
  lumbinfile = basedir+'mass2light/lumbin_rad260.fit'
  lumbin=mrdfits(lumbinfile, 1,/silent)
  fend = 'N1.fit'
  message=strarr(2)

  ytitle = '!S'+!tsym.sigma_cap+'!R!A'+!tsym.minus+'!N ('+!tsym.ltequal+'R) '+!tsym.minus+' !S'+$
    !tsym.sigma_cap+'!R!A'+!tsym.minus+$
    '!N (R) (h M!DSun!N pc!U'+!tsym.minus+'2!N)!X'
  xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)'

  IF keyword_set(loglog) THEN BEGIN
      addstr='_log' 
      xlog=1 & ylog=1
  ENDIF ELSE addstr=''

  IF keyword_set(bybin) THEN BEGIN
      psfile = basedir+'mass2light/lumbin_overplot_gmcf'+addstr+'_bybin.ps'
  ENDIF ELSE psfile = basedir+'mass2light/lumbin_overplot_gmcf'+addstr+'.ps'
  begplot,name=psfile,/color

  ;; plot gmcf in bins, bins in each color

  FOR iclr=0L, nclr-1 DO BEGIN 
      
      file=paramdir[iclr]+'lum1_zgal_gal_'+stripestr+'comb_corr_'+fend
      lum1=mrdfits(file,1,/silent)
      file=paramdir[iclr]+'lum2_zgal_gal_'+stripestr+'comb_corr_'+fend
      lum2=mrdfits(file,1,/silent)
      file=paramdir[iclr]+'lum3_zgal_gal_'+stripestr+'comb_corr_'+fend
      lum3=mrdfits(file,1,/silent)
      file=paramdir[iclr]+'lum4_zgal_gal_'+stripestr+'comb_corr_'+fend
      lum4=mrdfits(file,1,/silent)

      file=paramdir[iclr]+'lum1_zgal_gal_'+stripestr+'fitparam_'+fend
      lum1par=mrdfits(file,1,/silent)
      file=paramdir[iclr]+'lum2_zgal_gal_'+stripestr+'fitparam_'+fend
      lum2par=mrdfits(file,1,/silent)
      file=paramdir[iclr]+'lum3_zgal_gal_'+stripestr+'fitparam_'+fend
      lum3par=mrdfits(file,1,/silent)
      file=paramdir[iclr]+'lum4_zgal_gal_'+stripestr+'fitparam_'+fend
      lum4par=mrdfits(file,1,/silent)

      lum1.sigmaerr = lum1.sigmaerr*errfac
      lum2.sigmaerr = lum2.sigmaerr*errfac
      lum3.sigmaerr = lum3.sigmaerr*errfac
      lum4.sigmaerr = lum4.sigmaerr*errfac

      xx=arrscl(findgen(1000), 0.0, 10000.)
      w=lindgen(n_elements(lum1.meanr))
      binmax=2

      yrange=fltarr(2)
      yrange[0] = min([min(lum1.sigma[w])-max(lum1.sigmaerr[w]),$
                       min(lum2.sigma[w])-max(lum2.sigmaerr[w]),$
                       min(lum3.sigma[w])-max(lum3.sigmaerr[w]),$
                       min(lum4.sigma[w])-max(lum4.sigmaerr[w])])
      yrange[1] = max([max(lum1.sigma[w])+max(lum1.sigmaerr[w]),$
                       max(lum2.sigma[w])+max(lum2.sigmaerr[w]),$
                       max(lum3.sigma[w])+max(lum3.sigmaerr[w]),$
                       max(lum4.sigma[w])+max(lum4.sigmaerr[w])])
      IF keyword_set(loglog) THEN BEGIN
          yrange[0] = .1
          yrange[1] = yrange[1]*1.1
          xrange=[50,1100]
          xstyle=1
          ystyle=1
      ENDIF 

      color=!red
      aploterror,!gratio,lum1.meanr[w],lum1.sigma[w],lum1.sigmaerr[w],$
        psym=1,yrange=yrange,xrange=xrange,xlog=xlog,ylog=ylog,$
        xstyle=xstyle,ystyle=ystyle,xtitle=xtitle,ytitle=ytitle
      oploterror,lum1.meanr[w],lum1.sigma[w],lum1.sigmaerr[w],$
        color=color,errcolor=color,psym=1
      IF keyword_set(bybin) THEN sigmav = sqrt(lum1par.sissigma2[binmax]) $
      ELSE sigmav=lum1par.sigmav
      sigmaverr = sqrt(lum1par.sissigmaerr2[binmax])
      yfit=sigmasis(sigmav, xx,/core)
      print,sigmav,sqrt(lum1par.sissigma2[binmax])
      oplot,xx,yfit,color=color
      oplot,[0,10000],[0,0]      


      color=!magenta
      oploterror,lum2.meanr[w],lum2.sigma[w],lum2.sigmaerr[w],$
        psym=2,color=color,errcolor=color
      IF keyword_set(bybin) THEN sigmav = sqrt(lum2par.sissigma2[binmax]) $
      ELSE sigmav=lum2par.sigmav
      sigmaverr = sqrt(lum2par.sissigmaerr2[binmax])
      yfit=sigmasis(sigmav, xx,/core)
      print,sigmav,sqrt(lum2par.sissigma2[binmax])
      oplot,xx,yfit,color=color
      oplot,[0,10000],[0,0]


      color=!green
      oploterror,lum3.meanr[w],lum3.sigma[w],lum3.sigmaerr[w],$
        psym=4,color=color,errcolor=color
      IF keyword_set(bybin) THEN sigmav = sqrt(lum3par.sissigma2[binmax]) $
      ELSE sigmav=lum3par.sigmav
      sigmaverr = sqrt(lum3par.sissigmaerr2[binmax])
      yfit=sigmasis(sigmav, xx,/core)
      print,sigmav,sqrt(lum3par.sissigma2[binmax])
      oplot,xx,yfit,color=color
      oplot,[0,10000],[0,0]

      color=!blue
      oploterror,lum4.meanr[w],lum4.sigma[w],lum4.sigmaerr[w],$
        psym=5,color=color,errcolor=color
      IF keyword_set(bybin) THEN sigmav = sqrt(lum4par.sissigma2[binmax]) $
      ELSE sigmav=lum4par.sigmav
      sigmaverr = sqrt(lum4par.sissigmaerr2[binmax])
      yfit=sigmasis(sigmav, xx,/core)
      print,sigmav,sqrt(lum4par.sissigma2[binmax])
      oplot,xx,yfit,color=color
      oplot,[0,10000],[0,0]


      print

  ENDFOR 

  endplot

END 
