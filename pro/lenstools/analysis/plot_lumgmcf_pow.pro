PRO plot_lumgmcf_pow, stripes, paper=paper,diffyrange=diffyrange

  ;; plot combined gmcf for each luminosity bin
  ;; this is for the old mckay et al. paper

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

  IF keyword_set(loglog) THEN BEGIN
      addstr='_log' 
      xlog=1 & ylog=1
  ENDIF ELSE addstr=''

  IF keyword_set(paper) THEN addstr=addstr+'_paper'
  IF keyword_set(diffyrange) THEN addstr=addstr+'_diffyrange'

  psfile = basedir+'mass2light/lumbin_gmcf_pow'+addstr+'.ps'

  begplot,name=psfile

  ytitle = '!S'+!tsym.sigma_cap+'!R!A'+!tsym.minus+'!N ('+!tsym.ltequal+'R) '+!tsym.minus+' !S'+$
    !tsym.sigma_cap+'!R!A'+!tsym.minus+$
    '!N (R) (h M'+sunsymbol()+' pc!U'+!tsym.minus+'2!N)!X'
  xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)'

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

      w=lindgen(n_elements(lum1.meanr))
      binmax=2

      ;; Fit to R^{-0.8} power law
      alpha = -0.8
      xtt='A (h M'+sunsymbol()+' pc!U'+!tsym.minus+'2!N)'
      normrange = [-5,30]

      nnorm = 1000
      norm=arrscl( findgen(nnorm), normrange[0], normrange[1] )
      pow_chisq_conf_1par,lum1.meanr[w]/1000.,$
        lum1.sigma[w],lum1.sigmaerr[w],norm,alpha,chsqsf,norm1,normlow1,normhigh1, $
        xtitle=xtt
      pow_chisq_conf_1par,lum2.meanr[w]/1000.,$
        lum2.sigma[w],lum2.sigmaerr[w],norm,alpha,chsqsf,norm2,normlow2,normhigh2, $
        xtitle=xtt
      pow_chisq_conf_1par,lum3.meanr[w]/1000.,$
        lum3.sigma[w],lum3.sigmaerr[w],norm,alpha,chsqsf,norm3,normlow3,normhigh3, $
        xtitle=xtt
      pow_chisq_conf_1par,lum4.meanr[w]/1000.,$
        lum4.sigma[w],lum4.sigmaerr[w],norm,alpha,chsqsf,norm4,normlow4,normhigh4, $
        xtitle=xtt

      clevel=0
      errhigh1=normhigh1[clevel]-norm1 & errlow1 = norm1-normlow1[clevel]
      errhigh2=normhigh2[clevel]-norm2 & errlow2 = norm2-normlow2[clevel]    
      errhigh3=normhigh3[clevel]-norm3 & errlow3 = norm3-normlow3[clevel]
      errhigh4=normhigh4[clevel]-norm4 & errlow4 = norm4-normlow4[clevel]

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
          yrange[1] = 1.1*yrange[1]
          yrange[0] = .1
          xrange=[50,1100]
          xstyle=1
          ystyle=1
      ENDIF 
      erase & multiplot,[1,4]

      xx=arrscl(findgen(1000), .020, 1.)

      IF keyword_set(diffyrange) THEN delvarx,yrange,xrange

      ploterror,lum1.meanr[w],lum1.sigma[w],lum1.sigmaerr[w],$
        psym=1,yrange=yrange,xrange=xrange,xlog=xlog,ylog=ylog,$
        xstyle=xstyle,ystyle=ystyle

      norm = norm1
      err = errhigh1
      yfit = norm*xx^(alpha)
      oplot,xx*1000.,yfit
      oplot,[0,10000],[0,0]      
      message[0] = '<L> = '+ntostr(lumbin.lum[iclr,0],5)
      message[1] = 'A = '+$
        ntostr(rnd(norm,2),4)+' '+!tsym.plusminus+ntostr(rnd(err,2),4)
      IF NOT keyword_set(paper) THEN BEGIN 
          legend,message,/right
          legend,!colorsp[iclr],box=0,/center,/top,charsize=1.5
      ENDIF ELSE legend,!colorsp[iclr],box=0,/right,charsize=1.5
      

      IF keyword_set(diffyrange) THEN delvarx,yrange,xrange
      multiplot
      ploterror,lum2.meanr[w],lum2.sigma[w],lum2.sigmaerr[w],$
        psym=1,yrange=yrange,xrange=xrange,xlog=xlog,ylog=ylog,$
        xstyle=xstyle,ystyle=ystyle


      norm = norm2
      err = errhigh2
      yfit = norm*xx^(alpha)
      oplot,xx*1000.,yfit
      oplot,[0,10000],[0,0]      
      message[0] = '<L> = '+ntostr(lumbin.lum[iclr,1],5)
      message[1] = 'A = '+$
        ntostr(rnd(norm,2),4)+' '+!tsym.plusminus+ntostr(rnd(err,2),4)
      IF NOT keyword_set(paper) THEN legend,message,/right

      IF keyword_set(diffyrange) THEN delvarx,yrange,xrange
      multiplot
      ploterror,lum3.meanr[w],lum3.sigma[w],lum3.sigmaerr[w],$
        psym=1,yrange=yrange,xrange=xrange,xlog=xlog,ylog=ylog,$
        xstyle=xstyle,ystyle=ystyle

      norm = norm3
      err = errhigh3
      yfit = norm*xx^(alpha)
      oplot,xx*1000.,yfit
      oplot,[0,10000],[0,0]      
      message[0] = '<L> = '+ntostr(lumbin.lum[iclr,2],5)
      message[1] = 'A = '+$
        ntostr(rnd(norm,2),4)+' '+!tsym.plusminus+ntostr(rnd(err,2),4)
      IF NOT keyword_set(paper) THEN legend,message,/right

      IF keyword_set(diffyrange) THEN delvarx,yrange,xrange
      multiplot
      ploterror,lum4.meanr[w],lum4.sigma[w],lum4.sigmaerr[w],$
        psym=1,yrange=yrange,xrange=xrange,$
        xtitle=xtitle,ytitle=ytitle,xlog=xlog,ylog=ylog,$
        xstyle=xstyle,ystyle=ystyle

      norm = norm4
      err = errhigh4
      yfit = norm*xx^(alpha)
      oplot,xx*1000.,yfit
      oplot,[0,10000],[0,0]      
      message[0] = '<L> = '+ntostr(lumbin.lum[iclr,3],5)
      message[1] = 'A = '+$
        ntostr(rnd(norm,2),4)+' '+!tsym.plusminus+ntostr(rnd(err,2),4)
      IF NOT keyword_set(paper) THEN legend,message,/right

      multiplot,/reset

      print

  ENDFOR 

  endplot

END 
