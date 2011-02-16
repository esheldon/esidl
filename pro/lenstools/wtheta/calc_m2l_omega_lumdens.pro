PRO calc_m2l_omega_lumdens, powclr, m2l, m2lel, m2leh, $
                            omegahm2, omegahm2_el, omegahm2_eh,$
                            M0, M0el, M0eh, L0, L0el, L0eh,$
                            lpowers,lpowersel,lpowerseh,$
                            doplot=doplot, dofits=dofits, fourbin=fourbin, $
                            fixlumnorm=fixlumnorm,usefixlum=usefixlum,$
                            lumextno=lumextno, mextno=mextno,$
                            indir=indir,mindir=mindir

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: calc_m2l_omega_lumdens, powclr, m2l, m2lel, m2leh, '
      print,'                omegahm2, omegahm2_el, omegahm2_eh,'
      print,'                M0, M0el, M0eh, L0, L0el, L0eh,'
      print,'                lpowers,lpowersel,lpowerseh,'
      print,'                doplot=doplot, fourbin=fourbin, '
      print,'                fixlumnorm=fixlumnorm,usefixlum=usefixlum,'
      print,'                lumextno=lumextno, mextno=mextno'
      print,'                indir=indir,mindir=mindir'
      return
  ENDIF 

  delvarx,lpowers,lpowersel,lpowerseh

  IF NOT keyword_set(fourbin) THEN ttstr='twobin' ELSE ttstr=''

  lumfix = [1.314430862, 1.235023618, 1.194514446, 1.227677096, 1.220581657]
  meanlum = [0.784103, 0.889751, 1.51780, 2.07778, 2.58957]
  errfac=1.32
  
  ;; lum fits indir
  IF n_elements(indir) EQ 0 THEN indir = '/sdss5/data0/lensout/mass2light/'
  ;; mass indir
  IF n_elements(mindir) EQ 0 THEN mindir = '/sdss5/data0/lensout/stripe10/'

  IF keyword_set(usefixlum) AND keyword_set(fixlumnorm) THEN message,$
    'do not set /usefixlum and /fixlumnorm'

  ;; lumw extenstion
  ;; N3 is 0.1 L* stuff
  ;; N2 is magcut stuff
  ;; N1 is 0.1 L* stuff to 2 Mpc

  ;;default is ext = 'N3.fit'
  IF n_elements(lumextno) EQ 0 THEN lumextno = 3
  lextstr = 'N'+ntostr( long( round(lumextno) ) )
  lext = lextstr+'.fit'

  ;; mass extension
  ;; 1 is 1Mpc
  ;; 2 is 2Mpc
  IF n_elements(mextno) EQ 0 THEN mextno = 1
  mextstr = 'N'+ntostr( long( round(mextno) ) )
  mext = mextstr+'.fit'

  IF keyword_set(usefixlum) THEN BEGIN
      fstr = '_fixlum_' 
  ENDIF ELSE BEGIN
      fstr = '_'
  ENDELSE 
  psend = fstr+lextstr+'-'+mextstr+'.ps'
  fitend = fstr+lextstr+'-'+mextstr+'.fit'
  IF keyword_set(fixlumnorm) THEN BEGIN 
      psend = '_fixlumnorm'+psend
      fitend = '_fixlumnorm'+fitend
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output files, same extension as input
  ;; files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  outps = indir+'lumdense'+ttstr+'_m2l_omega_'+!colors[powclr]+psend
  outfits = indir+'lumdense'+ttstr+'_m2l_omega_'+!colors[powclr]+fitend

  print
  print,'Output ps file: ',outps
  print,'Output fits file: ',outfits
  print

  IF keyword_set(doplot) THEN begplot, name=outps, /color

  ;; !!! kludge to use new lum stuff E.S.S.

  allf = indir+'lumdens_wthetalumweq_stripe10_sum_rw_'+lext
  lall = mrdfits(allf,1,/silent)
print,allf
;  lall = mrdfits('/sdss6/data0/wtheta/lumdens_wthetalumweq_stripe10_sum_rw_N1.fit',1,/silent)
  lspiral = mrdfits(indir+'lumdense_spiral'+fstr+lext,1,/silent)
  lellip = mrdfits(indir+'lumdense_ellip'+fstr+lext,1,/silent)

  llum1 = mrdfits(indir+'lumdense_lumnum1'+ttstr+''+fstr+lext,1,/silent)
  llum2 = mrdfits(indir+'lumdense_lumnum2'+ttstr+''+fstr+lext,1,/silent)
  IF keyword_set(fourbin) THEN BEGIN 
      llum3 = mrdfits(indir+'lumdense_lumnum3'+ttstr+''+fstr+lext,1,/silent)
      llum4 = mrdfits(indir+'lumdense_lumnum4'+ttstr+''+fstr+lext,1,/silent)
  ENDIF 

  stripestr='stripe10_stripe36_stripe37_stripe42_stripe43_stripe82'
  stripestr='stripe10'

  all = mrdfits(mindir+'main_zgal_gal_'+stripestr+'_comb_corr_N6.fit',1,/silent)
  spiral=mrdfits(mindir+'spiral_zgal_gal_'+stripestr+'_comb_corr_'+mext,1,/silent)
  ellip=mrdfits(mindir+'ellip_zgal_gal_'+stripestr+'_comb_corr_'+mext,1,/silent)

  lum1=mrdfits(mindir+'sublum/i/lum1'+ttstr+'_zgal_gal_'+stripestr+'_comb_corr_'+mext,1,/silent)
  lum2=mrdfits(mindir+'sublum/i/lum2'+ttstr+'_zgal_gal_'+stripestr+'_comb_corr_'+mext,1,/silent)
  IF keyword_set(fourbin) THEN BEGIN 
      lum3=mrdfits(mindir+'sublum/i/lum3'+ttstr+'_zgal_gal_'+stripestr+'_comb_corr_'+mext,1,/silent)
      lum4=mrdfits(mindir+'sublum/i/lum4'+ttstr+'_zgal_gal_'+stripestr+'_comb_corr_'+mext,1,/silent)
  ENDIF 

  ;; account for 0.1 l_star cut 
  IF keyword_set(fixlumnorm) THEN BEGIN
      fixlum, lall, powclr
      fixlum, lspiral, powclr
      fixlum, lellip, powclr
      fixlum, llum1, powclr
      fixlum, llum2, powclr
      IF keyword_set(fourbin) THEN BEGIN 
          fixlum, llum3, powclr
          fixlum, llum4, powclr
      ENDIF 

  ENDIF 

  all.sigmaerr = all.sigmaerr*errfac
  spiral.sigmaerr = spiral.sigmaerr*errfac
  ellip.sigmaerr = ellip.sigmaerr*errfac
  lum1.sigmaerr = lum1.sigmaerr*errfac
  lum2.sigmaerr = lum2.sigmaerr*errfac
  IF keyword_set(fourbin) THEN BEGIN 
      lum3.sigmaerr = lum3.sigmaerr*errfac
      lum4.sigmaerr = lum4.sigmaerr*errfac
  ENDIF 
  normrange = [0, 25]
  powrange=[0.3, 1.4]
  norm=arrscl( findgen(400), normrange[0], normrange[1] )
  power=arrscl( findgen(400), -powrange[0], -powrange[1] )

  xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)'
  ytitle = '!S'+!tsym.sigma_cap+'!R!A'+!tsym.minus+'!N ('+!tsym.ltequal+'R) '+$
    !tsym.minus+' !S'+!tsym.sigma_cap+'!R!A'+!tsym.minus+$
    '!N (R) (h M'+sunsymbol()+' pc!U'+!tsym.minus+'2!N)'
  xx=arrscl( findgen(100), min(all.meanr), max(all.meanr) )
;  IF display_type() EQ 'X' THEN clear=0 ELSE clear=1


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Range in radius to fit
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nr=n_elements(all.meanr)
  minrbin = 4
  wr=( lindgen(nr) )[minrbin:nr-1]

  !p.multi=[0,1,2]
  aspect=1
  pow_chisq_conf_1par, all.meanr[wr]/1000., all.sigma[wr], all.sigmaerr[wr], $
    norm, lall.pow[powclr], chisq_surf, $
    allnorm, allnormlow, allnormhigh, chisq_per=chisq_per,$
    aspect=aspect,/center
  legend,['all power = '+ntostr(lall.pow[powclr]),$
          !tsym.chi+'!U2!N/'+!tsym.nu+' = '+ntostr(chisq_per)],/clear
  yfit = allnorm*(xx/1000.)^(lall.pow[powclr])
  ploterror, all.meanr, all.sigma, all.sigmaerr,psym=1,$
    xtitle=xtitle,ytitle=ytitle
  oplot, xx, yfit
  oplot, [0, 10000], [0,0]
  legend, !colorsp[powclr], /right, box=0
  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  pow_chisq_conf_1par, spiral.meanr[wr]/1000., spiral.sigma[wr], spiral.sigmaerr[wr], $
    norm, lspiral.pow[powclr], chisq_surf, $
    spiralnorm, spiralnormlow, spiralnormhigh, chisq_per=chisq_per,$
    aspect=aspect,/center
    legend,['spiral power = '+ntostr(lspiral.pow[powclr]),$
          !tsym.chi+'!U2!N/'+!tsym.nu+' = '+ntostr(chisq_per)],/clear
  yfit = spiralnorm*(xx/1000.)^(lspiral.pow[powclr])
  ploterror, spiral.meanr, spiral.sigma, spiral.sigmaerr,psym=1,$
    xtitle=xtitle,ytitle=ytitle
  oplot, xx, yfit
  oplot, [0, 10000], [0,0]
  legend, !colorsp[powclr], /right, box=0
  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  pow_chisq_conf_1par, ellip.meanr[wr]/1000., ellip.sigma[wr], ellip.sigmaerr[wr], $
    norm, lellip.pow[powclr], chisq_surf, $
    ellipnorm, ellipnormlow, ellipnormhigh, chisq_per=chisq_per,$
    aspect=aspect,/center
    legend,['ellip power = '+ntostr(lellip.pow[powclr]),$
          !tsym.chi+'!U2!N/'+!tsym.nu+' = '+ntostr(chisq_per)],/clear
  yfit = ellipnorm*(xx/1000.)^(lellip.pow[powclr])
  ploterror, ellip.meanr, ellip.sigma, ellip.sigmaerr,psym=1,$
    xtitle=xtitle,ytitle=ytitle
  oplot, xx, yfit
  oplot, [0, 10000], [0,0]
  legend, !colorsp[powclr], /right, box=0
  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  pow_chisq_conf_1par, lum1.meanr[wr]/1000., lum1.sigma[wr], lum1.sigmaerr[wr], $
    norm, llum1.pow[powclr], chisq_surf, $
    lum1norm, lum1normlow, lum1normhigh, chisq_per=chisq_per,$
    aspect=aspect,/center
    legend,['lum1 power = '+ntostr(llum1.pow[powclr]),$
          !tsym.chi+'!U2!N/'+!tsym.nu+' = '+ntostr(chisq_per)],/clear
  yfit = lum1norm*(xx/1000.)^(llum1.pow[powclr])
  ploterror, lum1.meanr, lum1.sigma, lum1.sigmaerr,psym=1,$
    xtitle=xtitle,ytitle=ytitle
  oplot, xx, yfit
  oplot, [0, 10000], [0,0]
  legend, !colorsp[powclr], /right, box=0
  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  pow_chisq_conf_1par, lum2.meanr[wr]/1000., lum2.sigma[wr], lum2.sigmaerr[wr], $
    norm, llum2.pow[powclr], chisq_surf, $
    lum2norm, lum2normlow, lum2normhigh, chisq_per=chisq_per,$
    aspect=aspect,/center
    legend,['lum2 power = '+ntostr(llum2.pow[powclr]),$
          !tsym.chi+'!U2!N/'+!tsym.nu+' = '+ntostr(chisq_per)],/clear
  yfit = lum2norm*(xx/1000.)^(llum2.pow[powclr])
  ploterror, lum2.meanr, lum2.sigma, lum2.sigmaerr,psym=1,$
    xtitle=xtitle,ytitle=ytitle
  oplot, xx, yfit
  oplot, [0, 10000], [0,0]
  legend, !colorsp[powclr], /right, box=0
  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  IF keyword_set(fourbin) THEN BEGIN 
      pow_chisq_conf_1par, lum3.meanr[wr]/1000., lum3.sigma[wr], lum3.sigmaerr[wr], $
        norm, llum3.pow[powclr], chisq_surf, $
        lum3norm, lum3normlow, lum3normhigh, chisq_per=chisq_per,$
        aspect=aspect,/center
      legend,['lum3 power = '+ntostr(llum3.pow[powclr]),$
              !tsym.chi+'!U2!N/'+!tsym.nu+' = '+ntostr(chisq_per)],/clear
      yfit = lum3norm*(xx/1000.)^(llum3.pow[powclr])
      ploterror, lum3.meanr, lum3.sigma, lum3.sigmaerr,psym=1,$
        xtitle=xtitle,ytitle=ytitle
      oplot, xx, yfit
      oplot, [0, 10000], [0,0]
      legend, !colorsp[powclr], /right, box=0
      IF display_type() EQ 'X' THEN key=get_kbrd(1)
      
      pow_chisq_conf_1par, lum4.meanr[wr]/1000., lum4.sigma[wr], lum4.sigmaerr[wr], $
        norm, llum4.pow[powclr], chisq_surf, $
        lum4norm, lum4normlow, lum4normhigh, chisq_per=chisq_per,$
        aspect=aspect,/center
      legend,['lum4 power = '+ntostr(llum4.pow[powclr]),$
              !tsym.chi+'!U2!N/'+!tsym.nu+' = '+ntostr(chisq_per)],/clear
      yfit = lum4norm*(xx/1000.)^(llum4.pow[powclr])
      ploterror, lum4.meanr, lum4.sigma, lum4.sigmaerr,psym=1,$
        xtitle=xtitle,ytitle=ytitle
      oplot, xx, yfit
      oplot, [0, 10000], [0,0]
      legend, !colorsp[powclr], /right, box=0
      IF display_type() EQ 'X' THEN key=get_kbrd(1)
  ENDIF 

  ;; now calculate global mass to light and error

  ;;;;;;;;;;;;;;;;;;;;
  ;; all galaxies
  ;;;;;;;;;;;;;;;;;;;;

  pow = -lall.pow[powclr] & powlow=-lall.powlow[powclr] & powhigh=-lall.powhigh[powclr]
  norm = lall.norm[powclr] & normlow=lall.normlow[powclr] & normhigh=lall.normhigh[powclr]
  fac = (2.0-pow)/pow
  range2error, powlow, pow, powhigh, lpeh, lpel ;power errors
  range2error, normlow, norm, normhigh, lnel, lneh ;lum norm errors
  range2error, allnormlow[0], allnorm, allnormhigh[0], nel, neh ;mass norm errors
  add_arrval, pow, lpowers & add_arrval, abs(lpel), lpowersel 
  add_arrval, abs(lpeh), lpowerseh

  allm2l = (allnorm*1.e12)/(lall.norm[powclr]*1.e10)*fac
  allel = allm2l*sqrt( (nel/allnorm)^2 + (lnel/norm)^2 + 4./(2.-pow)^2*(lpel/pow)^2 )
  alleh = allm2l*sqrt( (neh/allnorm)^2 + (lneh/norm)^2 + 4./(2.-pow)^2*(lpeh/pow)^2 )

  ;; Calculate M0 and L0
  allM0 = 2.0*!pi/pow*allnorm   ;10^12 M_{\sun}
  allM0el = allM0*sqrt( (nel/allnorm)^2 + (lpel/pow)^2 )
  allM0eh = allM0*sqrt( (neh/allnorm)^2 + (lpeh/pow)^2 )
  allL0 = 2.0*!pi/(2.0-pow)*norm ;10^10 L_{\sun}
  allL0el = allL0*sqrt( (lnel/norm)^2 + (lpel/(2.-pow))^2 )
  allL0eh = allL0*sqrt( (lneh/norm)^2 + (lpeh/(2.-pow))^2 )
  print,'relative errors'
  print,abs(neh/allnorm),abs(lneh/norm),abs(2./(2.-pow)*(lpeh/pow))

  ;;;;;;;;;;;;;;;;;;;;
  ;; spiral galaxies
  ;;;;;;;;;;;;;;;;;;;;

  pow = -lspiral.pow[powclr] & powlow=-lspiral.powlow[powclr] & powhigh=-lspiral.powhigh[powclr]
  norm = lspiral.norm[powclr] & normlow=lspiral.normlow[powclr] & normhigh=lspiral.normhigh[powclr]
  fac = (2.0-pow)/pow
  range2error, powlow, pow, powhigh, lpeh, lpel
  range2error, normlow, norm, normhigh, lnel, lneh
  range2error, spiralnormlow[0], spiralnorm, spiralnormhigh[0], nel, neh
  add_arrval, pow, lpowers & add_arrval, abs(lpel), lpowersel & add_arrval, abs(lpeh), lpowerseh

  spiralm2l = (spiralnorm*1.e12)/(lspiral.norm[powclr]*1.e10)*fac
  spiralel = spiralm2l*sqrt( (nel/spiralnorm)^2 + (lnel/norm)^2 + 4./(2.-pow)^2*(lpel/pow)^2 )
  spiraleh = spiralm2l*sqrt( (neh/spiralnorm)^2 + (lneh/norm)^2 + 4./(2.-pow)^2*(lpeh/pow)^2 )

  ;; Calculate M0 and L0
  spiralM0 = 2.0*!pi/pow*spiralnorm   ;10^12 M_{\sun}
  spiralM0el = spiralM0*sqrt( (nel/spiralnorm)^2 + (lpel/pow)^2 )
  spiralM0eh = spiralM0*sqrt( (neh/spiralnorm)^2 + (lpeh/pow)^2 )
  spiralL0 = 2.0*!pi/(2.0-pow)*norm ;10^10 L_{\sun}
  spiralL0el = spiralL0*sqrt( (lnel/norm)^2 + (lpel/(2.-pow))^2 )
  spiralL0eh = spiralL0*sqrt( (lneh/norm)^2 + (lpeh/(2.-pow))^2 )

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; elliptical galaxies
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  pow = -lellip.pow[powclr] & powlow=-lellip.powlow[powclr] & powhigh=-lellip.powhigh[powclr]
  norm = lellip.norm[powclr] & normlow=lellip.normlow[powclr] & normhigh=lellip.normhigh[powclr]
  fac = (2.0-pow)/pow
  range2error, powlow, pow, powhigh, lpeh, lpel
  range2error, normlow, norm, normhigh, lnel, lneh
  range2error, ellipnormlow[0], ellipnorm, ellipnormhigh[0], nel, neh
  add_arrval, pow, lpowers & add_arrval, abs(lpel), lpowersel & add_arrval, abs(lpeh), lpowerseh

  ellipm2l = (ellipnorm*1.e12)/(lellip.norm[powclr]*1.e10)*fac
  ellipel = ellipm2l*sqrt( (nel/ellipnorm)^2 + (lnel/norm)^2 + 4./(2.-pow)^2*(lpel/pow)^2 )
  ellipeh = ellipm2l*sqrt( (neh/ellipnorm)^2 + (lneh/norm)^2 + 4./(2.-pow)^2*(lpeh/pow)^2 )

  ;; Calculate M0 and L0
  ellipM0 = 2.0*!pi/pow*ellipnorm   ;10^12 M_{\sun}
  ellipM0el = ellipM0*sqrt( (nel/ellipnorm)^2 + (lpel/pow)^2 )
  ellipM0eh = ellipM0*sqrt( (neh/ellipnorm)^2 + (lpeh/pow)^2 )
  ellipL0 = 2.0*!pi/(2.0-pow)*norm ;10^10 L_{\sun}
  ellipL0el = ellipL0*sqrt( (lnel/norm)^2 + (lpel/(2.-pow))^2 )
  ellipL0eh = ellipL0*sqrt( (lneh/norm)^2 + (lpeh/(2.-pow))^2 )

  ;;;;;;;;;
  ;; lum1
  ;;;;;;;;;

  pow = -llum1.pow[powclr] & powlow=-llum1.powlow[powclr] & powhigh=-llum1.powhigh[powclr]
  norm = llum1.norm[powclr] & normlow=llum1.normlow[powclr] & normhigh=llum1.normhigh[powclr]
  fac = (2.0-pow)/pow
  range2error, powlow, pow, powhigh, lpeh, lpel
  range2error, normlow, norm, normhigh, lnel, lneh
  range2error, lum1normlow[0], lum1norm, lum1normhigh[0], nel, neh
  add_arrval, pow, lpowers & add_arrval, abs(lpel), lpowersel & add_arrval, abs(lpeh), lpowerseh

  lum1m2l = (lum1norm*1.e12)/(llum1.norm[powclr]*1.e10)*fac
  lum1el = lum1m2l*sqrt( (nel/lum1norm)^2 + (lnel/norm)^2 + 4./(2.-pow)^2*(lpel/pow)^2 )
  lum1eh = lum1m2l*sqrt( (neh/lum1norm)^2 + (lneh/norm)^2 + 4./(2.-pow)^2*(lpeh/pow)^2 )

  ;; Calculate M0 and L0
  lum1M0 = 2.0*!pi/pow*lum1norm   ;10^12 M_{\sun}
  lum1M0el = lum1M0*sqrt( (nel/lum1norm)^2 + (lpel/pow)^2 )
  lum1M0eh = lum1M0*sqrt( (neh/lum1norm)^2 + (lpeh/pow)^2 )
  lum1L0 = 2.0*!pi/(2.0-pow)*norm ;10^10 L_{\sun}
  lum1L0el = lum1L0*sqrt( (lnel/norm)^2 + (lpel/(2.-pow))^2 )
  lum1L0eh = lum1L0*sqrt( (lneh/norm)^2 + (lpeh/(2.-pow))^2 )

  ;;;;;;;;;
  ;; lum2
  ;;;;;;;;;

  pow = -llum2.pow[powclr] & powlow=-llum2.powlow[powclr] & powhigh=-llum2.powhigh[powclr]
  norm = llum2.norm[powclr] & normlow=llum2.normlow[powclr] & normhigh=llum2.normhigh[powclr]
  fac = (2.0-pow)/pow
  range2error, powlow, pow, powhigh, lpeh, lpel
  range2error, normlow, norm, normhigh, lnel, lneh
  range2error, lum2normlow[0], lum2norm, lum2normhigh[0], nel, neh
  add_arrval, pow, lpowers & add_arrval, abs(lpel), lpowersel & add_arrval, abs(lpeh), lpowerseh

  lum2m2l = (lum2norm*1.e12)/(llum2.norm[powclr]*1.e10)*fac
  lum2el = lum2m2l*sqrt( (nel/lum2norm)^2 + (lnel/norm)^2 + 4./(2.-pow)^2*(lpel/pow)^2 )
  lum2eh = lum2m2l*sqrt( (neh/lum2norm)^2 + (lneh/norm)^2 + 4./(2.-pow)^2*(lpeh/pow)^2 )

  ;; Calculate M0 and L0
  lum2M0 = 2.0*!pi/pow*lum2norm   ;10^12 M_{\sun}
  lum2M0el = lum2M0*sqrt( (nel/lum2norm)^2 + (lpel/pow)^2 )
  lum2M0eh = lum2M0*sqrt( (neh/lum2norm)^2 + (lpeh/pow)^2 )
  lum2L0 = 2.0*!pi/(2.0-pow)*norm ;10^10 L_{\sun}
  lum2L0el = lum2L0*sqrt( (lnel/norm)^2 + (lpel/(2.-pow))^2 )
  lum2L0eh = lum2L0*sqrt( (lneh/norm)^2 + (lpeh/(2.-pow))^2 )

  IF keyword_set(fourbin) THEN BEGIN 
      ;;;;;;;;;;;;
      ;; lum3
      ;;;;;;;;;;;;

      pow = -llum3.pow[powclr] & powlow=-llum3.powlow[powclr] & powhigh=-llum3.powhigh[powclr]
      norm = llum3.norm[powclr] & normlow=llum3.normlow[powclr] & normhigh=llum3.normhigh[powclr]
      fac = (2.0-pow)/pow
      range2error, powlow, pow, powhigh, lpeh, lpel
      range2error, normlow, norm, normhigh, lnel, lneh
      range2error, lum3normlow[0], lum3norm, lum3normhigh[0], nel, neh
      add_arrval, pow, lpowers & add_arrval, abs(lpel), lpowersel & add_arrval, abs(lpeh), lpowerseh

      lum3m2l = (lum3norm*1.e12)/(llum3.norm[powclr]*1.e10)*fac
      lum3el = lum3m2l*sqrt( (nel/lum3norm)^2 + (lnel/norm)^2 + 4./(2.-pow)^2*(lpel/pow)^2 )
      lum3eh = lum3m2l*sqrt( (neh/lum3norm)^2 + (lneh/norm)^2 + 4./(2.-pow)^2*(lpeh/pow)^2 )

      ;; Calculate M0 and L0
      lum3M0 = 2.0*!pi/pow*lum3norm ;10^12 M_{\sun}
      lum3M0el = lum3M0*sqrt( (nel/lum3norm)^2 + (lpel/pow)^2 )
      lum3M0eh = lum3M0*sqrt( (neh/lum3norm)^2 + (lpeh/pow)^2 )
      lum3L0 = 2.0*!pi/(2.0-pow)*norm ;10^10 L_{\sun}
      lum3L0el = lum3L0*sqrt( (lnel/norm)^2 + (lpel/(2.-pow))^2 )
      lum3L0eh = lum3L0*sqrt( (lneh/norm)^2 + (lpeh/(2.-pow))^2 )

      ;;;;;;;;;;;
      ;; lum4
      ;;;;;;;;;;;

      pow = -llum4.pow[powclr] & powlow=-llum4.powlow[powclr] & powhigh=-llum4.powhigh[powclr]
      norm = llum4.norm[powclr] & normlow=llum4.normlow[powclr] & normhigh=llum4.normhigh[powclr]
      fac = (2.0-pow)/pow
      range2error, powlow, pow, powhigh, lpeh, lpel
      range2error, normlow, norm, normhigh, lnel, lneh
      range2error, lum4normlow[0], lum4norm, lum4normhigh[0], nel, neh
      add_arrval, pow, lpowers & add_arrval, abs(lpel), lpowersel & add_arrval, abs(lpeh), lpowerseh

      lum4m2l = (lum4norm*1.e12)/(llum4.norm[powclr]*1.e10)*fac
      lum4el = lum4m2l*sqrt( (nel/lum4norm)^2 + (lnel/norm)^2 + 4./(2.-pow)^2*(lpel/pow)^2 )
      lum4eh = lum4m2l*sqrt( (neh/lum4norm)^2 + (lneh/norm)^2 + 4./(2.-pow)^2*(lpeh/pow)^2 )
 
      ;; Calculate M0 and L0
      lum4M0 = 2.0*!pi/pow*lum4norm ;10^12 M_{\sun}
      lum4M0el = lum4M0*sqrt( (nel/lum4norm)^2 + (lpel/pow)^2 )
      lum4M0eh = lum4M0*sqrt( (neh/lum4norm)^2 + (lpeh/pow)^2 )
      lum4L0 = 2.0*!pi/(2.0-pow)*norm ;10^10 L_{\sun}
      lum4L0el = lum4L0*sqrt( (lnel/norm)^2 + (lpel/(2.-pow))^2 )
      lum4L0eh = lum4L0*sqrt( (lneh/norm)^2 + (lpeh/(2.-pow))^2 )

  ENDIF 
  print
  print
  print,!colors[powclr]+'-band'
  print,'all M0 = ',ntostr(allM0),' + ',ntostr(allM0eh),' - ',ntostr(allM0el)
  print,'all L0 = ',ntostr(allL0),' + ',ntostr(allL0eh),' - ',ntostr(allL0el)
  print,'spiral M0 = ',ntostr(spiralM0),' + ',ntostr(spiralM0eh),' - ',ntostr(spiralM0el)
  print,'spiral L0 = ',ntostr(spiralL0),' + ',ntostr(spiralL0eh),' - ',ntostr(spiralL0el)
  print,'ellip M0 = ',ntostr(ellipM0),' + ',ntostr(ellipM0eh),' - ',ntostr(ellipM0el)
  print,'ellip L0 = ',ntostr(ellipL0),' + ',ntostr(ellipL0eh),' - ',ntostr(ellipL0el)
  print,'lum1 M0 = ',ntostr(lum1M0),' + ',ntostr(lum1M0eh),' - ',ntostr(lum1M0el)
  print,'lum1 L0 = ',ntostr(lum1L0),' + ',ntostr(lum1L0eh),' - ',ntostr(lum1L0el)
  print,'lum2 M0 = ',ntostr(lum2M0),' + ',ntostr(lum2M0eh),' - ',ntostr(lum2M0el)
  print,'lum2 L0 = ',ntostr(lum2L0),' + ',ntostr(lum2L0eh),' - ',ntostr(lum2L0el)
  IF keyword_set(fourbin) THEN BEGIN 
      print,'lum3 M0 = ',ntostr(lum3M0),' + ',ntostr(lum3M0eh),' - ',ntostr(lum3M0el)
      print,'lum3 L0 = ',ntostr(lum3L0),' + ',ntostr(lum3L0eh),' - ',ntostr(lum3L0el)
      print,'lum4 M0 = ',ntostr(lum4M0),' + ',ntostr(lum4M0eh),' - ',ntostr(lum4M0el)
      print,'lum4 L0 = ',ntostr(lum4L0),' + ',ntostr(lum4L0eh),' - ',ntostr(lum4L0el)
  ENDIF 
  print
  print,'All M/L: ',ntostr(allm2l),' + '+ntostr(alleh)+' - '+ntostr(allel)
  print,'Spiral M/L: ',ntostr(spiralm2l),' + '+ntostr(spiraleh)+' - '+ntostr(spiralel)
  print,'Ellip M/L: ',ntostr(ellipm2l),' + '+ntostr(ellipeh)+' - '+ntostr(ellipel)
  print,'Lum1 M/L: ',ntostr(lum1m2l),' + '+ntostr(lum1eh)+' - '+ntostr(lum1el)
  print,'Lum2 M/L: ',ntostr(lum2m2l),' + '+ntostr(lum2eh)+' - '+ntostr(lum2el)
  IF keyword_set(fourbin) THEN BEGIN 
      print,'Lum3 M/L: ',ntostr(lum3m2l),' + '+ntostr(lum3eh)+' - '+ntostr(lum3el)
      print,'Lum4 M/L: ',ntostr(lum4m2l),' + '+ntostr(lum4eh)+' - '+ntostr(lum4el)
  ENDIF 
  print

  ;; Characteristic Radius, where M/L falls to a fraction f of the
  ;; global value
  ;; Need to get mean luminosities


  colors = [!p.color, !blue, !red, !cyan, !green]
  psym = [1, 2, 4, 5, 6]
  m2l = [allm2l, spiralm2l, ellipm2l, lum1m2l, lum2m2l]
  m2leh = [alleh, spiraleh, ellipeh, lum1eh, lum2eh]
  m2lel = [allel, spiralel, ellipel, lum1el, lum2el]
  M0 = [allM0, spiralM0, ellipM0, lum1M0, lum2M0]
  M0eh = [allM0eh, spiralM0eh, ellipM0eh, lum1M0eh, lum2M0eh]
  M0el = [allM0el, spiralM0el, ellipM0el, lum1M0el, lum2M0el]
  L0 = [allL0, spiralL0, ellipL0, lum1L0, lum2L0]
  L0eh = [allL0eh, spiralL0eh, ellipL0eh, lum1L0eh, lum2L0eh]
  L0el = [allL0el, spiralL0el, ellipL0el, lum1L0el, lum2L0el]

  types = ['all','spiral','ellip','lum1','lum2']
  IF keyword_set(fourbin) THEN BEGIN 
      colors = [colors, !magenta, !yellow]
      psym = [psym, 7, 8]
      m2l = [m2l, lum3m2l, lum4m2l]
      m2leh = [m2leh, lum3eh, lum4eh]
      m2lel = [m2lel, lum3el, lum4el]

      M0 = [M0, lum3M0, lum4M0]
      M0eh = [M0eh, lum3M0eh, lum4M0eh]
      M0el = [M0el, lum3M0el, lum4M0el]
      L0 = [L0, lum3L0, lum4L0]
      L0eh = [L0eh, lum3L0eh, lum4L0eh]
      L0el = [L0el, lum3L0el, lum4L0el]

      types = [types, 'lum3','lum4']
  ENDIF 
  struct = create_struct('types',types,$
                         'm2l',m2l,$
                         'm2l_el',m2lel,$
                         'm2l_eh',m2leh)

  nnn = n_elements(types)
  myusersym, 'fill_circle'

  ;;;;;;;;;;;;;;;;;
  ;; M/L
  ;;;;;;;;;;;;;;;;;

  xrange = [0, nnn+1]
  nxx=n_elements(m2l)
  xx=findgen(nxx) + 1
  ytitle='M/L!D '+!colorsp[powclr]+'!N ( h M'+sunsymbol()+' /L'+sunsymbol()+' )'
  xtitle='Galaxy sub-sample'
  plotderror, xx, m2l, xx, xx, m2l-m2lel, m2l+m2leh, psym=3,$
    xtitle=xtitle,ytitle=ytitle,xrange=xrange,xstyle=1,yrange=[50,300]
  FOR i=0L, nxx-1 DO oplotderror, [xx[i]], [m2l[i]], [xx[i]], [xx[i]], [(m2l-m2lel)[i]], [(m2l+m2leh)[i]], $
    psym=psym[i], color=colors[i], errcolor=colors[i]
  legend, types, color=colors,$
    thick=replicate(!p.thick,nnn), psym=psym, /left
  legend,!colorsp[powclr],/right,box=0

  ;;;;;;;;;;;;;;;;
  ;; OMEGA
  ;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Blanton's luminosity densities
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  blanton_lumdensity, [0,1,2,3,4], ld, lderr;, omegamat=1.0

  ;; dan's formula
;;  rho_crit = double(2.767e-7)            ;Msolar/pc^3
  ;; my formula
  rho_crit = double(2.77545e-7)
  rho_crit = rho_crit*(10d)^18.


  ;; Here use h=1

  ;;;;;;;;;
  ;; I think there is no h in our measurement of omega, so
  ;; just use this one which has herr=0.

  h = 1.0
  herr = 0.0
  hstr='1'
  herrstr = '0'
  omegahm2 = m2l*ld[powclr]/rho_crit/h^2;*2./!dpi
  omegahm2_eh = omegahm2*sqrt( (lderr[powclr]/ld[powclr])^2 + $
                               (m2leh/m2l)^2 + 4.0*(herr/h)^2 )
  omegahm2_el = omegahm2*sqrt( (lderr[powclr]/ld[powclr])^2 + $
                               (m2lel/m2l)^2 + 4.0*(herr/h)^2 )
  struct = create_struct(struct, $
                         'omegahm2',omegahm2,$
                         'omegahm2_el',omegahm2_el,$
                         'omegahm2_eh',omegahm2_eh)

  print
  print,!colors[powclr]+'-band'
  print,'All Omega Matter : ',ntostr(omegahm2[0]),' + '+$
    ntostr(omegahm2_eh[0])+' - '+ntostr(omegahm2_el[0])
  print,'Spiral Omega Matter : ',ntostr(omegahm2[1]),' + '+$
    ntostr(omegahm2_eh[1])+' - '+ntostr(omegahm2_el[1])
  print,'Ellip Omega Matter : ',ntostr(omegahm2[2]),' + '+$
    ntostr(omegahm2_eh[2])+' - '+ntostr(omegahm2_el[2])
  print,'Lum1 Omega Matter : ',ntostr(omegahm2[3]),' + '+$
    ntostr(omegahm2_eh[3])+' - '+ntostr(omegahm2_el[3])
  print,'Lum2 Omega Matter : ',ntostr(omegahm2[4]),' + '+$
    ntostr(omegahm2_eh[4])+' - '+ntostr(omegahm2_el[4])
  IF keyword_set(fourbin) THEN BEGIN 
      print,'Lum3 Omega Matter : ',ntostr(omegahm2[5]),' + '+$
        ntostr(omegahm2_eh[5])+' - '+ntostr(omegahm2_el[5])
      print,'Lum4 Omega Matter : ',ntostr(omegahm2[6]),' + '+$
        ntostr(omegahm2_eh[6])+' - '+ntostr(omegahm2_el[6])
  ENDIF 
  print


  ;;ytitle=!tsym.omega_cap+'!DM!N h!U2!N'
  ytitle = !tsym.omega_cap+'!DM!N'
  xtitle='Galaxy sub-sample'
  plotderror, xx, omegahm2, xx, xx, omegahm2-omegahm2_el, omegahm2+omegahm2_eh, psym=3,xtitle=xtitle,ytitle=ytitle,xrange=xrange,xstyle=1,yrange=[0.1, 0.4]
  FOR i=0L, nxx-1 DO oplotderror, [xx[i]], [omegahm2[i]], [xx[i]], [xx[i]], [(omegahm2-omegahm2_el)[i]], [(omegahm2+omegahm2_eh)[i]], $
    psym=psym[i], color=colors[i], errcolor=colors[i]


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; up until now we used h=1. Now add in a prior
  ;; H0 = 72+/-8 (Freedman, et al. 2001
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  h = 0.72
  herr = 0.08
  pomega = 0.2/h
  pomegaerr = pomega*sqrt( (0.03/0.2)^2 + (herr/h)^2 )
  oplot,[0,100], [pomega,pomega],line=1
  oplot,[0,100], [pomega,pomega]+pomegaerr,line=2
  oplot,[0,100], [pomega,pomega]-pomegaerr,line=2
  print,'Peacock omega = ',ntostr(pomega),' +/- ',ntostr(pomegaerr)
  print

  legend, types, color=colors,$
    thick=replicate(!p.thick,nnn), psym=psym, /left,/clear
  ;;legend,[!colorsp[powclr], 'h = '+hstr+' '+!tsym.plusminus+' '+herrstr],/right,box=0
  legend,!colorsp[powclr],/right,box=0

  IF keyword_set(doplot) THEN endplot
  return

  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  hstr='0.72'
  herrstr = '0.08'
  omega = m2l*ld[powclr]/rho_crit/h^2;*2./!dpi
  eh = omega*sqrt( (lderr[powclr]/ld[powclr])^2 + (m2leh/m2l)^2 + 4.0*(herr/h)^2 )
  el = omega*sqrt( (lderr[powclr]/ld[powclr])^2 + (m2lel/m2l)^2 + 4.0*(herr/h)^2 )
  struct = create_struct(struct, $
                         'omega',omega,$
                         'omega_el',el,$
                         'omega_eh',eh)
  print
  print,!colors[powclr]+'-band'
  print,'All Omega Matter: ',ntostr(omega[0]),' + '+ntostr(eh[0])+' - '+ntostr(el[0])
  print,'Spiral Omega Matter: ',ntostr(omega[1]),' + '+ntostr(eh[1])+' - '+ntostr(el[1])
  print,'Ellip Omega Matter: ',ntostr(omega[2]),' + '+ntostr(eh[2])+' - '+ntostr(el[2])
  print,'Lum1 Omega Matter: ',ntostr(omega[3]),' + '+ntostr(eh[3])+' - '+ntostr(el[3])
  print,'Lum2 Omega Matter: ',ntostr(omega[4]),' + '+ntostr(eh[4])+' - '+ntostr(el[4])
  IF keyword_set(fourbin) THEN BEGIN 
      print,'Lum3 Omega Matter: ',ntostr(omega[5]),' + '+ntostr(eh[5])+' - '+ntostr(el[5])
      print,'Lum4 Omega Matter: ',ntostr(omega[6]),' + '+ntostr(eh[6])+' - '+ntostr(el[6])
  ENDIF 
  print

  ytitle=!tsym.omega_cap+'!DM!N'
  xtitle='Galaxy sub-sample'
  plotderror, xx, omega, xx, xx, omega-el, omega+eh, psym=3,xtitle=xtitle,ytitle=ytitle,xrange=xrange,xstyle=1,yrange=[0.1, 0.7]
  FOR i=0L, nxx-1 DO oplotderror, [xx[i]], [omega[i]], [xx[i]], [xx[i]], [(omega-el)[i]], [(omega+eh)[i]], $
    psym=psym[i], color=colors[i], errcolor=colors[i]
  legend, types, color=colors,$
    thick=replicate(!p.thick,nnn), psym=psym, /left
  legend,[!colorsp[powclr], 'h = '+hstr+' '+!tsym.plusminus+' '+herrstr],/right,box=0

  ;; overplot result from peacock et al.
  pomega = 0.2/h
  pomegaerr = pomega*sqrt( (0.03/0.2)^2 + (herr/h)^2 )
  print,'Peacock omega = ',ntostr(pomega),' +/- ',ntostr(pomegaerr)
  print

  oplot,[0,100], [pomega,pomega],line=1
  oplot,[0,100], [pomega,pomega]+pomegaerr,line=2
  oplot,[0,100], [pomega,pomega]-pomegaerr,line=2

  IF keyword_set(doplot) THEN endplot

  print,'Output fits file: ',outfits
  IF keyword_set(dofits) THEN mwrfits, struct, outfits, /create

  lpowers = [lall.pow[powclr], $
            lspiral.pow[powclr],$
            lellip.pow[powclr],$
            llum1.pow[powclr],$
            llum2.pow[powclr] ]
  lpowers_el = [lall.pow[powclr], $
                lspiral.pow[powclr],$
                lellip.pow[powclr],$
                llum1.pow[powclr],$
                llum2.pow[powclr] ]
  lpowers_eh = [lall.pow[powclr], $
                lspiral.pow[powclr],$
                lellip.pow[powclr],$
                llum1.pow[powclr],$
                llum2.pow[powclr] ]


END   
