PRO calc_m2l_omega_lumdens_onlyall, lumclr, $
                                    lallf, allf, m2l, m2lel, m2leh, $
                                    omegahm2, omegahm2_el, omegahm2_eh,$
                                    M0, M0el, M0eh, L0, L0el, L0eh,$
                                    lpowers,lpowersel,lpowerseh,$
                                    doplot=doplot, dofits=dofits, fourbin=fourbin, $
                                    fixlumnorm=fixlumnorm,usefixlum=usefixlum,$
                                    lumextno=lumextno, mextno=mextno,$
                                    indir=indir,mindir=mindir,$
                                    meanlum=meanlum, mclr=mclr

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: calc_m2l_omega_lumdens, lumclr, lallf, allf, m2l, m2lel, m2leh, '
      print,'                omegahm2, omegahm2_el, omegahm2_eh,'
      print,'                M0, M0el, M0eh, L0, L0el, L0eh,'
      print,'                lpowers,lpowersel,lpowerseh,'
      print,'                doplot=doplot, fourbin=fourbin, '
      print,'                fixlumnorm=fixlumnorm,usefixlum=usefixlum,'
      print,'                lumextno=lumextno, mextno=mextno'
      print,'                indir=indir,mindir=mindir'
      return
  ENDIF 

  ;; /masserror uses error in mass slope

  delvarx,lpowers,lpowersel,lpowerseh

  IF NOT keyword_set(fourbin) THEN ttstr='twobin' ELSE ttstr=''

  lumfix = [1.314430862, 1.235023618, 1.194514446, 1.227677096, 1.220581657]
  IF n_elements(meanlum) EQ 0 THEN $
    meanlum = [0.784103, 0.889751, 1.51780, 2.07778, 2.58957]
  ;;errfac=1.32
  errfac=1.0

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

;  outps = indir+'lumdense'+ttstr+'_m2l_omega_onlyall_'+!colors[lumclr]+psend
;  outfits = indir+'lumdense'+ttstr+'_m2l_omega_onlyall_'+!colors[lumclr]+fitend

  outfits = repstr(lallf, '.fit', '_m2l_omega_onlyall.fit')
  outps = repstr(lallf, '.fit', '_m2l_omega_onlyall.ps')
  outfits = repstr(outfits, 'rw', !colors[lumclr]+'w')
  outps = repstr(outps, 'rw', !colors[lumclr]+'w')
  print
  print,'Output ps file: ',outps
  print,'Output fits file: ',outfits
  print

  IF keyword_set(doplot) THEN begplot, name=outps, /color
  setup_mystuff

  ;; !!! kludge to use new lum stuff E.S.S.

  print,'Reading: ',lallf
  lall = mrdfits(lallf,1)

  print,'Reading: ',allf
  all = mrdfits(allf,1)

  ;lall.pow[lumclr] = -0.73
  ;lall.norm[lumclr] = 2.48


  ;; mass ranges
  range2error,all.powlow[0], all.power, all.powhigh[0],$
    mpeh, mpel
  
  print
  print,'Mass only slope: '+$
    ntostr(all.power)+' + '+ntostr(mpeh)+' - '+ntostr(mpel)
  print

  ;; account for 0.1 l_star cut 
  IF keyword_set(fixlumnorm) THEN BEGIN
      fixlum, lall, lumclr
  ENDIF 

  all.sigmaerr = all.sigmaerr*errfac

  normrange = [0, 25]
  powrange=[0.3, 1.4]
  norm=arrscl( findgen(400), normrange[0], normrange[1] )
  power=arrscl( findgen(400), -powrange[0], -powrange[1] )

  xtitle=!kpcxtitle2
  ytitle = !sigytitle
  xx=arrscl( findgen(100), min(all.meanr), max(all.meanr) )
;  IF display_type() EQ 'X' THEN clear=0 ELSE clear=1


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Range in radius to fit
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nr=n_elements(all.meanr)
  minrbin = 0
  wr=( lindgen(nr) )[minrbin:nr-1]

  ;!p.multi=[0,1,2]
  aspect=1
  pow_chisq_conf_1par, all.meanr[wr]/1000., all.sigma[wr], all.sigmaerr[wr], $
    norm, lall.pow[lumclr], chisq_surf, $
    allnorm, allnormlow, allnormhigh, chisq_per=chisq_per,$
    aspect=aspect,/center, yfit=yfit
  legend,['all power = '+ntostr(lall.pow[lumclr]),$
          'fit norm  = '+ntostr(allnorm),$
          !tsym.chi+'!U2!N/'+!tsym.nu+' = '+ntostr(chisq_per)],/clear

  IF display_type() EQ 'X' THEN key=get_kbrd(1)
  
;  plotboth_density_contrast, all, /logbin, /center, $
;                             xfit=all.meanr[wr], yfit=yfit, $
;                             xwindow1=xwindow1, ywindow1=ywindow1, clr=mclr,$
;                             yrange2=[-10,10];,yfit=all.orthosig[wr],fitsym=4
  plot_density_contrast, all, /logbin, /center, aspect=1
  oplot, all.meanr[wr], yfit
  legend, ['mass: '+!colorsp[mclr], 'lum: '+!colorsp[lumclr]], /right,box=0
  add_labels, xtickv=[10,20]

;  oploterror, lall.meanr, lall.ilumdense*(allnorm/lall.norm[lumclr]), $
;              lall.ilumdenserr*(allnorm/lall.norm[lumclr]), $
;              psym=4


  ;ploterror, all.meanr, all.sigma, all.sigmaerr,psym=1,$
  ;  xtitle=xtitle,ytitle=ytitle
  ;oplot, xx, yfit
  ;oplot, [0, 10000], [0,0]
  ;add_labels, xtickv=[10,20]

  IF display_type() EQ 'X' THEN key=get_kbrd(1)

  ;; now calculate global mass to light and error

  ;;;;;;;;;;;;;;;;;;;;
  ;; all galaxies
  ;;;;;;;;;;;;;;;;;;;;

  ;; make power positive number
  pow = -lall.pow[lumclr] & powlow=-lall.powlow[lumclr] & powhigh=-lall.powhigh[lumclr]
  norm = lall.norm[lumclr] & normlow=lall.normlow[lumclr] & normhigh=lall.normhigh[lumclr]

  ;; must switch low-high since we put pow=-pow
  range2error, powhigh, pow, powlow, lpeh, lpel ;power errors

  totpeh = sqrt(mpeh^2 + lpeh^2)
  totpel = sqrt(mpel^2 + lpel^2)

  range2error, normlow, norm, normhigh, lnel, lneh ;lum norm errors
  range2error, allnormlow[0], allnorm, allnormhigh[0], nel, neh ;mass norm errors
  add_arrval, pow, lpowers & add_arrval, abs(lpel), lpowersel 
  add_arrval, abs(lpeh), lpowerseh
  
  ;; fac takes care of conversions
  fac = (2.-pow)/pow
  m2l = (allnorm*1.e12)/(lall.norm[lumclr]*1.e10)*fac
  m2lel = m2l*sqrt( (nel/allnorm)^2 + (lnel/norm)^2 + 4./(2.-pow)^2*(lpel/pow)^2 )
  m2leh = m2l*sqrt( (neh/allnorm)^2 + (lneh/norm)^2 + 4./(2.-pow)^2*(lpeh/pow)^2 )
  m2lel_merr = $
    m2l*sqrt( (nel/allnorm)^2 + (lnel/norm)^2 + 4./(2.-pow)^2*(totpel/pow)^2 )
  m2leh_merr = $
    m2l*sqrt( (neh/allnorm)^2 + (lneh/norm)^2 + 4./(2.-pow)^2*(totpeh/pow)^2 )

  ;; Calculate M0 and L0
  M0 = 2.0*!pi/pow*allnorm   ;10^12 M_{\sun}
  M0el = M0*sqrt( (nel/allnorm)^2 + (lpel/pow)^2 )
  M0eh = M0*sqrt( (neh/allnorm)^2 + (lpeh/pow)^2 )
  M0el_merr = M0*sqrt( (nel/allnorm)^2 + (totpel/pow)^2 )
  M0eh_merr = M0*sqrt( (neh/allnorm)^2 + (totpeh/pow)^2 )

  L0 = 2.0*!pi/(2.0-pow)*norm ;10^10 L_{\sun}
  L0el = L0*sqrt( (lnel/norm)^2 + (lpel/(2.-pow))^2 )
  L0eh = L0*sqrt( (lneh/norm)^2 + (lpeh/(2.-pow))^2 )
  L0el_merr = L0*sqrt( (lnel/norm)^2 + (totpel/(2.-pow))^2 )
  L0eh_merr = L0*sqrt( (lneh/norm)^2 + (totpeh/(2.-pow))^2 )

  print,'relative errors'
  print,'norm        lum norm        power'
  print,abs(neh/allnorm),abs(lneh/norm),abs(2./(2.-pow)*(lpeh/pow))


  print
  print
  print,!colors[lumclr]+'-band'
  print,'Delta Sigma norm: '+ntostr(allnorm)
  print,'Lum norm: '+ntostr(norm)
  print,'all M0 = ',ntostr(M0),' + ',ntostr(M0eh),' - ',ntostr(M0el)
  print,'all L0 = ',ntostr(L0),' + ',ntostr(L0eh),' - ',ntostr(L0el)

  print
  print,'All M/L: ',ntostr(m2l),' + '+ntostr(m2leh)+' - '+ntostr(m2lel)
  print

  ;; Characteristic Radius, where M/L falls to a fraction f of the
  ;; global value
  ;; Need to get mean luminosities



  types = 'all'

  rmin = min(all.rmin_act)
  rmax = all.rmax_act

  mass_density = all.sigma*fac*1.e12 ;solar per Mpc^2
  mass_density_err = all.sigmaerr*fac*1.e12

  area = !pi*( (rmax/1000.)^2 - (rmin/1000.)^2 )
  massap = all.tsigma*fac*1.e12*area
  massaperr = all.tsigmaerr*fac*1.e12*area

  masstot = all.tsigma*fac*1.e12*!pi*(rmax/1000.)^2
  masstoterr = all.tsigmaerr*fac*1.e12*!pi*(rmax/1000.)^2

  lumdensity = norm*( rmax/1000.)^(-pow)

  ;; lum from model
  lumapmod = L0*( (rmax/1000.)^(2.-pow) - (rmin/1000.)^(2.-pow) )*1.e10
  lumtotmod = lumapmod + meanlum[lumclr]*1.e10
;  lumtotmod =  meanlum[lumclr]*1.e10 + $
;    L0*( (rmax/1000.)^(2.-pow) )*1.e10

  ;; lum right from data (missing inner 20 kpc)
  command = 'intlumdense = interpol(lall.'+$
    !colors[lumclr]+'tlumdense, lall.rmax, rmax)'
  IF NOT execute(command) THEN message,'Doh!'

  lumap = intlumdense*area*1.e10
  lumtot =  meanlum[lumclr]*1.e10 + lumap

  ;colprint,rmax,lumtotmod,lumtot


  struct = create_struct('types',types,$
                         'pow', pow,$
                         'm0',M0,$
                         'm0_eh',M0eh,$
                         'm0_el',M0el,$
                         'l0',L0,$
                         'l0_eh',L0eh,$
                         'l0_el',L0el,$
                         'm2l',m2l,$
                         'm2l_el',m2lel,$
                         'm2l_eh',m2leh,$
                         'meanr',all.meanr,$
                         'mass_density', mass_density,$ ;solar mass per Mpc^2
                         'mass_density_err', mass_density_err,$
                         'rmin',rmin,$
                         'rmax', rmax,$
                         'massap', massap,$
                         'massaperr',massaperr,$
                         'masstot',masstot,$
                         'masstoterr',masstoterr,$
                         'lumdensity',lumdensity,$
                         'lumap',lumap,$ ;solar lum per Mpc^2
                         'lumtot', lumtot,$
                         'meanlum', meanlum[lumclr])

  nnn = n_elements(types)

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
  rho_crit = double(2.77545e-7) ;Msolar/pc^3
  rho_crit = rho_crit*(10d)^18. ; Msolar/Mpc^3


  ;; Here use h=1

  ;;;;;;;;;
  ;; I think there is no h in our measurement of omega, so
  ;; just use this one which has herr=0.

  omega = m2l*ld[lumclr]/rho_crit
  omega_eh = omega*sqrt( (lderr[lumclr]/ld[lumclr])^2 + $
                               (m2leh/m2l)^2)
  omega_el = omega*sqrt( (lderr[lumclr]/ld[lumclr])^2 + $
                               (m2lel/m2l)^2)
  omega_eh_merr = omega*sqrt( (lderr[lumclr]/ld[lumclr])^2 + $
                               (m2leh_merr/m2l)^2)
  omega_el_merr = omega*sqrt( (lderr[lumclr]/ld[lumclr])^2 + $
                               (m2lel_merr/m2l)^2)

  struct = create_struct(struct, $
                         'omega',omega,$
                         'omega_el',omega_el,$
                         'omega_eh',omega_eh)

  print
  print,!colors[lumclr]+'-band'
  print,'All Omega Matter : ',ntostr(rnd(omega[0],3),5),' + '+$
    ntostr(rnd(omega_eh[0],3),5)+'('+ntostr(rnd(omega_eh_merr[0],3),5)+') - '+$
    ntostr(rnd(omega_el[0],3),5)+'('+ntostr(rnd(omega_el_merr[0],3),5)+')'

  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; H0 = 72+/-8 (Freedman, et al. 2001
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  h = 0.72
  herr = 0.08
  pomega = 0.2/h
  pomegaerr = pomega*sqrt( (0.03/0.2)^2 + (herr/h)^2 )
  print,'Peacock omega = ',ntostr(pomega),' +/- ',ntostr(pomegaerr)
  print

  ;;IF keyword_set(doplot) THEN endplot

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; calculate characteristic radius where
  ;; mass to light reaches fraction f of the
  ;; value at infinity
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  f = 0.5
  ml = meanlum[lumclr]

  Rs = ( f/(1.-f) * ml/L0 )^(1./(2.-pow))*1000.
  Rserr1 = Rs*sqrt( 1./(2.-pow)^2*(L0el/L0)^2 + $
                    (alog(ml/L0))^2/(2.-pow)^4 * lpel^2 )
  Rserr2 = Rs*sqrt( 1./(2.-pow)^2*(L0eh/L0)^2 + $
                    (alog(ml/L0))^2/(2.-pow)^4 * lpeh^2 )
  Rserr3 = Rs*sqrt( 1./(2.-pow)^2*(L0el/L0)^2 + $
                    (alog(ml/L0))^2/(2.-pow)^4 * lpeh^2 )
  Rserr4 = Rs*sqrt( 1./(2.-pow)^2*(L0eh/L0)^2 + $
                    (alog(ml/L0))^2/(2.-pow)^4 * lpel^2 )

  print,'Different errors: ',Rserr1, Rserr2, Rserr3, Rserr4
  Rserr = max([Rserr1,Rserr2,Rserr3,Rserr4])
  
  print,'Rs = '+ntostr(Rs) + ' +/- '+ntostr(Rserr)
  
  ;; mass2light function
  rratio = struct.rmax/Rs
  pp=2.-struct.pow
  m2lfunc = struct.M0/struct.L0*( rratio^pp/(1. + rratio^pp) ) 

  struct = create_struct(struct, $
                         'Rs', Rs,$
                         'Rserr', Rserr,$
                         'm2lfunc', m2lfunc)

  plotlumdensm2l, struct

  IF keyword_set(doplot) THEN endplot

  ;;;;;;;;;;;;;;;;;;;
  ;; make some tex
  ;;;;;;;;;;;;;;;;;;;

  IF lumclr EQ 2 THEN tstr='All' ELSE tstr = '   '
  m2lmes=tstr+' & '+!colors[lumclr]+'$^{\prime}$ & '+$
    ntostr(rnd(pow,3), 5)+$
    '$\pm '+ntostr(rnd(lpeh,3), 5)+'$' + $
    ' & '+$
    ntostr(rnd(L0,1), 4)+$
    '$^{+'+ntostr(rnd(L0eh,1), 3)+'}' + $
    '_{-'+ntostr(rnd(L0el,1), 3)+'}$' + $
    ' & '+$
    ntostr(long(rnd(Rs)))+$
    '$\pm '+ntostr(long(rnd(Rserr)))+'$' + $
    ' & '+$
    ntostr(rnd(M0,1), 4)+$
    '$^{+'+ntostr(rnd(M0eh,1), 3)+'}' + $
    '_{-'+ntostr(rnd(M0el,1), 3)+'}$' + $
    ' & '+$
    ntostr(long(rnd(m2l)), 5)+$
    '$^{+'+ntostr(long(rnd(m2leh)), 5)+'}' + $
    '_{-'+ntostr(long(rnd(m2lel)), 5)+'}$' + $
    ' & '+$
    ntostr(rnd(omega,2),4)+$
    '$ \pm '+ntostr(rnd(omega_eh,2), 4)+$
    ' ('+ntostr(rnd(omega_eh_merr,2),4)+') $'+' \\'

  print
  print,m2lmes


  IF keyword_set(dofits) THEN BEGIN
      print,'Output fits file: ',outfits  
      mwrfits, struct, outfits, /create
  ENDIF 
  return

  lpowers = [lall.pow[lumclr], $
            lspiral.pow[lumclr],$
            lellip.pow[lumclr],$
            llum1.pow[lumclr],$
            llum2.pow[lumclr] ]
  lpowers_el = [lall.pow[lumclr], $
                lspiral.pow[lumclr],$
                lellip.pow[lumclr],$
                llum1.pow[lumclr],$
                llum2.pow[lumclr] ]
  lpowers_eh = [lall.pow[lumclr], $
                lspiral.pow[lumclr],$
                lellip.pow[lumclr],$
                llum1.pow[lumclr],$
                llum2.pow[lumclr] ]


END   
