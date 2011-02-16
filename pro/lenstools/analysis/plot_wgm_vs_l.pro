PRO plot_wgm_vs_l_wuse, t1, t2, t3, t4, $
                        wuse1, wuse2, wuse3, wuse4, $
                        wusereb1, wusereb2, wusereb3, wusereb4

  COMMON wgmblock, lumclr, nbin, nodojack, cmess, bincolors,rbincolors, clr1, clr2, clr3, clr4, wgmfac, rt, sigt, siget, ytitle, yrange

  frac = 0.10
  wuse1 = where(t1.lumfrac LE frac AND t1.gmrfrac LE frac)
  wuse2 = where(t2.lumfrac LE frac AND t2.gmrfrac LE frac)

  minr = min(t1.meanr[wuse1])
  wusereb1 = where(t1.meanr_rebin GE minr)

  minr = min(t2.meanr[wuse2])
  wusereb2 = where(t2.meanr_rebin GE minr)

  IF nbin GT 2 THEN BEGIN 
      wuse3 = where(t3.lumfrac LE frac AND t3.gmrfrac LE frac)

      minr = min(t3.meanr[wuse3])
      wusereb3 = where(t3.meanr_rebin GE minr)
  ENDIF 
  
  IF nbin GT 3 THEN BEGIN 
      wuse4 = where(t4.lumfrac LE frac AND t4.gmrfrac LE frac)

      minr = min(t4.meanr[wuse4])
      wusereb4 = where(t4.meanr_rebin GE minr)
  ENDIF 

  help,wuse1,wuse2,wuse3,wuse4

return
  wuse1 = lindgen(18) & wuse2=wuse1 & wuse3=wuse1 & wuse4=wuse1
  wusereb1 = lindgen(9) & wusereb2=wusereb1 & $
    wusereb3=wusereb1 & wusereb4=wusereb1
return

  IF NOT keyword_set(rebin) THEN BEGIN 

      wuse1 = (lindgen(18))[2:17]
      wuse2 = (lindgen(18))[2:17]
      wuse3 = (lindgen(18))[2:17]
      IF lumclr EQ 4 THEN BEGIN 
          wuse4 = (lindgen(18))[2:13]
      ENDIF ELSE IF lumclr EQ 0 THEN BEGIN 
          wuse4 = (lindgen(18))[1:18]
      ENDIF ELSE BEGIN 
          wuse4 = (lindgen(18))[0:11]
      ENDELSE 

  ENDIF ELSE BEGIN 

      wuse1 = (lindgen(9))[1:8]
      wuse2 = (lindgen(9))[1:8]
      wuse3 = (lindgen(9))[1:8]
      IF lumclr EQ 4 THEN BEGIN 
          wuse4 = (lindgen(9))[1:6]
      ENDIF ELSE IF lumclr EQ 0 THEN BEGIN 
          wuse4 = (lindgen(9))[1:8]
      ENDIF ELSE BEGIN 
          wuse4 = (lindgen(9))[0:5]
      ENDELSE 

  ENDELSE 

END 

PRO plot_wgm_vs_l_readfiles,t1,t2,t3,t4,wgm=wgm,rebin=rebin

  COMMON wgmblock, lumclr, nbin, nodojack, cmess, bincolors,rbincolors, clr1, clr2, clr3, clr4, wgmfac, rt, sigt, siget, ytitle, yrange

  IF nodojack THEN jackstr='' ELSE jackstr = '_jack'

  CASE nbin OF 
      2: BEGIN 
          lumstr = 'twobin'
      END 
      3: BEGIN 
          lumstr = 'threebin'
      END 
      4: BEGIN 
          lumstr = 'fourbin'
      END 
      ELSE: message,'What!'
  ENDCASE 

  dir = '/net/cheops2/home/esheldon/lensout/combstripe/'
  pardir = dir + 'fits/'
  combdir = dir + 'comb/'

  f1 = combdir+'sublum/'+!colors[lumclr]+'/lum1'+lumstr+'_zgal_gal_stripe10_11_12_35_36_37_ri'+jackstr+'_comb_N4.fit'
  f2 = combdir+'sublum/'+!colors[lumclr]+'/lum2'+lumstr+'_zgal_gal_stripe10_11_12_35_36_37_ri'+jackstr+'_comb_N4.fit'
  f3 = combdir+'sublum/'+!colors[lumclr]+'/lum3'+lumstr+'_zgal_gal_stripe10_11_12_35_36_37_ri'+jackstr+'_comb_N4.fit'
  f4 = combdir+'sublum/'+!colors[lumclr]+'/lum4'+lumstr+'_zgal_gal_stripe10_11_12_35_36_37_ri'+jackstr+'_comb_N4.fit'

  print,f1
  t1 = mrdfits(f1, 1)
  print,f2
  t2 = mrdfits(f2, 1)

  IF nbin GT 2 THEN BEGIN
      print,f3
      t3 = mrdfits(f3, 1)
  ENDIF 
  IF nbin GT 3 THEN BEGIN
      print,f4
      t4 = mrdfits(f4, 1)
  ENDIF 

  tn = tag_names(t1)
  IF keyword_set(wgm) THEN BEGIN 
      ytitle = !wgmytitle
      delvarx, wgmfac
      add_arrval, t1.wgm[0]/t1.sigma[0], wgmfac
      t1.norm = t1.norm*wgmfac[0]

      add_arrval, t2.wgm[0]/t2.sigma[0], wgmfac
      t2.norm = t2.norm*wgmfac[1]

      IF nbin GT 2 THEN BEGIN 
          add_arrval, t3.wgm[0]/t3.sigma[0], wgmfac
          t3.norm = t3.norm*wgmfac[2]
      ENDIF 
      IF nbin GT 3 THEN BEGIN 
          add_arrval, t4.wgm[0]/t4.sigma[0], wgmfac
          t4.norm = t4.norm*wgmfac[3]
      ENDIF 

      IF keyword_set(rebin) THEN BEGIN 
          rt = where(tn EQ 'MEANR_REBIN')
          sigt = where(tn EQ 'WGM_REBIN')
          siget = where(tn EQ 'WGMERR_REBIN')
      ENDIF ELSE BEGIN 
          rt = where(tn EQ 'MEANR')
          sigt = where(tn EQ 'WGM')
          siget = where(tn EQ 'WGMERR')
      ENDELSE 

      ytitle=!wgmytitle
      yrange = [0.05,5.e3]*max(wgmfac)

  ENDIF ELSE BEGIN 
      ytitle = !deltaytitle
      wgmfac = replicate(1.0, nbin)

      IF keyword_set(rebin) THEN BEGIN 
          rt = where(tn EQ 'MEANR_REBIN')
          sigt = where(tn EQ 'SIGMA_REBIN')
          siget = where(tn EQ 'SIGMAERR_REBIN')
      ENDIF ELSE BEGIN 
          rt = where(tn EQ 'MEANR')
          sigt = where(tn EQ 'SIGMA')
          siget = where(tn EQ 'SIGMAERR')
      ENDELSE 

      ytitle=!deltaytitle
      yrange=[0.05,5.e3]

  ENDELSE 

END 

PRO plot_wgm_vs_l_getarr, t1, t2, t3, t4, $
                          norm, normerr, power, powerr, $
                          meanlum, meanlumerr, absmag

  COMMON wgmblock, lumclr, nbin, nodojack, cmess, bincolors,rbincolors, clr1, clr2, clr3, clr4, wgmfac, rt, sigt, siget, ytitle, yrange

  delvarx,cmess,norm, normerr, power, powerr, $
          meanlum, meanlumerr, absmag

  ;; for absolute magnitudes
  sun=[6.38,5.06,4.64,4.53,4.52]

  FOR i=1L, nbin DO BEGIN 

      tstr = 't'+ntostr(i)

      comstr = 'add_arrval, '+tstr+'.tmeanlum, meanlum'
      IF NOT execute(comstr) THEN message,'Doh'
      comstr = 'add_arrval, '+tstr+'.tmeanlumerr, meanlumerr'
      IF NOT execute(comstr) THEN message,'Doh'

      comstr = 'add_arrval, '+tstr+'.norm, norm'
      IF NOT execute(comstr) THEN message,'Doh'
      valstr = 'val = max(['+tstr+'.norm-'+tstr+'.normlow[0], '+$
                             tstr+'.normhigh[0]-'+tstr+'.norm])'
      comstr = valstr + ' & '+'add_arrval, val, normerr'
      IF NOT execute(comstr) THEN message,'Doh'


      comstr = 'add_arrval, '+tstr+'.power, power'
      IF NOT execute(comstr) THEN message,'Doh'
      valstr = 'val = max(['+tstr+'.power-'+tstr+'.powlow[0], '+$
                             tstr+'.powhigh[0]-'+tstr+'.power])'
      comstr = valstr + ' & '+'add_arrval, val, powerr'
      IF NOT execute(comstr) THEN message,'Doh'


      tabsmag = sun[lumclr] - 2.5*alog10(meanlum[i-1]*1.e10)

      val = 'L('+!colors[lumclr]+') = '+ntostr(meanlum[i-1],4,/round)+$
           ' M('+!colors[lumclr]+') = '+ntostr(tabsmag,6,/round)
      add_arrval, val, cmess, /front

  ENDFOR 

  ;; absolute magnitudes
  absmag = sun[lumclr] - 2.5*alog10(meanlum*1.e10)

  colprint,norm,$
           normerr, $
           power, $
           powerr, $
           meanlum,$
           meanlumerr, $
           absmag

END 

PRO plot_wgm_vs_l_plotstuff, t1, t2, t3, t4, $
                             wuse1=wuse1, wuse2=wuse2, $
                             wuse3=wuse3,wuse4=wuse4

  COMMON wgmblock, lumclr, nbin, nodojack, cmess, bincolors,rbincolors, clr1, clr2, clr3, clr4, wgmfac, rt, sigt, siget, ytitle, yrange

  nr = n_elements(t1.(rt))
  IF n_elements(wuse1) EQ 0 THEN BEGIN 
      wuse1 = lindgen(nr)
  ENDIF
  IF n_elements(wuse2) EQ 0 THEN BEGIN 
      wuse2 = lindgen(nr)
  ENDIF
  IF n_elements(wuse3) EQ 0 THEN BEGIN 
      wuse3 = lindgen(nr)
  ENDIF
  IF n_elements(wuse4) EQ 0 THEN BEGIN 
      wuse4 = lindgen(nr)
  ENDIF

  aploterror,1,t1.(rt)[wuse1],t1.(sigt)[wuse1],t1.(siget)[wuse1],$
             psym=8,/xlog,/ylog, $
             yrange=yrange,/yst, $
             xtitle=!kpcxtitle, ytitle=ytitle, /center

  oploterror, t1.(rt)[wuse1],t1.(sigt)[wuse1],t1.(siget)[wuse1],$
              psym=8, color=clr1, errc=clr1
  yfit1 = t1.norm*(t1.(rt)[wuse1]/1000.)^t1.power
  oplot, t1.(rt)[wuse1],yfit1,color=clr1

  oploterror, t2.(rt)[wuse2],t2.(sigt)[wuse2],t2.(siget)[wuse2],$
              psym=8, color=clr2, errc=clr2
  yfit2 = t2.norm*(t2.(rt)[wuse2]/1000.)^t2.power
  oplot, t2.(rt)[wuse2],yfit2,color=clr2

  IF nbin GT 2 THEN BEGIN 
      oploterror, t3.(rt)[wuse3],t3.(sigt)[wuse3],t3.(siget)[wuse3],$
                  psym=8, color=clr3, errc=clr3
      yfit3 = t3.norm*(t3.(rt)[wuse3]/1000.)^t3.power
      oplot, t3.(rt)[wuse3],yfit3,color=clr3
  ENDIF 
  IF nbin GT 3 THEN BEGIN 
      oploterror, t4.(rt)[wuse4],t4.(sigt)[wuse4],t4.(siget)[wuse4],$
                  psym=8, color=clr4, errc=clr4
      yfit4 = t4.norm*(t4.(rt)[wuse4]/1000.)^t4.power
      oplot, t4.(rt)[wuse4],yfit4,color=clr4

  ENDIF 

  legend, cmess,$
          color=rbincolors,textc=rbincolors,/right,box=0,$
          charsize=1, thick=replicate(!p.thick,nbin);, line=cline

END 

PRO plot_wgm_vs_l_plotfill, t1, t2, t3, t4, $
                            wuse1=wuse1, wuse2=wuse2, wuse3=wuse3,wuse4=wuse4

  COMMON wgmblock, lumclr, nbin, nodojack, cmess, bincolors,rbincolors, clr1, clr2, clr3, clr4, wgmfac, rt, sigt, siget, ytitle, yrange


  nr = n_elements(t1.meanr)
  IF n_elements(wuse1) EQ 0 THEN BEGIN 
      wuse1 = lindgen(nr)
  ENDIF
  IF n_elements(wuse2) EQ 0 THEN BEGIN 
      wuse2 = lindgen(nr)
  ENDIF
  IF n_elements(wuse3) EQ 0 THEN BEGIN 
      wuse3 = lindgen(nr)
  ENDIF
  IF n_elements(wuse4) EQ 0 THEN BEGIN 
      wuse4 = lindgen(nr)
  ENDIF

  aplot,1,t1.meanr[wuse1],t1.yfit[wuse1]*wgmfac[0], /xlog,/ylog, $
        yrange=yrange,/ystyle, $
        xtitle=!kpcxtitle, ytitle=ytitle, /center
  
  x = [t1.meanr[wuse1],reverse(t1.meanr[wuse1])]
  y = [t1.yallow_low[wuse1], reverse(t1.yallow_high[wuse1])]*wgmfac[0]
  polyfill, x, y, color=clr1,/data
  
  x = [t2.meanr[wuse1],reverse(t2.meanr[wuse2])]
  y = [t2.yallow_low[wuse2], reverse(t2.yallow_high[wuse2])]*wgmfac[1]
  polyfill, x, y, color=clr2,/data
  
  IF nbin GT 2 THEN BEGIN 
      x = [t3.meanr[wuse1],reverse(t3.meanr[wuse3])]
      y = [t3.yallow_low[wuse3], reverse(t3.yallow_high[wuse3])]*wgmfac[2]
      polyfill, x, y, color=clr3,/data
  ENDIF 

  IF nbin GT 3 THEN BEGIN 
      x = [t4.meanr[wuse4],reverse(t4.meanr[wuse4])]
      y = [t4.yallow_low[wuse4], reverse(t4.yallow_high[wuse4])]*wgmfac[3]
      polyfill, x, y, color=clr4,/data
  ENDIF 

  legend, cmess,$
          color=rbincolors,textc=rbincolors,/right,box=0,$
          charsize=1, thick=replicate(!p.thick,nbin) ;, line=cline


END 


PRO plot_wgm_vs_l_fitpower, t, usewgmfac, wuse=wuse

  COMMON wgmblock, lumclr, nbin, nodojack, cmess, bincolors,rbincolors, clr1, clr2, clr3, clr4, wgmfac, rt, sigt, siget, ytitle, yrange

  clevel = 0
  nr = n_elements(t.meanr)
  
  IF n_elements(wuse) EQ 0 THEN BEGIN 
      wuse = lindgen(nr)
  ENDIF

  charsize = 1

;  print
;  print,'Initial power law fitting'
  fitpower, t.meanr[wuse]/1000., t.sigma[wuse]*usewgmfac, $
            t.sigmaerr[wuse]*usewgmfac, $
            [5.*usewgmfac, -.8], tyfit, Aout, Asig,/silent

  rangefac = 10.0; ELSE rangefac=5.0
  Asig[1] = Asig[1]*rangefac > 0.3 < 2.0
  Asig[0] = Asig[0]*rangefac > 2.0 < 30.0

  spowrange = Aout[1] + [-Asig[1], +Asig[1]]
  normrange = Aout[0] + [Asig[0], -Asig[0]]

  print
  print,'Power law fits'
  
  nnorm = 200
  npow = 200
  
  pow_chisq_conf_gen, t.meanr/1000., t.sigma*usewgmfac, t.covariance*usewgmfac^2, $
                      spowrange, normrange, npow, nnorm, $
                      chisq_surf, $
                      bestpow, bestnorm, $
                      powlow, powhigh, $
                      normlow, normhigh, $
                      xtitle='Index', $
                      ytitle='Norm [h M'+sunsymbol()+$
                      ' pc!U'+!csym.minus+'2!N]', $
                      xstyle=1, ystyle=1, charsize=1.0, $
                      aspect=1, /center, $
                      yfit=combsigdiff_fit, $
                      minchisq=minchisq, degfree=degfree, $
                      yallow_low = sigallow_low, $
                      yallow_high=sigallow_high,wuse=wuse

  nhighdiff = normhigh[clevel]-bestnorm & nlowdiff = bestnorm - normlow[clevel]
  mess = 'Norm = '+ntostr(bestnorm,5,/round)+$
    '!S!U'+!csym.plus+ntostr(nhighdiff,5,/round)+'!R!D'+!csym.minus+ntostr(nlowdiff,5,/round)
  mess = [mess, '']
  phighdiff = powhigh[clevel]-bestpow & plowdiff = bestpow - powlow[clevel]
  mess = [mess, $
          'Index = '+ntostr(bestpow, 5,/round) + $
          '!S!U'+!csym.plus+ntostr(phighdiff,5,/round)+'!R!D'+!csym.minus+ntostr(plowdiff,5,/round)]
  legend, mess, /right, charsize=charsize, box=0
  
  mess2 = !csym.chi+'!U2!N/'+!csym.nu+' = '+ntostr(minchisq,4,/round)+'/'+ntostr(degfree)+$
    ' = '+ntostr(minchisq/degfree,4,/round)
  legend, mess2,charsize=charsize, box=0
  
  print,'Best fit Power Index: '+ntostr(bestpow)+'+'+ntostr(phighdiff)+'-'+ntostr(plowdiff)
  print,'Best fit Power Norm: '+ntostr(bestnorm)+'+'+ntostr(nhighdiff)+'-'+ntostr(nlowdiff)


;  IF !d.name eq 'X' THEN key=get_kbrd(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; <Sigma>-sigma plots: lower panel
  ;; is orthisig
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  to = prange(t.orthosig,t.orthosigerr,/symm,/slack)
;  omax = 2.0*sdev(t.orthosig) < abs(to[0])
;  obotyrange = [-omax,omax]
;  plotboth_density_contrast, t, /logbin, $
;                             xfit=t.meanr[wuse], $
;                             yfit=combsigdiff_fit[wuse],$
;                             botyrange=obotyrange, /center,charsize=1
  
  t.norm = bestnorm
  t.power = bestpow
  t.normlow = normlow
  t.normhigh = normhigh
  t.powlow = powlow
  t.powhigh = powhigh

  t.yfit = combsigdiff_fit
  t.yallow_low = sigallow_low
  t.yallow_high = sigallow_high

  !p.multi=0

END 

PRO plot_wgm_vs_l_fitnorm, t, usewgmfac, onorm, nerrhigh, nerrlow,mess=mess, wuse=wuse

  COMMON wgmblock, lumclr, nbin, nodojack, cmess, bincolors,rbincolors, clr1, clr2, clr3, clr4, wgmfac, rt, sigt, siget, ytitle, yrange

  leg_charsize = 0.7

  fac = 3.0
  nr = n_elements(t.meanr)

  IF n_elements(wuse) EQ 0 THEN BEGIN 
      wuse = lindgen(nr)
  ENDIF 

  norm = t.norm
  power = -0.85
  range2error, t.normlow, norm, t.normhigh, errlow, errhigh
  nrange = norm + [-errlow[0]*fac, errhigh[0]*fac]
  nnorm = 500

  xtitle='Norm [h M'+sunsymbol()+' pc!U'+!csym.minus+'2!N]'
  ytitle=!csym.delta_cap+!csym.chi+'!U2!N'
  pow_chisq_conf_1par_gen, t.meanr/1000., t.sigma*usewgmfac, t.covariance*usewgmfac^2, $
                           nrange, nnorm, power, $
                           chisq_surf, onorm, onormlow, onormhigh, $
                           xtitle=xtitle, ytitle=ytitle, charsize=1, $
                           wuse=wuse, minchisq=minchisq, degfree=degfree

  range2error, onormlow, onorm, onormhigh, nerrhigh, nerrlow

  IF n_elements(mess) NE 0 THEN BEGIN 
      legend, mess,/right, /clear, charsize=leg_charsize
  ENDIF 

  chmess = !csym.chi+'!U2!N/'+!csym.nu+' = '+$
    ntostr(minchisq,4,/round)+'/'+ntostr(degfree)+$
    ' = '+ntostr(minchisq/degfree,4,/round)
  legend, chmess,charsize=leg_charsize, /clear


END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Main
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO plot_wgm_vs_l, tlumclr, tnbin, t1, t2, t3, t4, nojack=nojack, dops=dops, rebin=rebin,$
                   wgm=wgm

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: plot_wgm_vs_l, lumclr, nbin, nojack=nojack, dops=dops, rebin=rebin,wgm=wgm'
      return
  ENDIF 

  COMMON wgmblock, lumclr, nbin, nodojack, cmess, bincolors,rbincolors, clr1, clr2, clr3, clr4, wgmfac, rt, sigt, siget, ytitle, yrange

  lumclr = tlumclr
  nbin = tnbin

  CASE nbin OF 
      2: BEGIN 
          lumstr = 'twobin'
      END 
      3: BEGIN 
          lumstr = 'threebin'
      END 
      4: BEGIN 
          lumstr = 'fourbin'
      END 
      ELSE: message,'What!'
  ENDCASE 

  IF keyword_set(nojack) THEN nodojack=1 ELSE nodojack=0
  IF nodojack THEN jackstr = '' ELSE jackstr = '_jack'

  IF keyword_set(rebin) THEN rebstr = '_rebin' ELSE rebstr = ''


  outdir = !lensout_dir + 'wgm_vs_l/'
  psfile = outdir + $
    'lum'+lumstr+'_wgm_vs_l_'+!colors[lumclr]+jackstr+rebstr+'.ps'
  IF keyword_set(dops) THEN begplot,name=psfile,/color

  !p.multi=0

  setup_mystuff
  IF !d.name EQ 'PS' THEN BEGIN 
      clr1 = !green
      clr2 = !red
      clr3 = !blue
      clr4 = !black
  ENDIF ELSE BEGIN 
      clr1 = !green
      clr2 = !red
      clr3 = !yellow
      clr4 = !cyan
  ENDELSE 

  bincolors = [clr1,clr2]

  IF nbin GT 2 THEN bincolors = [bincolors, clr3]
  IF nbin GT 3 THEN bincolors = [bincolors, clr4]
  rbincolors = reverse(bincolors)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; read the files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  plot_wgm_vs_l_readfiles,t1,t2,t3,t4,wgm=wgm,rebin=rebin

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; make the mean/err arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  plot_wgm_vs_l_getarr, t1, t2, t3, t4, $
                        norm, normerr, power, powerr, $
                        meanlum, meanlumerr, absmag

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plot stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  plot_wgm_vs_l_plotstuff, t1, t2, t3, t4


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; re-fit for norm, fixed index
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; skipping this
GOTO,jump

  plot_wgm_vs_l_wuse, t1, t2, t3, t4, wuse1, wuse2, wuse3, wuse4, $
                      wusereb1, wusereb2, wusereb3, wusereb4


  if !d.name eq 'X' then key=get_kbrd(1)
  plot_wgm_vs_l_fitpower, t1, wgmfac[0], wuse=wuse1
  if !d.name eq 'X' then key=get_kbrd(1)
  plot_wgm_vs_l_fitpower, t2, wgmfac[1], wuse=wuse2
  if !d.name eq 'X' then key=get_kbrd(1)

  IF nbin GT 2 THEN BEGIN 
      plot_wgm_vs_l_fitpower, t3, wgmfac[2], wuse=wuse3
      if !d.name eq 'X' then key=get_kbrd(1)
  ENDIF 
  IF nbin GT 3 THEN BEGIN 
      plot_wgm_vs_l_fitpower, t4, wgmfac[3], wuse=wuse4
      if !d.name eq 'X' then key=get_kbrd(1)
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; remake the arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  plot_wgm_vs_l_getarr, t1, t2, t3, t4, $
                        norm, normerr, power, powerr, $
                        meanlum, meanlumerr, absmag

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; replot stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(rebin) THEN BEGIN 
      plot_wgm_vs_l_plotstuff, t1, t2, t3, t4, $
                               wuse1=wusereb1, wuse2=wusereb2, $
                               wuse3=wusereb3, wuse4=wusereb4, $
                               /rebin
  ENDIF ELSE BEGIN 
      plot_wgm_vs_l_plotstuff, t1, t2, t3, t4, $
                               wuse1=wuse1, wuse2=wuse2, $
                               wuse3=wuse3, wuse4=wuse4
  ENDELSE 

jump:

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plot filled in allowed region
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  plot_wgm_vs_l_plotfill, t1, t2, t3, t4, $
                          wuse1=wuse1, wuse2=wuse2, wuse3=wuse3,wuse4=wuse4


  if !d.name eq 'X' then key=get_kbrd(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plot up the fits when we fix the norm
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  !p.multi=[0,2,2]
  xmold = !x.margin
  ymold = !y.margin

  !x.margin = [5.0, 1.5]
  !y.margin = [4.0, 1.0]


  mess = ['L('+!colors[lumclr]+') = '+ntostr(meanlum[0],4,/round), $
          'M('+!colors[lumclr]+') = '+ntostr(absmag[0],6,/round)]
  plot_wgm_vs_l_fitnorm, t1, wgmfac[0], norm1, nerrhigh1, nerrlow1, $
                         mess = mess

  add_arrval, norm1, lnorm
  add_arrval, max([nerrlow1, nerrhigh1]), lnormerr

  mess = ['L('+!colors[lumclr]+') = '+ntostr(meanlum[1],4,/round), $
          'M('+!colors[lumclr]+') = '+ntostr(absmag[1],6,/round)]
  plot_wgm_vs_l_fitnorm, t2, wgmfac[1], norm2, nerrhigh2, nerrlow2, $
                         mess = mess
  add_arrval, norm2, lnorm
  add_arrval, max([nerrlow2, nerrhigh2]), lnormerr

  IF nbin GT 2 THEN begin
      mess = ['L('+!colors[lumclr]+') = '+ntostr(meanlum[2],4,/round), $
              'M('+!colors[lumclr]+') = '+ntostr(absmag[2],6,/round)]

      plot_wgm_vs_l_fitnorm, t3, wgmfac[2], norm3, nerrhigh3, nerrlow3, $
                         mess = mess
      add_arrval, norm3, lnorm
      add_arrval, max([nerrlow3, nerrhigh3]), lnormerr

  ENDIF 
  IF nbin GT 3 THEN BEGIN 
      mess = ['L('+!colors[lumclr]+') = '+ntostr(meanlum[3],4,/round), $
              'M('+!colors[lumclr]+') = '+ntostr(absmag[3],6,/round)]
      plot_wgm_vs_l_fitnorm, t4, wgmfac[3], norm4, nerrhigh4, nerrlow4, $
                             mess = mess

      add_arrval, norm4, lnorm
      add_arrval, max([nerrlow4, nerrhigh4]), lnormerr
  ENDIF 

  !p.multi=0
  !x.margin = xmold
  !y.margin = ymold

  if !d.name eq 'X' then key=get_kbrd(1)

  IF nbin GT 2 THEN BEGIN 
      fitpower, meanlum, lnorm, lnormerr, [1.0, 1.0], tlyfit, Aout, Asig
      fac = 5.0
      nrange = Aout[0] + [-Asig[0]*fac, Asig[0]*fac]
      prange = Aout[1] + [-Asig[1]*fac, Asig[1]*fac]
      nn = 200

      !p.multi=[0,0,2]

      pow_chisq_conf_gen, meanlum, lnorm, lnormerr, prange, nrange, $
                          nn, nn, chi, pmin, nmin, $
                          perrlow=perrlow, nerrlow=nerrlow, $
                          perrhigh=perrhigh, nerrhigh=nerrhigh, $
                          yfit=lyfit, /dolegend, nkeep = [4, 4],$
                          xtitle = 'Index', $
                          ytitle='Norm', $
                          aspect=1, /center
      
      print,'Norm2  = '+ntostr(nmin)+'+'+ntostr(nerrhigh[0])+'-'+ntostr(nerrlow[0])
      print,'Index2 = '+ntostr(pmin)+'+'+ntostr(perrhigh[0])+'-'+ntostr(perrlow[0])
  ENDIF 

;  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  aploterror, 1, meanlum, lnorm, meanlumerr, lnormerr, psym=8, $
              xtitle = 'L('+!colors[lumclr]+') [10!U10!N L'+sunsymbol()+']',$
              ytitle = 'Norm [h M'+sunsymbol()+' pc!U'+!csym.minus+'2!N]',$
              /center
  oploterror, meanlum[0], lnorm[0], meanlumerr[0], lnormerr[0], $
              color=clr1, errc=clr1, psym=8
  oploterror, meanlum[1], lnorm[1], meanlumerr[1], lnormerr[1], $
              color=clr2, errc=clr2, psym=8
  IF nbin GT 2 THEN BEGIN 
      oploterror, meanlum[2], lnorm[2], meanlumerr[2], lnormerr[2], $
                  color=clr3, errc=clr3, psym=8
  ENDIF 
  IF nbin GT 3 THEN BEGIN 
      oploterror, meanlum[3], lnorm[3], meanlumerr[3], lnormerr[3], $
                  color=clr4, errc=clr4, psym=8
  ENDIF 

  IF nbin GT 2 THEN oplot, meanlum, lyfit
  legend, cmess,$
          color=rbincolors,textc=rbincolors,/left,box=0,$
          charsize=1, thick=replicate(!p.thick,nbin);, line=cline

  !p.multi=0

  IF keyword_set(dops) THEN endplot

END 
