PRO plot_typexi_getstrings, struct, clr, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, $
                              one_error=one_error
  
  nlstring = ntostr(long(struct.nlenses))

  meanabsmag = alog10(struct.tmeanlum*1.e10)*(-2.5) + !sunmag[clr]
  
  absmagstring = ntostr(rnd(meanabsmag,3),7)
  lumstring = ntostr(rnd(struct.tmeanlum,3),5)
  lumstring = lumstring + ' $\pm$ '+ntostr(rnd(struct.tmeanlumerr,3),5)

  gmrstring = ntostr(rnd(struct.tmeangmr,3),5)

  chistring = ntostr(rnd(struct.xi_chisq,2),4)+'/'+ntostr(struct.xi_degfree)

  r0 = ntostr(rnd(struct.r0,1),3)
  IF keyword_set(one_error) THEN BEGIN
      r0errlow = struct.r0 - struct.r0low[0]
      r0errhigh = struct.r0high[0] - struct.r0

      r0err = max([r0errlow, r0errhigh])
      r0err = ntostr(rnd( r0err, 1),3)
      r0string = r0+' $\pm$ '+r0err
  ENDIF ELSE BEGIN 
      r0errlow = ntostr(rnd(struct.r0 - struct.r0low[0],1),3)
      r0errhigh = ntostr(rnd(struct.r0high[0] - struct.r0,1),3)
  
      IF r0errlow EQ r0errhigh THEN BEGIN 
          r0string = r0+' $\pm$ '+r0errlow 
      ENDIF ELSE BEGIN 
          r0string = r0+'$^{+'+r0errhigh+'}_{-'+r0errlow+'}$' 
      ENDELSE 
  ENDELSE 
  gam = ntostr(rnd(struct.gamma,2),4)
  IF keyword_set(one_error) THEN BEGIN
      gamerrlow = struct.gamma-struct.gammalow[0]
      gamerrhigh = struct.gammahigh[0]-struct.gamma

      gamerr = max([gamerrlow,gamerrhigh])
      gamerr = ntostr(rnd(gamerr,2),4)
      gamstring = gam+' $\pm$ '+gamerr
  ENDIF ELSE BEGIN 
      gamerrlow = ntostr(rnd(struct.gamma-struct.gammalow[0],2),4)
      gamerrhigh = ntostr(rnd(struct.gammahigh[0]-struct.gamma,2),4)
      
      IF gamerrlow EQ gamerrhigh THEN BEGIN 
          gamstring = gam+' $\pm$ '+gamerrlow 
      ENDIF ELSE BEGIN 
          gamstring = gam+'$^{+'+gamerrhigh+'}_{-'+gamerrlow+'}$' 
      ENDELSE 
  ENDELSE 

END 

PRO plot_typexi_print_table, all, early, late, red, blue, vlim3

  one_error=1

  divider = '    & & & & & & \\'
  ;;divider = ' \hline '

  read_cuts, range_struct

  lun=-1
  
  printf,lun
  printf,lun,'\begin{deluxetable}{cccccccc}'
  printf,lun,'\tabletypesize{\small}'
  printf,lun,'\tablecaption{Correlation Functions for All, Early, and Late Type Galaxies \label{tab:allsamp}}'
      
  printf,lun,'\tablewidth{0pt}'
  printf,lun,'\tablecomments{Absolute magnitudes are \rmag-band Petrosian $M - 5 \log_{10} h$. '
  printf,lun,'               Values in parentheses are luminosity in units of $10^{10} h^{-2} L_{\sun}$.'
  printf,lun,'               The means are calculated using the same weights as the lensing measurement.'
  printf,lun,'               The value of $M_* (L_*)$ for the \rmag-band is -20.83 (1.54).'
  printf,lun,'               The $r_0$ and $\gamma$ are best fit parameters for $\xi_{gm} = '
  printf,lun,'               (r/r_0)^{-\gamma}$; $r_0$ is measured in $h^{-1}$ Mpc. A value of'
  printf,lun,"               $\Omega_m = 0.27$ was assumed. "
  printf,lun,"               The outer 12 bins were used for the ``early'' and ``red'' sample fits.  "
  printf,lun,"               For the ``late'' and ``blue'' samples, the data were rebinned from 17 to 8 radial bins."
  printf,lun,"}"
  
  printf,lun,'\tablehead{'
  printf,lun,'\colhead{Sample} &'
  printf,lun,'\colhead{Selection Criteria} &'
  printf,lun,'\colhead{Mean Abs. Mag.} &'
  printf,lun,'\colhead{Mean \gmr} &'
  printf,lun,'\colhead{N$_{Lenses}$} &'
  printf,lun,'\colhead{$r_0$} &'
  printf,lun,'\colhead{$\gamma$} &'
  printf,lun,'\colhead{$\chi^2/\nu$} '
  printf,lun,'}'
  
  printf,lun,'\'
  printf,lun,'\startdata'
  
  ;; All
  plot_typexi_getstrings, all, 2, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, $
                            one_error=one_error
  printf, lun, ' All & - & '+absmagstring+' ('+lumstring+') & '+gmrstring+' & '+nlstring+' & '+r0string+' & '+gamstring+' & '+chistring+' \\'


  ;; red
  plot_typexi_getstrings, red, 2, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, $
                            one_error=one_error
  printf, lun, ' Red & \gmr\ $ > $ \gmrcut\ & '+absmagstring+' ('+lumstring+') & '+gmrstring+' & '+nlstring+' & '+r0string+' & '+gamstring+' & '+chistring+' \\'
  
  ;; blue
  plot_typexi_getstrings, blue, 2, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, $
                            one_error=one_error
  printf, lun, ' Blue & \gmr\  $ < $ \gmrcut\ & '+absmagstring+' ('+lumstring+') & '+gmrstring+' & '+nlstring+' & '+r0string+' & '+gamstring+' & '+chistring+' \\'

  ;; Early
  plot_typexi_getstrings, early, 2, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, $
                            one_error=one_error
  printf, lun, ' Early & \eclass\ $ < $ \eclasscut\ & '+absmagstring+' ('+lumstring+') & '+gmrstring+' & '+nlstring+' & '+r0string+' & '+gamstring+' & '+chistring+' \\'

  ;; Late
  plot_typexi_getstrings, late, 2, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, $
                            one_error=one_error
  printf, lun, ' Late & \eclass\ $ > $ \eclasscut\ & '+absmagstring+' ('+lumstring+') & '+gmrstring+' & '+nlstring+' & '+r0string+' & '+gamstring+' & '+chistring+' \\'
  
;; Vlim
  plot_typexi_getstrings, vlim3, 2, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, $
                            one_error=one_error
  printf, lun, ' Vlim & $ -23.0 < M_r < -21.5 $ & '+absmagstring+' ('+lumstring+') & '+gmrstring+' & '+nlstring+' & '+r0string+' & '+gamstring+' & '+chistring+' \\'
  printf, lun, '   & $   0.1 < z   < 0.174 $ & & & & & & '

  printf,lun,'\enddata'  
  printf,lun,'\end{deluxetable}'
  
  print


END 

PRO plot_typexi_plotall_chisq, all, early, late, red, blue, $
                               color=color

  xt=!csym.gamma
  yt='r!D0!N [h!U'+!csym.minus+'1!N Mpc]'

  IF !d.name EQ 'X' THEN BEGIN
      bold = !p.background
      cold = !p.color
      !p.background = !white
      !p.color = !black

      thick1 = 2
      thick2 = !p.thick
  ENDIF ELSE BEGIN
      thick1 = 7
      thick2 = 1
  ENDELSE 

  IF keyword_set(color) THEN BEGIN 
      early_clr = !red
      late_clr = !blue

      line1 = 0
      line2 = 0
      line3 = 0
  ENDIF ELSE BEGIN 
      line1 = 0
      line2 = 1
      line3 = 2
      early_clr = !p.color
      late_clr = !p.color
  ENDELSE 


;  tot_gamrange = [0.9,2.5]
;  tot_r0range = [0,15]

  tot_gamrange = [1.3,2.4]
  tot_r0range = [1,12]

  erase & multiplot, [2,1], /square

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; early late
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  plot,  [0], /nodata, xrange=tot_gamrange, yrange=tot_r0range, $
         xstyle=1+2, ystyle=1+2, ytitle=yt, xtitle=xt, xticklen=0.04, yticklen=0.04, charsize=1.2

  contour, (all.xi_chisq_surf-all.xi_chisq), all.xi_gamvals, all.xi_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line1,line1,line1]
  contour, (early.xi_chisq_surf-early.xi_chisq), early.xi_gamvals, early.xi_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line1,line1,line1], color=early_clr
  ;;oplot, [early.gamma], [early.r0], psym=7
  contour, (late.xi_chisq_surf-late.xi_chisq), late.xi_gamvals, late.xi_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line2,line2,line2], color=late_clr
  ;;oplot, [late.gamma], [late.r0], psym=7, color=color2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; red/blue
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  multiplot
  plot,  [0], /nodata, xrange=tot_gamrange, yrange=tot_r0range, $
         xstyle=1+2, ystyle=1+2, xtitle=xt, xticklen=0.04, yticklen=0.04, charsize=1.2
  contour, (all.xi_chisq_surf-all.xi_chisq), all.xi_gamvals, all.xi_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line1,line1,line1]
  contour, (red.xi_chisq_surf-red.xi_chisq), red.xi_gamvals, red.xi_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line1,line1,line1], color=early_clr
  ;;oplot, [red.gamma], [red.r0], psym=7
  contour, (blue.xi_chisq_surf-blue.xi_chisq), blue.xi_gamvals, blue.xi_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line2,line2,line2], color=late_clr
  ;;oplot, [blue.gamma], [blue.r0], psym=7, color=color2

  multiplot,/reset

  IF !d.name EQ 'X' THEN BEGIN
      !p.background=bold
      !p.color=cold
  ENDIF 

END 
PRO plot_typexi_redo_fits, all, early, late, red, blue, vlim3

  prompt=0
  IF prompt THEN noprompt=0 ELSE noprompt=1

  ;; should we restrict the radii for fits?
  dowuse=0

  nr0=400
  ngam=400

  early_gamrange = [1.55, 2.15]
  early_r0range = [3.5,11]

  late_gamrange = [0.5, 3.0]
  late_r0range = [0,15]

  red_gamrange = early_gamrange
  red_r0range = early_r0range

  blue_r0range = late_r0range

  print,'--------------------------------'
  print,'Fitting All'
  fit_ximodel, all, nr0=nr0,ngam=ngam, $
               /dolegend,yfit=yfit,wuse=wuse, $
               gamrange=all_gamrange,r0range=all_r0range, /replace, $
               noprompt=noprompt

  delvarx, wuse
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)

  print,'--------------------------------'
  print,'Fitting vlim3'
  fit_ximodel, vlim3, nr0=nr0,ngam=ngam, /replace, $
               noprompt=noprompt

  IF prompt THEN key = prompt_kbrd()


  IF dowuse THEN wuse = where(early.r3 GT 0.1)
  print,'--------------------------------'
  print,'Fitting Early'
  fit_ximodel,early,nr0=nr0,ngam=ngam,$
              /dolegend,yfit=yfit,wuse=wuse, $
              gamrange=early_gamrange,r0range=early_r0range, /replace, $
              noprompt=noprompt

  delvarx, wuse
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)

  print,'--------------------------------'
  print,'Fitting Late'
  fit_ximodel,late,nr0=nr0,ngam=ngam,$
              /dolegend,yfit=yfit,wuse=wuse, $
              gamrange=late_gamrange,r0range=late_r0range, /replace, $
              noprompt=noprompt,/rebin,/nocov

  delvarx, wuse
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)

  print,'--------------------------------'
  print,'Fitting Red'
  IF dowuse THEN wuse = where(red.r3 GT 0.1)
  fit_ximodel,red,nr0=nr0,ngam=ngam,$
              /dolegend,yfit=yfit,wuse=wuse, $
              gamrange=red_gamrange,r0range=red_r0range, /replace, $
              noprompt=noprompt
  delvarx, wuse
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)

  print,'--------------------------------'
  print,'Fitting Blue'
  fit_ximodel,blue,nr0=nr0,ngam=ngam,$
              /dolegend,yfit=yfit,wuse=wuse, $
              gamrange=blue_gamrange,r0range=blue_r0range, /replace, $
              noprompt=noprompt,/rebin,/nocov
  print,'--------------------------------'

  delvarx, wuse
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)


END 

PRO plot_typexi, dops=dops, encap=encap, color=color, plotfits=plotfits

  ;; Maximum radius for plotting, since we have to extrapolate some
  ;; for the inversion

  max_radius = 7.0              ;Mpc

  lensout_dir = esheldon_config("lensout_dir")
  indir = lensout_dir+'combstripe/comb/'

  ;; comoving?
  comove_str = '_comoving'

  all = mrdfits(indir+'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi'+comove_str+'_N1.fit',1)

  early = mrdfits(indir+'eclass1_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi'+comove_str+'_N2.fit',1)
  late  = mrdfits(indir+'eclass2_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi'+comove_str+'_N2.fit',1)


  blue = mrdfits(indir+'gmr1_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi'+comove_str+'_N1.fit',1)
  red  = mrdfits(indir+'gmr2_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi'+comove_str+'_N1.fit',1)

  vlim3  = mrdfits(indir+'vlim3_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi'+comove_str+'_N1.fit',1)

  wa = where(all.r3   LT max_radius)
  we = where(early.r3 LT max_radius)
  wl = where(late.r3  LT max_radius)
  wr = where(red.r3   LT max_radius)
  wb = where(blue.r3  LT max_radius)
  wv = where(vlim3.r3 LT max_radius)

  fitstr=''
  IF keyword_set(plotfits) THEN BEGIN
      fitstr='_fits'
      plot_typexi_redo_fits, all, early, late, red, blue, vlim3
;return
      plot_typexi_print_table, all, early, late, red, blue, vlim3
;return
  ENDIF 
  IF keyword_set(color) THEN color_str = '_color' ELSE color_str=''

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Deltasig and xi plots for overall sample
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  psfile = '~/plots/xi_all'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+fitstr+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=7, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+fitstr+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 
  setup_mystuff

  xrange = [0.015, 10.0]
  yrange = [0.1, 3.e4]
  aploterror, 1, all.r3, all.xi, all.xierr, psym=8, /xlog, /ylog, $
              xrange=xrange, yrange=yrange, $
              xstyle=1+2, ystyle=1+2, xtickf='loglabels', ytickf='loglabels',$
              xtitle = !xigmxtitle,$
              ytitle = !xigmytitle, xticklen=0.04, yticklen=0.04

  IF keyword_set(plotfits) THEN BEGIN
      oplot, all.r3, (all.r0/all.r3)^all.gamma
      r0str = ntostr(all.r0, 3, /round)+' h!U'+!csym.minus+'1!N Mpc'
      gamstr = ntostr(all.gamma,4,/round)
      mess = '(r/'+r0str+')!U'+!csym.minus+gamstr+'!N'
      legend, mess, line=0, /right, box=0, thick=!p.thick, charsize=1.7
  END 

  key = prompt_kbrd()

  IF keyword_set(dops) THEN endplot, trim_bbox=encap
      
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; bias (all)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  idit = mrdfits(lensout_dir+'idit/idit_xi.fit',1)
  ;; correct idit's stuff for mean luminosity
  ;; function is (r0_idit/r)^1.8
  ;; therefor, the corrected one is xi*(r0_fixed/r0_idit)^1.8
  r0_idit = 5.77
  r0_fixed = 5.60
  idit_fix_factor = (r0_fixed/r0_idit)^1.8
  print,'Idit fix factor: ',idit_fix_factor
  xi_idit_fixed = idit.xi*idit_fix_factor

  ;; calculate the bias

  calc_xigg_xigm_bias, idit.r, xi_idit_fixed, idit.xierr, $
                       all.r3, all.xi, all.covxi, $
                       rcommon, $
                       bias, bias_cov, bias_err, $
                       cbias, cbias_err, $
                       bias_inv, bias_inv_cov, bias_inv_err, $
                       cbias_inv, cbias_inv_err

  IF NOT keyword_set(dops) THEN key=get_kbrd(1)

  ;; plot the bias
  print,'       radius     bias_inv bias_inv_err'
  colprint,rcommon,bias_inv, bias_inv_err
  print

  psfile = '~/plots/xi_all_bias_fit'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+fitstr+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=7, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+fitstr+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 
  setup_mystuff

  print,'b = '+ntostr(cbias)+' '+!plusminus+' '+ntostr(cbias_err)

  biasyt = 'b/r = '+!csym.xi+'!Dgg!N/'+!csym.xi+'!Dgm!N'
  bias_invyt = 'r/b = '+!csym.xi+'!Dgm!N/'+!csym.xi+'!Dgg!N'
  aploterror, 1.7, rcommon, bias, bias_err, /xlog, xstyle=2, $
              xtitle = !xigmxtitle, ytitle=biasyt, psym=8,/ynozero

  oplot, rcommon, bias, psym=8

  minr = min(rcommon, max=maxr)

  xx = [minr, maxr, maxr, minr]
  yy = [cbias+cbias_err,cbias+cbias_err,cbias-cbias_err,cbias-cbias_err]
  polyfill, xx, yy, /line_fill, orientation=45

  mess = '<b/r> = '+$
    ntostr(cbias,4,/round)+!csym.plusminus+$
    ntostr(cbias_err,4,/round)
  inv_mess = '<r/b> = '+$
    ntostr(cbias_inv,4,/round)+!csym.plusminus+$
    ntostr(cbias_inv_err,4,/round)
  legend, mess, /right, box=0

  IF NOT keyword_set(dops) THEN key=get_kbrd(1)

  ;; plot the inverse bias

  psfile = '~/plots/xi_all_bias_inv_fit'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+fitstr+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=7, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+fitstr+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 
  setup_mystuff

  aploterror, 1.7, rcommon, bias_inv, bias_inv_err, /xlog, xstyle=2, $
              xtitle = !xigmxtitle, ytitle=bias_invyt, psym=8,/ynozero

  oplot, rcommon, bias_inv, psym=8

  minr = min(rcommon, max=maxr)

  xx = [minr, maxr, maxr, minr]
  yy = [cbias_inv+cbias_inv_err,cbias_inv+cbias_inv_err,$
        cbias_inv-cbias_inv_err,cbias_inv-cbias_inv_err]
  polyfill, xx, yy, /line_fill, orientation=45


;  legend, inv_mess, /right, box=0



  key = prompt_kbrd()
  IF keyword_set(dops) THEN endplot, trim_bbox=encap

  ;; plot ours, idits, and bias

  psfile = '~/plots/xi_all_idit_bias'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+fitstr+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=10, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+fitstr+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 
  setup_mystuff

  xrange = [0.015, 10.0]
  plottwoerr, all.r3, all.xi, all.xierr, rcommon, bias, bias_err,$
              /xlog, /topylog, toppsym=8, botpsym=8, topaspect=1, $
              frac1 = 0.75, topytit=!xigmytitle, botyt = biasyt, $
              xtit=!xigmxtitle, $
              xoplot = idit.r, yoplot=xi_idit_fixed, $
              xrange=xrange, xstyle=1+2, $
              ytickformat1='loglabels',xtickformat='loglabels', $
              xticklen=0.04, yticklen=0.04
  polyfill, xx, yy, /line_fill, orientation=45

  ;; Plot the legend?
;  legend, mess, box=0,charsize=1.3;, /clear

  key = prompt_kbrd()
  IF keyword_set(dops) THEN endplot, trim_bbox=encap

  ;; plot ours, idits, and inverse bias

  psfile = '~/plots/xi_all_idit_bias_inv'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+fitstr+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=10, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+fitstr+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 
  setup_mystuff

  xrange = [0.015, 10.0]
  plottwoerr, all.r3, all.xi, all.xierr, rcommon, bias_inv, bias_inv_err,$
              /xlog, /topylog, toppsym=8, botpsym=8, topaspect=1, $
              frac1 = 0.75, topytit=!xigmytitle, botyt = bias_invyt, $
              xtit=!xigmxtitle, $
              xoplot = idit.r, yoplot=xi_idit_fixed, $
              xrange=xrange, xstyle=1+2, $
              ytickformat1='loglabels',xtickformat='loglabels', $
              xticklen=0.04, yticklen=0.04

  xx = [minr, maxr, maxr, minr]
  yy = [cbias_inv+cbias_inv_err,cbias_inv+cbias_inv_err,$
        cbias_inv-cbias_inv_err,cbias_inv-cbias_inv_err]
  polyfill, xx, yy, /line_fill, orientation=45

  ;; Plot the legend?
;  legend, inv_mess, box=0,charsize=1.3;, /clear

  key = prompt_kbrd()
  IF keyword_set(dops) THEN endplot, trim_bbox=encap
  


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; vlim sample
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  idit3 = mrdfits(lensout_dir+'idit/idit_xi_highlum.fit',1)

  w=where(vlim3.r3 GT 1.0,nw)
  nw = n_elements(vlim3.r3)
  w=lindgen(nw)
  covxi = vlim3.covxi[w[0]:w[nw-1], w[0]:w[nw-1] ]
  calc_xigg_xigm_bias, idit3.r, idit3.xi, idit3.xierr, $
                       vlim3.r3[w], vlim3.xi[w], covxi, $
                       rcommon, $
                       bias3, bias3_cov, bias3_err, $
                       cbias3, cbias3_err, $
                       bias3_inv, bias3_inv_cov, bias3_inv_err, $
                       cbias3_inv, cbias3_inv_err
  
  print,'       radius     bias_inv bias_inv_err'
  colprint,rcommon,bias3_inv, bias3_inv_err
  print

  psfile = '~/plots/xi_vlim3_bias_fit'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+fitstr+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=7, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+fitstr+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 
  setup_mystuff

  key=prompt_kbrd()

  biasyt = 'b/r = '+!csym.xi+'!Dgg!N/'+!csym.xi+'!Dgm!N'
  aploterror, 1.7, rcommon, bias3, bias3_err, /xlog, xstyle=2, $
              xtitle = !xigmxtitle, ytitle=biasyt, psym=8,/ynozero,$
              yrange=[-5,5]

  minr = min(rcommon, max=maxr)

  xx = [minr, maxr, maxr, minr]
  yy = [cbias3+cbias3_err,cbias3+cbias3_err,$
        cbias3-cbias3_err,cbias3-cbias3_err]
  polyfill, xx, yy, /line_fill, orientation=45

  mess = '<b/r> = '+$
    ntostr(cbias3,4,/round)+!csym.plusminus+$
    ntostr(cbias3_err,4,/round)
  inv_mess = '<r/b> = '+$
    ntostr(cbias3_inv,4,/round)+!csym.plusminus+$
    ntostr(cbias3_inv_err,4,/round)
  legend, mess, box=0,charsize=1.3;, /clear

  key = prompt_kbrd()
  IF keyword_set(dops) THEN endplot, trim_bbox=encap

  ;; plot ours, idits, and vlim bias

  psfile = '~/plots/xi_vlim3_idit_bias'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+fitstr+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=10, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+fitstr+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 
  setup_mystuff

  xrange = [0.015, 10.0]
  yrange = [0.1, 1.e5]
  botyrange = [-5,5]
  plottwoerr, vlim3.r3, vlim3.xi, vlim3.xierr, rcommon, bias3, bias3_err,$
              /xlog, /topylog, toppsym=8, botpsym=8, topaspect=1, $
              frac1 = 0.75, topytit=!xigmytitle, botyt = biasyt, $
              xtit=!xigmxtitle, $
              xoplot = idit3.r, yoplot=idit3.xi, $
              xrange=xrange, xstyle=1+2, topyrange=yrange, ystyle=1+2,$
              ytickformat1='loglabels',xtickformat='loglabels', $
              xticklen=0.04, yticklen=0.04, botyrange=botyrange, botystyle=3
  polyfill, xx, yy, /line_fill, orientation=45

  ;; Plot the legend?
;  legend, mess, box=0,charsize=1.3 ;, /clear

  key = prompt_kbrd()
  IF keyword_set(dops) THEN endplot, trim_bbox=encap

  ;; plot ours, idits, and inverse vlim bias

  psfile = '~/plots/xi_vlim3_idit_bias_inv'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+fitstr+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=10, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+fitstr+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 
  setup_mystuff

  xrange = [0.015, 10.0]
  yrange = [0.1, 1.e5]
  botyrange = [-1,2]
  plottwoerr, vlim3.r3, vlim3.xi, vlim3.xierr, rcommon, bias3_inv, bias3_inv_err,$
              /xlog, /topylog, toppsym=8, botpsym=8, topaspect=1, $
              frac1 = 0.75, topytit=!xigmytitle, botyt = bias_invyt, $
              xtit=!xigmxtitle, $
              xoplot = idit3.r, yoplot=idit3.xi, $
              xrange=xrange, xstyle=1+2, topyrange=yrange, ystyle=1+2,$
              ytickformat1='loglabels',xtickformat='loglabels', $
              xticklen=0.04, yticklen=0.04, botyrange=botyrange, botystyle=3

  xx = [minr, maxr, maxr, minr]
  yy = [cbias3_inv+cbias3_inv_err,cbias3_inv+cbias3_inv_err,$
        cbias3_inv-cbias3_inv_err,cbias3_inv-cbias3_inv_err]
  polyfill, xx, yy, /line_fill, orientation=45

  ;; Plot the legend?
;  legend, inv_mess, box=0,charsize=1.3 ;, /clear

  key = prompt_kbrd()
  IF keyword_set(dops) THEN endplot, trim_bbox=encap

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Deltasig Plots
  ;; early/late/red/blue all on same plot
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  psfile = '~/plots/deltasig_early_late'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+fitstr+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=7, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+fitstr+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 
  setup_mystuff

  psfile_chisq = repstr(psfile, 'deltasig','deltasig_chisq')

  IF !d.name EQ 'X' THEN BEGIN
      thick1 = 2
      IF keyword_set(color) THEN thick2 = thick1 ELSE thick2 = 1
  ENDIF ELSE BEGIN
      thick1 = 7
      IF keyword_set(color) THEN thick2 = thick1 ELSE thick2 = 1
  ENDELSE 

  IF keyword_set(color) THEN BEGIN 
      early_clr = !red
      late_clr = !blue
  ENDIF ELSE BEGIN 
      early_clr = !p.color
      late_clr = !p.color
  ENDELSE 

  ;; Red/blue are linestyle = line2 (2, dashed)
  ;; Early/Late are linestyle = line1 (0, solid)

  ;; Early/Red are thicker (thick1)
  ;; Late/Blue are less thick (thick2)

  line1 = 0
  line2 = 2

  yrange = [0.15, 200]
  xrange = [0.011, 15]
  aplot, 1, [0],[0],/nodata,line=0, /xlog, /ylog, $
         xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, $
         xtickf='loglabels', ytickf='loglabels', $
         xtitle=!mpcxtitle2, ytitle=!deltaytitle, xticklen=0.04, yticklen=0.04

  oploterror, early.meanr/1000, early.sigma, early.sigmaerr, $
              line=line1, color=early_clr, errc=early_clr, $
              thick=thick1
;  oploterror, early.meanr_rebin/1000, early.sigma_rebin, early.sigmaerr_rebin, $
;              line=2, color=early_clr, errc=early_clr, $
;              thick=thick1

  oploterror, late.meanr_rebin/1000, late.sigma_rebin, late.sigmaerr_rebin, $
              line=line1, color=late_clr, errc=late_clr, $
              thick=thick2

  oplot, red.meanr/1000,red.sigma,line=line2,color=early_clr,thick=thick1
  oplot, blue.meanr_rebin/1000,blue.sigma_rebin,line=line2,color=late_clr, thick=thick2

  key = prompt_kbrd()
  IF keyword_set(dops) THEN endplot, trim_bbox=encap

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Deltasig: now do them side-by-side
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT keyword_set(plotfits) THEN BEGIN 
      earlyline = 0
      lateline = 1
  ENDIF ELSE BEGIN 
      early_psym = 8
      late_psym = 1
  ENDELSE 
  hatlength = !D.X_VSIZE / 200
  symsizes = [0.7, 1]

  psfile = '~/plots/deltasig_early_late_sidebyside'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+fitstr+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=4, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+fitstr+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 
  setup_mystuff

  erase & multiplot,[2,1],/square
  plot, [0],[0],/nodata,line=0, /xlog, /ylog, $
        xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, xtickf='loglabels', ytickf='loglabels', $
        xtitle=!mpcxtitle2, ytitle=!deltaytitle, xticklen=0.04, yticklen=0.04

  oploterror, early.meanr/1000, early.sigma, early.sigmaerr, $
              color=early_clr, errc=early_clr, psym=early_psym, line=earlyline,symsize=symsizes[0]
  IF keyword_set(plotfits) THEN BEGIN 
      wuse=where(early.yfit GT 0) & oplot,early.meanr[wuse]/1000,early.yfit[wuse],color=early_clr
  ENDIF 
;  oploterror, early.meanr_rebin/1000, early.sigma_rebin, early.sigmaerr_rebin, $
;              line=2, color=early_clr, errc=early_clr
  oploterror, late.meanr_rebin/1000, late.sigma_rebin, late.sigmaerr_rebin, $
              color=late_clr, errc=late_clr, psym=late_psym, line=lateline,symsize=symsizes[1]
  IF keyword_set(plotfits) THEN BEGIN 
      wuse=where(late.yfit GT 0) & oplot,late.meanr_rebin[wuse]/1000,late.yfit[wuse],color=late_clr
;      wuse=where(late.yfit GT 0) & oplot,late.meanr[wuse]/1000,late.yfit[wuse],color=late_clr
  ENDIF 
  multiplot
  plot, [0],[0],/nodata,line=0, /xlog, /ylog, $
        xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, xtickf='loglabels', $
        xtitle=!mpcxtitle2, xticklen=0.04, yticklen=0.04

  oploterror, red.meanr/1000, red.sigma, red.sigmaerr, $
              color=early_clr, errc=early_clr, psym=early_psym, line=earlyline,symsize=symsizes[0]
  IF keyword_set(plotfits) THEN BEGIN 
      wuse=where(red.yfit GT 0) & oplot,red.meanr[wuse]/1000,red.yfit[wuse],color=early_clr
  ENDIF 
  oploterror, blue.meanr_rebin/1000, blue.sigma_rebin, blue.sigmaerr_rebin, $
              color=late_clr, errc=late_clr, psym=late_psym, line=lateline,symsize=symsizes[1]
  IF keyword_set(plotfits) THEN BEGIN 
      wuse=where(blue.yfit GT 0) & oplot,blue.meanr_rebin[wuse]/1000,blue.yfit[wuse],color=late_clr
;      wuse=where(blue.yfit GT 0) & oplot,blue.meanr[wuse]/1000,blue.yfit[wuse],color=late_clr
  ENDIF 
  multiplot,/reset
 
  key = prompt_kbrd()
  IF keyword_set(dops) THEN endplot, trim_bbox=encap


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; xigm Plots
  ;; early/late/red/blue all on same plot
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  psfile = '~/plots/xi_early_late'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+fitstr+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=7, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+fitstr+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 
  setup_mystuff

  psfile_chisq = repstr(psfile, 'xi','xi_chisq')

  yrange = [0.1, 4.e4]
  xrange = [0.015, 7]
  aplot, 1, [0],[0],/nodata,line=0, /xlog, /ylog, $
         xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, $
         xtickf='loglabels', ytickf='loglabels', $
         xtitle=!xigmxtitle, ytitle=!xigmytitle, xticklen=0.04, yticklen=0.04

  oploterror, early.r3[we], early.xi[we], early.xierr[we], $
              line=line1, color=early_clr, errc=early_clr, $
              thick=thick1


  oploterror, late.r3[wl], late.xi[wl], late.xierr[wl], $
              line=line1, color=late_clr, errc=late_clr, $
              thick=thick2

  oplot, red.r3[wr],red.xi[wr],line=line2,color=early_clr,thick=thick1
  oplot, blue.r3[wb],blue.xi[wb],$
         line=line2,color=late_clr, thick=thick2

  key = prompt_kbrd()
  IF keyword_set(dops) THEN endplot, trim_bbox=encap

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; xigm: now do them side-by-side
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT keyword_set(plotfits) THEN BEGIN 
      earlyline = 0
      lateline = 1
  ENDIF ELSE BEGIN 
      early_psym = 8
      late_psym = 1
  ENDELSE 
  hatlength = !D.X_VSIZE / 200
  symsizes = [0.7, 1]

  psfile = '~/plots/xi_early_late_sidebyside'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+fitstr+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=4, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+fitstr+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 
  setup_mystuff

  erase & multiplot,[2,1],/square
  plot, [0],[0],/nodata,line=0, /xlog, /ylog, $
        xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, xtickf='loglabels', ytickf='loglabels', $
        xtitle=!xigmxtitle, ytitle=!xigmytitle, xticklen=0.04, yticklen=0.04

  oploterror, early.r3[we], early.xi[we], early.xierr[we], $
              color=early_clr, errc=early_clr, psym=early_psym, line=earlyline,symsize=symsizes[0]
  IF keyword_set(plotfits) THEN BEGIN 
      oplot, early.r3[we], (early.r0/early.r3[we])^early.gamma
  ENDIF 

  oploterror, late.r3[wl], late.xi[wl], late.xierr[wl], $
              color=late_clr, errc=late_clr, psym=late_psym, line=lateline,symsize=symsizes[1]
  IF keyword_set(plotfits) THEN BEGIN 
      oplot, late.r3[wl], (late.r0/late.r3[wl])^late.gamma
  ENDIF 
  
  multiplot
  plot, [0],[0],/nodata,line=0, /xlog, /ylog, $
        xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, xtickf='loglabels', $
        xtitle=xtitle, xticklen=0.04, yticklen=0.04

  oploterror, red.r3[wr], red.xi[wr], red.xierr[wr], $
              color=early_clr, errc=early_clr, psym=early_psym, line=earlyline,symsize=symsizes[0]
  IF keyword_set(plotfits) THEN BEGIN 
      oplot, red.r3[wr], (red.r0/red.r3[wr])^red.gamma
  ENDIF 
  oploterror, blue.r3[wb], blue.xi[wb], blue.xierr[wb], $
              color=late_clr, errc=late_clr, psym=late_psym, line=lateline,symsize=symsizes[1]
  IF keyword_set(plotfits) THEN BEGIN 
      oplot, blue.r3[wb], (blue.r0/blue.r3[wb])^blue.gamma
  ENDIF 
  multiplot,/reset
 
  key=prompt_kbrd()
  IF keyword_set(dops) THEN endplot, trim_bbox=encap

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now put all on same plot
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  psfile = '~/plots/deltasig_xi_early_late'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+fitstr+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=10, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+fitstr+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 
  setup_mystuff

  ;;;;;;;;;;;;;;;;;;;
  ;; First, deltasig
  ;;;;;;;;;;;;;;;;;;;

  yrange = [0.1, 200]
  xrange = [0.011, 15]

  erase & multiplot,[1,2],/square

  plot, [0],[0],/nodata,line=0, /xlog, /ylog, $
        xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, $
        ytickf='loglabels', $
        ytitle=!deltaytitle, xticklen=0.04, yticklen=0.04

  oploterror, early.meanr/1000, early.sigma, early.sigmaerr, $
              line=line1, color=early_clr, errc=early_clr, $
              thick=thick1
;  oploterror, early.meanr_rebin/1000, early.sigma_rebin, early.sigmaerr_rebin, $
;              line=2, color=early_clr, errc=early_clr, $
;              thick=thick1

  oploterror, late.meanr_rebin/1000, late.sigma_rebin, late.sigmaerr_rebin, $
              line=line1, color=late_clr, errc=late_clr, $
              thick=thick2

  oplot, red.meanr/1000,red.sigma,line=line2,color=early_clr,thick=thick1
  oplot, blue.meanr_rebin/1000,blue.sigma_rebin,line=line2,color=late_clr, thick=thick2

  axis, xaxis=1, xtitle=!mpcxtitle, xrange=xrange, xstyle=1+2, xtickf='loglabels'

  message = ['Early', 'Red', 'Late', 'Blue']
  thick = [thick1, thick1, thick2, thick2]
  lines = [line1, line2, line1, line2]
  legend, message, line=lines, thick=thick, /right, box=0, charsize=1, $
          color=[early_clr, early_clr, late_clr, late_clr]

  ;;;;;;;;;;;;;;;
  ;; Now xi
  ;;;;;;;;;;;;;;;

  multiplot

  yrange = [0.1, 4.e4]
  ;;xrange = [0.015, 7]
  plot, [0],[0],/nodata,line=0, /xlog, /ylog, $
         xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, $
         xtickf='loglabels', ytickf='loglabels', $
         xtitle=!xigmxtitle, ytitle=!xigmytitle, xticklen=0.04, yticklen=0.04

  oploterror, early.r3[we], early.xi[we], early.xierr[we], $
              line=line1, color=early_clr, errc=early_clr, $
              thick=thick1

  oploterror, late.r3[wl], late.xi[wl], late.xierr[wl], $
              line=line1, color=late_clr, errc=late_clr, $
              thick=thick2

  oplot, red.r3[wr],red.xi[wr],line=line2,color=early_clr,thick=thick1
  oplot, blue.r3[wb],blue.xi[wb],$
         line=line2,color=late_clr, thick=thick2



  multiplot,/reset

  IF keyword_set(dops) THEN endplot, trim_bbox=encap

return
  IF keyword_set(plotfits) THEN BEGIN 

      IF keyword_set(dops) THEN BEGIN 

          IF keyword_set(encap) THEN BEGIN 
              begplot, name=psfile_chisq, color=color, xsize=7,ysize=4, /encap
          ENDIF ELSE BEGIN 
              begplot, name=psfile_chisq, color=color, /land
          ENDELSE 
      ENDIF 
      setup_mystuff

      IF NOT keyword_set(dops) THEN key = get_kbrd(1)
      plot_typexi_plotall_chisq, all, early, late, red, blue, $
                                   color=color

      IF keyword_set(dops) THEN BEGIN 
          endplot, trim_bbox=encap
          IF NOT keyword_set(encap) THEN pslandfix, psfile_chisq
          IF keyword_set(encap) THEN set_bbox, psfile_chisq, '%%BoundingBox: 50 20 480 255'
      ENDIF 

  ENDIF 

END 
