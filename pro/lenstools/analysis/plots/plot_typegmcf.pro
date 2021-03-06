PRO plot_typegmcf_getstrings, struct, clr, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, $
                              one_error=one_error
  
  nlstring = ntostr(long(struct.nlenses))

  meanabsmag = alog10(struct.tmeanlum*1.e10)*(-2.5) + !sunmag[clr]
  
  absmagstring = ntostr(rnd(meanabsmag,3),7)
  lumstring = ntostr(rnd(struct.tmeanlum,3),5)
  lumstring = lumstring + ' $\pm$ '+ntostr(rnd(struct.tmeanlumerr,3),5)

  gmrstring = ntostr(rnd(struct.tmeangmr,3),5)

  chistring = ntostr(rnd(struct.wgm_chisq,2),4)+'/'+ntostr(struct.wgm_degfree)

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

PRO plot_typegmcf_print_table, all, early, late, red, blue

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
  printf,lun,"               $\Omega_m = 0.27$ was assumed. The outer 13 bins were used for the ``early'' and ``red'' sample fits.  For the ``late'' and ``blue'' samples, the data were rebinned from 18 to 9 radial bins.}"
  
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
  plot_typegmcf_getstrings, all, 2, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, $
                            one_error=one_error
  printf, lun, ' All & - & '+absmagstring+' ('+lumstring+') & '+gmrstring+' & '+nlstring+' & '+r0string+' & '+gamstring+' & '+chistring+' \\'
  
  ;; red
  plot_typegmcf_getstrings, red, 2, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, $
                            one_error=one_error
  printf, lun, ' Red & \gmr\ $ > $ \gmrcut\ & '+absmagstring+' ('+lumstring+') & '+gmrstring+' & '+nlstring+' & '+r0string+' & '+gamstring+' & '+chistring+' \\'
  
  ;; blue
  plot_typegmcf_getstrings, blue, 2, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, $
                            one_error=one_error
  printf, lun, ' Blue & \gmr\  $ < $ \gmrcut\ & '+absmagstring+' ('+lumstring+') & '+gmrstring+' & '+nlstring+' & '+r0string+' & '+gamstring+' & '+chistring+' \\'

  ;; Early
  plot_typegmcf_getstrings, early, 2, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, $
                            one_error=one_error
  printf, lun, ' Early & \eclass\ $ < $ \eclasscut\ & '+absmagstring+' ('+lumstring+') & '+gmrstring+' & '+nlstring+' & '+r0string+' & '+gamstring+' & '+chistring+' \\'

  ;; Late
  plot_typegmcf_getstrings, late, 2, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, $
                            one_error=one_error
  printf, lun, ' Late & \eclass\ $ > $ \eclasscut\ & '+absmagstring+' ('+lumstring+') & '+gmrstring+' & '+nlstring+' & '+r0string+' & '+gamstring+' & '+chistring
  
  printf,lun,'\enddata'
  
;  printf,lun,'\tablenotetext{a}{Bandpass in which the luminosities were measured}'
;  printf,lun,'\tablenotetext{b}{Range in absolute magnitude for this bin and bandpass}'
;  printf,lun,'\tablenotetext{c}{Mean luminosity for this bin}'
;  printf,lun,'\tablenotetext{d}{Number of lenses used.}'
;  printf,lun,'\tablenotetext{d}{Best fit power law index}'
;  printf,lun,'\tablenotetext{d}{Best fit scale length}'
;  printf,lun,'\tablenotetext{d}{$\chi^2$ per degree of freedom. Only the outer 13 radial bins were used for the two highest luminosity bins}'
  
  printf,lun,'\end{deluxetable}'
  


END 

PRO plot_typegmcf_plotall_chisq, all, early, late, red, blue, $
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

  contour, (all.wgm_chisq_surf-all.wgm_chisq), all.wgm_gamvals, all.wgm_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line1,line1,line1]
  contour, (early.wgm_chisq_surf-early.wgm_chisq), early.wgm_gamvals, early.wgm_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line1,line1,line1], color=early_clr
  ;;oplot, [early.gamma], [early.r0], psym=7
  contour, (late.wgm_chisq_surf-late.wgm_chisq), late.wgm_gamvals, late.wgm_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line2,line2,line2], color=late_clr
  ;;oplot, [late.gamma], [late.r0], psym=7, color=color2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; red/blue
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  multiplot
  plot,  [0], /nodata, xrange=tot_gamrange, yrange=tot_r0range, $
         xstyle=1+2, ystyle=1+2, xtitle=xt, xticklen=0.04, yticklen=0.04, charsize=1.2
  contour, (all.wgm_chisq_surf-all.wgm_chisq), all.wgm_gamvals, all.wgm_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line1,line1,line1]
  contour, (red.wgm_chisq_surf-red.wgm_chisq), red.wgm_gamvals, red.wgm_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line1,line1,line1], color=early_clr
  ;;oplot, [red.gamma], [red.r0], psym=7
  contour, (blue.wgm_chisq_surf-blue.wgm_chisq), blue.wgm_gamvals, blue.wgm_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line2,line2,line2], color=late_clr
  ;;oplot, [blue.gamma], [blue.r0], psym=7, color=color2

  multiplot,/reset

  IF !d.name EQ 'X' THEN BEGIN
      !p.background=bold
      !p.color=cold
  ENDIF 

END 
PRO plot_typegmcf_redo_fits, all, early, late, red, blue

  prompt=0
  IF prompt THEN noprompt=0 ELSE noprompt=1

  nr0=400
  ngam=400

  early_gamrange = [1.55, 2.15]
  early_r0range = [3.5,11]

  late_gamrange = [0.5, 3.0]
  late_r0range = [0,15]

  red_gamrange = early_gamrange
  red_r0range = early_r0range

  blue_r0range = late_r0range

  calc_r0_gamma_wgm,all,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse, $
                    gamrange=all_gamrange,r0range=all_r0range, /replace, $
                    noprompt=noprompt
  delvarx, wuse
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)

  wuse = where(early.meanr GT 100)
  calc_r0_gamma_wgm,early,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse, $
                    gamrange=early_gamrange,r0range=early_r0range, /replace, $
                    noprompt=noprompt
  delvarx, wuse
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)

  ;;wuse = where(late.meanr GT 100)
  calc_r0_gamma_wgm,late,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse, $
                    gamrange=late_gamrange,r0range=late_r0range, /replace, $
                    noprompt=noprompt,/rebin
;  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
;  calc_r0_gamma_wgm,late,nr0=nr0,ngam=ngam,$
;                    /dolegend,yfit=yfit,wuse=wuse, $
;                    gamrange=late_gamrange,r0range=late_r0range, /replace, $
;                    noprompt=noprompt, input_best1=late.gamma, input_best2=late.r0

  delvarx, wuse
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)

  wuse = where(red.meanr GT 100)
  calc_r0_gamma_wgm,red,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse, $
                    gamrange=red_gamrange,r0range=red_r0range, /replace, $
                    noprompt=noprompt
  delvarx, wuse
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)

  ;;wuse = where(blue.meanr GT 100)
  calc_r0_gamma_wgm,blue,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse, $
                    gamrange=blue_gamrange,r0range=blue_r0range, /replace, $
                    noprompt=noprompt,/rebin
  delvarx, wuse
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)


END 

PRO plot_typegmcf, dops=dops, encap=encap, color=color, plotfits=plotfits

  indir = esheldon_config("lensout_dir")+'combstripe/comb/'

  comove_str = '_comoving'

  all = mrdfits(indir+'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb'+comove_str+'_N1.fit',1)


  early = mrdfits(indir+'eclass1_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb'+comove_str+'_N2.fit',1)
  late  = mrdfits(indir+'eclass2_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb'+comove_str+'_N2.fit',1)

  cov = corr2cov(all.corr, late.sigmaerr)
  late.covariance = cov
  blue = mrdfits(indir+'gmr1_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb'+comove_str+'_N1.fit',1)
  red  = mrdfits(indir+'gmr2_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb'+comove_str+'_N1.fit',1)
  cov = corr2cov(all.corr, blue.sigmaerr)
  blue.covariance = cov

  fitstr=''
  IF keyword_set(plotfits) THEN BEGIN
      fitstr='_fits'
      plot_typegmcf_redo_fits, all, early, late, red, blue

      plot_typegmcf_print_table, all, early, late, red, blue
  ENDIF 
  IF keyword_set(color) THEN color_str = '_color' ELSE color_str=''

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

  psfile_chisq = repstr(psfile, 'deltasig','deltasig_chisq')

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
  ENDIF ELSE BEGIN 
      early_clr = !p.color
      late_clr = !p.color
  ENDELSE 

  setup_mystuff

  yrange = [0.1, 200]
  xrange = [0.011, 15]
  aplot, 1, [0],[0],/nodata,line=0, /xlog, /ylog, $
         xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, xtickf='loglabels', ytickf='loglabels', $
         xtitle=!mpcxtitle2, ytitle=!deltaytitle, xticklen=0.04, yticklen=0.04

  oploterror, early.meanr/1000, early.sigma, early.sigmaerr, $
              line=2, color=early_clr, errc=early_clr, $
              thick=thick1
;  oploterror, early.meanr_rebin/1000, early.sigma_rebin, early.sigmaerr_rebin, $
;              line=2, color=early_clr, errc=early_clr, $
;              thick=thick1

  oploterror, late.meanr_rebin/1000, late.sigma_rebin, late.sigmaerr_rebin, $
              line=0, color=late_clr, errc=late_clr, $
              thick=thick1

  oplot, red.meanr/1000,red.sigma,line=2,color=early_clr,thick=thick2
  oplot, blue.meanr_rebin/1000,blue.sigma_rebin,line=0,color=late_clr, thick=thick2

  IF keyword_set(dops) THEN endplot
  IF keyword_set(encap) THEN set_bbox, psfile, '%%BoundingBox: 35 15 480 445'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now do them side-by-side
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT keyword_set(dops) THEN key=get_kbrd(1)

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
 
  IF !d.name EQ 'X' THEN BEGIN
      !p.background=bold
      !p.color=cold
  ENDIF 

  IF keyword_set(dops) THEN endplot
  IF keyword_set(encap) THEN set_bbox, psfile, '%%BoundingBox: 35 15 480 260'

  IF keyword_set(plotfits) THEN BEGIN 

      IF keyword_set(dops) THEN BEGIN 

          IF keyword_set(encap) THEN BEGIN 
              begplot, name=psfile_chisq, color=color, xsize=7,ysize=4, /encap
          ENDIF ELSE BEGIN 
              begplot, name=psfile_chisq, color=color, /land
          ENDELSE 
      ENDIF 

      IF NOT keyword_set(dops) THEN key = get_kbrd(1)
      plot_typegmcf_plotall_chisq, all, early, late, red, blue, $
                                   color=color

      IF keyword_set(dops) THEN BEGIN 
          endplot
          IF NOT keyword_set(encap) THEN pslandfix, psfile_chisq
          IF keyword_set(encap) THEN set_bbox, psfile_chisq, '%%BoundingBox: 50 20 480 255'
      ENDIF 

  ENDIF 

END 
