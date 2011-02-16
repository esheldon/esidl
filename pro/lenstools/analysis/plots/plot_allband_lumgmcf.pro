PRO plot_allband_lumgmcf_getrange, type, range_struct, clr, mmin, mmax

  CASE type OF
      'all': BEGIN 
          mmin = range_struct.main_minmag_threebin[clr,*]
          mmax = range_struct.main_maxmag_threebin[clr,*]
      END 
      'earlytwo': BEGIN 
          mmin = range_struct.early_minmag_twobin[clr,*]
          mmax = range_struct.early_maxmag_twobin[clr,*]
      END 
      'earlythree': BEGIN 
          mmin = range_struct.early_minmag_threebin[clr,*]
          mmax = range_struct.early_maxmag_threebin[clr,*]
      END 
      'redtwo': BEGIN 
          mmin = range_struct.red_minmag_twobin[clr,*]
          mmax = range_struct.red_maxmag_twobin[clr,*]
      END 
      'redthree': BEGIN 
          mmin = range_struct.red_minmag_threebin[clr,*]
          mmax = range_struct.red_maxmag_threebin[clr,*]
      END 
  ENDCASE 

  mmin = ntostr(rnd(mmin,1),5)
  mmax = ntostr(rnd(mmax,1),5)

END 

PRO plot_allband_lumgmcf_plotall_chisq, type, $
                        u1, u2, u3, $
                        g1, g2, g3, $
                        r1, r2, r3, $
                        i1, i2, i3, $
                        z1, z2, z3, color=color

  xt=!csym.gamma
  yt='r!D0!N [h!U'+!csym.minus+'1!N Mpc]'

  IF keyword_set(color) THEN BEGIN 
      IF !d.name EQ 'X' THEN BEGIN 
          color2 = !green
          color3 = !red
      ENDIF ELSE BEGIN 
          color2 = !blue
          color3 = !red
      ENDELSE 
      line1 = 0
      line2 = 0
      line3 = 0
  ENDIF ELSE BEGIN 
      line1 = 0
      line2 = 1
      line3 = 3
      color2 = !p.color
      color3 = !p.color
  ENDELSE 
line1=0
line2=1
line3=3
  CASE type OF
      'all': BEGIN 
          tot_gamrange = [1.31, 2.7]
          tot_r0range = [1,12]
      END 
      'earlytwo': BEGIN 
          tot_gamrange = [1.31, 2.5]
          tot_r0range = [2,12]
      END 
      'earlythree': BEGIN 
          tot_gamrange = [1.31, 2.6]
          tot_r0range = [2,12]
      END 
      'redtwo': BEGIN 
          tot_gamrange = [1.31, 2.5]
          tot_r0range = [2,12]
      END 
      'redthree': BEGIN 
          tot_gamrange = [1.31, 2.6]
          tot_r0range = [2,13]
      END 
  ENDCASE 

  tot_gamrange = [1.31,2.7]
  tot_r0range = [2,13]

  erase & multiplot, [3,2], /square

  ;; plot chisq together

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; u-band plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  plot,  [0], /nodata, xrange=tot_gamrange, yrange=tot_r0range, $
         xstyle=1+2, ystyle=1+2, ytitle=yt, xticklen=0.04, yticklen=0.04, charsize=1.2
  contour, (u1.wgm_chisq_surf-u1.wgm_chisq), u1.wgm_gamvals, u1.wgm_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line1,line1,line1]
  ;;oplot, [u1.gamma], [u1.r0], psym=7
  contour, (u2.wgm_chisq_surf-u2.wgm_chisq), u2.wgm_gamvals, u2.wgm_r0vals, $
           levels=!siglevels2, /overplot, color=color2, c_line=[line2,line2,line2]
  ;;oplot, [u2.gamma], [u2.r0], psym=7, color=color2
  IF n_elements(u3) NE 0 THEN BEGIN
      contour, (u3.wgm_chisq_surf-u3.wgm_chisq), u3.wgm_gamvals, u3.wgm_r0vals, $
               levels=!siglevels2, /overplot, color=color3, c_line=[line3,line3,line3]
      ;;oplot, [u3.gamma], [u3.ro], psym=7, color=color3
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; g-band plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  multiplot
  plot,  [0], /nodata, xrange=tot_gamrange, yrange=tot_r0range, $
         xstyle=1+2, ystyle=1+2, xticklen=0.04, yticklen=0.04, charsize=1.2
  contour, (g1.wgm_chisq_surf-g1.wgm_chisq), g1.wgm_gamvals, g1.wgm_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line1,line1,line1]
  contour, (g2.wgm_chisq_surf-g2.wgm_chisq), g2.wgm_gamvals, g2.wgm_r0vals, $
           levels=!siglevels2, /overplot, color=color2, c_line=[line2,line2,line2]
  IF n_elements(g3) NE 0 THEN contour, (g3.wgm_chisq_surf-g3.wgm_chisq), g3.wgm_gamvals, g3.wgm_r0vals, $
    levels=!siglevels2, /overplot, color=color3, c_line=[line3,line3,line3]

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; r-band plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  multiplot, /doxaxis
  plot,  [0], /nodata, xrange=tot_gamrange, yrange=tot_r0range, $
         xstyle=1+2, ystyle=1+2, xtitle=xt, xticklen=0.04, yticklen=0.04, charsize=1.2
  contour, (r1.wgm_chisq_surf-r1.wgm_chisq), r1.wgm_gamvals, r1.wgm_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line1,line1,line1]
  contour, (r2.wgm_chisq_surf-r2.wgm_chisq), r2.wgm_gamvals, r2.wgm_r0vals, $
           levels=!siglevels2, /overplot, color=color2, c_line=[line2,line2,line2]
  IF n_elements(r3) NE 0 THEN contour, (r3.wgm_chisq_surf-r3.wgm_chisq), r3.wgm_gamvals, r3.wgm_r0vals, $
    levels=!siglevels2, /overplot, color=color3, c_line=[line3,line3,line3]


  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; i-band plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  multiplot
  plot,  [0], /nodata, xrange=tot_gamrange, yrange=tot_r0range, $
         xstyle=1+2, ystyle=1+2, ytitle=yt, xtitle=xt, xticklen=0.04, yticklen=0.04, charsize=1.2
  contour, (i1.wgm_chisq_surf-i1.wgm_chisq), i1.wgm_gamvals, i1.wgm_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line1,line1,line1]
  contour, (i2.wgm_chisq_surf-i2.wgm_chisq), i2.wgm_gamvals, i2.wgm_r0vals, $
           levels=!siglevels2, /overplot, color=color2, c_line=[line2,line2,line2]
  IF n_elements(i3) NE 0 THEN contour, (i3.wgm_chisq_surf-i3.wgm_chisq), i3.wgm_gamvals, i3.wgm_r0vals, $
    levels=!siglevels2, /overplot, color=color3, c_line=[line3,line3,line3]

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; z-band plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  multiplot
  plot,  [0], /nodata, xrange=tot_gamrange, yrange=tot_r0range, $
         xstyle=1+2, ystyle=1+2, xtitle=xt, xticklen=0.04, yticklen=0.04, charsize=1.2
  contour, (z1.wgm_chisq_surf-z1.wgm_chisq), z1.wgm_gamvals, z1.wgm_r0vals, $
           levels=!siglevels2, /overplot, c_line=[line1,line1,line1]
  contour, (z2.wgm_chisq_surf-z2.wgm_chisq), z2.wgm_gamvals, z2.wgm_r0vals, $
           levels=!siglevels2, /overplot, color=color2, c_line=[line2,line2,line2]
  IF n_elements(z3) NE 0 THEN contour, (z3.wgm_chisq_surf-z3.wgm_chisq), z3.wgm_gamvals, z3.wgm_r0vals, $
    levels=!siglevels2, /overplot, color=color3, c_line=[line3,line3,line3]

  multiplot,/reset

END 

PRO plot_allband_lumgmcf_getstrings, struct, clr, nlstring, absmagstring, lumstring, chistring, r0string, gamstring, gmrstring, one_error=one_error

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

PRO plot_allband_lumgmcf_print_table, type, $
                      u1, u2, u3, $
                      g1, g2, g3, $
                      r1, r2, r3, $
                      i1, i2, i3, $
                      z1, z2, z3

  one_error=1

  divider = '    & & & & & & \\'
  ;;divider = ' \hline '

  read_cuts, range_struct

  lun=-1
  
  printf,lun
  printf,lun,'\begin{deluxetable}{cccccccc}'
  printf,lun,'\tabletypesize{\small}'
  IF type EQ 'all' THEN BEGIN 
      printf,lun,'\tablecaption{Luminosity Bins for All Galaxies \label{tab:lumbin}}'
  ENDIF ELSE IF type EQ 'earlytwo' OR type EQ 'earlythree' THEN BEGIN 
      printf,lun,'\tablecaption{Luminosity Bins for Early Types \label{tab:earlylumbin}}'
  ENDIF ELSE IF type EQ 'redtwo' OR type EQ 'redthree' THEN BEGIN 
      printf,lun,'\tablecaption{Luminosity Bins for Red Galaxies \label{tab:earlylumbin}}'
  ENDIF 
  printf,lun,'\tablewidth{0pt}'
  IF type EQ 'all' THEN BEGIN 
      printf,lun,'\tablecomments{See table \ref{tab:allsamp} for explanation of columns. The '
      printf,lun,'     value of $M_* (L_*)$ is -18.34(0.77), -20.04(1.10), '
      printf,lun,'     -20.83(1.54), -21.26(2.07), -21.55(2.68) for $u, g, r, i, z$ '
      printf,lun,'     respectively.}'
  ENDIF ELSE IF type EQ 'earlytwo' OR type EQ 'earlythree' THEN BEGIN 
      printf,lun,'\tablecomments{Same as table \ref{tab:lumbin} but for early types.  Early types are defined'
      printf,lun,'               as \eclass$ < -0.02$.}'
  ENDIF ELSE IF type EQ 'redtwo' OR type EQ 'redthree' THEN BEGIN 
      printf,lun,'\tablecomments{Same as table \ref{tab:lumbin} but for red galaxies, defined'
      printf,lun,'               as galaxies with \gmr\ $ > 0.7$.}'
  ENDIF 
  printf,lun,'\tablehead{'
  printf,lun,'\colhead{Bandpass} &'
  printf,lun,'\colhead{Abs. Mag. Range} &'
  printf,lun,'\colhead{Mean Abs. Mag.} &'
  printf,lun,'\colhead{Mean \gmr} &'
  printf,lun,'\colhead{N$_{Lenses}$} &'
  printf,lun,'\colhead{$r_0$} &'
  printf,lun,'\colhead{$\gamma$} &'
  printf,lun,'\colhead{$\chi^2/\nu$} '
  printf,lun,'}'
  
  printf,lun,'\'
  printf,lun,'\startdata'
  
  ;; u-band
  plot_allband_lumgmcf_getrange, type, range_struct, 0, mmin, mmax
  
  plot_allband_lumgmcf_getstrings, u1, 0, nl1, abs1, lum1, chi1, r0_1, gam1,gmr1,one_error=one_error
  plot_allband_lumgmcf_getstrings, u2, 0, nl2, abs2, lum2, chi2, r0_2, gam2,gmr2,one_error=one_error
  IF n_elements(u3) ne 0 THEN plot_allband_lumgmcf_getstrings, u3, 0, nl3, abs3, lum3, chi3, r0_3, gam3,gmr3,one_error=one_error
  
  printf,lun,'$u$ & '+mmin[0]+' $ < M_{u} < $ '+mmax[0]+' & '+abs1+' ('+lum1+') & '+gmr1+' & '+nl1+' & '+r0_1+' & '+gam1+' & '+chi1+' \\'
  printf,lun,' -  & '+mmin[1]+' $ < M_{u} < $ '+mmax[1]+' & '+abs2+' ('+lum2+') & '+gmr2+' & '+nl2+' & '+r0_2+' & '+gam2+' & '+chi2+' \\'
  IF n_elements(u3) ne 0 THEN printf,lun,' -  & '+mmin[2]+' $ < M_{u} < $ '+mmax[2]+' & '+abs3+' ('+lum3+') & '+gmr3+' & '+nl3+' & '+r0_3+' & '+gam3+' & '+chi3+' \\'
  
  ;; g-band
  printf,lun,divider
  plot_allband_lumgmcf_getrange, type, range_struct, 1, mmin, mmax
  
  plot_allband_lumgmcf_getstrings, g1, 1, nl1, abs1, lum1, chi1, r0_1, gam1,gmr1,one_error=one_error
  plot_allband_lumgmcf_getstrings, g2, 1, nl2, abs2, lum2, chi2, r0_2, gam2,gmr2,one_error=one_error
  IF n_elements(g3) ne 0 THEN plot_allband_lumgmcf_getstrings, g3, 1, nl3, abs3, lum3, chi3, r0_3, gam3,gmr3,one_error=one_error

  printf,lun,'$g$ & '+mmin[0]+' $ < M_{g} < $ '+mmax[0]+' & '+abs1+' ('+lum1+') & '+gmr1+' & '+nl1+' & '+r0_1+' & '+gam1+' & '+chi1+' \\'
  printf,lun,' -  & '+mmin[1]+' $ < M_{g} < $ '+mmax[1]+' & '+abs2+' ('+lum2+') & '+gmr2+' & '+nl2+' & '+r0_2+' & '+gam2+' & '+chi2+' \\'
  IF n_elements(g3) ne 0 THEN printf,lun,' -  & '+mmin[2]+' $ < M_{g} < $ '+mmax[2]+' & '+abs3+' ('+lum3+') & '+gmr3+' & '+nl3+' & '+r0_3+' & '+gam3+' & '+chi3+' \\'

  ;; r-band
  printf,lun,divider
  plot_allband_lumgmcf_getrange, type, range_struct, 2, mmin, mmax

  plot_allband_lumgmcf_getstrings, r1, 2, nl1, abs1, lum1, chi1, r0_1, gam1,gmr1,one_error=one_error
  plot_allband_lumgmcf_getstrings, r2, 2, nl2, abs2, lum2, chi2, r0_2, gam2,gmr2,one_error=one_error
  IF n_elements(r3) ne 0 THEN plot_allband_lumgmcf_getstrings, r3, 2, nl3, abs3, lum3, chi3, r0_3, gam3,gmr3,one_error=one_error
  
  printf,lun,'$r$ & '+mmin[0]+' $ < M_{r} < $ '+mmax[0]+' & '+abs1+' ('+lum1+') & '+gmr1+' & '+nl1+' & '+r0_1+' & '+gam1+' & '+chi1+' \\'
  printf,lun,' -  & '+mmin[1]+' $ < M_{r} < $ '+mmax[1]+' & '+abs2+' ('+lum2+') & '+gmr2+' & '+nl2+' & '+r0_2+' & '+gam2+' & '+chi2+' \\'
  IF n_elements(r3) ne 0 THEN printf,lun,' -  & '+mmin[2]+' $ < M_{r} < $ '+mmax[2]+' & '+abs3+' ('+lum3+') & '+gmr3+' & '+nl3+' & '+r0_3+' & '+gam3+' & '+chi3+' \\'
  
  ;; i-band
  printf,lun,divider
  plot_allband_lumgmcf_getrange, type, range_struct, 3, mmin, mmax
  
  plot_allband_lumgmcf_getstrings, i1, 3, nl1, abs1, lum1, chi1, r0_1, gam1,gmr1,one_error=one_error
  plot_allband_lumgmcf_getstrings, i2, 3, nl2, abs2, lum2, chi2, r0_2, gam2,gmr2,one_error=one_error
  IF n_elements(i3) ne 0 THEN plot_allband_lumgmcf_getstrings, i3, 3, nl3, abs3, lum3, chi3, r0_3, gam3,gmr3,one_error=one_error
  
  printf,lun,'$i$ & '+mmin[0]+' $ < M_{i} < $ '+mmax[0]+' & '+abs1+' ('+lum1+') & '+gmr1+' & '+nl1+' & '+r0_1+' & '+gam1+' & '+chi1+' \\'
  printf,lun,' -  & '+mmin[1]+' $ < M_{i} < $ '+mmax[1]+' & '+abs2+' ('+lum2+') & '+gmr2+' & '+nl2+' & '+r0_2+' & '+gam2+' & '+chi2+' \\'
  IF n_elements(i3) ne 0 THEN printf,lun,' -  & '+mmin[2]+' $ < M_{i} < $ '+mmax[2]+' & '+abs3+' ('+lum3+') & '+gmr3+' & '+nl3+' & '+r0_3+' & '+gam3+' & '+chi3+' \\'
  
  ;; z-band
  printf,lun,divider
  plot_allband_lumgmcf_getrange, type, range_struct, 4, mmin, mmax
  
  plot_allband_lumgmcf_getstrings, z1, 4, nl1, abs1, lum1, chi1, r0_1, gam1,gmr1,one_error=one_error
  plot_allband_lumgmcf_getstrings, z2, 4, nl2, abs2, lum2, chi2, r0_2, gam2,gmr2,one_error=one_error
  IF n_elements(z3) ne 0 THEN plot_allband_lumgmcf_getstrings, z3, 4, nl3, abs3, lum3, chi3, r0_3, gam3,gmr3,one_error=one_error
  
  printf,lun,'$z$ & '+mmin[0]+' $ < M_{z} < $ '+mmax[0]+' & '+abs1+' ('+lum1+') & '+gmr1+' & '+nl1+' & '+r0_1+' & '+gam1+' & '+chi1+' \\'
  IF n_elements(z3) ne 0 THEN BEGIN 
      printf,lun,' -  & '+mmin[1]+' $ < M_{z} < $ '+mmax[1]+' & '+abs2+' ('+lum2+') & '+gmr2+' & '+nl2+' & '+r0_2+' & '+gam2+' & '+chi2+' \\'
      printf,lun,' -  & '+mmin[2]+' $ < M_{z} < $ '+mmax[2]+' & '+abs3+' ('+lum3+') & '+gmr3+' & '+nl3+' & '+r0_3+' & '+gam3+' & '+chi3
  ENDIF ELSE BEGIN 
      printf,lun,' -  & '+mmin[1]+' $ < M_{z} < $ '+mmax[1]+' & '+abs2+' ('+lum2+') & '+gmr2+' & '+nl2+' & '+r0_2+' & '+gam2+' & '+chi2
  ENDELSE 
  
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

PRO plot_allband_lumgmcf_redo_fits, type, $
                    u1, u2, u3, $
                    g1, g2, g3, $
                    r1, r2, r3, $
                    i1, i2, i3, $
                    z1, z2, z3

  nr0 = 400
  ngam = 400

  prompt=0
  IF prompt THEN noprompt=0 ELSE noprompt=1
  
  CASE type of
      'all': BEGIN 

          ;;u1_gamrange=[1.4,2.0]
          ;;u1_r0range=[2,11]
          u2_gamrange=[1.6,3.0]
          u2_r0range = [1,9]
          u3_gamrange=[1.7,2.6]
          u3_r0range=[2,12]
          
          tot_gamrange = [1.4, 3.0]
          tot_r0range = [1,12]
          
          ;;g1_gamrange=[1.4,2.0]
          ;;g1_r0range=[2,11]
          g2_gamrange=[1.6,2.8]
          g2_r0range=[1,9]
          g3_gamrange=[1.8,2.6]
          g3_r0range=[2,11]
          
          ;;r1_gamrange=[1.4,2.0]
          ;;r1_r0range=[2,11]
          r2_gamrange=[1.6,2.8]
          r2_r0range=[1,9]
          r3_gamrange=[1.8,2.5]
          r3_r0range=[3,11]
          
          ;;i1_gamrange=[1.4,2.0]
          ;;i1_r0range=[2,11]
          i2_gamrange=[1.6,2.8]
          i2_r0range=[1,9]
          i3_gamrange=[1.8,2.5]
          i3_r0range=[3,11]
          
          ;;z1_gamrange=[1.4,2.0]
          ;;z1_r0range=[2,11]
          z2_gamrange=[1.6,2.8]
          z2_r0range=[1,9]
          z3_gamrange=[1.9,2.6]
          z3_r0range=[2,10]
      END 
      'earlytwo': BEGIN 

          tot_gamrange = [1.4, 3]
          tot_r0range =  [1,12]
          
          ;;u1_gamrange=[1.4,2.0]
          ;;u1_r0range=[2,11]
          u2_gamrange=[1.8,2.5]
          u2_r0range = [3,9.5]
          
          ;;g1_gamrange=[1.4,2.0]
          ;;g1_r0range=[2,11]
          g2_gamrange=[1.9,2.5]
          g2_r0range=[3,9.5]
          
          ;;r1_gamrange=[1.4,2.0]
          ;;r1_r0range=[2,11]
          r2_gamrange=[1.9,2.5]
          r2_r0range=[3,10]
          
          ;;i1_gamrange=[1.4,2.0]
          ;;i1_r0range=[2,11]
          i2_gamrange=[1.9,2.4]
          i2_r0range=[3,10]

          ;;z1_gamrange=[1.4,2.0]
          ;;z1_r0range=[2,11]
          z2_gamrange=[1.9,2.6]
          z2_r0range=[2,9]
      END 
      'earlythree': BEGIN 
          
          u1_gamrange = [1.51,1.95]
          u1_r0range  = [3.9,12]
          u2_gamrange = [1.8,2.6]
          u2_r0range  = [2.1,9.0]
          u3_gamrange = [1.8,2.5]
          u3_r0range  = [3,13]
          
          g1_gamrange = [1.51,1.95]
          g1_r0range  = [3.9,12.5] ;
          g2_gamrange = [1.8,2.6]
          g2_r0range  = [2.1,9.0]
          g3_gamrange = [1.85,2.5]
          g3_r0range  = [3,12]
              
          r1_gamrange = [1.51,1.95]
          r1_r0range  = [3.9,12]
          r2_gamrange = [1.8,2.5] ;
          r2_r0range  = [2.1,10.0]
          r3_gamrange = [1.9,2.5]
          r3_r0range  = [3,12]
              
          i1_gamrange = [1.51,1.95]
          i1_r0range  = [3.9,12]
          i2_gamrange = [1.65,2.45]
          i2_r0range  = [3.0,12.0]
          i3_gamrange = [1.9,2.5]
          i3_r0range  = [3,12]
              
          z1_gamrange = [1.51,1.95]
          z1_r0range  = [3.9,12]
          z2_gamrange = [1.6,2.75]
          z2_r0range  = [1.4,11.0]
          z3_gamrange = [1.95,2.6]
          z3_r0range  = [3,10]
            
      END 
      'redtwo': BEGIN 

          tot_gamrange = [1.4, 3]
          tot_r0range =  [1,12]
          
          ;;u1_gamrange=[1.4,2.0]
          ;;u1_r0range=[2,11]
          u2_gamrange=[1.8,2.5]
          u2_r0range = [3,9.5]
          
          ;;g1_gamrange=[1.4,2.0]
          ;;g1_r0range=[2,11]
          g2_gamrange=[1.9,2.5]
          g2_r0range=[3,9.5]
          
          ;;r1_gamrange=[1.4,2.0]
          ;;r1_r0range=[2,11]
          r2_gamrange=[1.9,2.5]
          r2_r0range=[3,10]
          
          ;;i1_gamrange=[1.4,2.0]
          ;;i1_r0range=[2,11]
          i2_gamrange=[1.9,2.4]
          i2_r0range=[3,10]

          ;;z1_gamrange=[1.4,2.0]
          ;;z1_r0range=[2,11]
          z2_gamrange=[1.9,2.6]
          z2_r0range=[2,9]
      END 
      'redthree': BEGIN 
          
          u1_gamrange = [1.51,1.95]
          u1_r0range  = [3.9,12]
          u2_gamrange = [1.8,2.6]
          u2_r0range  = [2.1,9.0]
          u3_gamrange = [1.8,2.5]
          u3_r0range  = [3.5,13.5]
          
          g1_gamrange = [1.51,1.95]
          g1_r0range  = [3.9,12.5] ;
          g2_gamrange = [1.8,2.6]
          g2_r0range  = [2.1,9.0]
          g3_gamrange = [1.85,2.5]
          g3_r0range  = [4,12]
              
          r1_gamrange = [1.51,1.95]
          r1_r0range  = [3.9,13]
          r2_gamrange = [1.8,2.6] ;
          r2_r0range  = [2.1,10.0]
          r3_gamrange = [1.9,2.5]
          r3_r0range  = [4,12]
              
          i1_gamrange = [1.51,1.95]
          i1_r0range  = [3.9,13]
          i2_gamrange = [1.8,2.6]
          i2_r0range  = [2.5,11.0]
          i3_gamrange = [1.9,2.5]
          i3_r0range  = [3,12]
              
          z1_gamrange = [1.51,1.95]
          z1_r0range  = [3.9,13]
          z2_gamrange = [1.6,2.75]
          z2_r0range  = [1.4,11.0]
          z3_gamrange = [1.95,2.55]
          z3_r0range  = [3,10]
            
      END 
  ENDCASE 

  ;;;;;;;;;;;;;;;;;;
  ;; u
  ;;;;;;;;;;;;;;;;;;

  print,'** Doing u-band fits'
  calc_r0_gamma_wgm,u1,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse, $
                    gamrange=u1_gamrange,r0range=u1_r0range, /replace, noprompt=noprompt
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  wuse = where(u2.meanr GT 100)
  calc_r0_gamma_wgm,u2,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse,$
                    gamrange=u2_gamrange,r0range=u2_r0range, /replace, noprompt=noprompt
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  IF n_elements(u3) NE 0 THEN BEGIN 
      calc_r0_gamma_wgm,u3,nr0=nr0,ngam=ngam,$
                        /dolegend,yfit=yfit,wuse=wuse,$
                        gamrange=u3_gamrange,r0range=u3_r0range, /replace, noprompt=noprompt
      IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  ENDIF 

  delvarx, wuse

  ;;;;;;;;;;;;;;;;;;
  ;; g
  ;;;;;;;;;;;;;;;;;;

  print,'** Doing g-band fits'
  calc_r0_gamma_wgm,g1,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse, /replace, $
                    gamrange=g1_gamrange,r0range=g1_r0range, noprompt=noprompt
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  wuse = where(g2.meanr GT 100)
  calc_r0_gamma_wgm,g2,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse,$
                    gamrange=g2_gamrange,r0range=g2_r0range, /replace, noprompt=noprompt
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  IF n_elements(g3) NE 0 THEN BEGIN 
      calc_r0_gamma_wgm,g3,nr0=nr0,ngam=ngam,$
                        /dolegend,yfit=yfit,wuse=wuse,$
                        gamrange=g3_gamrange,r0range=g3_r0range, /replace, noprompt=noprompt
      IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  ENDIF 

  delvarx, wuse

  ;;;;;;;;;;;;;;;;;;
  ;; r
  ;;;;;;;;;;;;;;;;;;

  print,'** Doing r-band fits'
  calc_r0_gamma_wgm,r1,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse, /replace, $
                    gamrange=r1_gamrange,r0range=r1_r0range, noprompt=noprompt
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  wuse = where(r2.meanr GT 100)
  calc_r0_gamma_wgm,r2,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse,$
                    gamrange=r2_gamrange,r0range=r2_r0range, /replace, noprompt=noprompt
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  IF n_elements(r3) NE 0 THEN BEGIN 
      calc_r0_gamma_wgm,r3,nr0=nr0,ngam=ngam,$
                        /dolegend,yfit=yfit,wuse=wuse,$
                        gamrange=r3_gamrange,r0range=r3_r0range, /replace, noprompt=noprompt
      IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  ENDIF 

  delvarx, wuse

  ;;;;;;;;;;;;;;;;;;
  ;; i
  ;;;;;;;;;;;;;;;;;;

  print,'** Doing i-band fits'
  calc_r0_gamma_wgm,i1,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse, /replace, $
                    gamrange=i1_gamrange,r0range=i1_r0range, noprompt=noprompt
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  wuse = where(i2.meanr GT 100)
  calc_r0_gamma_wgm,i2,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse,$
                    gamrange=i2_gamrange,r0range=i2_r0range, /replace, noprompt=noprompt
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  IF n_elements(i3) NE 0 THEN BEGIN 
      calc_r0_gamma_wgm,i3,nr0=nr0,ngam=ngam,$
                        /dolegend,yfit=yfit,wuse=wuse,$
                        gamrange=i3_gamrange,r0range=i3_r0range, /replace, noprompt=noprompt
      IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  ENDIF 

  delvarx, wuse

  ;;;;;;;;;;;;;;;;;;
  ;; z
  ;;;;;;;;;;;;;;;;;;

  print,'** Doing z-band fits'
  calc_r0_gamma_wgm,z1,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse, /replace, $
                    gamrange=z1_gamrange,r0range=z1_r0range, noprompt=noprompt
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  wuse = where(z2.meanr GT 100)
  calc_r0_gamma_wgm,z2,nr0=nr0,ngam=ngam,$
                    /dolegend,yfit=yfit,wuse=wuse,$
                    gamrange=z2_gamrange,r0range=z2_r0range, /replace, noprompt=noprompt
  IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  IF n_elements(z3) NE 0 THEN BEGIN 
      calc_r0_gamma_wgm,z3,nr0=nr0,ngam=ngam,$
                        /dolegend,yfit=yfit,wuse=wuse,$
                        gamrange=z3_gamrange,r0range=z3_r0range, /replace, noprompt=noprompt
      IF (!d.name EQ 'X') and prompt THEN key = get_kbrd(1)
  ENDIF 

  
END
 
PRO plot_allband_lumgmcf, type, wgm=wgm, dops=dops, encap=encap, oldcuts=oldcuts, plotfits=plotfits, color=color

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: plot_allband_lumgmcf, type, wgm=wgm, dops=dops'
      print
      print,'type="all" or "early"'
      return
  ENDIF 

  ;; plot the bins for each bandpass
  ;; this is new gmcf stuff (DR1+)
  ;; not wgm is based on power-law fits, so not as immediate as
  ;; deltasig

  dir = esheldon_config("lensout_dir")+'combstripe/comb/sublum/'

  IF keyword_set(plotfits) THEN pfitstr='_fits' ELSE pfitstr=''

  ;; NOTE: %%BoundingBox: 35 15 770 515
  IF keyword_set(dops) AND keyword_set(encap) THEN endstr='.eps' ELSE endstr='.ps'

  tstr = ''
  psstr = ''
  IF type EQ 'all' THEN BEGIN 
      IF keyword_set(oldcuts) THEN BEGIN 
          tstr='' 
          psstr = '_oldcuts'
      ENDIF ELSE BEGIN 
          tstr='num'
          psstr = ''
      ENDELSE 
  ENDIF

  IF keyword_set(color) THEN colorstr = '_color' ELSE colorstr = ''

  dothird=0
  fend = '_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_comoving_N1.fit'
  CASE type OF 
      'all': BEGIN 
          ffront = 'lum'      & fend = 'threebin'+tstr+fend & dothird=1
      END 
      'earlytwo': BEGIN 
          ffront = 'earlylum' & fend = 'twobin'+fend
      END  
      'earlythree': BEGIN 
          ffront = 'earlylum' & fend = 'threebin'+fend & dothird=1
      END 
      'redtwo': BEGIN 
          ffront = 'redlum'   & fend = 'twobin'+fend
      END 
      'redthree': BEGIN 
          ffront = 'redlum'   & fend = 'threebin'+fend & dothird=1
      END 
      ELSE: message,'No such type: '+ntostr(type)
  ENDCASE 

  psfile = '~/plots/deltasig_'+type+'_allband_bylum'+psstr+colorstr+pfitstr+endstr

  indir = dir +'u/'
  u1 = mrdfits(indir+ffront+'1'+fend,1)
  u2 = mrdfits(indir+ffront+'2'+fend,1)
  IF dothird THEN u3 = mrdfits(indir+ffront+'3'+fend,1)

  indir = dir +'g/'
  g1 = mrdfits(indir+ffront+'1'+fend,1)
  g2 = mrdfits(indir+ffront+'2'+fend,1)
  IF dothird THEN g3 = mrdfits(indir+ffront+'3'+fend,1)

  indir = dir +'r/'
  r1 = mrdfits(indir+ffront+'1'+fend,1)
  r2 = mrdfits(indir+ffront+'2'+fend,1)
  IF dothird THEN r3 = mrdfits(indir+ffront+'3'+fend,1)

  print,indir+ffront+'1'+fend
  print,indir+ffront+'2'+fend
  print,indir+ffront+'3'+fend

  indir = dir +'i/'
  i1 = mrdfits(indir+ffront+'1'+fend,1)
  i2 = mrdfits(indir+ffront+'2'+fend,1)
  IF dothird THEN i3 = mrdfits(indir+ffront+'3'+fend,1)

  indir = dir +'z/'
  z1 = mrdfits(indir+ffront+'1'+fend,1)
  z2 = mrdfits(indir+ffront+'2'+fend,1)
  IF dothird THEN z3 = mrdfits(indir+ffront+'3'+fend,1)

  IF keyword_set(plotfits) THEN BEGIN 
      plot_allband_lumgmcf_redo_fits, type, $
                      u1, u2, u3, $
                      g1, g2, g3, $
                      r1, r2, r3, $
                      i1, i2, i3, $
                      z1, z2, z3

      plot_allband_lumgmcf_print_table, type, $
                        u1, u2, u3, $
                        g1, g2, g3, $
                        r1, r2, r3, $
                        i1, i2, i3, $
                        z1, z2, z3
  ENDIF 

  IF keyword_set(dops) THEN BEGIN 
      
      IF keyword_set(encap) THEN BEGIN 
          begplot, name=psfile, color=color, xsize=11,ysize=8.5, /encap
      ENDIF ELSE BEGIN 
          begplot, name=psfile, color=color, /land
      ENDELSE 
  ENDIF 

  IF keyword_set(color) THEN BEGIN 
      IF !d.name EQ 'X' THEN BEGIN 
          color2 = !green
          color3 = !red
      ENDIF ELSE BEGIN 
          color2 = !blue
          color3 = !red
      ENDELSE 
  ENDIF ELSE BEGIN 
      color2 = !p.color
      color3 = !p.color
  ENDELSE 
  setup_mystuff

  line1 = 0
  line2 = 1
  line3 = 3
  IF NOT keyword_set(plotfits) THEN BEGIN 
      psym1 = 8
      psym2 = 1
      psym3 = 5
  ENDIF 

;  IF NOT keyword_set(plotfits) THEN BEGIN 
;      line1 = 0
;      line2 = 1
;      line3 = 3
;  ENDIF ELSE BEGIN 
;      psym1 = 8
;      psym2 = 1
;      psym3 = 5
;  ENDELSE 
  hatlength = !D.X_VSIZE / 200
  symsizes = [0.7, 1, 0.7]


  tags = strlowcase(tag_names(u1))
  IF keyword_set(wgm) THEN BEGIN 
      ytag = (where(tags EQ 'wgm_rebin'))[0]
      yerrtag = (where(tags EQ 'wgmerr_rebin'))[0]
      rtag = (where(tags EQ 'meanr_rebin'))[0]
      yrange=[2, 8000]
      ytitle = !wgmytitle      
  ENDIF ELSE BEGIN 

;      IF type EQ 'all' THEN BEGIN 
          ytag = (where(tags EQ 'sigma_rebin'))[0]
          yerrtag = (where(tags EQ 'sigmaerr_rebin'))[0]
          rtag = (where(tags EQ 'meanr_rebin'))[0]
;      ENDIF ELSE IF type EQ 'early' THEN BEGIN  
;          ytag = (where(tags EQ 'sigma'))[0]
;          yerrtag = (where(tags EQ 'sigmaerr'))[0]
;          rtag = (where(tags EQ 'meanr'))[0]
;      ENDIF 
      yrange=[0.1, 700]
      ytitle = !deltaytitle
  ENDELSE 

  xrange=[.015, 17]
  xtitle = !mpcxtitle2

  erase & multiplot, [3,2], /square

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; u-band plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  ploterror, u1.(rtag)/1000,u1.(ytag),u1.(yerrtag),/xlog,/ylog, $
             ytitle=ytitle, xrange=xrange, yrange=yrange,xstyle=2+1,ystyle=2+1, $
             ytickf='loglabels', xticklen=0.04, yticklen=0.04,$
             hatlength=hatlength, symsize=symsizes[0], psym=psym1, line=line1
  IF keyword_set(plotfits) THEN BEGIN 
      wuse=where(u1.yfit GT 0) & oplot,u1.meanr[wuse]/1000,u1.yfit[wuse]
  ENDIF ELSE BEGIN 
      oplot, u1.(rtag)/1000.,(u1.(ytag) > 0.01), line=line1
  ENDELSE 
             
  oploterror,u2.(rtag)/1000,u2.(ytag),u2.(yerrtag),color=color2,errc=color2,$
             hatlength=hatlength, symsize=symsizes[1], psym=psym2, line=line2
  IF keyword_set(plotfits) THEN BEGIN 
      wuse=where(u2.yfit GT 0) & oplot,u2.meanr[wuse]/1000,u2.yfit[wuse],color=color2, line=line2
  ENDIF ELSE BEGIN 
      oplot, u2.(rtag)/1000.,(u2.(ytag) > 0.01), line=line2, color=color2
  ENDELSE 

  IF n_elements(u3) ne 0 THEN BEGIN 
      oploterror,u3.(rtag)/1000,u3.(ytag),u3.(yerrtag),$
                 color=color3,errc=color3,$
                 hatlength=hatlength, symsize=symsizes[2], $
                 psym=psym3, line=line3
      IF keyword_set(plotfits) THEN BEGIN
          wuse=where(u3.yfit GT 0) & oplot,u3.meanr[wuse]/1000,u3.yfit[wuse],color=color3, line=line3
      ENDIF ELSE BEGIN 
          oplot, u3.(rtag)/1000.,(u3.(ytag) > 0.01), line=line3, color=color3
      ENDELSE  
  ENDIF 
  legend,'u',/right,charsize=2, box=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; g-band plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  multiplot
  ploterror, g1.(rtag)/1000,g1.(ytag),g1.(yerrtag),/xlog,/ylog, $
             xrange=xrange, yrange=yrange,xstyle=2+1,ystyle=2+1, $
             xticklen=0.04, yticklen=0.04,$
             hatlength=hatlength, symsize=symsizes[0], psym=psym1, line=line1
  IF keyword_set(plotfits) THEN BEGIN
      wuse=where(g1.yfit GT 0) & oplot,g1.meanr[wuse]/1000,g1.yfit[wuse]
  ENDIF ELSE BEGIN 
      oplot, g1.(rtag)/1000.,(g1.(ytag) > 0.01), line=line1
  ENDELSE 

  oploterror,g2.(rtag)/1000,g2.(ytag),g2.(yerrtag),color=color2,errc=color2,$
             hatlength=hatlength, symsize=symsizes[1], psym=psym2, line=line2
  IF keyword_set(plotfits) THEN BEGIN
      wuse=where(g2.yfit GT 0) & oplot,g2.meanr[wuse]/1000,g2.yfit[wuse],color=color2, line=oline2
  ENDIF ELSE BEGIN 
      oplot, g2.(rtag)/1000.,(g2.(ytag) > 0.01), line=line2, color=color2
  ENDELSE 

  IF n_elements(g3) ne 0 THEN BEGIN
      oploterror,g3.(rtag)/1000,g3.(ytag),g3.(yerrtag),color=color3,errc=color3,$
                 hatlength=hatlength, symsize=symsizes[2], psym=psym3, line=line3
      IF keyword_set(plotfits) THEN BEGIN
          wuse=where(g3.yfit GT 0) & oplot,g3.meanr[wuse]/1000,g3.yfit[wuse],color=color3, line=oline3
      ENDIF ELSE BEGIN 
      oplot, g3.(rtag)/1000.,(g3.(ytag) > 0.01), line=line3, color=color3
  ENDELSE 
  ENDIF 

  legend,'g',/right,charsize=2, box=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; r-band plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  multiplot, /doxaxis
  ploterror, r1.(rtag)/1000,r1.(ytag),r1.(yerrtag),/xlog,/ylog, $
             xrange=xrange, yrange=yrange,xstyle=2+1,ystyle=2+1, $
             xticklen=0.04, yticklen=0.04, xtitle=xtitle, xtickf='loglabels',$
             hatlength=hatlength, symsize=symsizes[0], psym=psym1, line=line1
  IF keyword_set(plotfits) THEN BEGIN
      wuse=where(r1.yfit GT 0) & oplot,r1.meanr[wuse]/1000,r1.yfit[wuse]
  ENDIF ELSE BEGIN 
      oplot, r1.(rtag)/1000.,(r1.(ytag) > 0.01), line=line1
  ENDELSE 

  oploterror,r2.(rtag)/1000,r2.(ytag),r2.(yerrtag),color=color2,errc=color2,$
             hatlength=hatlength, symsize=symsizes[1], psym=psym2, line=line2
  IF keyword_set(plotfits) THEN BEGIN
      wuse=where(r2.yfit GT 0) & oplot,r2.meanr[wuse]/1000,r2.yfit[wuse],color=color2, line=line2
  ENDIF ELSE BEGIN 
      oplot, r2.(rtag)/1000.,(r2.(ytag) > 0.01), line=line2, color=color2
  ENDELSE 

  IF n_elements(r3) ne 0 THEN BEGIN
      oploterror,r3.(rtag)/1000,r3.(ytag),r3.(yerrtag),color=color3,errc=color3,$
                 hatlength=hatlength, symsize=symsizes[2], psym=psym3, line=line3
      IF keyword_set(plotfits) THEN BEGIN
          wuse=where(r3.yfit GT 0) & oplot,r3.meanr[wuse]/1000,r3.yfit[wuse],color=color3, line=line3
      ENDIF ELSE BEGIN 
      oplot, r3.(rtag)/1000.,(r3.(ytag) > 0.01), line=line3, color=color3
  ENDELSE 
  ENDIF 
  legend,'r',/right,charsize=2, box=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; i-band plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  multiplot
  ploterror, i1.(rtag)/1000,i1.(ytag),i1.(yerrtag),/xlog,/ylog, $
             xrange=xrange, yrange=yrange,xstyle=2+1,ystyle=2+1, xtitle=xtitle, ytitle=ytitle, $
             ytickf='loglabels', xtickf='loglabels', xticklen=0.04, yticklen=0.04,$
             hatlength=hatlength, symsize=symsizes[0], psym=psym1, line=line1
  IF keyword_set(plotfits) THEN BEGIN
      wuse=where(i1.yfit GT 0) & oplot,i1.meanr[wuse]/1000,i1.yfit[wuse]
  ENDIF ELSE BEGIN 
      oplot, i1.(rtag)/1000.,(i1.(ytag) > 0.01), line=line1
  ENDELSE 

  oploterror,i2.(rtag)/1000,i2.(ytag),i2.(yerrtag),color=color2,errc=color2,$
             hatlength=hatlength, symsize=symsizes[1], psym=psym2, line=line2
  IF keyword_set(plotfits) THEN BEGIN
      wuse=where(i2.yfit GT 0) & oplot,i2.meanr[wuse]/1000,i2.yfit[wuse],color=color2, line=line2
  ENDIF ELSE BEGIN 
      oplot, i2.(rtag)/1000.,(i2.(ytag) > 0.01), line=line2, color=color2
  ENDELSE 

  IF n_elements(i3) ne 0 THEN BEGIN
      oploterror,i3.(rtag)/1000,i3.(ytag),i3.(yerrtag),color=color3,errc=color3,$
                 hatlength=hatlength, symsize=symsizes[2], psym=psym3, line=line3
      IF keyword_set(plotfits) THEN BEGIN
          wuse=where(i3.yfit GT 0) & oplot,i3.meanr[wuse]/1000,i3.yfit[wuse],color=color3, line=line3
      ENDIF ELSE BEGIN 
      oplot, i3.(rtag)/1000.,(i3.(ytag) > 0.01), line=line3, color=color3
  ENDELSE 
  ENDIF 
  legend,'i',/right,charsize=2, box=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; z-band plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  multiplot
  ploterror, z1.(rtag)/1000,z1.(ytag),z1.(yerrtag),/xlog,/ylog, $
             xrange=xrange, yrange=yrange,xstyle=2+1,ystyle=2+1, xtitle=xtitle, $
             xtickf='loglabels', xticklen=0.04, yticklen=0.04,$
             hatlength=hatlength, symsize=symsizes[0], psym=psym1, line=line1
  IF keyword_set(plotfits) THEN BEGIN
      wuse=where(z1.yfit GT 0) & oplot,z1.meanr[wuse]/1000,z1.yfit[wuse]
  ENDIF ELSE BEGIN 
      oplot, z1.(rtag)/1000.,(z1.(ytag) > 0.01), line=line1
  ENDELSE

  oploterror,z2.(rtag)/1000,z2.(ytag),z2.(yerrtag),color=color2,errc=color2,$
             hatlength=hatlength, symsize=symsizes[1], psym=psym2, line=line2
  IF keyword_set(plotfits) THEN BEGIN
      wuse=where(z2.yfit GT 0) & oplot,z2.meanr[wuse]/1000,z2.yfit[wuse],color=color2, line=line2
  ENDIF ELSE BEGIN 
      oplot, z2.(rtag)/1000.,(z2.(ytag) > 0.01), line=line2, color=color2
  ENDELSE 

  IF n_elements(z3) ne 0 THEN BEGIN
      oploterror,z3.(rtag)/1000,z3.(ytag),z3.(yerrtag),color=color3,errc=color3,$
                 hatlength=hatlength, symsize=symsizes[2], psym=psym3, line=line3
      IF keyword_set(plotfits) THEN BEGIN
          wuse=where(z3.yfit GT 0) & oplot,z3.meanr[wuse]/1000,z3.yfit[wuse],color=color3, line=line3
      ENDIF ELSE BEGIN 
      oplot, z3.(rtag)/1000.,(z3.(ytag) > 0.01), line=line3, color=color3
  ENDELSE 
  ENDIF 
  legend,'z',/right,charsize=2, box=0

  multiplot,/reset

  IF keyword_set(dops) THEN BEGIN 
      endplot
      IF NOT keyword_set(encap) THEN pslandfix, psfile
      IF keyword_set(encap) THEN set_bbox, psfile, '%%BoundingBox: 35 15 770 515'
  ENDIF 

  IF keyword_set(plotfits) THEN BEGIN 
      
      IF keyword_set(dops) THEN BEGIN 

          psfile = repstr(psfile, 'deltasig','deltasig_chisq')
          IF keyword_set(encap) THEN BEGIN 
              begplot, name=psfile, color=color, xsize=11,ysize=8.5, /encap
          ENDIF ELSE BEGIN 
              begplot, name=psfile, color=color, /land
          ENDELSE 
      ENDIF 

      IF !d.name EQ 'X' THEN key = get_kbrd(1)
      plot_allband_lumgmcf_plotall_chisq, type, $
                          u1, u2, u3, $
                          g1, g2, g3, $
                          r1, r2, r3, $
                          i1, i2, i3, $
                          z1, z2, z3, color=color

      IF keyword_set(dops) THEN BEGIN 
          endplot
          IF NOT keyword_set(encap) THEN pslandfix, psfile
          IF keyword_set(encap) THEN set_bbox, psfile, '%%BoundingBox: 35 15 770 515'
      ENDIF 

  ENDIF 

return
END 

