PRO plot_vdisxi, dops=dops, encap=encap, color=color

  indir = esheldon_config("lensout_dir")+'combstripe/comb/'

  vdis1 = mrdfits(indir+'vdis1_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit',1)
  vdis2  = mrdfits(indir+'vdis2_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit',1)

  IF keyword_set(color) THEN color_str = '_color' ELSE color_str=''

  psfile = '~/plots/xi_vdis'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=7, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 

  IF !d.name EQ 'X' THEN BEGIN
      bold = !p.background
      cold = !p.color
      !p.background = !white
      !p.color = !black
  ENDIF

  IF keyword_set(color) THEN BEGIN 
      vdis2_clr = !red
      vdis1_clr = !blue
  ENDIF ELSE BEGIN 
      vdis2_clr = !p.color
      vdis1_clr = !p.color
  ENDELSE 

  setup_mystuff

  yrange = [0.1, 7.e4]
  xrange = [0.015, 12]
  aplot, 1, [0],[0],/nodata,line=0, /xlog, /ylog, $
         xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, xtickf='loglabels', ytickf='loglabels', $
         xtitle=!xigmxtitle, ytitle=!xigmytitle, xticklen=0.04, yticklen=0.04

  wplot = where(vdis1.xi_rebin GT 0)
  wplot=lindgen(n_elements(vdis1.xi_rebin))
  oploterror, vdis2.r3_rebin, vdis2.xi_rebin, vdis2.xierr_rebin, $
              line=2, color=vdis2_clr, errc=vdis2_clr
  oploterror, vdis1.r3_rebin[wplot], (vdis1.xi_rebin[wplot] > 0.01), vdis1.xierr_rebin[wplot], $
              line=0, color=vdis1_clr, errc=vdis1_clr

;  IF !d.name EQ 'X' THEN key=get_kbrd(1)

;  aplot, 1, [0],[0],/nodata,line=0, /xlog, /ylog, $
;         xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, xtickf='loglabels', ytickf='loglabels', $
;         xtitle=!xigmxtitle, ytitle=!xigmytitle, xticklen=0.04, yticklen=0.04
;  wplot = where(vdis1.xi GT 0)
;  wplot=lindgen(n_elements(vdis1.xi))
;  oploterror, vdis2.r3, vdis2.xi, vdis2.xierr, $
;              line=2, color=vdis2_clr, errc=vdis2_clr
;  oploterror, vdis1.r3[wplot], (vdis1.xi[wplot] > 0.01), vdis1.xierr[wplot], $
;              line=0, color=vdis1_clr, errc=vdis1_clr

  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  IF keyword_set(dops) THEN endplot
  IF keyword_set(encap) THEN set_bbox, psfile, '%%BoundingBox: 35 15 480 445'

  

  ;; both together

  psfile = '~/plots/deltasig_xi_vdis'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+color_str+'.eps' 
          begplot, name=psfile, color=color, xsize=7, ysize=10, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+color_str+'.ps'
          begplot, name=psfile, color=color
      ENDELSE 
  ENDIF 

  erase & multiplot, [1,2], /square

  yrange = [0.3, 200]
  plot, [0],[0],/nodata,/xlog, /ylog, $
         xrange=xrange, yrange=yrange, $
         xstyle=1+2, ystyle=1+2, $
         ytickf='loglabels', $
         ytitle=!deltaytitle, $
         xticklen=0.04, yticklen=0.04

  oploterror, vdis2.meanr_rebin/1000,vdis2.sigma_rebin,vdis2.sigmaerr_rebin, $
              line=2, color=vdis2_clr, errc=vdis2_clr
  oploterror, vdis1.meanr_rebin/1000,vdis1.sigma_rebin,vdis1.sigmaerr_rebin, $
              line=0, color=vdis1_clr, errc=vdis1_clr
  axis, xaxis=1, xtitle=!mpcxtitle, xrange=xrange, xstyle=1+2, $
        xtickf='loglabels'

  multiplot

  yrange = [0.1, 7.e4]
  plot, [0],[0],/nodata,line=0, /xlog, /ylog, $
         xrange=xrange, yrange=yrange, xstyle=1+2, ystyle=1+2, $
         xtickf='loglabels', ytickf='loglabels', $
         xtitle=!xigmxtitle, ytitle=!xigmytitle, xticklen=0.04, yticklen=0.04

  wplot = where(vdis1.xi_rebin GT 0)
  wplot=lindgen(n_elements(vdis1.xi_rebin))
  oploterror, vdis2.r3_rebin, vdis2.xi_rebin, vdis2.xierr_rebin, $
              line=2, color=vdis2_clr, errc=vdis2_clr
  oploterror,vdis1.r3_rebin[wplot],$
             (vdis1.xi_rebin[wplot] > 0.01),$
             vdis1.xierr_rebin[wplot], $
             line=0, color=vdis1_clr, errc=vdis1_clr

  multiplot,/reset

  IF !d.name EQ 'X' THEN BEGIN
      !p.background=bold
      !p.color=cold
  ENDIF 

  IF keyword_set(dops) THEN endplot
  IF keyword_set(encap) THEN set_bbox, psfile,'%%BoundingBox: 35 15 420 730'
 
END 
