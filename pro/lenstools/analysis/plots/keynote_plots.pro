PRO keynote_plots

  dir = '~/lensout/combstripe/comb/sublum/r/'

  flum1 = dir + 'lum1threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
  flum2 = dir + 'lum2threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
  flum3 = dir + 'lum3threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'

  lum1 = mrdfits(flum1,1)
  lum2 = mrdfits(flum2,1)
  lum3 = mrdfits(flum3,1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; r-band three bin deltasig
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  file = '~/plots/deltasig_rband_bylum_color.eps'
  begplot, name=file, /encap, /color
  !p.charsize=2.0

  setup_mystuff

  colors = [!p.color, !blue, !red]
  lines = [0, 1, 3]

  xtitle = !mpcxtitle2
  ytitle = !deltaytitle
  xrange = [0.025, 15.0]
  yrange = [0.3, 6.e2]
  aploterror, 1, $
              lum1.meanr_rebin/1000, lum1.sigma_rebin, lum1.sigmaerr_rebin, $
              /xlog, /ylog, xrange=xrange, yrange=yrange, $
              xstyle=1+2, ystyle=1+2, $
              xtitle = xtitle, ytitle=ytitle, $
              xtickf='loglabels', ytickf='loglabels', $
              xticklen=0.04,yticklen=0.04
  oploterror, lum2.meanr_rebin/1000, lum2.sigma_rebin, lum2.sigmaerr_rebin, $
              color=colors[1], errc=colors[1], $
              line = lines[1]
  oploterror, lum3.meanr_rebin/1000, lum3.sigma_rebin, lum3.sigmaerr_rebin, $
              color=colors[2], errc=colors[2], $
              line = lines[2]

  legend, 'r',/right,box=0, charsize=2

  endplot, /trim_bbox

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; r-band three bin xi
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  file = '~/plots/xi_rband_bylum_color.eps'
  begplot, name=file, /encap, /color
  !p.charsize=2.0
  setup_mystuff

  ytitle = !xigmytitle
  xtitle = !xigmxtitle
  xrange = [0.015, 15.0]
  yrange=[0.11, 1e5]

  aploterror, 1, $
              lum1.r3_rebin, lum1.xi_rebin, lum1.xierr_rebin, $
              /xlog, /ylog, xrange=xrange, yrange=yrange, $
              xstyle=1+2, ystyle=1+2, $
              xtitle = xtitle, ytitle=ytitle, $
              xtickf='loglabels', ytickf='loglabels', $
              xticklen=0.04,yticklen=0.04
  oploterror, lum2.r3_rebin, lum2.xi_rebin, lum2.xierr_rebin, $
              color=colors[1], errc=colors[1], $
              line = lines[1]
  oploterror, lum3.r3_rebin, lum3.xi_rebin, lum3.xierr_rebin, $
              color=colors[2], errc=colors[2], $
              line = lines[2]

  legend, 'r',/right,box=0, charsize=2

  endplot,/trim_bbox

END 
