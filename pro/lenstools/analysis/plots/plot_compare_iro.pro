PRO plot_compare_iro

  dir = '~/lensout/combstripe/comb/sublum/r/'
  lum1=mrdfits(dir + 'lum1threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_N1.fit',1)
  lum2=mrdfits(dir + 'lum2threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_N1.fit',1)
  lum3=mrdfits(dir + 'lum3threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_N1.fit',1)

  idir = '~/iro/'

  ifile1dsig = idir+'gmx_17.0_to_21.6_0.1_withss_sdss.dat'
  ifile2dsig = idir+'gmx_21.6_to_22.2_0.1_withss_sdss.dat'

  ifile1xi = idir+'threed_17.0_to_21.6_0.1_withss_sdss.dat'
  ifile2xi = idir+'threed_21.6_to_22.2_0.1_withss_sdss.dat'
  readcol,ifile1xi,r3_1,junk,xi_1,xi_poisson_1,xi_jack_1,xi_errbest_1
  readcol,ifile2xi,r3_2,junk,xi_2,xi_poisson_2,xi_jack_2,xi_errbest_2

  readcol,ifile1dsig,meanr_1,junk,dsig_1,sigma_1,tsigma_1,dsig_poisson_1,dsig_jack_1,dsig_errbest_1
  readcol,ifile2dsig,meanr_2,junk,dsig_2,sigma_2,tsigma_2,dsig_poisson_2,dsig_jack_2,dsig_errbest_2

  fac = 8.28375e10
  xi_1 = xi_1/fac - 1.0
  xi_poisson_1 = xi_poisson_1/fac
  xi_jack_1 = xi_jack_1/fac
  xi_errbest_1 = xi_errbest_1/fac

  xi_2 = xi_2/fac - 1.0
  xi_poisson_2 = xi_poisson_2/fac
  xi_jack_2 = xi_jack_2/fac
  xi_errbest_2 = xi_errbest_2/fac

  forprint,r3_1,xi_1
  setup_mystuff

  !x.ticklen = 0.04
  !y.ticklen = 0.04

  xrange = [0.015,10.0]


  erase & multiplot, [1,2],/square
  yrange = [0.2, 100.0]
  ploterror, lum1.meanr/1000., lum1.sigma, lum1.sigmaerr,/xlog,/ylog,psym=8,$
             xrange=xrange,xstyle=1+2,yrange=yrange,ystyle=1+2,$
             ytickf='loglabels', ytitle=!deltaytitle
  axis, xaxis=1, xtitle=!mpcxtitle,xrange=xrange,xstyle=1+2, xtickf='loglabels'
  oplot, meanr_1, dsig_1, color=!green

  multiplot

  yrange = [0.2,6.e4]
  ploterror, lum1.r3, lum1.xi, lum1.xierr,/xlog,/ylog,psym=8,$
             xrange=xrange,xstyle=1+2,yrange=yrange,ystyle=1+2,$
             xtickf='loglabels',ytickf='loglabels',$
             ytitle=!xigmytitle,xtitle=!xigmxtitle
  oplot, r3_1, xi_1, color=!green
  multiplot,/reset


  key=get_kbrd(1)

  

  erase & multiplot, [1,2],/square
  yrange = [0.2, 300.0]
  ploterror, lum2.meanr/1000., lum2.sigma, lum2.sigmaerr,/xlog,/ylog,psym=8,$
             xrange=xrange,xstyle=1+2,yrange=yrange,ystyle=1+2,$
             ytickf='loglabels', ytitle=!deltaytitle
  axis, xaxis=1, xtitle=!mpcxtitle,xrange=xrange,xstyle=1+2, xtickf='loglabels'
  oplot, meanr_2, dsig_2, color=!green
  oplot,lum3.meanr_rebin/1000.,lum3.sigma_rebin,color=!red

  multiplot

  yrange = [0.5, 1.e5]
  ploterror, lum2.r3, lum2.xi, lum2.xierr,/xlog,/ylog,psym=8,$
             xrange=xrange,xstyle=1+2,yrange=yrange,ystyle=1+2,$
             xtickf='loglabels',ytickf='loglabels',$
             ytitle=!xigmytitle,xtitle=!xigmxtitle
  oplot, r3_2, xi_2, color=!green
  oplot, lum3.r3_rebin, lum3.xi_rebin, color=!red
  multiplot,/reset

  key=get_kbrd(1)

  

;  aploterror, 1, lum1.r3, lum1.xi, lum1.xierr,/xlog,/ylog,psym=8,$
;              yrange=yrange, ystyle=1+2,$
;             ytitle=!xigmytitle,xtitle=!xigmxtitle
;  oploterror, lum2.r3, lum2.xi, lum2.xierr, color=!green, psym=8,errc=!green

  aploterror, 1, lum1.r3_rebin, lum1.xi_rebin, lum1.xierr_rebin,$
              /xlog,/ylog,psym=8,$
              yrange=yrange, ystyle=1+2,$
             ytitle=!xigmytitle,xtitle=!xigmxtitle
  oploterror, lum2.r3_rebin, lum2.xi_rebin, lum2.xierr_rebin, $
              color=!green, psym=8,errc=!green

  
  oplot, r3_1, xi_1
  oplot, r3_2, xi_2, color=!green
  oplot, lum3.r3_rebin, lum3.xi_rebin, color=!red

  !x.ticklen = 0.02
  !y.ticklen = 0.02


  ;; print out ascii files
  outdir = '/net/cheops1/data3/lensout/combstripe/comb/sublum/r/xigm_ascii/'
  colprint,lum1.r3,lum1.xi,lum1.xierr,$
           file=outdir+'lum1_xi_r.dat'
  colprint,lum2.r3,lum2.xi,lum2.xierr,$
           file=outdir+'lum2_xi_r.dat'
  colprint,lum3.r3,lum3.xi,lum3.xierr,$
           file=outdir+'lum3_xi_r.dat'

  colprint,lum1.meanr/1000.,lum1.sigma,lum1.sigmaerr,$
           file=outdir+'lum1_deltasig_r.dat'
  colprint,lum2.meanr/1000.,lum2.sigma,lum2.sigmaerr,$
           file=outdir+'lum2_deltasig_r.dat'
  colprint,lum3.meanr/1000.,lum3.sigma,lum3.sigmaerr,$
           file=outdir+'lum3_deltasig_r.dat'


END 
