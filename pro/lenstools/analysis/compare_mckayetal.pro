PRO compare_mckayetal, dops=dops

  IF !d.name EQ 'PS' THEN color=!blue ELSE color=!green

  indir = '~/lensout/combstripe/comb/sublum/r/'

  f1=indir + 'lum1threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_N1.fit'
  f2=indir + 'lum2threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_N1.fit'
  f3=indir + 'lum3threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_N1.fit'

;  f1=indir + 'lum1threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
;  f2=indir + 'lum2threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
;  f3=indir + 'lum3threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'


  l1 = mrdfits(f1,1)
  l2 = mrdfits(f2,1)
  l3 = mrdfits(f3,1)


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read info about inverse critical density 
  ;; make a fix factor from it and Uros, star-galaxy
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  scf = '/net/cheops1/data0/corrected/sigmacrit/sigcritinv_func_stripe27_28_29_30_31_32_33_34_35_36_37_gri_h_lambda.fit'
  sc = mrdfits(scf,1)

  wf = '~/lensout/weights/weight_vs_z.fit'
  wst = mrdfits(wf,1)

  scinv = interpol(sc.sigcritinv, sc.zlens, wst.z)
  scinv_old2 = interpol(sc.sigcritinv2, sc.zlens, wst.z)

  aplot, !gratio, wst.z, scinv_old2
  oplot, wst.z, scinv, line=2

  mscinv = $
    qgauss(scinv*wst.num, wst.z,1000)/qgauss(wst.num,wst.z,1000)
  mscinv_old2 = $
    qgauss(scinv_old2*wst.num, wst.z,1000)/qgauss(wst.num,wst.z,1000)

  fixfac = mscinv_old2/mscinv
  print,fixfac

  ;; also include Hirara & Seljak fixes
  fixfac = fixfac*1.1
  print,fixfac
  ;; Assume there were some stars
  fixfac = fixfac*1.1
  print,fixfac


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Mckay et al. numbers
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  mk = [1.8, 7.4, 6.6, 16.1]
  mkerr = [0.4, 1.0, 1.7, 3.7]
  lk = [0.991, 3.070, 4.95, 8.23]
  mlk = [166., 228., 126., 187.]
  mlkerr = [32., 31., 33., 43.]

  ;; Mass is in units of solar masses. Convert to units of 10^10 
  m1 = interpol(l1.sismass, l1.meanr, 260.0)/1.e10
  m1err = interpol(l1.sismasserr, l1.meanr, 260.0)/1.e10
  ;; luminosity is in units of 10^10
  lum1 = mean(l1.meanlum[0:7])
  ml1 = m1/lum1
  ml1err = m1err/lum1

  m2 = interpol(l2.sismass, l2.meanr, 260.0)/1.e10
  m2err = interpol(l2.sismasserr, l2.meanr, 260.0)/1.e10
  lum2 = mean(l2.meanlum[0:7])
  ml2 = m2/lum2
  ml2err = m2err/lum2

  m3 = interpol(l3.sismass, l3.meanr, 260.0)/1.e10
  m3err = interpol(l3.sismasserr, l3.meanr, 260.0)/1.e10
  lum3 = mean(l3.meanlum[0:7])
  ml3 = m3/lum3
  ml3err = m3err/lum3

  m = [m1, m2, m3]
  merr = [m1err, m2err, m3err]
  lum = [lum1, lum2, lum3]
  ml = [ml1, ml2, ml3]
  mlerr = [ml1err, ml2err, ml3err]

  print,'          lum         mass      masserr          m/l      m/l_err'
  colprint, lum, m/100., merr/100., ml, mlerr

  ;;;;;;;;;;;;;;;;;
  ;; Plots
  ;;;;;;;;;;;;;;;;;

  yt1 = 'M [10!U12!N M'+sunsymbol()+' ]'
  yt2 = 'M/L [solar units]'
  xt = 'Luminosity [10!U10!N L'+sunsymbol()+' ]'

  !p.multi=[0,0,2]

  ploterror, lk, mk, mkerr, psym=8, yrange=[0,40],$
    xtitle=xt, ytitle=yt1
  oplot, lk, mk
  oploterror, $
    lum, m/100, merr/100, $
    psym=8, color=color, errc=color
  oplot, lum, m/100, color=color

  ploterror, lk, mlk, mlkerr, psym=8, yrange=[0,600], $
    xtitle=xt, ytitle=yt2
  oplot, lk, mlk
  oplot, lum, ml, color=color
  oploterror, $
    lum, ml, mlerr,$
    psym=8, color=color,errc=color

  !p.multi=0

  key = prompt_kbrd()

  ;; Now see if we can account for difference in the inverse critical
  ;; densities.
  ;; integrate inverse critical density over the lens redshift
  ;; distribution: interpolate to lens redshift histogram z's



  key = prompt_kbrd()

  !p.multi=[0,0,2]

  ploterror, lk, mk*fixfac, mkerr*fixfac, psym=8, yrange=[0,40],$
    xtitle=xt, ytitle=yt1
  oplot, lk, mk*fixfac
  oploterror, $
    lum, m/100, merr/100, $
    psym=8, color=color, errc=color
  oplot, lum, m/100, color=color

  ploterror, lk, mlk*fixfac, mlkerr*fixfac, psym=8, yrange=[0,600], $
    xtitle=xt, ytitle=yt2
  oplot, lk, mlk*fixfac
  oplot, lum, ml, color=color
  oploterror, $
    lum, ml, mlerr,$
    psym=8, color=color,errc=color

  legend,'Fix factor = '+ntostr(fixfac,4,/round),/right,box=0

  !p.multi=0

END 
