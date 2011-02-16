PRO wtheta_test, lens10, rand10, lens82, rand82


  indir = '/sdss4/data1/esheldon/GAL_GAL/spectra/'

  nend = 'N2.fit'

  lensfiles = indir+['wtheta_stripe10_r_sum_'+nend,'wtheta_stripe82_r_sum_'+nend]
  randfiles = indir+['wtheta_zrand_stripe10_r_sum_'+nend,'wtheta_zrand_stripe82_r_sum_'+nend]

  combine_zshear, lensfiles[0], lens10
  combine_zshear, randfiles[0], rand10

  frac_overdense, lens10.npair, lens10.nlenses, rand10.npair, rand10.nlenses, $
                  frac10, frac10err

  ploterror, lens10.meanr, frac10, frac10err, psym=1

  key=get_kbrd(1)

  combine_zshear, lensfiles[1], lens82
  combine_zshear, randfiles[1], rand82

  frac_overdense, lens82.npair, lens82.nlenses, rand82.npair, rand82.nlenses, $
                  frac82, frac82err

  ploterror, lens82.meanr, frac82, frac82err, psym=1

END 
