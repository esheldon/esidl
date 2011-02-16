PRO print_for_iro

  outdir = '~/lensout/iro/xigm_ascii/'

  ;; All
  indir = '~/lensout/combstripe/comb/'
  infile = indir + 'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
  t=mrdfits(infile,1)

  outfile = outdir + 'deltasig_all.dat'
  print,'Outfile: ',outfile
  colprint, t.meanr/1000., t.sigma, t.sigmaerr, file=outfile

  outfile = outdir + 'xigm_all.dat'
  print,'Outfile: ',outfile
  colprint, t.r3, t.xi, t.xierr, file=outfile

  outfile = outdir + 'xigm_all_cov.fit'
  print,'Outfile: ',outfile
  

  ;; Vlim
  indir = '~/lensout/combstripe/comb/'
  infile = indir + 'vlim3_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
  t=mrdfits(infile,1)

  outfile = outdir + 'deltasig_vlim3.dat'
  print,'Outfile: ',outfile
  colprint, t.meanr/1000., t.sigma, t.sigmaerr, file=outfile

  outfile = outdir + 'xigm_vlim3.dat'
  print,'Outfile: ',outfile
  colprint, t.r3, t.xi, t.xierr, file=outfile

  ;; Red Galaxies
  indir = '~/lensout/combstripe/comb/'
  infile = indir + 'gmr2_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
  t=mrdfits(infile,1)

  outfile = outdir + 'deltasig_red.dat'
  print,'Outfile: ',outfile
  colprint, t.meanr/1000., t.sigma, t.sigmaerr, file=outfile

  outfile = outdir + 'xigm_red.dat'
  print,'Outfile: ',outfile
  colprint, t.r3, t.xi, t.xierr, file=outfile

  ;; Blue Galaxies
  indir = '~/lensout/combstripe/comb/'
  infile = indir + 'gmr1_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
  t=mrdfits(infile,1)

  outfile = outdir + 'deltasig_blue.dat'
  print,'Outfile: ',outfile
  colprint, t.meanr/1000., t.sigma, t.sigmaerr, file=outfile

  outfile = outdir + 'xigm_blue.dat'
  print,'Outfile: ',outfile
  colprint, t.r3, t.xi, t.xierr, file=outfile


  ;; Luminosity Bins
  FOR clr=0,4 DO BEGIN 

      ;; Lum1
      clrstr = !colors[clr]

      indir = '~/lensout/combstripe/comb/sublum/'+clrstr+'/'
      infile = indir + 'lum1threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
      t=mrdfits(infile,1)
      
      outfile = outdir + 'deltasig_lum1_'+clrstr+'.dat'
      print,'Outfile: ',outfile
      colprint, t.meanr/1000., t.sigma, t.sigmaerr, file=outfile
      
      outfile = outdir + 'xigm_lum1_'+clrstr+'.dat'
      print,'Outfile: ',outfile
      colprint, t.r3, t.xi, t.xierr, file=outfile
      
      ;; Lum2
      indir = '~/lensout/combstripe/comb/sublum/'+clrstr+'/'
      infile = indir + 'lum2threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
      t=mrdfits(infile,1)
      
      outfile = outdir + 'deltasig_lum2_'+clrstr+'.dat'
      print,'Outfile: ',outfile
      colprint, t.meanr/1000., t.sigma, t.sigmaerr, file=outfile
      
      outfile = outdir + 'xigm_lum2_'+clrstr+'.dat'
      print,'Outfile: ',outfile
      print,'Outfile: ',outfile
      colprint, t.r3, t.xi, t.xierr, file=outfile
      
      ;; Lum3
      indir = '~/lensout/combstripe/comb/sublum/'+clrstr+'/'
      infile = indir + 'lum3threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
      t=mrdfits(infile,1)
      
      outfile = outdir + 'deltasig_lum3_'+clrstr+'.dat'
      print,'Outfile: ',outfile
      colprint, t.meanr/1000., t.sigma, t.sigmaerr, file=outfile
      
      outfile = outdir + 'xigm_lum3_'+clrstr+'.dat'
      print,'Outfile: ',outfile
      colprint, t.r3, t.xi, t.xierr, file=outfile
  ENDFOR 

  ;; Red Luminosity Bins
  FOR clr=0,4 DO BEGIN 

      ;; Lum1
      clrstr = !colors[clr]

      indir = '~/lensout/combstripe/comb/sublum/'+clrstr+'/'
      infile = indir + 'redlum1threebin_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
      t=mrdfits(infile,1)
      
      outfile = outdir + 'deltasig_redlum1_'+clrstr+'.dat'
      print,'Outfile: ',outfile
      colprint, t.meanr/1000., t.sigma, t.sigmaerr, file=outfile
      
      outfile = outdir + 'xigm_redlum1_'+clrstr+'.dat'
      print,'Outfile: ',outfile
      colprint, t.r3, t.xi, t.xierr, file=outfile
      
      ;; Lum2
      indir = '~/lensout/combstripe/comb/sublum/'+clrstr+'/'
      infile = indir + 'redlum2threebin_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
      t=mrdfits(infile,1)
      
      outfile = outdir + 'deltasig_redlum2_'+clrstr+'.dat'
      print,'Outfile: ',outfile
      colprint, t.meanr/1000., t.sigma, t.sigmaerr, file=outfile
      
      outfile = outdir + 'xigm_redlum2_'+clrstr+'.dat'
      print,'Outfile: ',outfile
      print,'Outfile: ',outfile
      colprint, t.r3, t.xi, t.xierr, file=outfile
      
      ;; Lum3
      indir = '~/lensout/combstripe/comb/sublum/'+clrstr+'/'
      infile = indir + 'redlum3threebin_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit'
      t=mrdfits(infile,1)
      
      outfile = outdir + 'deltasig_redlum3_'+clrstr+'.dat'
      print,'Outfile: ',outfile
      colprint, t.meanr/1000., t.sigma, t.sigmaerr, file=outfile
      
      outfile = outdir + 'xigm_redlum3_'+clrstr+'.dat'
      print,'Outfile: ',outfile
      colprint, t.r3, t.xi, t.xierr, file=outfile
  ENDFOR 


END 
