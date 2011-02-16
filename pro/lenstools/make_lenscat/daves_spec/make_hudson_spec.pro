PRO make_hudson_spec, lcat, scat

  ;; Read the lens catalog, match to corresponding source catalog ascii file
  ;; and output the required info to an ascii file

  stripes = [9,10,11,12,13,14,15]
  stripestr = stripearr2string(stripes)

  outdir = sdssidl_config('shapecorr_dir')+'combined/'
  lcatfile = outdir + 'stripe'+stripestr+'_spec_hudson.dat'
  srcfits = outdir + 'stripe'+stripestr+'_hudson.fit'

  nlcat = n_elements(lcat)
  nscat = n_elements(scat)
  IF nlcat EQ 0 THEN get_spectra_lcat, stripes, lcat, /spec
  IF nscat EQ 0 THEN scat = mrdfits(srcfits, 1)

  nlcat = n_elements(lcat)
  nscat = n_elements(scat)


  ;; first match by id
  sphoto_match, lcat, scat, mlcat, mscat

  help,mscat,mlcat

  l_unmatched = lindgen(nlcat)
  s_unmatched = lindgen(nscat)

  remove, mlcat, l_unmatched
  remove, mscat, s_unmatched

  csurvey2eq, lcat.clambda, lcat.ceta, lra, ldec
  csurvey2eq, scat.clambda, scat.ceta, sra, sdec

  close_match_radec, lra[l_unmatched], ldec[l_unmatched], $
                     sra[s_unmatched], sdec[s_unmatched], $
                     ml, ms, 2.0/3600., 1
  
  help,ml,ms
  
  match, mlcat, l_unmatched[ml], m1, m2
  help,m1,m2
  match, mscat, s_unmatched[ms], m1, m2
  help,m1,m2

  source_matches = [mscat, s_unmatched[ms]]
  lens_matches   = [mlcat, l_unmatched[ml]]

  help, source_matches, rem_dup(source_matches)
  help, lens_matches,   rem_dup(lens_matches)

  rmd = rem_dup( source_matches )

  source_matches = source_matches[rmd]
  lens_matches = lens_matches[rmd]

  help, source_matches, lens_matches

  matches = replicate(-1L, nlcat)
  matches[lens_matches] = source_matches

  ;; Re-correct

  print,'Writing lenses to file: ',lcatfile
  openw, lun, lcatfile, /get_lun

  format = $
    '('+$
    '2(D17.11,:,1X),' + $       ; lambda, eta
    '25(F0,:,1X),'+$
    '(I0,:,1X)' + $
    ')'
  
  ;; Rotate shapes
  e1 = lcat.m_e1_corr_h[2]
  e2 = lcat.m_e2_corr_h[2]
  e1e1err = lcat.m_e1e1err[2]
  e1e2err = lcat.m_e1e2err[2]
  e2e2err = lcat.m_e2e2err[2]
  rotation = lcat.rotation[2]

  rotate_e1e2error, rotation, $
                    e1, e2, $
                    e1e1err, e1e2err, e2e2err,$
                    e1out, e2out, $
                    e1e1errout, e1e2errout, e2e2errout

  petro = lcat.petrocounts - lcat.reddening

  tm = systime(1)
  FOR i=0L, nlcat-1 DO BEGIN 
      
      printf, lun, $
              lcat[i].clambda, lcat[i].ceta, $
              e1out[i], e2out[i], $
              e1e1errout[i], e1e2errout[i], e2e2errout[i], $
              lcat[i].photoz_z, lcat[i].photoz_zerr, $
              lcat[i].photoz_abscounts[0], $
              lcat[i].photoz_abscounts[1], $
              lcat[i].photoz_abscounts[2], $
              lcat[i].photoz_abscounts[3], $
              lcat[i].photoz_abscounts[4], $
              petro[0,i], $
              petro[1,i], $
              petro[2,i], $
              petro[3,i], $
              petro[4,i], $
              lcat[i].petrorad[2], $
              lcat[i].z, $
              lcat[i].z_err, $
              lcat[i].absmag[0], $
              lcat[i].absmag[1], $
              lcat[i].absmag[2], $
              lcat[i].absmag[3], $
              lcat[i].absmag[4], $
              matches[i], $
              format=format

  ENDFOR 
  free_lun, lun
  ptime,systime(1)-tm


END 
