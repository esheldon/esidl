PRO read_hudson_spec, lcat

  stripes = [9,10,11,12,13,14,15]
  stripestr = stripearr2string(stripes)

  outdir = sdssidl_config('shapecorr_dir')+'combined/'
  lcatfile = outdir + 'stripe'+stripestr+'_spec_hudson.dat'

  instruct = create_struct('clambda', 0d, $
                           'ceta',    0d, $
                           'e1', 0.0, $
                           'e2', 0.0, $
                           'e1e1err',   0.0, $
                           'e1e2err',   0.0, $
                           'e2e2err',   0.0, $
                           'photoz_z',  0.0, $
                           'photoz_zerr', 0.0, $
                           'photoz_abscounts_u', 0.0, $
                           'photoz_abscounts_g', 0.0, $
                           'photoz_abscounts_r', 0.0, $
                           'photoz_abscounts_i', 0.0, $
                           'photoz_abscounts_z', 0.0, $
                           'upetro', 0.0, $
                           'gpetro', 0.0, $
                           'rpetro', 0.0, $
                           'ipetro', 0.0, $
                           'zpetro', 0.0, $
                           'rpetrorad', 0.0, $
                           'z', 0.0, $
                           'z_err', 0.0, $
                           'uabsmag', 0.0, $
                           'gabsmag', 0.0, $
                           'rabsmag', 0.0, $
                           'iabsmag', 0.0, $
                           'zabsmag', 0.0, $
                           'source_match', 0L)

  print,'Counting lines in: ',lcatfile
  nl = numlines(lcatfile)
  lcat = replicate(instruct, nl)

  print,'Reading file'
  openr, lun, lcatfile, /get_lun

  tt = systime(1)
  readf, lun, lcat
  ptime,systime(1)-tt
return
  
  ;; would have to do the following if there were strings perhaps

  tt = systime(1)
  source_match = 0L
  clambda = 0d
  ceta = 0d
  FOR i=0L, nl-1 DO BEGIN 
      readf, lun, $
             clambda, ceta, $
             e1, e2, e1e1err, e1e2err, e2e2err, $
             phz, phzerr, $
             phz_absu, phz_absg, phz_absr, phz_absi, phz_absz, $
             up, gp, rp, ip, zp, $
             rpetrorad, $
             z, z_err, $
             uabs, gabs, rabs, iabs, zabs, $
             source_match

      lcat[i].clambda = clambda
      lcat[i].ceta = ceta
      lcat[i].e1 = e1
      lcat[i].e2 = e2
      lcat[i].e1e1err = e1e1err
      lcat[i].e1e2err = e1e2err
      lcat[i].e2e2err = e2e2err
      lcat[i].photoz_z = phz
      lcat[i].photoz_zerr = phzerr
      lcat[i].photoz_abscounts_u = phz_absu
      lcat[i].photoz_abscounts_g = phz_absg
      lcat[i].photoz_abscounts_r = phz_absr
      lcat[i].photoz_abscounts_i = phz_absi
      lcat[i].photoz_abscounts_z = phz_absz
      lcat[i].upetro = up
      lcat[i].gpetro = gp
      lcat[i].rpetro = rp
      lcat[i].ipetro = ip
      lcat[i].zpetro = zp
      lcat[i].rpetrorad = rpetrorad
      lcat[i].z = z
      lcat[i].z_err = z_err
      lcat[i].uabsmag = uabs
      lcat[i].gabsmag = gabs
      lcat[i].rabsmag = rabs
      lcat[i].iabsmag = iabs
      lcat[i].zabsmag = zabs
      lcat[i].source_match = source_match

  ENDFOR 
  ptime,systime(1) - tt

  return

  tt = systime(1)
  rdfloatstr, lcatfile, instruct, lcat, /double
  ptime,systime(1) - tt

END 
	
