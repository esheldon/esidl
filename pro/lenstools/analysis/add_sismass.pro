PRO add_sismass, lensumfile, stripe, clr, outstruct

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: add_sismass, lensumfile, stripe, clr, [outstruct]'
      return
  ENDIF 

  nstripe = 2
  dens = create_struct('stripe', 0L, 'density', fltarr(5))
  dens = replicate(dens, nstripe)
  ;; stripe 10
  dens[0].stripe = 10
  dens[0].density[1] = 1.2312680 ;g # per square arcminute
  dens[0].density[2] = 1.7684189 ;r
  dens[0].density[3] = 1.5075363 ;i

  ;; stripe 82
  dens[1].stripe = 82
  dens[1].density[1] = 0.80527492 ;g
  dens[1].density[2] = 1.3028618 ;r
  dens[1].density[3] = 1.0964445 ;i

  w=where(stripe EQ dens.stripe,nw)
  IF nw EQ 0 THEN BEGIN 
      print,'No density info on stripe '+ntostr(stripe)
      return
  ENDIF 

  outfile = repstr(lensumfile, 'lensum', 'lensum_massadd')
  print
  print,'Output file: ',outfile
  print

  lensum = mrdfits(lensumfile, 1, hdr)

  hval = sxpar(hdr, 'H')
  rminkpc = sxpar(hdr, 'RMINKPC')
  rmaxkpc = sxpar(hdr, 'RMAXKPC')
  binsize = sxpar(hdr, 'BINWIDTH')
  
  IF (hval EQ 0) OR (rminkpc EQ 0) OR (rmaxkpc EQ 0) OR (binsize EQ 0) THEN BEGIN 
      print,'Bad header info'
      return
  ENDIF 

  fxhclean, hdr

  lensum_sis_fit, temporary(lensum), dens[w].density[clr], $
    binsize, rminkpc, rmaxkpc, hval, $
    outstruct

  mwrfits, outstruct, outfile, hdr, /create

  return
END 
