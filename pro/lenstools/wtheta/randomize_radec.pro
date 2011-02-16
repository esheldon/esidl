PRO randomize_radec,t

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: randomize_radec,t'
      return
  ENDIF 

  print,'Randomizing ra/dec'

;  COMMON wthetarand_block, seed, etarange, etarange_spectra, $
;    mlambda, meta, FIRSTLAMBDA, LASTLAMBDA, Mlimit

  IF NOT tag_exist(t, 'ra') THEN BEGIN 

      add_tags, t, ['ra','dec'], ['0d','0d'], newt
      delvarx,t
      t=temporary(newt)

      survey2eq, t.lambda, t.eta, tra, tdec
      t.ra = tra
      t.dec = tdec

  ENDIF 

  seed=long(systime(1))

  nn=n_elements(t)

  maxra=max(t.ra, min=minra)
  maxdec=max(t.dec, min=mindec)

  ra = double( arrscl(randomu(seed, nn), minra, maxra, $
                      arrmin=0., arrmax=1.) )
  dec = double( arrscl(randomu(seed,nn), mindec, maxdec, $
                       arrmin=0., arrmax=1.) )

  print,'Sorting ras'
  s=sort(ra)

  t.ra = ra
  t.dec = dec

  t=temporary(t[s])

END 
