PRO compare_frac, stripe, clr, shstruct

  tsh = shstruct

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
  dens[1].density[2] = 1.3028618  ;r
  dens[1].density[3] = 1.0964445  ;i

  w=where(dens.stripe EQ 82)

  ;; density is in #/kpc^2, convert to #/arcmin^2
  tsh.density = tsh.density*(angdist_lambda(tsh.zmean)*1000.)^2 ;#/radians^2
  tsh.density = tsh.density*(!pi/180.)^2 ;#/deg^2
  tsh.density = tsh.density*(1./60.)^2   ;#/arcmin^2

;  tsh.area = tsh.area/(angdist_lambda(tsh.zmean)*1000.)^2
;  tsh.area = tsh.area*(180./!pi)^2
;  tsh.area = tsh.area*(60.)^2

  frac = tsh.density/dens[w].density[clr]
  corr = frac-1.

  print,corr

END 
