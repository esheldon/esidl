PRO runquadshear_model, ell, q, etan, etan1, etan2

  ;; run tanshear_mom for various ellipticities

  ell = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,$
         0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95 ]

  nell = n_elements(ell)

  maxr = 100
  bin = 10
  grid=1.0
  
  delvarx,q,etan,etan1,etan2
  FOR i=0L, nell-1 DO BEGIN 

      ee = ell[i]
      print,'e = ',ee
      
      quadshear_model, ee, maxr, bin, grid, rr, tetan, tq, tetan1, tetan2

      nr=n_elements(rr)
      add_arrval, tq[nr-1], q
      add_arrval, tetan[nr-1], etan
      add_arrval, tetan1[nr-1], etan1
      add_arrval, tetan2[nr-1], etan2

  ENDFOR 

  outfile = '/sdss5/data0/lensout/ellipmodel_vse.fit'
  ss=create_struct('ellip', ell, $
                   'q', q, $
                   'etan', etan, $
                   'etan1', etan1, $
                   'etan2', etan2)

  mwrfits, ss,outfile, /create

END 
