PRO make_hudson_scat, ndup, scat=scat

  ;; output good sources to an ascii catalog with the required tags
  ;; Also, make a random catalog with the same density, cutting out
  ;; the unmasked regions: hudson needs this

  stripes = [9,10,11,12,13,14,15]
  stripestr = stripearr2string(stripes)

  sdss_shapecorr_dir = sdssidl_config('shapecorr_dir')
  outdir = sdss_shapecorr_dir+'combined/'
  srcfile = outdir + 'stripe'+stripestr+'_hudson.dat'
  srcfits = outdir + 'stripe'+stripestr+'_hudson.fit'
  rsrcfile = outdir + 'stripe'+stripestr+'_rand_hudson.dat'
  pngfile = '~/plots/stripe'+stripestr+'_hudson.png'


  ;; Area
;  area = 1365.512015            ;normal
  area = 1366.025729            ;/nomissf
;  area = 1199.095569            ;combined
  nomissf = 0
  combined = 0

  print
  print,'Stripes: ',stripes
  print,'Area = ',area

  clrs = [1,2,3]
  get_nz_photoz, stripes, clrs, nz, /hirata, /hudson
  primary_bound_multi, stripes, pb, overall_bound

  IF n_elements(scat) EQ 0 THEN BEGIN 
      get_scat, stripes, clrs, scat, /hirata, /hudson
      scat = temporary(scat[nz.useind])

      apply_pixel_mask, scat.clambda, scat.ceta, masked, unmasked, $
                        combined=combined, nomissf=nomissf

      scat = temporary(scat[unmasked])

  ENDIF 

  title = 'Stripes  '+repstr(stripestr, '_',' ')
  plot, scat.clambda, scat.ceta, psym=3, /iso, $
        xtitle=!csym.lambda, ytitle=!csym.eta, $
        title=title

  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; generate random points
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; density
  nkeep = n_elements(scat)
  density = nkeep/area

  lammin = overall_bound.lammin
  lammax = overall_bound.lammax
  etamin = overall_bound.etamin
  etamax = overall_bound.etamax

  strad2Deg = 360.0*360.0/(4.0*!pi*!pi)
  overall_area = sin(!d2r*lammax) - sin(!d2r*lammin);
  overall_area = overall_area*strad2Deg*(!d2r*(etamax - etamin));

  nrand = ndup*round(overall_area*density)
  
  ;; generate uniformly in sin(lambda) and eta
  print,'Generating '+ntostr(nrand)+' Randoms at density '+ntostr(density)+' per square degree'
  sinflam = sin( lammin*!d2r )
  sinllam = sin( lammax*!d2r )

  slambda = arrscl( randomu(seed, nrand, /double), $
                    sinflam, sinllam, $
                    arrmin=0., arrmax=1.)
  
  lambda = asin(temporary(slambda))*!r2d
  
  eta = arrscl(randomu(seed, nrand, /double), $
               etamin, etamax, $
               arrmin=0., arrmax=1.)

  apply_pixel_mask, lambda, eta, rmasked, runmasked, $
                    combined=combined, nomissf=nomissf


  nrand = n_elements(rmasked)
  lambda = lambda[rmasked]
  eta = eta[rmasked]

  print,'Of these, '+ntostr(nrand)+' kept outside mask'

  ;; now associate each random with a real galaxy.  Need to make sure the
  ;; spectro galaxies are correctly represented here
  index = long( arrscl( randomu(seed, nrand), 0, nkeep, arrmin=0, arrmax=1) )

  oplot, lambda, eta, color=!blue, psym=3
  print,'Writing png file: ',pngfile
  
  legend,['Real', 'Random'], psym=[8,8], color=[!p.color, !blue], $
         /right,box=0, /clear

  write_png, pngfile, tvrd(/true)

; GOTO,jump
  ;; write source catalog
  print,'Writing sources to file: ',srcfile
  openw, lun, srcfile, /get_lun

  format = $
    '('+$
    '2(D17.11,:,1X),' + $       ; lambda, eta
    '18(F0,:,1X)'+$
    ')'

  FOR i=0L, nkeep-1 DO BEGIN 

      printf, lun, $
              scat[i].clambda, scat[i].ceta, $
              scat[i].e1_recorr, scat[i].e2_recorr, $
              scat[i].e1e1err, scat[i].e1e2err, scat[i].e2e2err, $
              scat[i].photoz_z, scat[i].photoz_zerr, $
              scat[i].photoz_abscounts[0], $
              scat[i].photoz_abscounts[1], $
              scat[i].photoz_abscounts[2], $
              scat[i].photoz_abscounts[3], $
              scat[i].photoz_abscounts[4], $
              scat[i].upetro, $
              scat[i].gpetro, $
              scat[i].rpetro, $
              scat[i].ipetro, $
              scat[i].zpetro, $
              scat[i].rpetrorad, $
              format=format
              
  ENDFOR 

  free_lun, lun


;return

  print,'Writing Fits file',srcfits
  mwrfits2, scat, srcfits, /create, /destroy
jump:

  ;; write random catalog
  print,'Writing random sources to file: ',rsrcfile
  openw, lun, rsrcfile, /get_lun

  format = $
    '('+$
    '2(D17.11,:,1X),' + $       ; lambda, eta
    '(I0,:,1X)'+$
    ')'

  FOR ii=0L, nrand-1 DO BEGIN 
      
      printf, lun, lambda[ii], eta[ii], index[ii], format=format

  ENDFOR 

  free_lun, lun


END 
