
PRO wthetarand_lensgal_lcut_logbin, nbin, wclr, stripe, allcat, $
                                    rmin, rmax, zrand_in, ndup, $
                                    use_lambda=use_lambda, $
                                    step=step, addstr=addstr, $
                                    outdir=outdir, $
                                    wgood=wgood, $
                                    usecat=usecat, $
                                    datfile=datfile, sumfile=sumfile, zfile=zfile, $
                                    lensumfile=lensumfile, $
                                    maxe=maxe,fraclstar=fraclstar,numfile=numfile,$
                                    nbeta=nbeta, tsgals=tsgals

  IF n_params() LT 8 THEN BEGIN
      print,'-Syntax: wthetarand_lensgal_lumweight, nbin, wclr, stripe, allcat, rmin, rmax, binsize, zrand_in, ndup,'
      print,'  use_lambda=use_lambda,'
      print,'  step=step, addstr=addstr, '
      print,'  outdir=outdir, '
      print,'  wgood=wgood, '
      print,'  usecat=usecat, '
      print,'  datfile=datfile, sumfile=sumfile, zfile=zfile, '
      print,'  lensumfile=lensumfile, '
      print,'  maxe=maxe , tsgals=tsgals'
      print,'-> lenscat may have ra,dec or lambda, eta but allcat must contain lambda,eta'
      print,'-> allcat must be sorted by lambda'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Common blocks
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  COMMON wthetarand_block, seed, etarange, etarange_spectra, $
    mlambda, meta, FIRSTLAMBDA, LASTLAMBDA, Mlimit

  Mstar = [-18.34,-20.04,-20.83,-21.26,-21.55]
  Msun = [6.39,5.07,4.62,4.52,4.48]
  Lstar = 10.0^((Mstar-Msun)/(-2.5))
  Lstar = Lstar/1.e10           ;in units of 10^10 solar lum

  IF n_elements(fraclstar) EQ 0 THEN fraclstar = 0.1
  lcut = lstar[wclr]*fraclstar

;  Mlimit = Mstar[2]-2.5*alog10(fraclstar)

  print
  print,'Going to '+ntostr(fraclstar)+' of L*'
  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  colors = ['u','g','r','i','z']

  time = systime(1)

  IF NOT keyword_set(use_lambda) THEN omegamat=1.0 ELSE omegamat=0.3

  ;; Size by which to group things  300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = 'wthetarandlumw'
  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss4/data1/esheldon/TMP/'

  IF wclr LT 0 OR wclr GT 4 THEN message,'wclr must be in [0,4]!'

  IF n_elements(numfile) EQ 0 THEN fend=colors[wclr]+'w_N1' $
  ELSE fend=colors[wclr]+'w_N'+ntostr(long(numfile))

  stripestr = '_stripe'+ntostr(stripe)

  d2r = !dpi/180.0d0            ;Change degrees to radians.
  r2d = 180./!dpi               ;radians to degrees

  ;; Now nbin is input
  ;;nbin = long( (rmax - rmin)/binsize ) + 1 
  pbin = 5 < (nbin-1)
  
  print
  print,'-------------------------------------------------------------'
  IF NOT keyword_set(use_lambda) THEN print,'Using lambda = 0' $
  ELSE print,'Using lambda = 0.3'
  print,'Rmin: ',rmin
  print,'Rmax: ',rmax
  print,'Step: ',step
  print,'Using ',ntostr(nbin),' bins between ',ntostr(rmin), $
        ' and ',ntostr(rmax),' kpc'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; duplicate the input z distribution
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Duplicating: ',ntostr(ndup)

  sdssidl_setup,/silent
  nrand_in = n_elements(zrand_in)
  nrand = nrand_in*ndup
  ninit = nrand

  FOR dd=0L, ndup-1 DO BEGIN 

      add_arrval, zrand_in, zrand
      add_arrval, lindgen(nrand_in), zindex

  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read in etarange info
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  read_etarange, stripe, etarange
  read_etarange, stripe, etarange_spectra, /spectra

  mlambda = [etarange.lambda*d2r, etarange.lambda*d2r]
  meta = [etarange.mineta*d2r, etarange.maxeta*d2r]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read Mask info for this stripe
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  read_stripe_mask, stripe, mask, tsgals=tsgals

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Beta values
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;IF n_elements(nbeta) EQ 0 THEN nbeta = 41 ;this does beta in [0,2], steps of 0.05
  IF n_elements(nbeta) EQ 0 THEN nbeta = 1 
  betamin = 0.0
  betastep = 0.05
  beta = fltarr(nbeta)
  FOR i=0L, nbeta-1 DO beta[i] = betamin+i*betastep
  print
  print,'n_beta = ',nbeta
  print,'max beta = ',max(beta)
  print,'min beta = ',betamin
  print
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; declare some arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; use nbin+1 since histogram often has one extra small bin
  ubin = nbin+1

  arrval = fltarr(ubin)
  arrval2 = fltarr(ubin, nbeta)

  wthetastr=wthetaLumwStruct(arrval, beta)
  wthetasumstr=wthetaLumwSumstruct(arrval, beta)

  wthetalensum=wthetaLumwLensumStruct(arrval, beta)
  wthetalensum = create_struct('z', 0., wthetalensum)
  wthetalensum=replicate(wthetalensum, nrand)
  wthetalensum.zindex = zindex
  wthetalensum.z = zrand

  wsum       = arrval
  rsum       = arrval
  lsum       = arrval
  lwsum      = arrval
  npsum      = arrval2
  npair      = arrval

  rmax_act   = arrval
  rmax_act_tmp = arrval
  rmin_act = replicate(1.e7, ubin)
  rmin_act_tmp = arrval

  ;; temporary
  tmpwsum = arrval
  tmprsum = arrval
  tmplsum = arrval
  tmplwsum = arrval
  tmpnpsum = arrval2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up output files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  prename = outdir+addstr+stripestr
  
  datfile = prename + '_'+fend+'.fit'
  sumfile = prename + '_sum_'+fend+'.fit'
  zfile   = prename + '_z_'+fend+'.fit'
  psfile  = prename + '_points_'+fend+'.ps'
  lensumfile = prename + '_lumlensum_'+fend+'.fit'

  WHILE fexist(datfile) DO BEGIN
      datfile = newname(datfile)
      sumfile = newname(sumfile)
      zfile = newname(zfile)
      psfile = newname(psfile)
      lensumfile = newname(lensumfile)
  ENDWHILE 
  print
  print,'Dat file: ',datfile
  print,'Sum file: ',sumfile
  print,'Lens sum file: ',lensumfile
  print,'Points file: ',psfile
  print,'Redshift file: ',zfile
  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; source cat
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; they are already sorted
  nall = n_elements(allcat)
  wall = lindgen(nall)

  FIRSTLAMBDA = allcat[0].lambda
  LASTLAMBDA  = allcat[nall-1].lambda 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Lens cat: Set up sigma crit and Dlens. 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wlens = lindgen(nrand)

  h=1.
  print,'Using h = ',h
  print,'Calculating 1/sigma_crit'
  print

  clr=2
  sigcritinv = sdss_sigma_crit(stripe, clr, zrand, $
                               wgood=wgood, use_lambda=use_lambda)
  sigmacrit = 1./sigcritinv

  ;;Convert from Mpc to kpc
  DL = angdist_lambda( zrand, h=h, omegamat=omegamat)*1000. 
  angmax = rmax/DL*180./!pi     ;angles in degrees

  ;; don't want angle to be larger than 1/4 of a stripe (until we use 
  ;; multiple stripes)
  max_allowed_angle = 2.5/4.0
  wgood2 = where(angmax[wgood] LE max_allowed_angle)
  wgood = wgood[wgood2]
  wlens = wlens[wgood]

  ;; don't want angle to be larger than 1/4 of a stripe (until we use 
  ;; multiple stripes)
  max_allowed_angle = 2.5/4.0
  wgood2 = where(angmax[wgood] LE max_allowed_angle)
  wgood = wgood[wgood2]
  wlens = wlens[wgood]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; If it is a southern stripe just rotate it so it for the
  ;; initial creation of lambda's. Then rotate back to use
  ;; the mask
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF stripe GT 45 THEN BEGIN
      print
      print,'This is a Southern Stripe.'
      print
      issouth = 1
      rotate_lambda, FIRSTLAMBDA
      rotate_lambda, LASTLAMBDA
  ENDIF ELSE issouth=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; create lambda's and eta's
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ngood = 0L
  llambda = dblarr(nrand)
  leta = llambda

  bad = wlens
  nwlens = n_elements(wlens)
  nbad = nwlens
  print
  print,"Generating lambda's and eta's: ",ntostr(nbad)
  WHILE ngood LT nwlens DO BEGIN 

      IF ngood EQ 0 THEN print,'Generating ',ntostr(nbad) $
      ELSE print,'Re-Generating ',ntostr(nbad)
      tmplam = randomu(seed, nbad)

      ;; generate lambdas
      wthetarand_generate_lambdas, nbad, angmax[bad], tmplam

      ;; rotate the lambda'a if southern stripe
      IF issouth THEN rotate_lambda, tmplam

      ;; now apply mask
      apply_mask, mask, tmplam, tbad, tgood

      IF tgood[0] NE -1 THEN BEGIN 
          ;; we found some good ones
          good = bad[tgood]
          IF tbad[0] NE -1 THEN BEGIN 
              bad = bad[tbad]
              nbad = n_elements(bad)
          ENDIF ELSE nbad=0L
          llambda[good] = tmplam[tgood]

          ntgood = n_elements(tgood)
          ngood = ngood + ntgood

          ;; make generate etas
          wthetarand_generate_etas, llambda[good], ntgood, angmax[good], tmpeta
          leta[good] = tmpeta
          
      ENDIF 

  ENDWHILE 

  w=where(leta[wlens] EQ 1.e10,nw)
  IF nw NE 0 THEN remove, w, wlens

  ;; sort the subscript array
  IF issouth THEN BEGIN 
      ;; rotate and sort
      rotate_lambda, llambda
      s = sort(llambda[wlens])
      wlens = wlens[s]
      ;; rotate back
      rotate_lambda, llambda

      rotate_lambda, FIRSTLAMBDA
      rotate_lambda, LASTLAMBDA
  ENDIF ELSE BEGIN 
      s = sort(llambda[wlens])
      wlens = wlens[s]
  ENDELSE 

  begplot,name=psfile,/landscape,/color
  setupplot

  plot, [etarange_spectra.lambda, etarange_spectra.lambda], $
    [etarange_spectra.maxeta,etarange_spectra.mineta],psym=4
  plot, llambda[wlens], leta[wlens], psym=3, /ynozero, ystyle=1
  oplot, etarange_spectra.lambda,etarange_spectra.maxeta,psym=4,color=!red
  oplot, etarange_spectra.lambda,etarange_spectra.mineta,psym=4,color=!red
  ep
  pslandfix, psfile
  setupplot

  nlens = n_elements(wlens)
  print,'Threw out ',ntostr(ninit - nlens),' lenses because too deep or close'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; print some stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'FIRSTLAMBDA = ',FIRSTLAMBDA
  print,'LASTLAMBDA = ',LASTLAMBDA
  print,'Using '+ntostr(nlens)+'/'+ntostr(ninit)+' lenses'
  print
  IF nlens LT 100 THEN sym = 1 ELSE sym = 3
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; count neighbors
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nstepOld = nlens/step
  nstep=nstepOld
  left = nlens MOD step
  
  ;; To account for leftover stuff
  IF nstepOld EQ 0 THEN BEGIN
      nstepOld = -10
      step = left
      nstep = 1
  ENDIF ELSE BEGIN
      IF left NE 0 THEN nstep = nstepOld + 1
  ENDELSE 

  indices = lindgen(nlens)
  FOR group = 0L, nstep-1 DO BEGIN
      IF group EQ nstepOld THEN BEGIN
          ind = indices[ group*step: group*step+left-1  ]
          ii = wlens[ind]
          step = left
      ENDIF ELSE BEGIN
          ind = indices[ group*step : (group+1)*step -1 ]
          ii = wlens[ind]
      ENDELSE 

      ;; Choose sources around this lens group
      ;; they are sorted by lambda, but it could go over lambda=180      
      angmax2 = max(angmax[ii])
      
      maxii = wlens[ max(ind) ] & minii = wlens[ min(ind) ]

      maxlam1 = llambda[maxii]+angmax2
      minlam1 = llambda[minii]-angmax2

      wtheta_get_neighbors, allcat.lambda, minlam1, maxlam1, wneigh, nneigh, $
        issouth=issouth

      IF nneigh NE 0 THEN BEGIN 
          IF lwsum[pbin] NE 0 THEN $
            print,'Mean Lum '+ntostr(pbin)+' : '+ntostr(lsum[pbin]/lwsum[pbin])

          print,'Lens Group = ',ntostr(group+1)+'/'+ntostr(nstep)
          FOR gi=0L, step-1 DO BEGIN ;Loop over lenses in group
              IF (gi MOD 10) EQ 0 THEN  print,'.',format='(a,$)'
              index = ii[gi]
              ceneta = leta[index] ;lens center
              cenlam  = llambda[index]

              sig_inv = sigcritinv[index]
              lensw = sig_inv^2/DL[index]^2*1.e19 ;1.e19 for order 1 numbers

              angmax_i = angmax[index]

              ;; choose sources around this lens
              maxlam2 = cenlam+angmax_i
              minlam2 = cenlam-angmax_i

              wtheta_get_neighbors, allcat[wneigh].lambda, minlam2, maxlam2, $
                twneigh, nneigh2, issouth=issouth

              ;; If there are any sources left, measure the shear
              IF nneigh2 GT 0 THEN BEGIN 

                  wneigh2 = wneigh[twneigh]

                  ;; find luminosity, assuming they are at same redshift
                  ;; as the lens galaxy.
                  ;; returns lumsolar in units of 10^10
                  wtheta_absmag, zrand[index], wclr, $
                    allcat[wneigh2].petrocounts[wclr], $
                    allcat[wneigh2].gr, absmag, lumsolar

                  wneigh3 = where(lumsolar GE lcut, nneigh2)
                  IF nneigh2 EQ 0 THEN GOTO,jump
                  wneigh2 = wneigh2[wneigh3]
                  lumsolar=lumsolar[wneigh3]

                  ;; luminosity weight will be (lum/lstar)^beta
                  ;; must both be in units of 10^10 solar
                  lumw=lumsolar/lstar[wclr]

                  ;; calculate angular distance (dec->lambda in distance calc)
                  tscateta = allcat[wneigh2].eta
                  etadiff = (ceneta-tscateta)*d2r
                  cosetadiff = cos(etadiff)
                  
                  tcenlam = cenlam*d2r ;lens center in radians
                  sincenlam = sin(tcenlam)
                  coscenlam = cos(tcenlam)
                  
                  tscatlam = allcat[wneigh2].lambda*d2r
                  sinscatlam = sin(tscatlam)
                  cosscatlam = cos(tscatlam)
;                  tscatlam = tscatlam*r2d
                  ;; Find distance in kpc
                  args = sincenlam*sinscatlam + coscenlam*cosscatlam*cosetadiff
                  warg = where(args LT -1., narg)
                  IF narg NE 0 THEN args[warg] = -1.
                  warg = where(args GT 1., narg)
                  IF narg NE 0 THEN args[warg] = 1.

                  R=acos( args )*DL[index]
                  IF nneigh2 EQ 1 THEN R=[R]

;                  hist=histogram(R, binsize=binsize, min=rmin, $
;                                 max=rmax,rever=rev_ind)

                  ;; binsizes increase exponentially
                  logbin, R, rmin, rmax, nbin, hist, rev_ind

                  numbin=n_elements(hist)
                  IF numbin GT ubin THEN BEGIN
                      message,'What!!!   #hist = '+ntostr(numbin)+$
                        '   nbin = '+ntostr(ubin)
                  ENDIF 
                  whist = where(hist NE 0, nhist)
                  ng=0
                  ;; Check if there are any in this annulus rmin-rmax
                  IF nhist NE 0 THEN BEGIN 
                          
                      ;;nbin is predefined. Last bin will hold few
                      npair[*] = 0L
                      rmax_act_tmp[*] = 0.
                      rmin_act_tmp[*] = 0.
                      FOR i=0L, nhist-1 DO BEGIN 
                          binnum = whist[i]
                          w=rev_ind( rev_ind(binnum):rev_ind(binnum+1)-1 )
                          
                          rmax_act_tmp[binnum] = max(R[w], min=minrr)
                          rmin_act_tmp[binnum] = minrr

                          rmax_act[binnum] = max( [rmax_act[binnum], rmax_act_tmp[binnum]] )
                          rmin_act[binnum] = min( [rmin_act[binnum], rmin_act_tmp[binnum]] )
                          npair[binnum] = n_elements(w)

                          ;; loop over the different beta values
                          tmprsum[binnum]=total(R[w]*lumw[w]^beta[0])*lensw
                          FOR ib=0L, nbeta-1 DO BEGIN 
                              ;; luminosity weighted counts
                              tmpnpsum[binnum,ib]=total(lumw[w]^beta[ib])*lensw
                          ENDFOR 
                          tmplsum[binnum] = total(lumsolar[w])*lensw
                          tmplwsum[binnum] = lensw
                          tmpwsum[binnum] = lensw

                      ENDFOR 

                      rsum[whist] = rsum[whist] + tmprsum[whist]
                      npsum[whist,*] = npsum[whist,*] + tmpnpsum[whist,*]
                      wsum[whist] = wsum[whist] + tmpwsum[whist]
                      lsum[whist] = lsum[whist] + tmplsum[whist]
                      lwsum[whist] = lwsum[whist] + tmplwsum[whist]

                      ;; copy in individual lens stuff
                      ;; Note totpairs is unweighted
                      wthetalensum[index].totpairs = total(npair[whist])
                      wthetalensum[index].npsum[whist,*] = tmpnpsum[whist,*]
                      wthetalensum[index].rmax_act = R2
                      wthetalensum[index].rmin_act = R1
                      wthetalensum[index].rsum[whist] = tmprsum[whist]
                      wthetalensum[index].lsum[whist] = tmplsum[whist]
                      wthetalensum[index].lwsum[whist] = tmplwsum[whist]
                      wthetalensum[index].wsum[whist] = tmpwsum[whist]

                      ;; reinitialize the arrays
                      tmprsum[*] = 0.
                      tmpnpsum[*] = 0
                      tmpwsum[*] = 0.
                      tmplsum[*] = 0.
                      tmplwsum[*] = 0.
                      npair[*] = 0.

                  ENDIF ELSE BEGIN 
                      print,'/',format='(a,$)'
                      ;indices[ ind[gi] ] = -1
                  ENDELSE 

                  R = 0
                  etadiff=0      ;Free memory
              ENDIF ELSE BEGIN
jump:
                  print,'|',format='(a,$)'
                  ;indices[ ind[gi] ] = -1
              ENDELSE 

          ENDFOR 
          print
          ;xrel = 0 & yrel = 0
      ENDIF ELSE BEGIN 
          print,'Two-'
          ;indices[ ind ] = -1
      ENDELSE 
      metadiff = 0               ;Free memory
      mdist = 0
  ENDFOR 
  ;; Remove unused lenses
  wbad = where(indices EQ -1, nwbad)
  IF nwbad NE 0 THEN remove, wbad, wlens
  lensused = n_elements(wlens)
  
  wgood = wlens
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Find averages from sums
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Finally used ',ntostr(lensused),'/',ntostr(nlens),' Lenses'

  IF npsum[ubin-1,0] LT npsum[ubin-2,0] THEN BEGIN 
      ;; This means last bin was incomplete
      ubin = ubin-1
  ENDIF 

  wthetastr.rmax_act = rmax_act
  wthetastr.rmin_act = rmin_act
  FOR i=0L, ubin-1 DO BEGIN 

      ;; calculate area, density of background galaxies
;      R1 = rmin + i*binsize
;      R2 = rmin + (i+1)*binsize
      R1 = wthetastr.rmin_act[i]
      R2 = wthetastr.rmax_act[i]
      wthetastr.area[i] = !pi*(R2^2 - R1^2)

      ;; shear struct for this set of lenses
      wthetastr.meanr[i] = rsum[i]/npsum[i,0]
      ;; does not include central galaxy
      wthetastr.meanlum[i] = lsum[i]/lwsum[i]

      FOR ib=0L, nbeta-1 DO BEGIN 
          ;; now an average per lens
          wthetastr.npair[i,ib] = npsum[i,ib]/wsum[i]

          ;; Now totals within each rmax_act[i]
          wthetastr.tnpair[i,ib] = total(npsum[0:i,ib])/total(wsum[0:i])

          ;; npair is average per lens, no need to divide by lensused
          wthetastr.density[i,ib]=wthetastr.npair[i,ib]/wthetastr.area[i] ;/lensused
      ENDFOR 
      ; sum struct for adding to other sets of lenses
      wthetasumstr.rsum[i] = rsum[i]
      wthetasumstr.npsum[i,*] = npsum[i,*]
      wthetasumstr.lsum[i] = lsum[i]
      wthetasumstr.lwsum[i] = lwsum[i]
      wthetasumstr.wsum[i] = wsum[i]

  ENDFOR 

  totpairs = total(wthetalensum.totpairs)

  wthetastr.nlenses = lensused
  wthetastr.totpairs = totpairs
  wthetastr.binsize = -1.
  wthetastr.rmin = rmin
  wthetastr.rmax = rmax
  wthetastr.h = h

  wthetasumstr.nlenses = lensused
  wthetasumstr.totpairs = totpairs
  wthetasumstr.binsize = -1.
  wthetasumstr.rmin = rmin
  wthetasumstr.rmax = rmax
  wthetasumstr.rmax_act = rmax_act
  wthetasumstr.rmin_act = rmin_act
  wthetasumstr.h = h

  wthetalensum.binsize = -1.
  wthetalensum.rmin = rmin
  wthetalensum.rmax = rmax
  wthetalensum.h = h

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Ouput the data
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; headers
  hdr = whdr(wthetastr)      
  sumhdr=hdr
  lhdr=hdr

  ;; Shear file for these lenses    
  mwrfits, wthetastr, datfile, hdr, /create

  ;; Sum file for these lenses
  mwrfits, wthetasumstr, sumfile, sumhdr, /create

  ;; File with structure for each lens used
  wthetalensum = temporary(wthetalensum[wlens])
  mwrfits2, wthetalensum, lensumfile, lhdr, /create, /destroy

  ;; File containing ra,dec,redshift for each used lens
  zs = create_struct('z', 0., $
                     'lambda', double(0.), $
                     'eta', double(0.) )
  zstruct = replicate(zs, lensused)
  zstruct.z = zrand[wlens]
  zstruct.lambda = llambda[wlens]
  zstruct.eta = leta[wlens]

  mwrfits, zstruct, zfile, /create

  step=oldstep

  ptime, systime(1)-time
  return
END 
