PRO sigmacrit_zcuts, zL, nsig, nlens, doplot=doplot, usetrue=usetrue

  ;; test various ways of getting the mean shear: when we have a photoz
  ;; distribution, using the overall redshift distribution, making cuts
  ;; on the photoz.

  ;;
  ;; Assume an overall source redshift distribution P(zs). Then draw from
  ;; this and make some gaussians P(zs)_i.  
  ;; The mean zs is drawn from the redshift distribution P(zs).  These P_i are 
  ;; the redshift distributions of the individual source galaxies.  Find the 
  ;; stacked P(zs) and the fraction f of sources in front of this lens. Then
  ;; compute the mean density contrast.

  ;; Then look at adding up all gaussians over all lenses and look at the
  ;; overall correction factor, does this give the same result?

  ;; do not add shot noise

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: sigmacrit_zcuts, zL, nsig, nlens, doplot=doplot, usetrue=usetrue'
      return
  ENDIF 

  tt = systime(1)

  ;; where to put the results
  outdir = '/net/cheops1/data0/esheldon/test_sigmacriterr/sim/'
  outfile = outdir + 'test_sigmacrit_'+ntostr(rnd(zL,3),5)+'_N1.fit'

  i=1
  WHILE fexist(outfile) DO BEGIN
      ;;outfile = newname(outfile)
      istr = ntostr(i)
      outfile = outdir + 'test_sigmacrit_'+ntostr(rnd(zL,3),5)+'_N'+istr+'.fit'
      i=i+1
  ENDWHILE 
  print,'Output file: ',outfile

  IF n_elements(nzstruct) EQ 0 THEN BEGIN 
      zfile = sdssidl_config('shapecorr_dir')+'sigmacrit/nzstruct_stripe10.fit'
      nzstruct = mrdfits(zfile,1,/silent)
  ENDIF 

  zs = nzstruct.z
  deltazs = zs[1]-zs[0]
  pofzs = nzstruct.rnz/total(nzstruct.rnz*deltazs)
  
  ;; mean density contrast
  meandenscont = 10.0
;  denscont = meandenscont + randomn(seed, nlens)

  ;; use same for all
  denscont = meandenscont + fltarr(nlens)

  ;; min denscont
  mindenscont = 4.0
  w=where(denscont GT mindenscont, nlens)
  denscont = denscont[w]

  ;; For now all lenses at same redshift

  zL = replicate(zL, nlens)
  print,'reading in SDSS sigma crit'
  stripe = 10
  clr = 2
  msigma_crit = 1./sdss_sigma_crit(stripe, clr, zL, $
                                   /use_lambda, /silent)

  ;; source width in redshift
  sigzs = 0.03

  ;; mean number of sources behind lens. Nsource will be drawn from
  ;; a poisson distribution with this mean
  mean_nsource = 10

  ;; number of z values in each P(zs)_i
  npzs = 500L

  ;; P(z) for each lens
  pofzi = fltarr(npzs)
  zsi = arrscl( findgen(npzs), min(zs), max(zs) )

  ;; used distribution
  pofzitot = pofzi

  ;; the same,but this time the true distributions
  pofzitrue = pofzi
  pofzitottrue = pofzi

  ;; different ways:
  ;; 1) measure denscont for each lens, correcting for contamination as you go,
  ;;    then average over the lenses
  ;; 2) use mean density contrast msigma_crit, just summing over all lenses
  mdenscont = fltarr(Nlens)
  frac = fltarr(Nlens)
  fractrue = frac               ;the true frac

  mdenscontsum = mdenscont
  nsourcesum  = lonarr(Nlens)
  nsourcei = lonarr(Nlens)
  nstotal = 0L

  indices = lindgen(nlens)

  delvarx, zstrue, zsest

  ;; generate the cumulative prob function for redshifts
  print,'Generating cumulative p(z)'
  genrand, pofzs, zs, 10, tzs, cumul=cumulpofzs

  ;; generate values from a gaussian with width 1 and mean 0 
  ;; from which we will draw values for our photoz gaussians, scaled
  ;; properly.  Generate 2*Nlens*mean_nsource
  print,"Generating source z's from normal distribution"
  ngen = 2L*Nlens*mean_nsource
  gaussrand = randomu(seed, ngen, /normal)

  ;; We draw from gaussrand in order. This index tells us where
  ;; we are in the array
  gausind = 0L

  ;; arrays to store values from individual sources
  dval = -9999.
  zstrue = replicate(dval, ngen)
  zsest  = zstrue
  sigzs = replicate(sigzs, ngen)

  FOR Li=0L, Nlens-1 DO BEGIN 

      ;; how many sources?
      nsource = randomn(seed, poisson=mean_nsource)

      ;;print,'Nsource['+ntostr(Li)+'] = ',nsource
      IF nsource GT 0 THEN BEGIN 
          
          nstotal = nstotal + nsource
          nsourcesum[Li] = nsource

          ;; generate true source redshifts
          genrand, cumulpofzs, zs, nsource, zsource, /cumul

          ;; True shear, defined as denscont*Dls/Ds
          DLs = angdist_lambda(zsource, zL[Li], $
                               h=1.0, omegamat=0.3)
          Ds  = angdist_lambda(zsource, $
                               h=1.0, omegamat=0.3)
          DL = angdist_lambda(zL[Li], $
                               h=1.0, omegamat=0.3)
          D = Dls*DL/Ds/1000.   ;Gpc
          sigmac = 1.663e3/D

          ;;shear = denscont[Li]*DLs/Ds
          shear = denscont[Li]/sigmac

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Thus, if we chose an object that is actually in front
          ;; we average in zero
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          w=where(zsource LT zL[Li],  nw)
          IF nw NE 0 THEN shear[w] = 0.0

          ;; now density contrast from mean sigmacrit
          ;; will later use to get mean/sdev
          mdenscontsum[Li] = total(msigma_crit[Li]*shear)

          ;; Now do this with photoz's

          Ngood = 0L
          pofzi[*] = 0.0
          FOR si=0L, nsource-1 DO BEGIN 
                  
              ;; draw from the normal distribution and scale appropriately
              ;; for mean and width of this photoz gaussian
              zstmp = gaussrand[gausind]*sigzs[gausind] + zsource[si]

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Redshift Cut:
              ;; Now, is our _estimate_ of the zs Nsig*sig behind
              ;; the lens redshift?
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              IF (zstmp - zL[Li]) GT nsig*sigzs[gausind] THEN BEGIN 

                  ;;add_arrval, zstmp, zsest
                  ;;add_arrval, zsource[si], zstrue

                  zsest[gausind] = zstmp
                  zstrue[gausind] = zsource[si]

                  Ngood = Ngood + 1L

                  ;; now generate new p(zs) from this
                  ;; to get our estimated distribution
                  vals = (zsi - zstmp)^2/2./sigzs[gausind]^2 < 9.0
                  pzs = exp(-vals)/sqrt(2.*!pi)/sigzs[gausind]

                  ;; The true distribution
                  vals = (zsi - zsource[si])^2/2./sigzs[gausind]^2 < 9.0
                  pzstrue = exp(-vals)/sqrt(2.*!pi)/sigzs[gausind]

                  ;; add to total estimated distribution around this lens

                  pofzi[*] = pofzi[*] + pzs[*]
                  pofzitot[*] = pofzitot[*] + pzs[*]

                  pofzitrue[*] = pofzitrue[*] + pzstrue[*]
                  pofzitottrue[*] = pofzitottrue[*] + pzstrue[*]

                  ;; estimate of dls/Ds from estimated redshift
                  Dlsi = angdist_lambda(zstmp, zL[Li], $
                                       h=1.0, omegamat=0.3)
                  Dsi  = angdist_lambda(zstmp, $
                                       h=1.0, omegamat=0.3)
                  
                  Di = DLsi*DL/Dsi/1000.
                  sigmaci = 1.663e3/Di
                  ;; now our estimate of denscont for this source
                  ;;mdenscont[Li] = mdenscont[Li] + shear[si]*Dsi/DLsi
                  mdenscont[Li] = mdenscont[Li] + shear[si]*sigmaci
                  nsourcei[Li] = nsourcei[Li] + 1

                  ;;print,'Added: ',shear[si]*sigmaci
              ENDIF 

              ;; increment the index
              gausind = gausind + 1L
          ENDFOR

          IF Ngood GT 0 THEN BEGIN 

              ;; mean denscont for this lens
              mdenscont[Li] = mdenscont[Li]/Ngood
          
              IF keyword_set(doplot) THEN BEGIN 
                  plot,zsi, pofzi
                  key=get_kbrd(1)
              ENDIF 
              
              ws = where(zsi GT zL[Li], nws)
              IF nws EQ 0 THEN BEGIN
                  print,'no zs gt zL'
                  mdenscont[Li] = 0.0 
                  indices[Li] = -1
              ENDIF ELSE BEGIN 
                  fracgood = total(pofzi[ws],/double)/total(pofzi,/double)
                  fracgoodtrue = total(pofzitrue[ws],/double)/total(pofzitrue,/double)
                  fractrue[Li] = fracgoodtrue
                  IF fracgood EQ 0 THEN BEGIN
                      print,'bad frac'
                      mdenscont[Li] = 0.0
                      frac[Li] = -1.
                      indices[Li] = -1
                  ENDIF ELSE BEGIN 
                      ;;print,'frac = ',fracgood ;now printing below
                      frac[Li] = fracgood
                      mdenscont[Li] = mdenscont[Li]/fracgood
                  ENDELSE 
              ENDELSE 
          ENDIF ELSE BEGIN 
              print,'Found no sources far enough behind: Nsource['+ntostr(Li)+'] = ',nsource
              indices[Li] = -1
          ENDELSE 

          IF (Li MOD 100) EQ 0 THEN BEGIN
              print,'Nsource['+ntostr(Li)+'] = ',nsource
              print,'frac = ',fracgood
          ENDIF 
      ENDIF ELSE BEGIN ;; nsource > 0
          print,'No sources' 
      ENDELSE 
  ENDFOR 

  ;; usable lenses
  w=where(indices NE -1, nw)
  IF nw NE 0 THEN indices = indices[w]
  
  ;; sources
  ws=where(zsest NE dval, nw)
  zsest = zsest[ws]
  zstrue = zstrue[ws]

  ;; normalize the distributions: we added total(nsourcei) normalized
  ;; gaussians, so divide by total(nsourcei)
  nsourceitot = total(nsourcei[w])
  ;;deltazsi = zsi[1]-zsi[0]
  pofzitot = pofzitot/nsourceitot

  pofzitottrue = pofzitottrue/nsourceitot

  plot, zsi, pofzitot

  tmp=create_struct('nlens', nlens, $
                    'zs', zs, $
                    'pofzs', pofzs, $
                    'meandenscont', meandenscont, $
                    'denscont', denscont[w], $
                    'zl', zL[w], $
                    'msigma_crit', msigma_crit[w], $
                    'sigzs', sigzs, $
                    'mean_nsource', mean_nsource, $
                    'npzs', npzs, $
                    'nsig', nsig, $
                    'zsi', zsi, $
                    'pofzitot', pofzitot, $ ; reconstructed pofz after cuts
                    'pofzitottrue', pofzitottrue, $; true distribution
                    'fracgood', frac[w], $ ; fraction of objects expected 
                    $                   ; to be in front of the lens
                    'fracgoodtrue',fractrue[w], $ ;true fraction
                    'mdenscont', mdenscont[w], $ ;zcuts
                    'nsourcei', nsourcei[w], $   ;# in each after zcuts
                    'nsourceitot', nsourceitot, $
                    'mdenscontsum', mdenscontsum[w], $ ; using overall dist
                    'nsourcesum', nsourcesum[w],$ ;number sources used here
                    'nstotal', nstotal, $; total sources used in sums
                    'indices', indices, $; good objects
                    'zstrue', zstrue, $  ; input redshift for sources
                    'zsest', zsest)      ; estimated (drawn from dist)
  print
  print,'Output file: ',outfile 
  mwrfits, tmp, outfile, /create

  ptime,systime(1)-tt

END 
