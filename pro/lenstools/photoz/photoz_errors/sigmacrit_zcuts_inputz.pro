PRO sigmacrit_zcuts_inputz_sigmacrit, str, zL, msigma_crit, z, pofz

  photoz_dist, str, z, pofz, good
  deltaz = z[1]-z[0]
  ;; normalize
  pofz = pofz/total(pofz*deltaz)

  IF display_exists() THEN BEGIN 
      plothist, str[good].photoz_z, bin=0.01,/norm
      oplot, z, pofz,color=!red
  ENDIF

  ;; normalize 
  norm = total(pofz,/double)

  ;; generate average sigmacrit
  nL = 1000L
  minzL = 0.01
  maxzL = 1.5
  zL = arrscl( findgen(nL), minzL, maxzL )

  msigma_crit = fltarr(nL)
return
  FOR i=0L, nL-1 DO BEGIN 

      w=where( (z - zL[i]) GT 0.001, nw)
      IF nw NE 0 THEN BEGIN 
          DLs = angdist_lambda(z[w], zL[i], $
                               h=1.0, omegamat=0.3)
          Ds  = angdist_lambda(z[w], $
                               h=1.0, omegamat=0.3)
          DL = angdist_lambda(zL[i], $
                              h=1.0, omegamat=0.3)
          D = Dls*DL/Ds/1000.   ;Gpc
          sc = 1.663e3/D

          ;; now the mean
          msigma_crit[i] = total(pofz[w]*sc,/double)

      ENDIF 

  ENDFOR 

  msigma_crit = msigma_crit/norm
  w=where(msigma_crit GT 0)
  zL = zL[w]
  msigma_crit = msigma_crit[w]


END 

PRO sigmacrit_zcuts_inputz, zLin, nsig, nlens, stripe, str, $
                            doplot=doplot, usetrue=usetrue, $
                            weight_type=weight_type

  ;; this one we input the photoz struct

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

  IF n_params() LT 5 THEN BEGIN 
      print,'-Syntax: sigmacrit_zcuts, zL, nsig, nlens, str, doplot=doplot, usetrue=usetrue'
      return
  ENDIF 

  IF n_elements(weight_type) EQ 0 THEN weight_type = -1

  tt = systime(1)
  setup_mystuff

  ;; norm of power law (sigcriterr/sigcrit)/sigzs vs. zs
  ;; fracerr = enorm*zs^epow * sigzs
  ;; sigcriterr = fracerr*sigcrit
  enorm = !sigerr.ap_norm*zLin[0]^!sigerr.bp_norm
  epow  = !sigerr.a_pow  + !sigerr.b_pow *zLin[0]

  ;; where to put the results
  CASE sdssidl_config('hostname') OF
      'fatso': outdir = '/home/bruno1/users/esheldon/pzcuts_sigmacrit/inputz/'
      'cheops1': outdir = '/net/cheops1/data0/esheldon/pzcuts_sigmacrit/inputz/'
      'cheops2': outdir = '/net/cheops1/data0/esheldon/pzcuts_sigmacrit/inputz/'
      'cheops3': outdir = '/net/cheops1/data0/esheldon/pzcuts_sigmacrit/inputz/'
      'cheops4': outdir = '/net/cheops1/data0/esheldon/pzcuts_sigmacrit/inputz/'
      'treebeard':outdir = '/home/esheldon/pzcuts_sigmacrit/inputz/'
      ELSE: stop
  ENDCASE 
  
  i=1
  outfile = outdir + ntostr(stripe)+'_sigmacrit_'+ntostr(rnd(zLin,3),5)+'_N1.fit'
  WHILE fexist(outfile) DO BEGIN
      ;;outfile = newname(outfile)
      istr = ntostr(i)
      outfile = outdir + 'test_sigmacrit_'+ntostr(rnd(zLin,3),5)+'_N'+istr+'.fit'
      i=i+1
  ENDWHILE 
  print,'Output file: ',outfile

  ;; make cuts on this distribution
  rmag = str.petrocounts[2] - str.reddening[2]
  zerrmin = 0.01
  goodz=where( rmag lt 21.0 and rmag gt 18.0 and $
              str.objc_type eq 3 and $
              str.photoz_z gt 0.02 and str.photoz_z lt 0.9 and $
              str.photoz_quality gt 0 and str.photoz_quality lt 12,ngoodz)

  zsource = str[goodz].photoz_z
  sigzs   = str[goodz].photoz_zerr > zerrmin

  ;; mean density contrast: for reference
  meandenscont = 10.0

  ;; use same for all
  denscont = meandenscont + fltarr(nlens)

  ;; min denscont
  mindenscont = 4.0
  w=where(denscont GT mindenscont, nlens)
  denscont = denscont[w]

  ;; mean sigmacrit for lens redshift, given
  ;; pofz
  ;; !!!!!!DONT USE THE SIGMACRIT, IT ISN'T RIGHT
  sigmacrit_zcuts_inputz_sigmacrit, str[goodz], tzL, tmsigma_crit, zs, pofzs
  zL = replicate(zLin, nlens)

  msigma_crit = interpol(tmsigma_crit, tzL, zL)

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
  nsources = lonarr(Nlens)
  wsum = fltarr(Nlens)
  nstotal = 0L

  indices = lindgen(nlens)

  delvarx, zstrue, zsest

  ;; generate values from a gaussian with width 1 and mean 0 
  ;; from which we will draw values for our photoz gaussians, scaled
  ;; properly.  Generate 2*Nlens*mean_nsource
  print,"Generating source z's from normal distribution"

  ngen = 2L*Nlens*mean_nsource
  gausrand = randomu(seed, ngen, /normal)

  ;; We draw from gausrand in order. This index tells us where
  ;; we are in the array
  gausind = 0L

  ;; arrays to store values from individual sources
  dval = -9999.
  zstrue = replicate(dval, ngen)
  zsest  = zstrue
  sigzsout = zstrue

  FOR Li=0L, Nlens-1 DO BEGIN 

      ;; how many sources?
      nsource = randomn(seed, poisson=mean_nsource)

      ;;print,'Nsource['+ntostr(Li)+'] = ',nsource
      IF nsource GT 0 THEN BEGIN 
          
          nstotal = nstotal + nsource
          ;nsourcesum[Li] = nsource

          ;; get random true source redshifts
          sind = long( arrscl(randomu(seed, nsource), 0, ngoodz, $
                              arrmin=0.0, arrmax=1.0) )

          ;; True shear, defined as denscont*Dls/Ds
          DLs = angdist_lambda(zsource[sind], zL[Li], $
                               h=1.0, omegamat=0.3)
          Ds  = angdist_lambda(zsource[sind], $
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

          w=where(zsource[sind] LT zL[Li],  nw)
          IF nw NE 0 THEN shear[w] = 0.0

          ;; now density contrast from mean sigmacrit
          ;; will later use to get mean/sdev
;          mdenscontsum[Li] = total(msigma_crit[Li]*shear)

          ;; Now do this with photoz's

          Ngood = 0L
          pofzi[*] = 0.0
          FOR si=0L, nsource-1 DO BEGIN 

              ind = sind[si]

              ;; draw from the normal distribution and scale appropriately
              ;; for mean and width of this photoz gaussian
              zstmp = gausrand[gausind]*sigzs[ind] + zsource[ind]

              ;; Now, is our _estimate_ of the zs Nsig*sig behind
              ;; the lens redshift?
              IF (zstmp - zL[Li]) GT nsig*sigzs[ind] THEN BEGIN 

                  zsest[gausind] = zstmp
                  zstrue[gausind] = zsource[ind]
                  sigzsout[gausind] = sigzs[ind]

                  Ngood = Ngood + 1L

                  ;; now generate new p(zs) from this
                  ;; to get our estimated distribution
                  vals = (zsi - zstmp)^2/2./sigzs[ind]^2 < 9.0
                  pzs = exp(-vals)/sqrt(2.*!pi)/sigzs[ind]

                  ;; The true distribution
                  vals = (zsi - zsource[ind])^2/2./sigzs[ind]^2 < 9.0
                  pzstrue = exp(-vals)/sqrt(2.*!pi)/sigzs[ind]

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
                  CASE weight_type OF
                      -1: weight = 1
                      1:  BEGIN
                          sfracerr = ( enorm*zstmp^epow )*sigzs[ind]
                          sigcriterr2 = (sigmaci*sfracerr)^2
                          e_err2 = str[ind].e1e1err^2 + str[ind].e2e2err^2
                          e2 = str[ind].e1^2 + str[ind].e2^2

                          err2 = 0.25*( sigmaci^2*e_err2 + sigcriterr2*e2 )

                          weight = 1./err2
                          
                      END 
                      ELSE: message,'WHAT!'
                  ENDCASE 

                  wsum[Li] = wsum[Li] + weight
                  mdenscont[Li] = mdenscont[Li] + shear[si]*sigmaci*weight
                  nsources[Li] = nsources[Li] + 1

                  ;;print,'Added: ',shear[si]*sigmaci
              ENDIF 

              ;; increment the index
              gausind = gausind + 1L
          ENDFOR

          IF Ngood GT 0 THEN BEGIN 

              ;; mean denscont for this lens
              mdenscont[Li] = mdenscont[Li]/wsum[Li]
          
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
  sigzsout = sigzsout[ws]

  ;; normalize the distributions: we added total(nsources) normalized
  ;; gaussians, so divide by total(nsources)
  nsourcestot = total(nsources[w])
  ;;deltazsi = zsi[1]-zsi[0]
  pofzitot = pofzitot/nsourcestot

  pofzitottrue = pofzitottrue/nsourcestot

  plot, zsi, pofzitot

  tmp=create_struct('nlens', nlens, $
                    'weight_type', weight_type, $
                    'zs', zs, $
                    'pofzs', pofzs, $
                    'meandenscont', meandenscont, $
                    'denscont', denscont[w], $
                    'zl', zL[w], $
                    'msigma_crit', msigma_crit[w], $
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
                    'wsum', wsum[w], $ ;sum of weights after zcuts
                    'nsources', nsources[w], $   ;# in each after zcuts
                    'nsourcestot', nsourcestot, $
                    'mdenscontsum', mdenscontsum[w], $ ; using overall dist
                    'nsourcesum', nsourcesum[w],$ ;number sources used here
                    'nstotal', nstotal, $; total sources used in sums
                    'indices', indices, $; good objects
                    'sigzs', sigzsout, $ ; zerr
                    'zstrue', zstrue, $  ; input redshift for sources
                    'zsest', zsest)      ; estimated (drawn from dist)
  print
  print,'Output file: ',outfile 
  mwrfits, tmp, outfile, /create

  ptime,systime(1)-tt

END 
