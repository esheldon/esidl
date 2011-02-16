PRO shearjackknife_bystripe, lenses, random, $
                             sigma, sigmaerr, covariance, $
                             rsigma, rsigmaerr, rcovariance, $
                             wuse=wuse, Nsubin=Nsubin, ortho=ortho, $
                             rebin=rebin, sampval=sampval, rsampval=rsampval

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: shearjackknife_bystripe, lenses, random, sigma, sigmaerr, covariance, wuse=wuse, Nsub=Nsub'
      return
  ENDIF 

;  plot,lenses.clambda,lenses.ceta,psym=3,/ynozero
;  w=where(lenses.sigmaerr[12] NE 0)
;  oplot, lenses[w].clambda, lenses[w].ceta,psym=3,color=!green

  IF keyword_set(ortho) THEN BEGIN 
      tt = tag_exist(lenses[0], 'orthosig', index=sigtag)
      tt = tag_exist(lenses[0], 'orthosigerr', index=sigerrtag)
  ENDIF ELSE BEGIN 
      tt = tag_exist(lenses[0], 'sigma', index=sigtag)
      tt = tag_exist(lenses[0], 'sigmaerr', index=sigerrtag)
  ENDELSE 

  Nlenses = n_elements(lenses)
  Nrand = n_elements(random)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; get the stripe for each lens
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  runs = lenses[rem_dup(lenses.run)].run
  nrun = n_elements(runs)

  ;; get the stripes for these runs
  read_stripe_index, si
  match, runs, si.run, mr, msi
  IF n_elements(mr) NE n_elements(runs) THEN message,'Not all runs matched'
  mstripes = si[msi].stripe
  runs = runs[mr]

  ;; now associate a stripe with each lens
  stripes = intarr(Nlenses)
  FOR i=0L, nrun-1 DO BEGIN
      w=where(lenses.run EQ runs[i])
      stripes[w] = mstripes[i]
  ENDFOR 
  rmds = rem_dup(mstripes)
  ustripes = mstripes[rmds]

  wu = where(ustripes NE 30)
  ustripes = ustripes[wu]
  Nstripe = n_elements(ustripes)

  ;; random points stripes
  rstripes = eta2stripenum(random.ceta)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; the bins to use
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(wuse) NE 0 THEN BEGIN 
      Nbin = n_elements(wuse)
  ENDIF ELSE BEGIN 
      Nbin = n_elements(lenses[0].rsum)
      wuse = lindgen(nbin)
  ENDELSE  

  sigma  = dblarr(Nbin)
  sigmaerr = sigma
  covariance = dblarr(Nbin, Nbin)

  rsigma = sigma
  rsigmaerr = sigmaerr
  rcovariance = covariance

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; number of subsamples to use per stripe
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(Nsubin) EQ 0 THEN Nsubin = 200

  Nsub = Nsubin*Nstripe
  factor = (Nsub - 1.)/Nsub

  sampval = dblarr(Nsub, Nbin)
  jval = dblarr(Nbin)

  rsampval = sampval
  rjval = jval

  print
  IF keyword_set(ortho) THEN print,'Doing orthosig'
  print,'Stripes: ',ustripes
  print,'Number of stripes: ',Nstripe
  print,'Number of subsamples per stripe: ',Nsubin
  print,'Number of subsamples: ',Nsub
  print,'Number of objects: ',Nlenses
  print,'Number of random: ',Nrand
  print,'Number of bins: ',Nbin
  print,'Jackknifing'

  ;; Use histogram to get random for each lens
  ;; using histogram this way, with min=0, forces a bin for every
  ;; integer from 0 to the max. That way, we can use
  ;; arr1[i] as subscript for reverse_indices! 
  
  ;; match multi gets all matches, but also returnes reverse_indices
  ;; for our use.  Assumes each lens zindex is represented in random
;  match_multi, lenses.zindex, random.zindex, rmatch, reverse_indices=revind
  
  ;; this gets all matches

  ssh = total(lenses.sshsum)/total(lenses.wsum_ssh)

  jsubarr = lonarr(nsub,nbin)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; loop over the radial bins first
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR bin=0L, Nbin-1 DO BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; lenses which contributed to this bin, over
      ;; all stripes
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; Lenses
      IF keyword_set(ortho) THEN BEGIN 
          corr = 1.0 
      ENDIF ELSE BEGIN 
          corr = total(lenses.wsum[wuse[bin]], /double)/total(random.wsum[wuse[bin]], /double)*float(Nrand)/float(Nlenses)/ssh
      ENDELSE 

      weights = lenses.wsum[wuse[bin]]
      sigsum_tot = total( lenses.(sigtag)[wuse[bin]]*weights, /double)
      wtot = total( weights, /double)
      sigmean = sigsum_tot/wtot*corr
;      print,'reg mean: ',sigmean


      ;; Random points


      rweights = random.wsum[wuse[bin]]
      rsigsum_tot = total( random.sigma[wuse[bin]]*rweights, /double)
      rwtot = total( rweights, /double)
      rsigmean = rsigsum_tot/rwtot

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; loop over the stripes: this is first division for
      ;; subsampling
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      FOR ist=0L, Nstripe-1 DO BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; lenses that contribute from this stripe
          ;; wst subscripts the overall array "lenses"
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          wst = where(stripes EQ ustripes[ist], Nst)
          ;;wst = wgood[wst]

          ;; random points from this stripe
          rwst = where(rstripes EQ ustripes[ist], rNst)
          ;;rwst = rwgood[rwst]

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; correction factor for this radial bin.
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          ;; Now divide this stripe into chunks by lambda
          ;; s will also subscript overall arrays
          s = sort(lenses[wst].clambda)
          s = wst[s]
          
          mini = 0L
          maxi = Nst-1
          binsize = float(maxi-mini)/Nsubin
          
          rs = sort(random[rwst].clambda)
          rs = rwst[rs]

          rmini = 0L
          rmaxi = rNst-1
          rbinsize = float(rmaxi-rmini)/Nsubin

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; will use rev_ind for subscripting the sort array
          ;; (and thus the weights array, which is already
          ;;  sorted)
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          hist = histogram( lindgen(Nst), binsize = binsize, $
                            rever=rev_ind )
          rhist = histogram( lindgen(rNst), binsize = rbinsize, $
                             rever=rrev_ind)

          ;; weights already sorted now
          weights_st = weights[s]
          rweights_st = rweights[rs]

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; loop over the subsamples in this stripe
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          FOR j=0L, Nsubin-1 DO BEGIN 

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; index for this subsample
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              jsub = ist*Nsubin + j

              ;; Lenses
              ;; any in this subsample? the answer should always be yes
              IF rev_ind[j] NE rev_ind[j+1] THEN BEGIN 
                  ;; subscript lenses with s[wt] and weights with wt
                  wt = rev_ind[ rev_ind[j]:rev_ind[j+1]-1 ]

                  ;; subscript for lenses
                  wt2 = s[wt]

                  sigmod = sigsum_tot - $
                    total(lenses[wt2].(sigtag)[wuse[bin]]*weights_st[wt], /double)
                  wmod = wtot - total(weights_st[wt], /double)
                  
                  sampval[jsub, bin] = sigmod/wmod*corr
              ENDIF ELSE message,'Why empty bin?'
              
              ;; Random
              ;; any in this subsample? the answer should always be yes
              IF rrev_ind[j] NE rrev_ind[j+1] THEN BEGIN 
                  ;; subscript lenses with s[wt] and weights with wt
                  rwt = rrev_ind[ rrev_ind[j]:rrev_ind[j+1]-1 ]

                  ;; subscript for lenses
                  rwt2 = rs[rwt]

                  rsigmod = rsigsum_tot - $
                    total(random[rwt2].sigma[wuse[bin]]*rweights_st[rwt], /double)
                  rwmod = rwtot - total(rweights_st[rwt], /double)
                  
                  rsampval[jsub, bin] = rsigmod/rwmod
              ENDIF ELSE message,'Why empty bin?'

          ENDFOR ;; lambda bins
                    
      ENDFOR ;; stripes

      jval[bin] = mean_check(sampval[*, bin], /double)
      sigma[bin] = sigmean + (Nsub-1)*(sigmean - jval[bin])

      rjval[bin] = mean_check(rsampval[*, bin], /double)
      rsigma[bin] = rsigmean + (Nsub-1)*(rsigmean - rjval[bin])

;      print,'Jack mean: ',sigma[bin]
  ENDFOR 

  ;; Now jackknife covariance between variables
  FOR jbin=0L, Nbin-1 DO BEGIN 
      wj=where(sampval[*,jbin] NE 0)
      rwj=where(rsampval[*,jbin] NE 0)

      FOR ibin=jbin, Nbin-1 DO BEGIN 

          wi=where(sampval[*,ibin] NE 0)
          rwi=where(rsampval[*,ibin] NE 0)

          tmp = total( (sampval[wi,ibin]-jval[ibin])*(sampval[wj,jbin]-jval[jbin]), /double)
          covariance[jbin, ibin] = factor*tmp

          tmp = total( (rsampval[rwi,ibin]-rjval[ibin])*(rsampval[rwj,jbin]-rjval[jbin]), /double)
          rcovariance[jbin, ibin] = factor*tmp

          IF jbin NE ibin THEN BEGIN
              covariance[ibin, jbin] = covariance[jbin, ibin]
              rcovariance[ibin, jbin] = rcovariance[jbin, ibin]
          ENDIF 
          IF jbin EQ ibin THEN BEGIN
              sigmaerr[ibin] = sqrt( covariance[ibin, ibin] )
              rsigmaerr[ibin] = sqrt( rcovariance[ibin, ibin] )
          ENDIF 

      ENDFOR 
  ENDFOR 

  sampval=0
  jval=0

  rsampval=0
  rjval=0

return
END 

