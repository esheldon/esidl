FUNCTION create_shearjackknife_samples, lensin, lensout, $
  jackknife_ids=jackknife_ids, $
  jackknife_file=jackknife_file, $
  wuse=wuse

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: jackstruct = create_shearjackknife_samples(lensin, lensout, jackknife_ids=, jackknife_file=, wuse=)'
      return, -1
  ENDIF 

  
  nLenses = n_elements(lensout)
  IF n_elements(lensin) NE nlenses THEN $
    message,'lensin and lensout must be same size'

  IF n_elements(jackknife_ids) EQ 0 THEN BEGIN 
      jackknife_ids = csurvey2jack(lensin.clambda, lensin.ceta, $
                                   file=jackknife_file)
  ENDIF 

;  pixnums = csurvey2pix(lensin.clambda, lensin.ceta)
  simpctable, colorlist=colorlist
  nc = n_elements(colorlist)
;  rmd = rem_dup(pixnums)
;  display_pixel, pixnums[rmd], /iso, resolution=256

  minjack = min(jackknife_ids, max=maxjack)
  IF minjack LT 0 THEN message,'Some not found in a jackknife region'

  h = histogram(jackknife_ids-minjack, min=0, rev=rev)

  ;; How many unique jackknife regions do we have?
  wh = where(h NE 0, Nsub)


  IF n_elements(wuse) NE 0 THEN BEGIN 
      Nbin = n_elements(wuse)
  ENDIF ELSE BEGIN 
      Nbin = n_elements(lensout[0].rsum)
      wuse = lindgen(nbin)
  ENDELSE  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Define the arrays we will use
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sigsum_tot = dblarr(Nbin)
  wsum_tot   = dblarr(Nbin)

  sigsum_sub = dblarr(Nsub, Nbin)
  wsum_sub   = dblarr(Nsub, Nbin)


  osigsum_tot = dblarr(Nbin)
  owsum_tot   = dblarr(Nbin)

  osigsum_sub = dblarr(Nsub, Nbin)
  owsum_sub   = dblarr(Nsub, Nbin)


  jackKnifeID = lonarr(Nsub)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;
  ;; Create the jackknife sub-samples
  ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nLensesUsed = 0L

  FOR sub=0L, Nsub-1 DO BEGIN 

      isub = wh[sub]

      IF rev[isub] EQ rev[isub+1] THEN message,'error in histogramming?'

      wsub = rev[ rev[isub]:rev[isub+1] -1 ]
      nwsub = n_elements(wsub)

      ;; plot the pixels.  This takes longer than the jackknifing
;      rmd = rem_dup(pixnums[wsub])
;      display_pixel, pixnums[wsub[rmd]], $
;        color=colorlist[ sub MOD nc ], /over_plot, resolution=256

      nLensesUsed = nLensesUsed + nwsub

      jackknifeid[sub] = jackknife_ids[wsub[0]]

      FOR bin=0L, Nbin-1 DO BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Generate sums for this sub-region and bin
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          weights = lensout[wsub].wsum[wuse[bin]]
          oweights = lensout[wsub].owsum[wuse[bin]]

          sigsum_sub[sub, bin] = $
            total( lensout[wsub].sigma[wuse[bin]]*weights, /double)
          wsum_sub[sub, bin] = total(weights)
          
          osigsum_sub[sub, bin] = $
            total( lensout[wsub].orthosig[wuse[bin]]*oweights, /double)
          owsum_sub[sub, bin] = total(weights)



          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Add these to the totals
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          sigsum_tot[bin] = sigsum_tot[bin] + sigsum_sub[sub, bin]
          wsum_tot[bin] = wsum_tot[bin] + wsum_sub[sub, bin]

          osigsum_tot[bin] = osigsum_tot[bin] + osigsum_sub[sub, bin]
          owsum_tot[bin] = owsum_tot[bin] + owsum_sub[sub, bin]
          
      ENDFOR ;; Loop over bins

  ENDFOR 

  print,'Used '+ntostr(nLensesUsed)+'/'+ntostr(nLenses)

  jackStruct = create_struct('nLenses',     nLenses, $
                             'nLensesUsed', nLensesUsed, $
                             $
                             'Nsub',        Nsub, $
                             'jackKnifeID', jackKnifeID, $
                             $
                             'sigsum_tot',  sigsum_tot, $
                             'wsum_tot',    wsum_tot, $
                             'sigsum_sub',  sigsum_sub, $
                             'wsum_sub',    wsum_sub, $
                             $
                             'osigsum_tot',  osigsum_tot, $
                             'owsum_tot',    owsum_tot, $
                             'osigsum_sub',  osigsum_sub, $
                             'owsum_sub',    owsum_sub)

  return,jackStruct

END 
