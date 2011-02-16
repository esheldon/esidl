PRO find_jackknife_ids, jackKnifeRegionsPtr, clambda,  ceta, jackKnifeIDs

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: find_jackknife_ids, jackKnifeRegionsPtr, clambda,  ceta, jackKnifeIDs'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The lenses MUST correspond to the stripes
  ;; that are input
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Get the stripes from the jackKnifeRegionsPtr
  ;; Also count the regions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nLenses = n_elements(clambda)

  nstripe = n_elements(jackKnifeRegionsPtr)
  stripes = intarr(nstripe)

  print
  print,'There are '+ntostr(nstripe)+' stripes'
  print,'Counting sub-regions'
  Nsub = 0L
  FOR ist=0L, nstripe-1 DO BEGIN
      
      stripe_jRegions = *jackKnifeRegionsPtr[ist]
 
      stripes[ist] = stripe_jRegions.stripe

      Nsub = Nsub + n_elements(stripe_jRegions.jackKnifeID)

  ENDFOR 
  print,'Found '+ntostr(Nsub)+' sub regions'

  lstripes  = eta2stripenum(ceta)
  rmds = rem_dup(lstripes)
  ustripes = lstripes[rmds]

  match, stripes, ustripes, mst, must
  IF (mst[0] EQ -1) OR (n_elements(mst) NE n_elements(stripes)) THEN BEGIN 
      message,'coords contain stripes other than the input stripes',/inf
      print
      forprint,'Point stripes: ',ustripes
      print
      forprint,'Input stripes: ',stripes
      message,'Aborting'
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Define the arrays we will use
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  jackKnifeIDs = replicate(-1L, nLenses)

  nLensesUsed = 0L
  FOR ist=0L, nstripe-1 DO BEGIN 

      stripe = stripes[ist]
      wst = where(lstripes EQ stripe, nwst)
      print
      print,'Found '+strn(nwst,len=5)+$
        ' objects in stripe '+stripe2string(stripe)
      ;;wait,2

      stripe_jregions = *jackKnifeRegionsPtr[ist]

      nStripeSub = n_elements(stripe_jregions.jackKnifeID)

      FOR isub=0L, nStripeSub-1 DO BEGIN 
          
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Get objects in this subregion
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          lamMin = stripe_jregions.chunkLamMin[isub]
          lamMax = stripe_jregions.chunkLamMax[isub]

          wsub = $
            where(clambda[wst] GT lamMin AND clambda[wst] LE lamMax, niSub)

          print,'  Found '+strn(niSub,len=4)+' objects in sub-region '+$
            ntostr(stripe_jregions.jackKnifeID[isub])

;          IF niSub EQ 0 THEN stop

          nLensesUsed = nLensesUsed + niSub
          IF niSub NE 0 THEN BEGIN 

              wsub = wst[wsub]
              jackKnifeIds[wsub] = stripe_jregions.jackKnifeID[isub]


          ENDIF ;; any in sub-region?

      ENDFOR ;; loop over sub-regions in this stripe

  ENDFOR 

  print,'Used '+ntostr(nLensesUsed)+'/'+ntostr(nLenses)

END 
