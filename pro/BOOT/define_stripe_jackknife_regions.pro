PRO define_stripe_jackknife_regions, stripe, nPointsPerChunk, clusters=clusters

  maskDir = sdssidl_config('shapecorr_dir') +'masks/'
  maskFile = maskDir + 'sphpoly_mask_sample14.dat'
  outDir = maskDir

  IF keyword_set(clusters) THEN BEGIN 
      outfile = $
        outdir + 'jackknife_regions_stripe'+stripe2string(stripe)+'_MaxBCG.fit'
      psfile = $
        outdir + 'jackknife_regions_stripe'+stripe2string(stripe)+'_MaxBCG.ps'
  ENDIF ELSE BEGIN 
      outfile = $
        outdir+'jackknife_regions_stripe'+stripe2string(stripe)+'_sphpoly.fit'
      psfile = $
        outdir+'jackknife_regions_stripe'+stripe2string(stripe)+'_sphpoly.ps'
  ENDELSE 

  print
  print,'Writing to file: ',outfile
  
  begplot, name=psfile, /color, /landscape

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; length of the regions in lambda
  ;; There will be many empty chunks for small samples
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  lamChunkLength = 0.6                  ; degrees  

  IF !d.name EQ 'PS' THEN opcolor=!blue ELSE opcolor=!blue
  length_frac = 0.9

  primary_bound, stripe, bound_array

  lamMin = bound_array.lammin
  lamMax = bound_array.lammax
  etaMin = bound_array.etaMin
  etaMax = bound_array.etaMax
  eta = (etaMax + etaMin)/2.0

  ;; number of points to place within each 
  ;; lambda chunk for testing

  densityOfPoints = lamChunkLength/npointsPerChunk

  maxNchunk = round( (lamMax - lamMin)/lamChunkLength ) + 1

  ;; total number of test points
  nPoints = round( (lamMax-lamMin)/densityOfPoints )
  print
  print,'Number of check points: '+ntostr(nPoints)

  lam = arrscl(findgen(nPoints), lamMin, lamMax)
  eta = replicate(eta, nPoints)

  plot_eta = arrscl(randomu(seed, nPoints), $
                    etaMin, etaMax, $
                    arrmin=0.0, arrmax=1.0)
  plot, lam, plot_eta, psym=3, /ynozero

  ;; Now apply mask.
  IF NOT keyword_set(clusters) THEN BEGIN 
      compcut = 0.9
      print,'sphPoly Applying masks'
      print
      ;; First apply spectroscopic mask
      apply_sphpoly_mask, lam, eta, bad, good, $
        compcut=compcut, /lameta, $
        maskFile=maskFile

  ENDIF ELSE BEGIN 
      print,'Applying MaxBCG masks'
      print
      MaxBCG_apply_mask, lam, eta, bad, good
  ENDELSE 
  oplot, lam[bad], plot_eta[bad], psym=3, color=!red

  IF good[0] EQ -1 THEN BEGIN 
      print,'No unmasked regions (sphPoly)'
      endplot,/landfix
      return
  ENDIF 
  lam = lam[good]
  eta = eta[good]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; And the pixel mask too
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  apply_pixel_mask, lam, eta, bad, good, /nomissf
  IF good[0] EQ -1 THEN BEGIN 
      print,'No unmasked regions (pixel)'
      endplot,/landfix
      return
  ENDIF 
  lam = lam[good]
  eta = eta[good]

  ;; Now step by npoitsPerChunk to define boundaries
  ngood = n_elements(good)
  indarr = lindgen(ngood)
  h = histogram( indarr, binsize = npointsPerChunk, rev=rev)
  nh = n_elements(h)
  
  print,'maxNchunk = '+ntostr(maxNchunk)
  print,'number of bins = '+ntostr(nh)
  chunkLamMin = replicate(-1d, maxNchunk)
  chunkLamMax = chunkLamMin
  jackKnifeID = lonarr(maxNchunk)

  FOR i=0L, nh-1 DO BEGIN 

      jackKnifeID[i] = stripe*10L^4 + i
      IF rev[i] NE rev[i+1] THEN BEGIN 
          w = rev[ rev[i]:rev[i+1]-1 ]
          nw = n_elements(w)

          minw = min(w, max=maxw)

          IF i EQ 0 THEN BEGIN 
              lamMin = lam[minw]
              lamMax = lam[maxw]
          ENDIF ELSE BEGIN 
              ;; chunkLamMin[i] = lam[minw]
              lamMin = chunkLamMax[i-1]
              lamMin = lam[minw] - densityOfPoints ; the step
              lamMax = lam[maxw]
          ENDELSE 

          lamDiff = lamMax - lamMin
          ;;print,'lamDiff = ',lamDiff
          ;; make sure its at least 90% the size we want
          ;; and contiguous
          IF ( (float(nw)/npointsPerChunk GE length_frac) AND $
               lamDiff LE (lamChunkLength+densityOfPoints)) THEN BEGIN 

              chunkLamMin[i] = lamMin
              chunkLamMax[i] = lamMax

              print, i+1, chunkLamMin[i], chunkLamMax[i], $
                lamDiff

              plot_box, $
                lam[minw], lam[maxw], $
                etaMin, etaMax, color=opcolor

              
          ENDIF 

      ENDIF 
  ENDFOR 

  endplot,/landfix

  w=where(chunkLamMin NE -1, nChunk)

  IF nChunk EQ 0 THEN BEGIN 
      print,'No chunks were defined'
      return
  ENDIF 

  help,maxNchunk,nChunk
  chunkLamMin = chunkLamMin[w]
  chunkLamMax = chunkLamMax[w]
  jackKnifeID = jackKnifeID[w]

  s=create_struct('jackKnifeId', jackKnifeID, $
                  'stripe', stripe, $
                  'etaMin', etaMin, $
                  'etaMax', etaMax, $
                  'chunkLamMin',chunkLamMin, $
                  'chunkLamMax',chunkLamMax)

  print
  print,'Writing to file: ',outfile
  mwrfits, s, outfile, /create

END 
