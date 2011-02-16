PRO shear_jackknife, type, ext=ext, $
                     stripes=stripes, clr=clr, randNum=randNum, $
                     subLumClr=subLumClr, MaxBCG=MaxBCG

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; A wrapper for wjackknife that reads in the appropriate files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(type) EQ 0 THEN type=''
  IF n_elements(ext) EQ 0 THEN ext='N1.fit'

  IF n_elements(stripes) EQ 0 THEN BEGIN 
      stripes = [9,10,11,12,13,14,15,$
                 27,28,29,30,31,32,33,34,35,36,37,$
                 76,82,86]
  ENDIF 
  IF n_elements(clr) EQ 0 THEN BEGIN 
      clr = [1,2,3]
  ENDIF 

  nRand = n_elements(randNum)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; What type of files are we reading here?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF nRand EQ 0 THEN BEGIN 

      IF keyword_set(MaxBCG) THEN front='MaxBCG_' ELSE front = 'zgal_gal_'
      front = type + front
      outFront = front

  ENDIF ELSE BEGIN 

      IF keyword_set(MaxBCG) THEN radd='MaxBCG_' ELSE radd=''

      front    = type + radd + 'zrand'+ntostr(randNum)+'_'
      outFront = type + radd + 'zrand_'

  ENDELSE 

  lensumfile_name, stripes, outFront, clr, dir, name, $
    /hirata, /recorr, ext=ext, subLumClr=subLumClr
  sampDir = dir + 'jackknife_samples/'
  outFile = sampDir + 'jackknife_'+name
  print
  print,'Output jackknife File: ',outFile

  nFile = n_elements(front)

  ;; Read in the samples.  There can be multiple files

  FOR i=0L, nFile-1 DO BEGIN 

      lensumfile_name, stripes, front[i], clr, dir, name, $
        /hirata, /recorr, ext=ext, subLumClr=subLumClr

      sampFile = sampDir + 'jackknife_samples_'+name

      print
      print,'Input Sample File: ',sampFile
      tsamples      = mrdfits(sampFile, 1)
      tOrthoSamples = mrdfits(sampFile, 2)

      IF i EQ 0 THEN BEGIN 
          meanFile = dir + repstr(name, '_lensum', '')
          meanStruct = mrdfits(meanFile, 1)

          samples = tsamples
          orthoSamples = tOrthoSamples
      ENDIF ELSE BEGIN 

          samples.nLenses     = samples.nLenses + tsamples.nLenses
          samples.nLensesUsed = samples.nLensesUsed + tsamples.nLensesUsed
 
          samples.sigsum_tot = samples.sigsum_tot + tsamples.sigsum_tot
          samples.sigsum_sub = samples.sigsum_sub + tsamples.sigsum_sub
          samples.wsum_tot   = samples.wsum_tot + tsamples.wsum_tot
          samples.wsum_sub   = samples.wsum_sub + tsamples.wsum_sub

          orthoSamples.sigsum_tot = $
            orthoSamples.sigsum_tot + tOrthoSamples.sigsum_tot
          orthoSamples.sigsum_sub = $
            orthoSamples.sigsum_sub + tOrthoSamples.sigsum_sub
          orthoSamples.wsum_tot   = $
            orthoSamples.wsum_tot + tOrthoSamples.wsum_tot
          orthoSamples.wsum_sub   = $
            orthoSamples.wsum_sub + tOrthoSamples.wsum_sub

      ENDELSE 
  ENDFOR 

  print
  print,'Jackknifing'
  wjackknife, $
    samples.sigsum_sub, samples.wsum_sub, $
    sigma, sigmaerr, covariance
  wjackknife, $
    orthoSamples.sigsum_sub, orthoSamples.wsum_sub, $
    orthosig, orthosigerr, orthocov

  calc_correlation_matrix, covariance, corr, status
  IF status NE 0 THEN message,'Failed to calculate corr matrix for sigma'
  calc_correlation_matrix, orthocov, orthocorr, status
  IF status NE 0 THEN message,'Failed to calculate corr matrix for ortho'


  ;; We didn't use all the lenses. Scale the errors appropriately
;  IF nRand EQ 0 THEN BEGIN 
;      fac = sqrt(float(samples.nLensesUsed)/samples.nLenses)

;      sigmaerr = sigmaerr*fac
;      orthosigerr = orthosigerr*fac

;      covariance = covariance*fac^2
;      orthocov = orthocov*fac^2

;  ENDIF 


  s = create_struct('nLenses', samples.nLenses, $
                    'nLensesUsed', samples.nLensesUsed, $
                    'meanr', meanStruct.meanr, $
                    'sigma', sigma, $
                    'sigmaerr', sigmaerr, $
                    'covariance', covariance, $
                    'corr', corr, $
                    'orthosig', orthosig, $
                    'orthosigerr', orthosigerr, $
                    'orthocov', orthocov, $
                    'orthocorr', orthocorr)


  print
  print,'Writing to file: ',outfile

  mwrfits, s, outFile, /create

END 
