PRO compare_rlrg_sources, lrgLensum=lrgLensum, allLensum=allLensum, $
                          correct=correct

  IF !d.name EQ 'PS' THEN BEGIN 
      lrgColor = !blue
  ENDIF ELSE BEGIN 
      lrgColor = !green
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; First, we compare the case where we use the new lrg sources
  ;; There is the case with their new photozs and the old photozs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  pngFile = '~/plots/compare_rlrg_same_sources.png'
  window, xsize=640, ysize=450
  
  stripeString = $
    'stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_76_82_86'
  inDir = '~/lensout/'+stripeString+'/'

  rlrgOldPhotozFile = inDir + $
    'zgal_gal_'+stripeString+'_gri_rlrg_recorr_h_N1.fit'
  
  rlrgNewPhotozFile = inDir + $
    'zgal_gal_'+stripeString+'_gri_rlrg_recorr_h_N2.fit'

  rlrgOld = mrdfits(rlrgOldPhotozFile,1)
  rlrgNew = mrdfits(rlrgNewPhotozFile,1)

  newColor = lrgColor

  plot_density_contrast, rlrgOld, title='LRG Sources', aspect=!gratio, /center

  oploterror, rlrgNew.meanr, rlrgNew.sigma, rlrgNew.sigmaerr, $
    psym=4, color=newColor, errc=newColor

  legend, $
    ['Old LRG photoz', 'New LRG photoz'], $
    psym=[8, 4], /right, box=0, $
    color=[!p.color, newColor]

  print
  print,'Writing png file: ',pngFile
  write_png, pngFile, tvrd(/true)

  key = prompt_kbrd('Hit a key')
  GOTO, jump

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now we match up to the lenses run through with the "all" photoz
  ;; catalog
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  stripeString = $
    'stripe10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_76_82_86'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; First compare mean luminosities
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  inDir = '~/lensout/'+stripeString+'/'
  
  IF n_elements(allLensum) EQ 0 THEN BEGIN 
      lrgFile = inDir + $
        'matchAll_zgal_gal_'+stripeString+'_gri_recorr_h_lensum_N2.fit'

      lrgLensum = mrdfits(lrgFile,1)
  ENDIF 

  IF n_elements(allLensum) EQ 0 THEN BEGIN 
      allFile = inDir + $
        'matchLRGNew_zgal_gal_'+stripeString+'_gri_recorr_h_lensum_N1.fit'

      allLensum = mrdfits(allLensum, 1)
  ENDIF 

  clr = 2
  lrgLumSolar = calc_sdss_lumsolar(lrgLensum.absPetroMag[2],clr=clr)
  allLumSolar = calc_sdss_lumsolar(allLensum.absPetroMag[2],clr=clr)

  minmag = -24.
  maxmag = -17.
  wlrg = where(lrgLensum.absPetroMag[clr] GT minMag AND $
               lrgLensum.absPetroMag[clr] LT maxMag)
  wall = where(allLensum.absPetroMag[clr] GT minMag AND $
               allLensum.absPetroMag[clr] LT maxMag)

  wmom, lrgLumSolar[wlrg], blah, lrgMeanLum, blah, lrgMeanLumErr, $
    inputWeights = lrgLensum[wlrg].weight
  wmom, allLumSolar[wall], blah, allMeanLum, blah, allMeanLumErr, $
    inputWeights = allLensum[wall].weight

;  lrgwsum = total(lrgLensum[wlrg].weight)
;  lrgMeanLum = total(lrgLensum[wlrg].weight*lrgLumSolar[wlrg])/lrgwsum

;  allwsum = total(allLensum[wall].weight)
;  allMeanLum = total(allLensum[wall].weight*allLumSolar[wall])/allwsum

  print
  print,'Mean lum[2] lrg: '+$
    ntostr(lrgMeanLum)+' '+!plusminus+' '+ntostr(lrgMeanLumErr)
  print,'Mean lum[2] all: '+$
    ntostr(allMeanLum)+' '+!plusminus+' '+ntostr(allMeanLumErr)
  print

  factor = allMeanLum/lrgMeanLum

  key = prompt_kbrd('Hit a key')

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now plot the results: all lenses
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  stripeString = $
    'stripe10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_76_82_86'

  inDir = '~/lensout/combstripe/comb/'

  allFile = inDir + 'matchLRGNew_zgal_gal_'+$
    stripeString+'_gri_recorr_h_comb_N1.fit'
  lrgFile = inDir + 'matchAll_zgal_gal_'+$
    stripeString+'_gri_rlrg_recorr_h_comb_N2.fit'

  all = mrdfits(allFile, 1)
  lrg = mrdfits(lrgFile, 1)

  IF keyword_set(correct) THEN BEGIN
      rtitle = 'Corrected for weighting'
      lrg.sigma = lrg.sigma*factor
      lrg.sigmaerr = lrg.sigmaerr*factor
  ENDIF 

  plot_density_contrast, all, aspect=!gratio

  oploterror, lrg.meanr+10, lrg.sigma, lrg.sigmaerr, $
    psym=4, color=lrgColor, errc=lrgColor

  legend, $
    ['All photoz', 'New LRG photoz'], $
    psym=[8, 4], /right, box=0, $
    color=[!p.color, lrgColor]

  key = prompt_kbrd('Hit a key')

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now the ratio
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ratio = lrg.sigma/all.sigma
  ratioErr = ratio*sqrt( (lrg.sigmaerr/lrg.sigma)^2 + $
                         (all.sigmaerr/all.sigma)^2 )

  aploterror, !gratio, all.meanr, ratio, ratioerr, psym=8, $
    xtitle=!mpcxtitle2, ytitle='LRG_photoz/ All photoz'

  wmom, ratio, ratioerr, wmean, wsig, werr

  print
  print,'Mean ratio = '+ntostr(wmean)+' '+!plusminus+' '+ntostr(werr)
  print
  plot_box, 0.0, 1000.0, wmean-werr, wmean+werr, $
    /polyfill, /line_fill, orient=45

  key = prompt_kbrd('Hit a key')

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Luminosity bins
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

jump:

  stripeString = $
    'stripe10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_76_82_86'

  !p.multi=[0,0,2]
  inDir = '~/lensout/combstripe/comb/sublum/r/'

  window, xsize=640, ysize=800

  FOR lumNum=1,3 DO BEGIN 

      lnstr=ntostr(lumNum)

      pngFile = '~/plots/compare_rlrg_lum'+ntostr(lnstr)+'.png'

      allFront = 'lum'+lnstr+'threebinnum_matchLRGNew_'
      lrgFront = 'lum'+lnstr+'threebinnum_matchAll_'

      allFile = inDir + allFront + 'zgal_gal_'+$
        stripeString+'_gri_recorr_h_comb_N1.fit'
      lrgFile = inDir + lrgFront + 'zgal_gal_'+$
        stripeString+'_gri_rlrg_recorr_h_comb_N2.fit'
      
      all = mrdfits(allFile, 1)
      lrg = mrdfits(lrgFile, 1)

      print
      print,'All mean r: ',all.tmeanlum
      print,'LRG mean r: ',lrg.tmeanlum
      print,'LRG/All:    ',lrg.tmeanlum/all.tmeanlum

      IF keyword_set(correct) THEN BEGIN
          rtitle = 'Corrected for weighting'
          lrg.sigma = lrg.sigma*factor
          lrg.sigmaerr = lrg.sigmaerr*factor
      ENDIF 

      yrange = prange(all.sigma,lrg.sigma,all.sigmaerr,lrg.sigmaerr,/slack)
            
      plot_density_contrast, all, aspect=!gratio, yrange=yrange, /center
      
      oploterror, lrg.meanr+10, lrg.sigma, lrg.sigmaerr, $
        psym=4, color=lrgColor, errc=lrgColor
      
      legend, $
        ['All photoz', 'New LRG photoz'], $
        psym=[8, 4], /right, box=0, $
        color=[!p.color, lrgColor]
      
;      key = prompt_kbrd('Hit a key')
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Now the ratio
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ratio = lrg.sigma/all.sigma
      ratioErr = ratio*sqrt( (lrg.sigmaerr/lrg.sigma)^2 + $
                             (all.sigmaerr/all.sigma)^2 )

      aploterror, !gratio, all.meanr, ratio, ratioerr, psym=8, $
        xtitle=!mpcxtitle2, ytitle='LRG_photoz/ All photoz', /center

      wmom, ratio, ratioerr, wmean, wsig, werr
  
      print
      print,'Mean ratio = '+ntostr(wmean)+' '+!plusminus+' '+ntostr(werr)
      plot_box, 0.0, 1000.0, wmean-werr, wmean+werr, $
        /polyfill, /line_fill, orient=45

      mess = '<ratio> = '+$
        ntostr(wmean, 4, /round)+!csym.plusminus+ntostr(werr, 4, /round)
      legend, mess, box=0

      print
      print,'Writing png file: ',pngFile
      write_png, pngFile, tvrd(/true)

      IF lumnum NE 4 THEN key = prompt_kbrd('Hit a key')

  ENDFOR 

  !p.multi=0

END 
