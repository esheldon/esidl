PRO test_newphotoz_formula, meangauss, siggauss, modmean, modsig

  nmean = 20
  minmean = 0.01
  maxmean = 0.5
  meangauss = arrscl( findgen(nmean), minmean, maxmean )

  nsig = 20
  minsig = 0.01
  maxsig = 0.5
  siggauss = arrscl( findgen(nsig), minsig, maxsig )
  ;;siggauss = 0.05

  modmean = fltarr(nmean, nsig)
  modsig  = fltarr(nmean, nsig)

  FOR imean=0L, nmean-1 DO BEGIN 
      FOR isig=0L, nsig-1 DO BEGIN 
          
          newphotoz_formula, meangauss[imean], siggauss[isig], $
                             tmodmean, tmodsig

          modmean[imean, isig] = tmodmean
          modsig[imean, isig] = tmodsig

      ENDFOR 
  ENDFOR 

END 
