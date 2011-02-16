PRO run_spectra_wrapper, depth=depth, $
                         newPhotoZ=newPhotoZ, randNum=randNum, $
                         rlrg_sources=rlrg_sources, $
                         pixelMaskFlags=pixelMaskFlags, $
                         rlrgMask=rlrgMask, $
                         noTestQuad=noTestQuad, $
                         lenscat=lenscat,$
                         scat=scat, revind=revind,$
                         extno=extno, outdir=outdir

;  pixelMaskFlags = 1
;  rlrgMask = 1

  ;; no stripe 10 in rlrgMask catalog
;  stripes = [10,11,12,13,14,15,$
;             27,28,29,30,31,32,33,34,35,36,37,$
;             76,82,86]
;  nbinOrBinsize=100.

  stripes = [9,10,11,12,13,14,15,$
             27,28,29,30,31,32,33,34,35,36,37, $
             76,86]
  nbinOrBinsize = 18
  logbin = 1
  callCPP=1
  clr = [1,2,3]
  rmin = 20.0
  rmax = 10000.0

  outdir='~/tmp/'

  nrand = n_elements(randnum)
  IF nrand NE 0 THEN BEGIN 

      FOR i=0L, nrand-1 DO BEGIN 

          ;; We re-use the scat but not the random lens catalog

          run_spectra, stripes, clr, nbinOrBinsize, $
            depth=depth, $
            rmin=rmin, rmax=rmax, $
            rlrg_sources=rlrg_sources, noTestQuad=noTestQuad, $
            newPhotoZ=newPhotoZ, $
            scat=scat, revind=revind, $
            pixelMaskFlags=pixelMaskFlags, rlrgMask=rlrgMask, $
            randNum=randNum[i], $
            logbin=logbin, $
            callCPP=callCPP, $
            extno=extno, outdir=outdir

      ENDFOR 

  ENDIF ELSE BEGIN 

      run_spectra, stripes, clr, nbinOrBinsize, $
        depth=depth, $
        rmin=rmin, rmax=rmax, $
        rlrg_sources=rlrg_sources, noTestQuad=noTestQuad, $
        newPhotoZ=newPhotoZ, $
        scat=scat, lenscat=lenscat, revind=revind, $
        pixelMaskFlags=pixelMaskFlags, rlrgMask=rlrgMask, $
        logbin=logbin, $
        callCPP=callCPP, $
        extno=extno, outdir=outdir

  ENDELSE 


END 
