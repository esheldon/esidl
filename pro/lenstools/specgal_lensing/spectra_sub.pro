PRO spectra_sub, type, randNum, $
                 stripes=stripes, $
                 ext=ext, $
                 nooverwrite=nooverwrite, len=len

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: spectra_sub, type, randNum, stripes, $'
      print,'    ext=ext, $'
      print,'    basedir=basedir, $'
      print,'    nooverwrite=nooverwrite'
      return
  ENDIF 

  IF keyword_set(nooverwrite) THEN overwrite=0 ELSE overwrite=1
  IF n_elements(stripes) EQ 0 THEN BEGIN 
      stripes = [9,10,11,12,13,14,15,$
                 27,28,29,30,31,32,33,34,35,36,37,$
                 76,86]
  ENDIF 
  clr = [1,2,3]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Lenssum file name
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  front = 'zgal_gal_'
  lensumfile_name, stripes, front, clr, lensumDir, lensumFileName, $
    /hirata, /recorr, ext=ext

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Rand lensum file names
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  lensumfile_name, stripes, '', clr, lensumDir, rLensumFileName, $
    /hirata, /recorr, ext=ext
  
  rfronts = 'zrand'+ntostr(randNum)+'_'
  rLensumFileNames = rfronts + rLensumFileName

  CASE type OF
      'lumfourbinnum': BEGIN 

          lensumFile = lensumDir + lensumFileName
          rLensumFiles = lensumDir + rlensumFileNames

          binLumByNum_getPerc, 4, perc

          print
          print,'Reading: ',lensumFile
          IF n_elements(len) EQ 0 THEN len = mrdfits(lensumFile, 1)

          FOR iclr = 0,4 DO BEGIN 
              outdir = dir + 'sublum/'+!colors[iclr]+'/'

              ;;!!!! BINBYNUM HAS BUILT-IN Z CUTS!!!! 0.6,0.02
              binLumByNum, len, iclr, perc, wStruct

              addstr = 'lum1fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, lensumFile, rLensumFiles, $
                           addstr=addstr,uselens=len, $
                           overwrite=overwrite, indices=wStruct.w1, $
                           outdir=outdir
              
              addstr = 'lum2fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, lensumFile, rLensumFiles, $
                           addstr=addstr,uselens=len, $
                           overwrite=overwrite, indices=wStruct.w2, $
                           outdir=outdir
              
              addstr = 'lum3fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, lensumFile, rLensumFiles, $
                           addstr=addstr,uselens=len, $
                           overwrite=overwrite, indices=wStruct.w3, $
                           outdir=outdir
              
              addstr = 'lum4fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, lensumFile, rLensumFiles, $
                           addstr=addstr,uselens=len, $
                           overwrite=overwrite, indices=wStruct.w4, $
                           outdir=outdir

          ENDFOR 

          delvarx, len

      END 
      'matchLRGNew': BEGIN 

          ;; match our outputs from new lrgs to the old non-lrg photozs
          ;; Here lrg refers to lrg sources and all refers to all sources

          stripeString = $
            'stripe10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_76_82_86'
          inDir = '~/lensout/'+stripeString+'/'

          allFile = inDir + 'zgal_gal_' + $
            stripeString+'_gri_recorr_h_lensum_N1.fit'
          lrgFile = inDir + 'zgal_gal_' + $
            stripeString+'_gri_rlrg_recorr_h_lensum_N2.fit'

          randAllFiles = inDir + 'zrand' + ntostr(randNum)+'_'+$
            stripeString+'_gri_recorr_h_lensum_N1.fit'
          randlrgFiles = inDir + 'zrand' + ntostr(randNum)+'_'+$
            stripeString+'_gri_rlrg_recorr_h_lensum_N2.fit'

          print
          print,'Reading all photoz outputs: ',allFile
          allStruct = mrdfits(allFile,1)

          print
          print,'Reading new lrg outputs: ',lrgFile
          lrgStruct = mrdfits(lrgFile, 1)

          print
          print,'Matching'
          sphoto_match, lrgStruct, allStruct, matchesLRG, matchesAll

          help,matchesLRG,matchesAll

          addstr = 'matchLRGNew'
          print,'Sub All sample: ',addstr
          zsub_sample_multirand, allFile, randAllFiles, $
            addstr=addstr,uselens=allStruct, $
            overwrite=overwrite, indices=matchesAll, $
            outdir=outdir

          delvarx, allStruct

          addstr = 'matchAll'
          print,'Sub lrg sample: ',addstr
          zsub_sample_multirand, lrgFile, randlrgFiles, $
            addstr=addstr,uselens=lrgStruct, $
            overwrite=overwrite, indices=matchesLRG, $
            outdir=outdir

          
      END 
      'lumthreebinnum_matchLRGNew': BEGIN 

          stripeString = $
            'stripe10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_76_82_86'
          inDir = '~/lensout/'+stripeString+'/'

          allFront = 'matchLRGNew_'
          lrgFront = 'matchAll_'

          allFile = inDir + allFront + 'zgal_gal_' + $
            stripeString+'_gri_recorr_h_lensum_N1.fit'
          lrgFile = inDir + lrgFront + 'zgal_gal_' + $
            stripeString+'_gri_rlrg_recorr_h_lensum_N2.fit'

          randAllFiles = inDir + allFront + 'zrand' + ntostr(randNum)+'_'+$
            stripeString+'_gri_recorr_h_lensum_N1.fit'
          randlrgFiles = inDir + lrgFront + 'zrand' + ntostr(randNum)+'_'+$
            stripeString+'_gri_rlrg_recorr_h_lensum_N2.fit'


          print
          print,'Reading all file: ',allFile
          allStruct = mrdfits(allFile, 1)

          print
          print,'Reading lrg file: ',lrgFile
          lrgStruct = mrdfits(lrgFile, 1)

          ;; percentages
          binLumByNum_getPerc, 3, perc

          ;; only r-band for now
          FOR iclr = 2,2 DO BEGIN 
              outdir = inDir + 'sublum/'+!colors[iclr]+'/'

              ;;!!!! BINBYNUM HAS BUILT-IN Z CUTS!!!! 0.6,0.02
              binLumByNum, allStruct, iClr, perc, wAllStruct
              binLumByNum, lrgStruct, iClr, perc, wLRGStruct

              addstr = 'lum1threebinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, allFile, randAllFiles, $
                           addstr=addstr,uselens=allStruct, $
                           overwrite=overwrite, indices=wAllStruct.w1, $
                           outdir=outdir
              
              addstr = 'lum2threebinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, allFile, randAllFiles, $
                           addstr=addstr,uselens=allStruct, $
                           overwrite=overwrite, indices=wAllStruct.w2, $
                           outdir=outdir
              
              addstr = 'lum3threebinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, allFile, randAllFiles, $
                           addstr=addstr,uselens=allStruct, $
                           overwrite=overwrite, indices=wAllStruct.w3, $
                           outdir=outdir
              
              outdir = inDir + 'sublum/'+!colors[iclr]+'/'

              addstr = 'lum1threebinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, lrgFile, randlrgFiles, $
                           addstr=addstr,uselens=lrgStruct, $
                           overwrite=overwrite, indices=wLrgStruct.w1, $
                           outdir=outdir
              
              addstr = 'lum2threebinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, lrgFile, randlrgFiles, $
                           addstr=addstr,uselens=lrgStruct, $
                           overwrite=overwrite, indices=wLrgStruct.w2, $
                           outdir=outdir
              
              addstr = 'lum3threebinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, lrgFile, randlrgFiles, $
                           addstr=addstr,uselens=lrgStruct, $
                           overwrite=overwrite, indices=wLrgStruct.w3, $
                           outdir=outdir
              
          ENDFOR 

          delvarx, allStruct, lrgStruct


      END 
      'lumfourbinnum_matchLRGNew': BEGIN 

          stripeString = $
            'stripe10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_76_82_86'
          inDir = '~/lensout/'+stripeString+'/'

          allFront = 'matchLRGNew_'
          lrgFront = 'matchAll_'

          allFile = inDir + allFront + 'zgal_gal_' + $
            stripeString+'_gri_recorr_h_lensum_N1.fit'
          lrgFile = inDir + lrgFront + 'zgal_gal_' + $
            stripeString+'_gri_rlrg_recorr_h_lensum_N2.fit'

          randAllFiles = inDir + allFront + 'zrand' + ntostr(randNum)+'_'+$
            stripeString+'_gri_recorr_h_lensum_N1.fit'
          randlrgFiles = inDir + lrgFront + 'zrand' + ntostr(randNum)+'_'+$
            stripeString+'_gri_rlrg_recorr_h_lensum_N2.fit'


          print
          print,'Reading all file: ',allFile
          allStruct = mrdfits(allFile, 1)

          print
          print,'Reading lrg file: ',lrgFile
          lrgStruct = mrdfits(lrgFile, 1)

          ;; percentages
          binLumByNum_getPerc, 4, perc

          ;; only r-band for now
          FOR iclr = 2,2 DO BEGIN 
              outdir = inDir + 'sublum/'+!colors[iclr]+'/'

              ;;!!!! BINBYNUM HAS BUILT-IN Z CUTS!!!! 0.6,0.02
              binLumByNum, allStruct, iClr, perc, wAllStruct
              binLumByNum, lrgStruct, iClr, perc, wLRGStruct

              addstr = 'lum1fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, allFile, randAllFiles, $
                           addstr=addstr,uselens=allStruct, $
                           overwrite=overwrite, indices=wAllStruct.w1, $
                           outdir=outdir
              
              addstr = 'lum2fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, allFile, randAllFiles, $
                           addstr=addstr,uselens=allStruct, $
                           overwrite=overwrite, indices=wAllStruct.w2, $
                           outdir=outdir
              
              addstr = 'lum3fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, allFile, randAllFiles, $
                           addstr=addstr,uselens=allStruct, $
                           overwrite=overwrite, indices=wAllStruct.w3, $
                           outdir=outdir
              
              addstr = 'lum4fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, allFile, randAllFiles, $
                           addstr=addstr,uselens=allStruct, $
                           overwrite=overwrite, indices=wAllStruct.w4, $
                           outdir=outdir

              outdir = inDir + 'sublum/'+!colors[iclr]+'/'

              addstr = 'lum1fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, lrgFile, randlrgFiles, $
                           addstr=addstr,uselens=lrgStruct, $
                           overwrite=overwrite, indices=wLrgStruct.w1, $
                           outdir=outdir
              
              addstr = 'lum2fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, lrgFile, randlrgFiles, $
                           addstr=addstr,uselens=lrgStruct, $
                           overwrite=overwrite, indices=wLrgStruct.w2, $
                           outdir=outdir
              
              addstr = 'lum3fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, lrgFile, randlrgFiles, $
                           addstr=addstr,uselens=lrgStruct, $
                           overwrite=overwrite, indices=wLrgStruct.w3, $
                           outdir=outdir
              
              addstr = 'lum4fourbinnum'
              print,'Sub sample: ',addstr
              zsub_sample_multirand, lrgFile, randlrgFiles, $
                           addstr=addstr,uselens=lrgStruct, $
                           overwrite=overwrite, indices=wLrgStruct.w4, $
                           outdir=outdir


          ENDFOR 

          delvarx, allStruct, lrgStruct


      END 
      ELSE: print,'Unknown type: '+type
  ENDCASE 

END 
