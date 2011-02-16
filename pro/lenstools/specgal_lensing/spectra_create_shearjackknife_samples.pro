PRO spectra_create_shearjackknife_samples, type, ext=ext, $
                                           stripes=stripes, clr=clr, $
                                           subLumClr=subLumClr, $
                                           randNum=randNum, lenses=lenses

  IF n_elements(type) EQ 0 THEN type=''
  IF n_elements(ext) EQ 0 THEN ext='N1.fit'

  IF n_elements(stripes) EQ 0 THEN BEGIN 
      stripes = [9,10,11,12,13,14,15,$
                 27,28,29,30,31,32,33,34,35,36,37,$
                 76,86]
  ENDIF 
  IF n_elements(clr) EQ 0 THEN BEGIN 
      clr = [1,2,3]
  ENDIF 

  nRand = n_elements(randNum)

  IF nRand EQ 0 THEN BEGIN 
      front = type + 'zgal_gal_'
  ENDIF ELSE BEGIN 
      front = type + 'zrand'+ntostr(randNum)+'_'
  ENDELSE 

  ;; Read jackknife regions
  read_stripe_jackknife_regions, stripes, jackKnifeRegionsPtr

  nfile = n_elements(front)
  FOR i=0L, nfile-1 DO BEGIN 

      lensumfile_name, stripes, front[i], clr, dir, name, $
        /hirata, /recorr, ext=ext, subLumClr=subLumClr

      lensumFile = dir + name

      outDir = dir + 'jackknife_samples/'
      outFile = outDir + 'jackknife_samples_'+name

      print
      print,'Input file: ',lensumFile
      print,'Output file: ',outfile
      IF n_elements(lenses) EQ 0 OR nRand NE 0 THEN BEGIN 
          lenses = mrdfits(lensumFile,1)
      ENDIF 

      print
      print,'Jackknifing sigma'
      create_shearjackknife_samples, $
        jackKnifeRegionsPtr, lenses, jackStruct, $
        wuse=wuse

      print
      print,'Writing to file: ',outfile
      mwrfits, jackStruct, outfile, /create

      print
      print,'Jackknifing ortho'
      create_shearjackknife_samples, $
        jackKnifeRegionsPtr, lenses, orthoJackStruct, $
        wuse=wuse, /ortho

      ;; This will append the ortho
      print
      print,'Appending ortho to file: ',outfile
      mwrfits, orthoJackStruct, outfile

  ENDFOR 

  ptr_free, jackKnifeRegionsPtr

END 
