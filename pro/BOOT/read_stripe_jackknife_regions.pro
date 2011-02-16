PRO read_stripe_jackknife_regions, stripes, jackKnifeRegionsPtr, clusters=clusters

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: read_stripe_jackknife_regions, stripes, jackKnifeRegionsPtr, '
      print,'                  /clusters'
      return
  ENDIF 

  nstripe = n_elements(stripes)
  jackKnifeRegionsPtr = ptrarr(nstripe)

  maskDir = sdssidl_config('shapecorr_dir') +'masks/'
  IF keyword_set(clusters) THEN BEGIN 
      jackKnifeRegionFiles = $
        maskDir+$
        'jackknife_regions_stripe'+stripe2string(stripes)+'_MaxBCG.fit'

  ENDIF ELSE BEGIN 
      jackKnifeRegionFiles = $
        maskDir+$
        'jackknife_regions_stripe'+stripe2string(stripes)+'_sphpoly.fit'
  ENDELSE 

  FOR ist=0L, nstripe-1 DO BEGIN 
      IF NOT fexist(jackKnifeRegionFiles[ist]) THEN BEGIN 
          ptr_free, jackKnifeRegionsPtr
          message,'File does not exist: '+jackKnifeRegionFiles[ist]
      ENDIF 
      print,'Reading region file: ',jackKnifeRegionFiles[ist]
      jackKnifeRegionsPtr[ist] = $
        ptr_new(mrdfits(jackKnifeRegionFiles[ist],1,/silent) , /no_copy)
  ENDFOR 

END 
