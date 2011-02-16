PRO logbin_spacing, minx, maxx, nbin, $
                    bin_min, bin_max, $
                    multiple=multiple, $
                    bin_multiple=bin_multiple, $
                    bin_mean=bin_mean, dim=dim

  IF n_params() LT 5 THEN BEGIN 
      print,'-Syntax: logbin_spacing, minx, maxx, nbin, '
      print,'              bin_min, bin_max, '
      print,'              bin_mean=, '
      print,'              /multiple, dim='
      print,'If /multiple, then maxx is actually the multiple to apply'
      print,'Mean is the geometric mean for the input dimension'
      return
  ENDIF 

  IF keyword_set(multiple) THEN BEGIN 

      bin_multiple = maxx

  ENDIF ELSE BEGIN 

      logminx = alog10(minx)
      logmaxx = alog10(maxx)

      ;; this is the maximum binsize to get this nbin
      binsize = ( logmaxx - logminx )/nbin
      bin_multiple = 10.^binsize

  ENDELSE 

  bin_min  = dblarr(nbin)
  bin_max  = bin_min

  bin_min[0] = minx
  bin_max[0] = minx*bin_multiple

  FOR i=1, nbin-1 DO BEGIN 
      bin_min[i] = bin_max[i-1]
      bin_max[i] = bin_min[i]*bin_multiple
  ENDFOR 

  IF arg_present(bin_mean) THEN BEGIN 
      IF n_elements(dim) EQ 0 THEN BEGIN 
          print,'  You must enter dim= (the dimension) to calculate the mean'
          print,'  Ignoring'
          return
      ENDIF 
      bin_mean = dblarr(nbin)
      FOR i=0, nbin-1 DO bin_mean[i] = gmean(bin_min[i], bin_max[i], dim[0])
  ENDIF 

END 
