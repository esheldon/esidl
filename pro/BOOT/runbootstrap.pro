PRO runbootstrap, data, nbin, nresamp, datamean, dataerr, reverse=reverse

  ;; assumes data is of form, fltarr(nbin, #independent measures)
  ;; unless nbin = 1 then fltarr(# independent measures)
  ;; no weights

  IF keyword_set(reverse) THEN nm = n_elements(data[0,*]) $
  ELSE nm = n_elements(data[*,0])
  
  IF nbin EQ 1 THEN senddata = reform( replicate(data[0], nresamp ), nresamp, 1)$
  ELSE senddata = replicate(data[0], nresamp, nbin)
  FOR samp=0L, nresamp-1 DO BEGIN 
      
      ;; random set of indices
      ii = round( (nm-1)*randomu(seed, nm) )
      FOR bin=0L, nbin-1 DO BEGIN
          ;; what is order of input array?
          IF keyword_set(reverse) THEN BEGIN 
              senddata[samp, bin] = mean_check(data[bin, ii])
          ENDIF ELSE BEGIN 
              senddata[samp, bin] = mean_check(data[ii, bin])
          ENDELSE 
      ENDFOR 
  ENDFOR 
 
  bootstrap, senddata, datamean, dataerr

END 
