FUNCTION match_hist, x1, x2, binsize, min=min, max=max, status=status, $
                     x1h=x1h, x2h=x2h, new_x2h=new_x2h, $
                     indices1=indices1
    
  status = 1
  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: indices1 = match_hist(x1, x2, binsize, min=, max=, indices1=, status=)'
      print,'Chooses values of x2 such that the histograms of '
      print,'   x1 and x2 have the same shape'
      print,'Defaults:  min=min(x1), max=max(x1)'
      return,-1
  ENDIF 

  IF n_elements(min) EQ 0 THEN min = min(x1)
  IF n_elements(max) EQ 0 THEN max = max(x1)

  ;; Histogram the first array
  if arg_present(indices1) then begin
      x1h = histogram(x1, min=min, max=max, binsize=binsize, rev=rev1)
  endif else begin
      x1h = histogram(x1, min=min, max=max, binsize=binsize)
  endelse

  ;; second histogram
  x2h = histogram(x2, min=min, max=max, binsize=binsize, rev=rev2)

  nbin = n_elements(x1h)

  ;; The percentage of x2 to keep in each bin
  perc_keep = fltarr(nbin)
  ratio = fltarr(nbin)
  new_x2h = lonarr(nbin)
  
  wh = where(x1h GT 0 AND x2h GT 0, ngood)

  IF ngood EQ 0 THEN BEGIN 
      print,'No values from x2 found in x1 range ' + $
        '['+ntostr(min)+', '+ntostr(max)+']'
      return,-1
  ENDIF 

  ;; normalize by the maximum in x1 histogram
  ratio[wh] = float(x1h[wh])/x2h[wh]
  maxratio = max(ratio, maxi)

  ;colprint, x1h[wh], x2h[wh]

  ;pplot, ratio[wh], psym=10
  ;pplot, ratio[wh], psym=8,/over,color=c2i('red')

  ;; The percentage to keep
  perc_keep[wh] = ratio[wh]/maxratio

  ;; pointer array to keep track of the ones from x2 we want
  ptrlist2 = ptrarr(nbin)
  if arg_present(indices1) then ptrlist1 = ptrarr(nbin)

  FOR j=0L, nbin-1 DO BEGIN 

      IF x1h[j] GT 0 THEN BEGIN 

          IF rev2[j] NE rev2[j+1] THEN BEGIN 
              w = rev2[ rev2[j]:rev2[j+1]-1 ]
              
              nrtot = n_elements(w)
              
              IF perc_keep[j] LT 1 THEN BEGIN 
                  ;; randomly choose perc_keep of them
                  nkeep = round(perc_keep[j]*nrtot) > 1
                  ind = lindgen(nrtot)
                  ;; now randomize the order of the indices and keep the
                  ;; first nkeep of them
                  s = sort( randomu(seed, nrtot) )
                  ind = ind[s]
                  ind = ind[0:nkeep-1]
                  
                  w = w[ind]
                  
              ENDIF ELSE BEGIN
                  nkeep = nrtot
              ENDELSE 

              ptrlist2[j] = ptr_new(w, /no_copy)
              new_x2h[j] = nkeep
             
              if arg_present(indices1) then begin
                  w1 = rev1[ rev1[j]:rev1[j+1]-1]
                  ptrlist1[j] = ptr_new(w1,/no_copy)
              endif
              
          ENDIF 
      ENDIF 

      ;;print,j,x1h[j],x2h[j],ratio[j],perc_keep[j],new_x2h[j]

  ENDFOR 

  indices2 = combine_ptrlist(ptrlist2)
  if arg_present(indices1) then indices1 = combine_ptrlist(ptrlist1)

  status = 0
  return, indices2

                 
END 
