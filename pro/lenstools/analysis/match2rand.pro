PRO match2rand, lensum, randsum, wrand, numlist, hist, rev_ind

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: match2rand, lensum, randsum, wrand, numlist, histogram, rev_ind'
       print,'This works even for subsamples of subsamples'
      return
  ENDIF 

  match_multi, lensum.zindex, randsum.zindex, wrand, $
               numlist=numlist, hist=hist, reverse_indices=rev_ind

  return

  nlens = n_elements(lensum)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; using histogram this way, with min=0, forces a bin for every
  ;; integer from 0 to the max. That way, we can use
  ;; lensum[i].zindex as subscript for rev_ind! Then it
  ;; works for sub-samples
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  hist=histogram(randsum.zindex, reverse_indices=rev_ind, $
              min=0, max=max(lensum.zindex))

  ptrlist = ptrarr(nlens)
  numlist = lonarr(nlens)
  nrand = 0L
  FOR i=0L, nlens-1 DO BEGIN 
      
      zi = lensum[i].zindex

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Use of histogram with min=0 is why this works for 
      ;; subsamples of subsamples: when you do histogram
      ;; on integers with bin=1, min=0, bins are created
      ;; for each integer from 0 to max(zindex), and 
      ;; we can just subscript rev_ind with zi. Note that
      ;; if the wrong random points are sent then we may
      ;; accidentally subscript rev_ind with too large a
      ;; number and crash the program.
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF rev_ind[zi] NE rev_ind[zi+1] THEN BEGIN 
          
          srand = rev_ind[ rev_ind[zi]:rev_ind[zi+1] -1 ]
          nsrand = n_elements(srand)
          ptrlist[i] = ptr_new(srand)
          numlist[i] = nsrand
          nrand = nrand + nsrand

      ENDIF 

  ENDFOR 

  wrand = lonarr(nrand)
  beg=0L
  FOR i=0L, nlens-1 DO BEGIN 
      IF numlist[i] NE 0 THEN BEGIN 
          wrand[beg:beg+numlist[i]-1] = *ptrlist[i]
      ENDIF 
      ptr_free, ptrlist[i]
      beg = beg + numlist[i]
  ENDFOR 

  return
END 
