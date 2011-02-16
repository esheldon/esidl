PRO logbin, x, minx, maxx, nbin, hist, rev_ind, $
            xbinned=xbinned, $
            bin_multiple=bin_multiple

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;
  ;; Finds the rev_ind for logarithmically spaced
  ;; bins. 
  ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: logbin, x, minx, maxx, nbin, hist, rev_ind, '
      print,'            xbinned=, '
      print,'            bin_multiple='
      return
  ENDIF 

  ww=where(x LT 0,nww)
  IF nww NE 0 THEN message,'X values must be greater than zero!!'

  logx = alog10(x)
  logminx = alog10(minx)
  logmaxx = alog10(maxx)

  ;; this is the maximum binsize to get this nbin
  binsize = ( logmaxx - logminx )/nbin
  bin_multiple = 10.^binsize

  hist = histogram( logx, binsize=binsize, $
                    min=logminx, max=logmaxx, $
                    rever=rev_ind )
  h_nbin = n_elements(hist)

  ;; We only get the requested nbin
  IF h_nbin GT nbin THEN hist = hist[0:nbin-1]

  IF arg_present(xbinned) THEN BEGIN 
      xbinned = dblarr(nbin)
      FOR i=0L, nbin-1 DO BEGIN 
          IF hist[i] NE 0 THEN BEGIN 
              w = rev_ind[ rev_ind[i]:rev_ind[i+1] -1 ]
              xbinned[i] = mean( x[w] )
          ENDIF 
      ENDFOR 
  ENDIF 

  return

  ;; THis was for test purposes
  ;; to use put y back in input list
  xbinned = fltarr(nbin+1)
  ybinned= xbinned
  logxbinned = xbinned
  num = lonarr(nbin+1)
  whist=where(hist NE 0, nhist)
  IF nhist NE 0 THEN BEGIN 
      
      FOR bin=0, nhist-1 DO BEGIN 

          binnum = whist[bin]

          w=rev_ind[ rev_ind[binnum]:rev_ind[binnum+1] -1 ]

          num[binnum] = n_elements(w)
          xbinned[binnum] = total(x[w])
          ybinned[binnum] = total(y[w])
          logxbinned[binnum] = total( logx[w] )

      ENDFOR 

  ENDIF 

  xbinned[whist] = xbinned[whist]/num[whist]
  ybinned[whist] = ybinned[whist]/num[whist]


END 
