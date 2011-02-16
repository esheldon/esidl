PRO boot_confidence, datax, data, covariance, data_resamp, Nresamp, $
          pmin, nmin, pminr, nminr, level, chisq_surf, powvals, normvals, $
          diagonal=diagonal

  ;; fit a power law to the data
  ;; Then use bootstrap samples and refit, plotting the points
  ;; on the contour plot and recording where the point lies
  ;; with respect to the 1-, 2-, and 3-sigma regions.

  ndata=n_elements(datax)
  dataerr = fltarr(ndata)
  FOR i=0L, ndata-1 DO dataerr[i] = sqrt(covariance[i,i])

  IF keyword_set(diagonal) THEN BEGIN 
      powrange=[-0.85, -0.7]
      normrange = [1.1, 1.3]

      esend = dataerr
  ENDIF ELSE BEGIN 
      powrange=[-0.8, -0.65]
      normrange = [1.1, 1.4]

      esend = covariance
  ENDELSE 
  nnorm=200L
  npow=200L

  pminr = fltarr(Nresamp)
  nminr = fltarr(Nresamp)
  level = intarr(Nresamp)

  ;;plot, powrange, normrange, psym=3, /iso, /ynozero

  color=!p.color
  ;; delta-chisq levels
  levels2 = [2.30, 6.17, 9.21]

  FOR ir=0L, Nresamp-1 DO BEGIN 
      
      IF ir EQ 0 THEN BEGIN 
          pow_chisq_conf_gen,datax,data, esend,$
                             powrange,normrange,npow,nnorm,$
                             chisq_surf,pmin,nmin,pl,ph,nl,nh,/iso,$
                             powvals=powvals, normvals=normvals

          deltachi = chisq_surf - min(chisq_surf)

      ENDIF ELSE BEGIN 
          send = reform(data_resamp[ir,*])

          IF keyword_set(diagonal) THEN BEGIN 
              ;; might as well use comfit
              fitpower, datax, send, dataerr, [norm, pow], $
                        yfit, a, /silent
              oplot, [a[1]], [a[0]], psym=3, color=color

              pminr[ir] = a[1]
              nminr[ir] = a[0]
          ENDIF ELSE BEGIN 
              ;; only my function uses full covariance matrix
              pow_chisq_conf_gen,datax,send,esend,$
                                 powrange,normrange,npow,nnorm,$
                                 ch,tpmin,tnmin,pl,ph,nl,nh,/nodisplay
              oplot, [tpmin], [tnmin], psym=3, color=color
              pminr[ir] = tpmin
              nminr[ir] = tnmin
          ENDELSE 

          ;; where does the new point live?
          ;; Is it within powrange and normrange?
          IF (pminr[ir] GE powrange[0]) AND (pminr[ir] LE powrange[1]) THEN BEGIN 
              
              ;; its inside the region
              ;; find closest powvals and normvals
              pdiff = abs(powvals - pminr[ir])
              ndiff = abs(normvals - nminr[ir])

              wp = ( where(pdiff EQ min(pdiff)) )[0]
              wn = ( where(ndiff EQ min(ndiff)) )[0]

              ;; what is deltachi at this location?
              td = deltachi[wp, wn]
              IF (td LE levels2[0]) THEN BEGIN
                  level[ir] = 0
              ENDIF ELSE IF (td LE levels2[1]) THEN BEGIN 
                  level[ir] = 1
              ENDIF ELSE IF (td LE levels2[2]) THEN BEGIN 
                  level[ir] = 2
              ENDIF ELSE level[ir] = 3

          ENDIF ELSE level[ir] = 3

      ENDELSE 
      
      print,ir

  ENDFOR 
END 
