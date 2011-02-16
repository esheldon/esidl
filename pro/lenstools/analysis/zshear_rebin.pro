PRO zshear_rebin, t, meanr, sigma, sigmaerr, wsum, orthosig, orthosigerr

  IF n_params() LT 5 THEN BEGIN 
      print,'-Syntax: zshear_rebin, struct, meanr, sigma, sigmaerr, wsum, orthosig, orthosigerr'
      return
  ENDIF 

  IF tag_exist(t[0],'meanr') THEN BEGIN 
      domeanr = 1
  ENDIF ELSE BEGIN 
      domeanr = 0
  ENDELSE 

  nst = n_elements(t)
  nbin = n_elements(t[0].rsum)
  
  new_nbin = nbin/2
  
  nbinmod = nbin MOD 2

  IF nst EQ 1 THEN BEGIN 
      meanr = fltarr(new_nbin)
  ENDIF ELSE BEGIN 
      meanr = fltarr(new_nbin, nst)
  ENDELSE 
  sigma = meanr
  sigmaerr = meanr
  orthosig = meanr
  orthosigerr = meanr
  wsum = meanr

  FOR j=0L, nst-1 DO BEGIN 
      IF (nst GT 1) AND ((j MOD 100) EQ 0) THEN print,'.',format='(a,$)'
      FOR i=0L, new_nbin-1 DO BEGIN 

          IF nbinmod NE 0 AND i EQ new_nbin-1 THEN BEGIN 
              bin1 = i*2
              bin2 = i*2+2
          ENDIF ELSE BEGIN 
              bin1 = i*2
              bin2 = i*2 + 1
          ENDELSE 
          

          tsigma       = t[j].sigma[bin1:bin2]
          tsigmaerr    = t[j].sigmaerr[bin1:bin2]
          torthosig    = t[j].orthosig[bin1:bin2]
          torthosigerr = t[j].orthosigerr[bin1:bin2]

          trsum  = t[j].rsum[bin1:bin2]
          tnpair = t[j].npair[bin1:bin2]
          twsum  = t[j].rsum[bin1:bin2]

          w=where(tsigmaerr NE 0.0, nw)
          IF nw NE 0 THEN BEGIN 

              wmom, tsigma[w], tsigmaerr[w], $
                    wmean, wsig, werr
              sigma[i,j] = wmean
              sigmaerr[i,j] = werr

              wmom, torthosig[w], torthosigerr[w], $
                    wmean, wsig, werr
              orthosig[i,j] = wmean
              orthosigerr[i,j] = werr

              IF domeanr THEN BEGIN 
                  meanr[i,j] = $
                    total(trsum[w])/total(tnpair[w])
              ENDIF ELSE BEGIN 
                  meanr[i,j] = $
                    total(trsum[w])
              ENDELSE 
              
              wsum[i,j] = total(twsum[w])
          ENDIF 
      ENDFOR 
  ENDFOR 

END 
