PRO rebin_data, x, y, newx, newy, yerr=yerr, newyerr=newyerr, cov=cov, newcov=newcov

  ;; Rebin the data to half as many bins

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: rebin_data, x, y, newx, newy, yerr=yerr, newyerr=newyerr'
      return
  ENDIF 

  nx = n_elements(x)
  ny = n_elements(y)
  IF nx NE ny THEN message,'X must be same size as Y'
  
  nyerr = n_elements(yerr)
  IF nyerr NE 0 THEN BEGIN 
      doerr = 1
      IF nyerr NE nx THEN message,'Yerr must be same size as Y'
  ENDIF ELSE doerr = 0

  ncov = n_elements(cov)
  IF ncov NE 0 THEN BEGIN 
      docov = 1
      scov = size(cov)
      IF scov[0] NE 2 THEN $
        message,'cov must be a NX by NX matrix'
      IF scov[1] NE nx OR scov[2] NE nx THEN $
        message,'cov must be a NX by NX matrix'
  ENDIF ELSE docov = 0

  new_nx = nx/2

  nxmod = nx MOD 2

  newx = fltarr(new_nx)
  newy = newx

  delvarx, newyerr
  IF doerr THEN newyerr = newx

  delvarx, newcov
  IF docov THEN BEGIN
      cinv = invert(cov)
      newcinv = fltarr(new_nx, new_nx)
  ENDIF 
  FOR i=0L, new_nx-1 DO BEGIN 

      IF nxmod NE 0 AND i EQ new_nx-1 THEN BEGIN 
          bin1 = i*2
          bin2 = i*2+2
      ENDIF ELSE BEGIN 
          bin1 = i*2
          bin2 = i*2 + 1
      ENDELSE 

      tx = x[bin1:bin2]
      ty = y[bin1:bin2]
      nn = n_elements(tx)

      IF doerr THEN BEGIN 

          tyerr = yerr[bin1:bin2]

          ;; mean y and err
          wmom, ty, tyerr, wmean, wsig, werr
          newy[i] = wmean
          newyerr[i] = werr
          
          ;; mean x
          wmom, tx, tyerr, wmean, wsig, werr
          newx[i] = wmean

      ENDIF ELSE BEGIN 
          
          newy[i] = total(ty)/nn
          newx[i] = total(tx)/nn

      ENDELSE 

  ENDFOR 

  IF docov THEN BEGIN 
      FOR i=0L, new_nx-1 DO BEGIN 
          IF nxmod NE 0 AND i EQ new_nx-1 THEN BEGIN 
              xbin1 = i*2
              xbin2 = i*2+2
          ENDIF ELSE BEGIN 
              xbin1 = i*2
              xbin2 = i*2 + 1
          ENDELSE 
          FOR j=0L, new_nx-1 DO BEGIN 
              
              IF nxmod NE 0 AND j EQ new_nx-1 THEN BEGIN 
                  ybin1 = j*2
                  ybin2 = j*2+2
              ENDIF ELSE BEGIN 
                  ybin1 = j*2
                  ybin2 = j*2 + 1
              ENDELSE 
              
              IF docov THEN BEGIN 
                  newcinv[i,j] = total( cinv[xbin1:xbin2, ybin1:ybin2] )
              ENDIF 
          ENDFOR
      ENDFOR 
      
      newcov = invert(newcinv)
  ENDIF 

END 

