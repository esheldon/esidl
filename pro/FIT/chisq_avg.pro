PRO chisq_avg, data, covar, datamean, dataerr, status=status

  status = 0
  delvarx, datamean, dataerr, datavar

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: chisq_avg, data, covar [, datamean, dataerr, status=status]'
      return
  ENDIF 

  Ndata = n_elements(data)
  str = size(covar,/struct)
  IF str.n_dimensions NE 2 THEN BEGIN 
      IF str.n_elements NE Ndata THEN BEGIN
          IF str.n_elements EQ 1 THEN BEGIN
              ;; just replicate this n times
              dataerr2 = replicate(covar[0], Ndata)
              wmom, data, dataerr2, datamean, datavar, dataerr
              return
          ENDIF ELSE BEGIN
              message,'covar must be same size as data, an NdataXNdata array, or a scalar',/inf
              return
          ENDELSE 
      ENDIF ELSE BEGIN 
          wmom, data, covar, datamean, datavar, dataerr
          return
      ENDELSE 
  ENDIF ELSE BEGIN 
      IF str.dimensions[0] NE Ndata THEN BEGIN 
          message,'covar must be same size as data, an NdataXNdata array, or a scalar',/inf
          return
      ENDIF ELSE BEGIN 
          S = invert(covar, status, /double)
          IF status NE 0 THEN BEGIN
              print,covar
              stop
          ENDIF 
          
          ;;IF status NE 0 THEN return
          ;;print,'Doing full inversion'
          M = replicate(1., Ndata)

          dataerr2 = 1./total(S)

          datamean = dataerr2*( transpose(M)#reform(S##data) )
          datamean = datamean[0]
          dataerr = sqrt(dataerr2)

      ENDELSE 
  ENDELSE 
END 
