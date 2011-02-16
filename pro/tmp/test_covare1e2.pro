PRO test_covare1e2, cat

  ;; this is to test covariance between 
  ;; e1 and e2

  n=n_elements(cat)

  w=where(tag_names(cat) EQ 'E1')
  IF n_elements(cat[0].e1) EQ 1 THEN BEGIN 

      me1 = mean(cat.e1)
      me2 = mean(cat.e2)
      
      sdev = sdev(cat.e1)
      covar = total( (cat.e1 - me1)*(cat.e2-me2) )/(n-1)
      codev = sqrt( abs(covar) )

  ENDIF ELSE BEGIN 
      
      size2 = cat.ixx[2] + cat.iyy[2]
      e1 = (cat.ixx[2] - cat.iyy[2])/size2
      e2 = 2.*cat.ixy[2]/size2
      me1 = mean(e1)
      me2 = mean(e2)
      sdev = sdev(e1)
      covar = total( (e1-me1)*(e2-me2) )/(n-1)
      codev = sqrt(abs(covar))
      print,sdev,covar,codev,sdev/sqrt(n)

      me1 = mean(cat.e1[2])
      me2 = mean(cat.e2[2])
      
      sdev = sdev(cat.e1[2])
      covar = total( (cat.e1[2] - me1)*(cat.e2[2]-me2) )/(n-1)
      codev = sqrt(abs(covar))
  ENDELSE 

  print,sdev,covar,codev,sdev/sqrt(n)
return
END 
