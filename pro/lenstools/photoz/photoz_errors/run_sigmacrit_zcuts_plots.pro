PRO run_sigmacrit_zcuts_plots

  dir = '/net/cheops1/data0/esheldon/pzcuts_sigmacrit/inputz/'
  cd,dir
  files = findfile("test_sigmacrit*.fit")
  nf = n_elements(files)

  FOR i=0L, nf-1 DO BEGIN 
 
      sigmacrit_zcuts_plots,files[i],/doplot
      sigmacrit_zcuts_plots,files[i],/doplot,/refix

  ENDFOR 


END 
