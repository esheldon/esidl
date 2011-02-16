PRO  profmean_rad, pixrad, arcsec, pixels

  IF n_params() EQ 0 THEN BEGIN
      print,'-Syntax: profmean_rad, pixrad, arcsec, pixels'
      return
  ENDIF 

  pixrad = [ 0.56, 1.69, 2.58, 4.41, 7.51, 11.58, 18.58, 28.55, 45.50, $
            70.51, 110.5, 172.5, 269.5, 420.5, 657.5 ]
  
  arcsec = [ 0.23, 0.68, 1.03, 1.76, 3.00, 4.63, 7.43, 11.42, 18.20, $
             28.20, 44.21, 69.00, 107.81, 168.20, 263.00 ]
      
  pixels = [ 1, 9, 21, 61, 177, 421, 1085, 2561, 6505, 15619, $
             38381, 93475, 228207, 555525, 1358149 ]

  return
END 
