PRO sdss_specid_extract, specID, $
                         plate=plate, $
                         mjd=mjd, $
                         fiberID=fiberID, $
                         type=type, $
                         lineOrIndexOrZ=lineOrIndexOrZ

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: sdss_specid_extract, specID, '
      print,'            plate=, mjd=, fiberID=, type=, '
      print,'            lienOrIndexOrZ='
      return
  ENDIF 

  IF arg_present(plate) THEN $
    plate = long( ishft(specID,-48) )
  
  IF arg_present(mjd) THEN $
    mjd = long( ishft(specID,-32) AND '0000FFFF'XL )

  IF arg_present(fiberID) THEN $
    fiberID = fix( ishft(specID,-22) AND '000003FF'XL )

  IF arg_present(type) THEN $
    type = byte( ishft(specID,-16) AND '00000003F'XL )

  IF arg_present(lineOrIndexOrZ) THEN $
    lineOrIndexOrZ = long( specID AND '0000FFFF'XL )

END 
