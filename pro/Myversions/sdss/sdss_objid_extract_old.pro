PRO sdss_objid_extract, objID, $
                        run=run, $
                        rerun=rerun, $
                        camcol=camcol, $
                        field=field, $
                        id=id, $
                        sky_version=sky_version, $
                        first_field=first_field

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: sdss_objid_extract, objID, '
      print,'             run=, rerun=, camcol=, field=, id=, '
      print,'             sky_version=, first_field='
      return
  ENDIF 

  IF arg_present(run) THEN $
    run = long( ishft(objID,-32) AND '0000FFFF'XL )

  IF arg_present(rerun) THEN $
    rerun = fix( ishft(objID,-48) AND '000007FF'XL )

  IF arg_present(camcol) THEN $
    camcol = byte( ishft(objID,-29) AND '00000007'XL )

  IF arg_present(field) THEN $
    field = fix( ishft(objID,-16) AND '000000FFF'XL )

  IF arg_present(id) THEN $
    id = long( objID AND '0000FFFF'XL )

  IF arg_present(sky_version) THEN $
    sky_version = byte( ishft(objID,-59) AND '0000000F'XL )

  IF arg_present(first_field) THEN $
    first_field = byte( ishft(objID,-28) AND '000000001'XL )

END 
