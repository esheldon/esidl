PRO read_corshape, run, rerun, camcol, clr, struct, allcols=allcols, hirata=hirata

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: read_corshape, run, rerun, camcol, clr, struct, allcols=allcols, hirata=hirata'
      return
  ENDIF 

  corshape_names, run, rerun, camcol, clr, psfile, fitsfile, $
                    extno=extno, allcols=allcols, hirata=hirata

  print,'Reading corshape file for slope correction: ',fitsfile
  struct = mrdfits(fitsfile, 1, /silent)

END 
