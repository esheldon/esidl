FUNCTION sdss_photoid_range, run, rerun, camcol, field, id

  minfield = 0
  maxfield = 5000
  mincamcol = 1
  maxcamcol = 6
  minid = 0
  maxid = 5000

  IF n_params() EQ 5 THEN BEGIN 
      id = sdss_photoid(run,rerun,camcol,field,id)
      return,[id,id]
  ENDIF ELSE IF n_params() EQ 4 THEN BEGIN 

      min_photoid = sdss_photoid(run,rerun,camcol,field,minid)
      max_photoid = sdss_photoid(run,rerun,camcol,field,maxid)
      return,[min_photoid,max_photoid]

  ENDIF ELSE IF n_params() EQ 3 THEN BEGIN 

      min_photoid = sdss_photoid(run,rerun,camcol,minfield,minid)
      max_photoid = sdss_photoid(run,rerun,camcol,maxfield,maxid)
      return,[min_photoid,max_photoid]

  ENDIF ELSE IF n_params() EQ 2 THEN BEGIN 

      min_photoid = sdss_photoid(run,rerun,mincamcol,minfield,minid)
      max_photoid = sdss_photoid(run,rerun,maxcamcol,maxfield,maxid)
      return,[min_photoid,max_photoid]


  ENDIF ELSE IF n_params() EQ 1 THEN BEGIN 

      rerun = sdss_rerun(run)
      IF rerun EQ -1 THEN message,'Unknown run: you must enter rerun'
      min_photoid = sdss_photoid(run,rerun,mincamcol,minfield,minid)
      max_photoid = sdss_photoid(run,rerun,maxcamcol,maxfield,maxid)
      return,[min_photoid, max_photoid]

  ENDIF ELSE BEGIN 
      message,'-Syntax: idrange = sdss_photoid_range(run, [rerun,camcol,field,id])'
  ENDELSE 

END 
