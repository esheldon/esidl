PRO read_runlist, runstruct, silent=silent

  on_error,2

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: read_runlist, runstruct'
      return
  ENDIF 

;; runstat 94 51075 82 N 336 56 151 -129 11 546 good

  file = sdssidl_config('data_dir')+'dbm/steves_runlist/run.par'

  format = 'A,L,L,L,A,L,L,L,L,L,L,A'
  readcol, file, $
    stname, run, mjd, stripe, strip, muStart, muEnd, slongStart, slongEnd, $
    frameStart, frameEnd, quality, format=format, silent=silent

  nrun = n_elements(run)
  runstruct = create_struct('run', 0L, $
                            'mjd', 0L, $
                            'stripe', 0L, $
                            'strip', "", $
                            'muStart', 0L, $
                            'muEnd', 0L, $
                            'slongStart', 0L, $
                            'slongEnd', 0L, $
                            'frameStart', 0L, $
                            'frameEnd', 0L, $
                            'quality', "")

  runstruct = replicate(runstruct, nrun)
  runstruct.run = run
  runstruct.mjd = mjd
  runstruct.stripe = stripe
  runstruct.strip = strip
  runstruct.muStart = muStart
  runstruct.muEnd = muEnd
  runstruct.slongStart = slongStart
  runstruct.slongEnd = slongEnd
  runstruct.frameStart = frameStart
  runstruct.frameEnd = frameEnd
  runstruct.quality = quality

return

  yanny_read, file, pdata

  runstruct = *pdata
  ptr_free, pdata

END 
