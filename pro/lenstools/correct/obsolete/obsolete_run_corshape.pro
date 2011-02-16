PRO run_corshape,run,rerun, onefile=onefile
  

  FOR camcol=1,6 DO BEGIN 
      clr=1
      corshape2, run, rerun, camcol, clr, onefile=onefile
      clr=2
      corshape2, run, rerun, camcol, clr, onefile=onefile
      clr=3
      corshape2, run, rerun, camcol, clr, onefile=onefile
  ENDFOR 

  clr=1
  corshape_allcols2, run, rerun, clr, onefile=onefile
  clr=2
  corshape_allcols2, run, rerun, clr, onefile=onefile
  clr=3
  corshape_allcols2, run, rerun, clr, onefile=onefile

  bell,10

END 
