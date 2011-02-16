PRO vagc_survey_rot, runs, rerun

  sdss_spec_dir = sdssidl_config('spec_dir')
  indir = sdss_spec_dir + 'blanton/astrom/'
  outdir = sdss_spec_dir + 'blanton/survey_rot/'

  nrun = n_elements(runs)

  FOR i=0L, nrun-1 DO BEGIN 
      
      FOR camcol=1,6 DO BEGIN 
          
          sdss_survey_rot, runs[i], rerun, camcol, indir=indir, outdir=outdir

      ENDFOR 

  ENDFOR 

END 
