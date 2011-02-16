PRO runlumdens

  lumnum = 'N5'
  mnum = 'N8'

  indir = '/sdss5/data0/lensout/stripe10/'

;  file = indir+'matchg'+mnum+'_wthetalumweq_stripe10_sum_rw_'+lumnum+'.fit'
;  rfile = indir+'matchg'+mnum+'_wthetarandlumweq_stripe10_sum_rw_'+lumnum+'.fit'

;  lumdens, file, rfile, nrad=12, /doplot

  ;; only did r last time
  file = indir+'matchr'+mnum+'_wthetalumweq_stripe10_sum_rw_'+lumnum+'.fit'
  rfile = indir+'matchr'+mnum+'_wthetarandlumweq_stripe10_sum_rw_'+lumnum+'.fit'

  lumdens, file, rfile, nrad=12, /doplot

;  file = indir+'matchi'+mnum+'_wthetalumweq_stripe10_sum_rw_'+lumnum+'.fit'
;  rfile = indir+'matchi'+mnum+'_wthetarandlumweq_stripe10_sum_rw_'+lumnum+'.fit'

;  lumdens, file, rfile, nrad=12, /doplot



END 
