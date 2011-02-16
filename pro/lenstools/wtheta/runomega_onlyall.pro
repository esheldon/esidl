PRO runomega_onlyall, mclr, lumclr,doplot=doplot

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: runomega_onlyall, mclr, lumclr,doplot=doplot'
      return
  ENDIF 

  ;; wrapper for calc_m2l....

  ldir='/sdss5/data0/lensout/mass2light/'
  mdir='/sdss5/data0/lensout/stripe10/'

  lumnum = 'N5'
  mnum = 'N8'
  ;; always send rw, replaces rw with correct color from lumw within program
  lallf=ldir+'lumdens_match'+!colors[mclr]+''+mnum+'_wthetalumweq_stripe10_sum_rw_'+lumnum+'.fit'

  allf=mdir+'match'+$
    !colors[lumclr]+''+lumnum+'_main_zgal_gal_stripe10_'+!colors[mclr]+'_corr_'+mnum+'.fit'
;  allf=mdir+'match'+$
;    !colors[lumclr]+''+lumnum+'_main_zgal_gal_stripe10_comb_corr_'+mnum+'.fit'
  lumf=mdir+'match'+$
    !colors[lumclr]+''+lumnum+'_main_zgal_gal_stripe10_'+!colors[mclr]+'_meanlum_'+mnum+'.fit'
  IF NOT fexist(lumf) THEN BEGIN
      ;; This only works for single stripes
      print
      print,'Creating meanlum file'
      infile=repstr(allf, 'corr', 'lensum')
      make_meanlum_struct, infile
  ENDIF  
  lum=mrdfits(lumf,1)

  calc_m2l_omega_lumdens_onlyall, lumclr, lallf, allf, meanlum=lum.meanlum,/fixlum,doplot=doplot, /dofits, mclr=mclr

END 
