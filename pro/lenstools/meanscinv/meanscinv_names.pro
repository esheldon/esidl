PRO meanscinv_names, stripe, clr, fitfile, psfile, use_lambda=use_lambda, hirata=hirata, hudson=hudson, lrg_sources=lrg_sources, rlrg_sources=rlrg_sources

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: meanscinv_names, stripe(s), clr(s), fitfile, psfile, /use_lambda, /hirata, /lrg_sources, /rlrg_sources'
      print,'/hirata, /use_lambda are default'
      return
  ENDIF 

  ;; use_lambda is now the default
  IF n_elements(use_lambda) EQ 0 THEN use_lambda=1
  ;; hirata now default
  IF n_elements(hirata) EQ 0 THEN hirata=1

  IF keyword_set(use_lambda) THEN lamstr = '_lambda' ELSE lamstr = ''

  dir = sdssidl_config('SHAPECORR_DIR')+'sigmacrit/'
  psdir = dir + 'figures/'

  nclr = n_elements(clr)

  clr_string = clrarr2string(clr)
  stripe_str = stripearr2string(stripe)


  IF keyword_set(hirata) THEN hirstr = '_h' ELSE hirstr = ''
  IF keyword_set(hudson) THEN hudstr = '_hudson' ELSE hudstr=''

  IF keyword_set(lrg_sources) THEN BEGIN 
      lrgstr = '_lrg' 
  ENDIF ELSE BEGIN 
      IF keyword_set(rlrg_sources) THEN lrgstr = '_rlrg' ELSE lrgstr=''
  ENDELSE 

  fitfile = 'mean_scinv_interpolate_stripe'+stripe_str+'_'+clr_string+lrgstr+hirstr+hudstr+lamstr+'.fit'
  psfile = repstr(fitfile, '.fit', '.ps')

  fitfile = dir + fitfile
  psfile = psdir + psfile

END 
