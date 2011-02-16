PRO run_vagc_combine_stripes_spec, lss=lss, letter=letter, post=post, $
                                   outdir=outdir

  vagc_getstripes, allstripes, nstripe, lss=lss, letter=letter, post=post

  ;; Standard samples
  stripes = [9,10,11,12,13,14,15]
  vagc_combine_stripes_spec, stripes, $
    lss=lss, letter=letter, post=post
  IF NOT keyword_set(lss) THEN BEGIN 
      vagc_combine_stripes_spec, stripes, /lrg, $
        lss=lss, letter=letter, post=post
  ENDIF 

  stripes = [27,28,29,30,31,32,33,34,35,36,37]
  vagc_combine_stripes_spec, stripes, $
    lss=lss, letter=letter, post=post
  IF NOT keyword_set(lss) THEN BEGIN 
      vagc_combine_stripes_spec, stripes, /lrg, $
        lss=lss, letter=letter, post=post
  ENDIF 
return
  stripes = [9,10,11,12,13,14,15,27,28,29,30,31,32,33,34,35,36,37]
  vagc_combine_stripes_spec, stripes, $
    lss=lss, letter=letter, post=post
  IF NOT keyword_set(lss) THEN BEGIN 
      vagc_combine_stripes_spec, stripes, /lrg, $
        lss=lss, letter=letter, post=post
  ENDIF 

  vagc_combine_stripes_spec, allstripes, $
    lss=lss, letter=letter, post=post
  IF NOT keyword_set(lss) THEN BEGIN 
      vagc_combine_stripes_spec, allstripes, /lrg, $
        lss=lss, letter=letter, post=post
  ENDIF 

END 
