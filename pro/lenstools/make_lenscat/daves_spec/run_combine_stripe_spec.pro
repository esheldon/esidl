PRO run_combine_stripe_spec

  read_stripe_index,sid

  stripes = sid[rem_dup(sid.stripe)].stripe
  nst = n_elements(stripes)

  FOR i=0L, nst-1 DO BEGIN 

      combine_stripe_spec, stripes[i]
      combine_stripe_spec, stripes[i], /lrg

  ENDFOR 

  stripes = [10,11,12]
  combine_stripes_spec, stripes
  combine_stripes_spec, stripes, /lrg

  stripes = [9,10,11,12,13,14]
  combine_stripes_spec, stripes
  combine_stripes_spec, stripes, /lrg

  stripes = [9,10,11,12,13,14,15]
  combine_stripes_spec, stripes
  combine_stripes_spec, stripes, /lrg

  stripes = [28,29,30,31,32,33,34,35,36,37]
  combine_stripes_spec, stripes
  combine_stripes_spec, stripes, /lrg

  stripes = [35,36,37]
  combine_stripes_spec, stripes
  combine_stripes_spec, stripes, /lrg

END 
