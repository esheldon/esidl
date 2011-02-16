PRO run_vagc_combine_stripe_spec, lss=lss, letter=letter, post=post, $
                                  sample=sample, $
                                  outdir=outdir

  vagc_getstripes, stripes, nstripe, $
    lss=lss, letter=letter, post=post, sample=sample


  nstripe = n_elements(stripes)
  FOR i=0L, nstripe-1 DO BEGIN 

      stripe = stripes[i]
      print,'Doing stripe '+ntostr(stripe)

      vagc_combine_stripe_spec, stripe, $
        outdir=outdir, lrg=lrg, $
        lss=lss, letter=letter, post=post, sample=sample

  ENDFOR 

END 
