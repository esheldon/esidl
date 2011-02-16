PRO run_generate_meanscinv, lrg_sources=lrg_sources, rlrg_sources=rlrg_sources

  npts = 200
  stripes = [9,10,11,12,13,14,15,$
             27,28,29,30,31,32,33,34,35,36,37,$
             76, 86]

;  generate_meanscinv, stripes, [1,2,3], npts


  stripes = [9,10,11,12,13,14,15,$
             27,28,29,30,31,32,33,34,35,36,37,$
             76, 82, 86]

;  generate_meanscinv, stripes, [1,2,3], npts
;  generate_meanscinv, stripes, [1,2,3], npts, /lrg_sources

  ;; they don't include stripe 9
  stripes = [10,11,12,13,14,15,$
             27,28,29,30,31,32,33,34,35,36,37,$
             76, 82, 86]

  generate_meanscinv, stripes, [1,2,3], npts, $
    /rlrg_sources

END 
