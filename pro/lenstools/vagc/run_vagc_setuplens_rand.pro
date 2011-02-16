PRO run_vagc_setuplens_rand

  stripes = [9,10,11,12,13,14,15,$
             27,28,29,30,31,32,33,34,35,36,37, $
             76,86]

  rmin = 20.0
  rmax = 10000.0

  nbin = 18
  logbin = 1

;  vagc_setuplens, rmin, rmax, nbin, stripes, /logbin

  randFileNums = 10+lindgen(10)

  vagc_setuprand, randFileNums, rmin, rmax, nbin, stripes, /logbin

return
  ;; this is for rlrg catalog masks
  
  stripes = [10,11,12,13,14,15,$
             27,28,29,30,31,32,33,34,35,36,37, $
             76,82,86]

  rmin = 20.0
  rmax = 1020.0
  binsize = 100.0

  vagc_setuplens, rmin, rmax, binsize, stripes, $
    /rlrgMask

  randFileNums = 20+lindgen(10)

  vagc_setuprand, randFileNums, rmin, rmax, binsize, stripes, $
                    /rlrgMask

END 
