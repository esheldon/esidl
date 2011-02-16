PRO make_compare_wthetalens, stripe, ext=ext, overwrite=overwrite

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: sub_lumwtheta, stripe'
      return
  ENDIF 
  
  ;; /lumw for files with total luminosity and
  ;; npair weighted by lum^beta

  ;; N3 is 0.1 L* stuff
  ;; N2 is magcut stuff
  ;; N1 is 0.1 L* stuff to 2 Mpc

  colors=['u','g','r','i','z']
  msun = [6.39,5.07,4.62,4.52,4.48]
  mstar = [-18.34,-20.04,-20.83,-21.26,-21.55]
  lstar = 10.0^((mstar-msun)/(-2.5))

  IF n_elements(ext) EQ 0 THEN ext = 'N3.fit'

  runstr = 'stripe'+ntostr(long(stripe))+'_'

  indir = '/sdss6/data0/wtheta/'

  IF keyword_set(logbin) THEN lgstr = 'lg_' ELSE lgstr='_'
  front = 'wthetalumw_'
  rfront = 'wthetarandlumw_'

  ff = indir + front+runstr+'lumlensum_rw_'+ext
  rff = indir + rfront+runstr+'lumlensum_rw_'+ext
  
  print,ff,rff
  lensum=mrdfits(ff,1)
  randsum=mrdfits(rff,1)

  compare = mrdfits('/sdss5/data0/lensout/stripe10/wthetalumw_stripe10_lumlensum_rw_N5.fit',1)
  rcompare = mrdfits('/sdss5/data0/lensout/stripe10/wthetarandlumw_stripe10_lumlensum_rw_N5.fit',1)

  sphoto_match, lensum, compare, ml, mc
  help,lensum,compare,ml,mc
  addstr='compare'
  sub_sample_wtheta, ff, rff, addstr=addstr, outdir=outdir, $
    uselens = lensum, userand = randsum, indices=ml,$
    overwrite=overwrite


END   
