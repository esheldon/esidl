PRO sub_lowhighwtheta, stripe, massadd=massadd
  
  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: sub_lowhighwtheta, stripe'
      return
  ENDIF 

  msun = [6.39,5.07,4.62,4.52,4.48]
  mstar = [-18.34,-20.04,-20.83,-21.26,-21.55]
  lstar = 10.0^((mstar-msun)/(-2.5))

  ext = 'N1.fit'

  runstr = 'stripe'+ntostr(long(stripe))+'_'

  indir = '/sdss5/data0/lensout/'+'stripe'+ntostr(long(stripe))+'/'

  maxz = '0.6'
  minz = '0.02'

  ;; median voronoi_dens from stripes 10/82
  med_vor = '105.8'

  IF keyword_set(massadd) THEN mstr = 'massadd_' ELSE mstr=''

  ff = indir + 'wtheta_'+runstr+'lensum_'+mstr+ext
  rff = indir + 'wthetarand_'+runstr+'lensum_'+ext

  lensum=mrdfits(ff,1)
  randsum=mrdfits(rff,1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; cut on z, voronoi_dens defined
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wstr = '(lensum.z1d le '+maxz+') and (lensum.z1d gt '+minz+')'
  wstr = wstr + ' and ( lensum.voronoi_dens gt 0.0 )'

  addstr = 'tot'
  sub_sample_wtheta, ff, rff, wstr, addstr=addstr,uselens=lensum,userand=randsum

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; cut on z, high density
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wstr = '(lensum.z1d le '+maxz+') and (lensum.z1d gt '+minz+')'
  wstr = wstr + ' and ( lensum.voronoi_dens gt 0.0 )'
  wstr = wstr + ' and ( lensum.voronoi_dens ge '+med_vor+' )'

  addstr = 'high'
  sub_sample_wtheta, ff, rff, wstr, addstr=addstr,uselens=lensum,userand=randsum

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; cut on z, low density
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wstr = '(lensum.z1d le '+maxz+') and (lensum.z1d gt '+minz+')'
  wstr = wstr + ' and ( lensum.voronoi_dens gt 0.0 )'
  wstr = wstr + ' and ( lensum.voronoi_dens lt '+med_vor+' )'

  addstr = 'low'
  sub_sample_wtheta, ff, rff, wstr, addstr=addstr,uselens=lensum,userand=randsum


END 
