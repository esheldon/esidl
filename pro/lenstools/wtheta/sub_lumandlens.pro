PRO sub_lumandlens, stripe, lumclr, mclr, lumext=lumext, mext=mext, radec=radec

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: sub_lumandlens, stripe, lumclr, mclr'
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

  IF n_elements(lumext) EQ 0 THEN lumext = 'N1.fit'
  IF n_elements(mext) EQ 0 THEN mext = 'N1.fit'

  runstr = 'stripe'+ntostr(long(stripe))+'_'
;  lumindir = '/sdss6/data0/wtheta/'
;  mindir=lumindir
  mindir = '/sdss5/data0/lensout/stripe'+ntostr(long(stripe))+'/'
  lumindir=mindir

  IF keyword_set(logbin) THEN lgstr = 'lg_' ELSE lgstr='_'

  IF keyword_set(radec) THEN BEGIN 
      front = 'wthetalumweq_'
      rfront = 'wthetarandlumweq_'
  ENDIF ELSE BEGIN 
      front = 'wthetalumw_'
      rfront = 'wthetarandlumw_'
  ENDELSE 

  lff = lumindir + front+runstr+'lumlensum_'+!colors[lumclr]+'w_'+lumext
  lrff = lumindir + rfront+runstr+'lumlensum_'+!colors[lumclr]+'w_'+lumext

  mff = mindir + 'main_zgal_gal_'+runstr+!colors[mclr]+'_lensum_'+mext
  mrff = mindir + 'main_zrand_'+runstr+!colors[mclr]+'_lensum_'+mext

  print,'lum files: ',lff,lrff
  llensum=mrdfits(lff,1)

  print,'lensing files: ',mff,mrff
  mlensum=mrdfits(mff,1)

  sphoto_match, llensum, mlensum, ml, mm
  help,llensum,mlensum,ml,mm

  delvarx,llensum,mlensum

  ;; use overwrite, because matchN6 will be screwed up by newname()
  addstr='match'+!colors[mclr]+repstr(mext,'.fit','')
  sub_sample_wtheta, lff, lrff, addstr=addstr, outdir=outdir, $
    indices=ml,/overwrite

  addstr='match'+!colors[lumclr]+repstr(lumext,'.fit','')
  zsub_sample, mff, mrff, addstr=addstr, outdir=outdir, $
    indices=mm,/overwrite


END 
