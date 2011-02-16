PRO sub_type_lumwtheta, stripe, lumw=lumw, ext=ext, logbin=logbin,overwrite=overwrite
  
  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: sub_type_lumwtheta, stripe'
      return
  ENDIF 

  ;; N3 is 0.1 L* stuff
  ;; N2 is magcut stuff
  ;; N1 is 0.1 L* stuff to 2 Mpc

  colors=['u','g','r','i','z']
  msun = [6.39,5.07,4.62,4.52,4.48]
  mstar = [-18.34,-20.04,-20.83,-21.26,-21.55]
  lstar = 10.0^((mstar-msun)/(-2.5))

  IF n_elements(ext) EQ 0 THEN ext = 'N3.fit'

  runstr = 'stripe'+ntostr(long(stripe))+'_'

  indir = '/sdss5/data0/lensout/'+'stripe'+ntostr(long(stripe))+'/'

  IF keyword_set(logbin) THEN lgstr = 'lg_' ELSE lgstr='_'

  front = 'wtheta'+lgstr
  rfront = 'wthetarand'+lgstr

  IF NOT keyword_set(lumw) THEN BEGIN 

      ff = indir + front+runstr+'lensum_'+ext
      rff = indir + rfront+runstr+'lensum_'+ext
      
      lensum=mrdfits(ff,1)
      randsum=mrdfits(rff,1)
  ENDIF 

  maxz = '0.6'
  minz = '0.02'

  FOR iclr=0,4 DO BEGIN 

      outdir = indir + 'sublum/'+colors[iclr]+'/'

      ;; for sub-sampling the weighted stuff
      IF keyword_set(lumw) THEN BEGIN 
          setzero, lensum, rlensum
          astr=!colors[iclr]+'w_' 

          front = 'wthetalumw'+lgstr
          rfront = 'wthetarandlumw'+lgstr

          ff = indir + front+runstr+'lumlensum_'+astr+ext
          rff = indir + rfront+runstr+'lumlensum_'+astr+ext

          lensum=mrdfits(ff,1)
          randsum=mrdfits(rff,1)
      ENDIF 


      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; cut on z, spiral
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      wstr = '(lensum.z1d le '+maxz+') and (lensum.z1d gt '+minz+')'
      wstr = wstr + ' and ( (lensum.primtarget and 2L^6) ne 0)'
      wstr = wstr + ' and ( (lensum.classification and 2L^3) ne 0)'
      
      addstr = 'spiral'
      sub_sample_wtheta, ff, rff, wstr, addstr=addstr,uselens=lensum,$
        userand=randsum,lumw=lumw, outdir=outdir,$
        overwrite=overwrite

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; cut on z, ellip
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      wstr = '(lensum.z1d le '+maxz+') and (lensum.z1d gt '+minz+')'
      wstr = wstr + ' and ( (lensum.primtarget and 2L^6) ne 0)'
      wstr = wstr + ' and ( (lensum.classification and 2L^0) ne 0)'
      wstr = wstr + ' and ( (lensum.classification and 2L^8) eq 0)'
      
      addstr = 'ellip'
      sub_sample_wtheta, ff, rff, wstr, addstr=addstr,uselens=lensum,$
        userand=randsum,lumw=lumw, outdir=outdir,$
        overwrite=overwrite


  ENDFOR 


END 
