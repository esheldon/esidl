PRO sub_lumwtheta, stripe, lumw=lumw, ext=ext, logbin=logbin, overwrite=overwrite
  
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

  ;; skip to 2bin stuff
  GOTO,jump1

  ;;;;;; Bin by number
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

          print,'reading files: ',ff,rff

          lensum=mrdfits(ff,1)
          randsum=mrdfits(rff,1)
      ENDIF 

      ;;!!!! BINBYNUM HAS BUILT-IN Z CUTS!!!! 0.6,0.02
      binbynum, lensum, iclr, w1, w2, w3, w4
      addstr = 'lum1'
      sub_sample_wtheta, ff, rff, addstr=addstr, outdir=outdir, $
        uselens = lensum, userand = randsum, indices=w1, lumw=lumw,$
        overwrite=overwrite
      addstr = 'lum2'
      sub_sample_wtheta, ff, rff, addstr=addstr, outdir=outdir, $
        uselens = lensum, userand = randsum, indices=w2, lumw=lumw,$
        overwrite=overwrite
      addstr = 'lum3'
      sub_sample_wtheta, ff, rff, addstr=addstr, outdir=outdir, $
        uselens = lensum, userand = randsum, indices=w3, lumw=lumw,$
        overwrite=overwrite
      addstr = 'lum4'
      sub_sample_wtheta, ff, rff, addstr=addstr, outdir=outdir, $
        uselens = lensum, userand = randsum, indices=w4, lumw=lumw,$
        overwrite=overwrite

  ENDFOR 

jump1:
  ;;;;;; 2 Bins by number
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

          print,'reading files: ',ff,rff

          lensum=mrdfits(ff,1)
          randsum=mrdfits(rff,1)
      ENDIF 

      ;;!!!! BINBYNUM HAS BUILT-IN Z CUTS!!!! 0.6,0.02
      binbynum_2bin, lensum, iclr, w1, w2
      addstr = 'lum1twobin'
      sub_sample_wtheta, ff, rff, addstr=addstr, outdir=outdir, $
        uselens = lensum, userand = randsum, indices=w1, lumw=lumw,$
        overwrite=overwrite
      addstr = 'lum2twobin'
      sub_sample_wtheta, ff, rff, addstr=addstr, outdir=outdir, $
        uselens = lensum, userand = randsum, indices=w2, lumw=lumw,$
        overwrite=overwrite

  ENDFOR 
END 
