PRO make_wthetaout_newweight, noweight=noweight

  ;; make new lens outputs using only the 1/Sigma_crit^2 weighting
  ;; or no weighting

;  stripes = [10,82,36,37,42,43]
  stripes = [10]
;  stripestr = ['stripe10','stripe82','stripe36','stripe37','stripe42','stripe43']
  stripestr = ['stripe10']
  nstripe = n_elements(stripestr)
  nclr=n_elements(!colors)

  basedir = '/sdss5/data0/lensout/'+stripestr+'/'

  IF keyword_set(noweight) THEN addstr = '_noweight' ELSE addstr='_sigonly'

  FOR st=0L, nstripe-1 DO BEGIN 

      print,basedir[st]

      ;; loop over bandpasses luminosities were measured in
      FOR clr=0L, nclr-1 DO BEGIN 

          indir = basedir[st]
          print,indir
          
          file = indir+'wthetalumw_'+stripestr[st]+'_'+!colors[clr]+'w_N2.fit'
          rfile = indir+'wthetarandlumw_'+stripestr[st]+'_'+!colors[clr]+'w_N2.fit'
          sumfile=indir+'wthetalumw_'+stripestr[st]+'_sum_'+!colors[clr]+'w_N2.fit'
          rsumfile=indir+'wthetarandlumw_'+stripestr[st]+'_sum_'+!colors[clr]+'w_N2.fit'
          zfile = indir+'wthetalumw_'+stripestr[st]+'_z_'+!colors[clr]+'w_N2.fit'
          rzfile = indir+'wthetarandlumw_'+stripestr[st]+'_z_'+!colors[clr]+'w_N2.fit'
          lumlensumfile=indir+'wthetalumw_'+stripestr[st]+'_lumlensum_'+!colors[clr]+'w_N2.fit'
          rlumlensumfile=indir+'wthetarandlumw_'+stripestr[st]+'_lumlensum_'+!colors[clr]+'w_N2.fit'
          lensumfile=indir+'wthetalumw_'+stripestr[st]+'_lensum_'+!colors[clr]+'w_N2.fit'
          rlensumfile=indir+'wthetarandlumw_'+stripestr[st]+'_lensum_'+!colors[clr]+'w_N2.fit'
;(
          print,file,rfile
          sh=mrdfits(file, 1, shhdr,/silent)
          rsh=mrdfits(rfile, 1, rshhdr,/silent)
          sum=mrdfits(sumfile, 1, sumhdr,/silent)
          rsum=mrdfits(rsumfile, 1, rsumhdr,/silent)
          zst=mrdfits(zfile, 1, zhdr,/silent)
          rzst=mrdfits(rzfile, 1, rzhdr,/silent)

          lumlensum=mrdfits(lumlensumfile, 1, lumlensumhdr, /silent)
          rlumlensum=mrdfits(rlumlensumfile, 1, rlumlensumhdr, /silent)

          print,lensumfile
          lensum=mrdfits(lensumfile, 1, lensumhdr,/silent)
          print,rlensumfile
          rlensum=mrdfits(rlensumfile, 1, rlensumhdr,/silent)

          nbeta = n_elements(lensum[0].npsum[0,*])
 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; produce new weighting
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              
          ;; we know these all passed check
          scritinv=sdss_sigma_crit(stripes[st], 2, lensum.z1d, /use_lambda, /silent)
          DL = angdist_lambda(lensum.z1d)*1000.

          rscritinv=sdss_sigma_crit(stripes[st], 2, rlensum.z, /use_lambda, /silent)
          rDL = angdist_lambda(rlensum.z)*1000.

          nbin=n_elements(lensum[0].rsum)
          
          IF keyword_set(noweight) THEN BEGIN 
              wdiv = DL^2/scritinv^2/1.e19
              rwdiv = rDL^2/rscritinv^2/1.e19
          ENDIF ELSE BEGIN 
              wdiv = DL^2/2.e11
              rwdiv = rDL^2/2.e11
          ENDELSE 
          FOR i=0L, nbin-1 DO BEGIN 

              ;; lenses
                  
              print,'.',format='(a,$)'
              
              lensum.wsum[i] = lensum.wsum[i]*wdiv
              
              lensum.rsum[i] = lensum.rsum[i]*wdiv
              FOR ib=0L, nbeta-1 DO BEGIN 
                  lensum.npsum[i,ib] = lensum.npsum[i,ib]*wdiv
              ENDFOR 
              lensum.lsum[i] = lensum.lsum[i]*wdiv
              lensum.lwsum[i] = lensum.lwsum[i]*wdiv
              
              lumlensum.lsum[i] = lumlensum.lsum[i]*wdiv
              lumlensum.lwsum[i] = lumlensum.lwsum[i]*wdiv
              
              ;; random points
                  
              print,'.',format='(a,$)'
             
              rlensum.wsum[i] = rlensum.wsum[i]*rwdiv
              
              rlensum.rsum[i] = rlensum.rsum[i]*rwdiv
              FOR ib=0L, nbeta-1 DO BEGIN 
                  rlensum.npsum[i,ib] = rlensum.npsum[i,ib]*rwdiv
              ENDFOR 
              rlensum.lsum[i] = rlensum.lsum[i]*rwdiv
              rlensum.lwsum[i] = rlensum.lwsum[i]*rwdiv
              
              rlumlensum.lsum[i] = rlumlensum.lsum[i]*rwdiv
              rlumlensum.lwsum[i] = rlumlensum.lwsum[i]*rwdiv
                            
          ENDFOR 

          combine_wthetalumw_lensum, lensum, sh.binsize, sh.rmin, sh.rmax, sh.h, newsum, newsh
          combine_wthetalumw_lensum, rlensum, rsh.binsize, rsh.rmin, rsh.rmax, rsh.h, newrsum, newrsh

          FXHCLEAN, shhdr & FXHCLEAN, rshhdr
          FXHCLEAN, sumhdr & FXHCLEAN, rsumhdr
          FXHCLEAN, zhdr & FXHCLEAN, rzhdr
          FXHCLEAN, lumlensumhdr & FXHCLEAN, rlumlensumhdr

          FXHCLEAN, lensumhdr & FXHCLEAN, rlensumhdr
          
          newfile = indir+'wthetalumw'+addstr+'_'+stripestr[st]+'_'+!colors[clr]+'w_N2.fit'
          newrfile = indir+'wthetarandlumw'+addstr+'_'+stripestr[st]+'_'+!colors[clr]+'w_N2.fit'
          newsumfile=indir+'wthetalumw'+addstr+'_'+stripestr[st]+'_sum_'+!colors[clr]+'w_N2.fit'
          newrsumfile=indir+'wthetarandlumw'+addstr+'_'+stripestr[st]+'_sum_'+!colors[clr]+'w_N2.fit'
          newzfile = indir+'wthetalumw'+addstr+'_'+stripestr[st]+'_z_'+!colors[clr]+'w_N2.fit'
          newrzfile = indir+'wthetarandlumw'+addstr+'_'+stripestr[st]+'_z_'+!colors[clr]+'w_N2.fit'
          newlumlensumfile=indir+'wthetalumw'+addstr+'_'+stripestr[st]+'_lumlensum_'+!colors[clr]+'w_N2.fit'
          newrlumlensumfile=indir+'wthetarandlumw'+addstr+'_'+stripestr[st]+'_lumlensum_'+!colors[clr]+'w_N2.fit'
          newlensumfile=indir+'wthetalumw'+addstr+'_'+stripestr[st]+'_lensum_'+!colors[clr]+'w_N2.fit'
          newrlensumfile=indir+'wthetarandlumw'+addstr+'_'+stripestr[st]+'_lensum_'+!colors[clr]+'w_N2.fit'
          print
          print,newfile,newrfile
          print
;stop
;(
;aploterror,!gratio,newsh.meanr,newsh.sigma,newsh.sigmaerr,psym=1
;help,newsum,/str
          
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; output new files
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          mwrfits2, lensum, newlensumfile, lensumhdr, /create, /destroy
          mwrfits2, rlensum, newrlensumfile, rlensumhdr, /create, /destroy

          mwrfits2, lumlensum, newlumlensumfile, lumlensumhdr, /create, /destroy
          mwrfits2, rlumlensum, newrlumlensumfile, rlumlensumhdr, /create, /destroy

          mwrfits2, newsh, newfile, shhdr, /create, /destroy
          mwrfits2, newrsh, newrfile, rshhdr, /create, /destroy
          mwrfits2, newsum, newsumfile, sumhdr, /create, /destroy
          mwrfits2, newrsum, newrsumfile, rsumhdr, /create, /destroy
          mwrfits2, zst, newzfile, zhdr, /create, /destroy
          mwrfits2, rzst, newrzfile, rzhdr, /create, /destroy
          

      ENDFOR 
  ENDFOR 

END 
