PRO make_lensout_newweight, noweight=noweight

  ;; make new lens outputs using only the 1/Sigma_crit^2 weighting
  ;; or no weighting

  stripestr = ['stripe10','stripe82','stripe36','stripe37','stripe42','stripe43']
;  stripestr = ['stripe82','stripe36','stripe37','stripe42','stripe43']
  nstripe = n_elements(stripestr)
  nclr=n_elements(!colors)

  basedir = '/sdss5/data0/lensout/'+stripestr+'/sublum/'
  type = 'lum'

  IF keyword_set(noweight) THEN addstr = '_noweight' ELSE addstr='_sigonly'

  FOR st=0L, nstripe-1 DO BEGIN 

      print,basedir[st]

      ;; loop over bandpasses luminosities were measured in
      FOR clr=0L, nclr-1 DO BEGIN 

          indir = basedir[st] + !colors[clr]+'/'
          print,indir

          ;; loop over bandpasses of lensing measurements
          FOR lclr=1,3 DO BEGIN 
          
              file = indir+type+'_zgal_gal_'+stripestr[st]+'_'+!colors[lclr]+'_N1.fit'
              rfile = indir+type+'_zrand_'+stripestr[st]+'_'+!colors[lclr]+'_N1.fit'
              sumfile=indir+type+'_zgal_gal_'+stripestr[st]+'_'+!colors[lclr]+'_sum_N1.fit'
              rsumfile=indir+type+'_zrand_'+stripestr[st]+'_'+!colors[lclr]+'_sum_N1.fit'
              zfile = indir+type+'_zgal_gal_'+stripestr[st]+'_'+!colors[lclr]+'_z_N1.fit'
              rzfile = indir+type+'_zrand_'+stripestr[st]+'_'+!colors[lclr]+'_z_N1.fit'
              lensumfile=indir+type+'_zgal_gal_'+stripestr[st]+'_'+!colors[lclr]+'_lensum_N1.fit'
              rlensumfile=indir+type+'_zrand_'+stripestr[st]+'_'+!colors[lclr]+'_lensum_N1.fit'

              print,file,rfile
              
              sh=mrdfits(file, 1, shhdr,/silent)
              rsh=mrdfits(rfile, 1, rshhdr,/silent)
              sum=mrdfits(sumfile, 1, sumhdr,/silent)
              rsum=mrdfits(rsumfile, 1, rsumhdr,/silent)
              zst=mrdfits(zfile, 1, zhdr,/silent)
              rzst=mrdfits(rzfile, 1, rzhdr,/silent)
              
              lensum=mrdfits(lensumfile, 1, lensumhdr,/silent)
              rlensum=mrdfits(rlensumfile, 1, rlensumhdr,/silent)
              
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; produce new weighting
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              
              nbin=n_elements(lensum[0].rsum)
              FOR i=0L, nbin-1 DO BEGIN 

                  ;; lenses
                  ww=where(lensum.npair[i] GT 0,nww)
                  IF nww NE 0 THEN BEGIN 
                      
                      print,'.',format='(a,$)'
                      IF keyword_set(noweight) THEN $
                        wdiv = 1./lensum[ww].npair[i]/lensum[ww].scritinv^2 $
                      ELSE wdiv = 1./lensum[ww].npair[i]
                      
                      lensum[ww].wsum[i] = lensum[ww].wsum[i]*wdiv
                      
                      lensum[ww].tansigsum[i] = lensum[ww].tansigsum[i]*wdiv
                      lensum[ww].tansigerrsum[i] = lensum[ww].tansigerrsum[i]*wdiv^2
                      lensum[ww].radsigsum[i] = lensum[ww].radsigsum[i]*wdiv
                      lensum[ww].radsigerrsum[i] = lensum[ww].radsigerrsum[i]*wdiv^2
                      
                      lensum[ww].etansum[i] = lensum[ww].etansum[i]*wdiv
                      lensum[ww].etanerrsum[i] = lensum[ww].etanerrsum[i]*wdiv^2
                      lensum[ww].eradsum[i] = lensum[ww].eradsum[i]*wdiv
                      lensum[ww].eraderrsum[i] = lensum[ww].eraderrsum[i]*wdiv^2
                  ENDIF 
                  
                  ;; random points
                  ww=where(rlensum.npair[i] GT 0,nww)
                  IF nww NE 0 THEN BEGIN 

                      print,'.',format='(a,$)'
                      IF keyword_set(noweight) THEN $
                        wdiv = 1./rlensum[ww].npair[i]/rlensum[ww].scritinv^2 $
                      ELSE wdiv = 1./rlensum[ww].npair[i]
                      
                      rlensum[ww].wsum[i] = rlensum[ww].wsum[i]*wdiv
                      
                      rlensum[ww].tansigsum[i] = rlensum[ww].tansigsum[i]*wdiv
                      rlensum[ww].tansigerrsum[i] = rlensum[ww].tansigerrsum[i]*wdiv^2
                      rlensum[ww].radsigsum[i] = rlensum[ww].radsigsum[i]*wdiv
                      rlensum[ww].radsigerrsum[i] = rlensum[ww].radsigerrsum[i]*wdiv^2
                      
                      rlensum[ww].etansum[i] = rlensum[ww].etansum[i]*wdiv
                      rlensum[ww].etanerrsum[i] = rlensum[ww].etanerrsum[i]*wdiv^2
                      rlensum[ww].eradsum[i] = rlensum[ww].eradsum[i]*wdiv
                      rlensum[ww].eraderrsum[i] = rlensum[ww].eraderrsum[i]*wdiv^2
                  ENDIF 
                  
              ENDFOR 
              combine_zlensum, lensum, sh.binsize, sh.rmin, sh.rmax, sh.h, newsum, newsh
              combine_zlensum, rlensum, rsh.binsize, rsh.rmin, rsh.rmax, rsh.h, newrsum, newrsh
              
              FXHCLEAN, shhdr & FXHCLEAN, rshhdr
              FXHCLEAN, sumhdr & FXHCLEAN, rsumhdr
              FXHCLEAN, zhdr & FXHCLEAN, rzhdr

              FXHCLEAN, lensumhdr & FXHCLEAN, rlensumhdr
              
              newfile = indir+type+'_zgal_gal'+addstr+'_'+stripestr[st]+'_'+!colors[lclr]+'_N1.fit'
              newrfile = indir+type+'_zrand'+addstr+'_'+stripestr[st]+'_'+!colors[lclr]+'_N1.fit'
              newsumfile=indir+type+'_zgal_gal'+addstr+'_'+stripestr[st]+'_'+!colors[lclr]+'_sum_N1.fit'
              newrsumfile=indir+type+'_zrand'+addstr+'_'+stripestr[st]+'_'+!colors[lclr]+'_sum_N1.fit'
              newzfile = indir+type+'_zgal_gal'+addstr+'_'+stripestr[st]+'_'+!colors[lclr]+'_z_N1.fit'
              newrzfile = indir+type+'_zrand'+addstr+'_'+stripestr[st]+'_'+!colors[lclr]+'_z_N1.fit'
              newlensumfile=indir+type+'_zgal_gal'+addstr+'_'+stripestr[st]+'_'+!colors[lclr]+'_lensum_N1.fit'
              newrlensumfile=indir+type+'_zrand'+addstr+'_'+stripestr[st]+'_'+!colors[lclr]+'_lensum_N1.fit'
              print
              print,newfile,newrfile
              print
;aploterror,!gratio,newsh.meanr,newsh.sigma,newsh.sigmaerr,psym=1
;help,newsum,/str

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; output new files
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              
              mwrfits2, newsh, newfile, shhdr, /create, /destroy
              mwrfits2, newrsh, newrfile, rshhdr, /create, /destroy
              mwrfits2, newsum, newsumfile, sumhdr, /create, /destroy
              mwrfits2, newrsum, newrsumfile, rsumhdr, /create, /destroy
              mwrfits2, zst, newzfile, zhdr, /create, /destroy
              mwrfits2, rzst, newrzfile, rzhdr, /create, /destroy
              
              mwrfits2, lensum, newlensumfile, lensumhdr, /create, /destroy
              mwrfits2, rlensum, newrlensumfile, rlensumhdr, /create, /destroy


          ENDFOR 
      ENDFOR 
  ENDFOR 

END 
