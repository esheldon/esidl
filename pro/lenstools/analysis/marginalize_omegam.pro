PRO marginalize_omegam, gam, r0, likelihood, nocalc=nocalc, domean=domean

  ;; generate chi^2 contours for gamma, r0, marginalizing over
  ;; a gaussian in omega matter

  dir = '~/lensout/combstripe/comb/'
;  file = dir + 'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_N1.fit'
  file = dir + 'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_N1.fit'

  t = mrdfits(file,1)

  gamrange = [1.5, 2.2]
  r0range = [2.0, 10.0]

  r0range = [0.0, 12.0]

  ngam = 200
  nr0 = 200

  rprange = [0.01, 10.0]
  nrp = 200

  omegam_range = !omegam + 4.0*[-!omegamerr, !omegamerr]

  nomegam = 50
  omegam_vals = arrscl( dindgen(nomegam), omegam_range[0], omegam_range[1])

  omegam_like = gaussprob( omegam_vals, !omegam, !omegamerr, /double )

  plot, omegam_vals, omegam_like, psym=8

  likelihood_omega = dblarr(nomegam, ngam, nr0)

  IF NOT keyword_set(nocalc) THEN BEGIN 
      likelihood = dblarr(ngam, nr0)
      
      wpmodel, gamrange, ngam, r0range, nr0, rprange, nrp, !omegam, $
               gam, r0, rp, wp
      
      FOR im = 0L, nomegam-1 DO BEGIN 
          
          omegam = omegam_vals[im]
          
          wpsend = wp*(omegam/!omegam)
          print,'Omegam = '+ntostr(omegam)
          
          chisq_conf,t.meanr/1000.,t.sigma,t.covariance,rp,wpsend,gam,r0,$
                     ch,bgam,br0,$
                     errlow1=gerrlow,errhigh1=gerrhigh,$
                     errlow2=r0errlow,errhigh2=r0errhigh,yfit=yfit, $
                     xtit=!csym.gamma,$
                     ytit='r!D0!N [h!U'+!csym.minus+'1!N Mpc]',$
                     /dolegend,names=[!csym.gamma,'r!D0!N'],nkeep=[4,4], $
                     likelihood = tlike, /noprompt
          
          likelihood_omega[im, *, *] = tlike
          
      ENDFOR 

      ;; integrate over omegam likelihood
      x1 = omegam_vals[0]
      x2 = omegam_vals[nomegam-1]
      
      npts = nomegam*2L
      gauleg, x1, x2, npts, XXi, WWi
      
      FOR ig=0L, ngam-1 DO BEGIN 
          FOR ir0=0L, nr0-1 DO BEGIN 
              ;; the likelihood
              funcvals = $
                interpol( likelihood_omega[*, ig, ir0], omegam_vals, XXi )

              ;; the prior on omega
              tomegam_like = interpol( omegam_like, omegam_vals, XXi )

              ;; multiply the likelihood by the prior
              funcvals = funcvals*tomegam_like
              
              ;; integrate over omega matter
              likelihood[ig, ir0] = total( funcvals*WWi )
          ENDFOR 
      ENDFOR 
      likelihood = likelihood/max(likelihood)

  ENDIF 
  
  
  chisquared = -2.0*alog(likelihood)
  minchisq = min(chisquared)
  chisq_diff = chisquared - minchisq
  
  acontour, 1.0, chisq_diff, gam, r0, levels=!siglevels2, /center, $
            xtitle=!csym.gamma, ytitle='r!D0!N [h!U'+!csym.minus+'1!N Mpc]'
  
  index = lindgen(long(ngam)*long(nr0))
  x = index MOD ngam
  y = index/nr0
  w=where(chisquared EQ minchisq, nw)
  bestgam = (gam[x[w]])[0]
  bestr0 = (r0[y[w]])[0]

  oplot,[bestgam],[bestr0], psym=7

  print,'Hit a key'
  key=get_kbrd(1)

  marginalize_like2d, gam, r0, likelihood, gam_like, r0_like
  ;; find 1,2,3-sigma regions
  
  get_like_conf,gam, gam_like, $
                bestgam, gamlow, gamhigh, $
                gamerrlow, gamerrhigh, domean=domean
  get_like_conf, r0, r0_like, $
                 bestr0, r0low, r0high, $
                 r0errlow, r0errhigh, domean=domean

  oplot,[bestgam],[bestr0], psym=7,color=!green


END 
