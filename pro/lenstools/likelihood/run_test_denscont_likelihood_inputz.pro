PRO run_test_denscont_likelihood_inputz, scat, useind, $
                                         lcat=lcat, $
                                         doplot=doplot,$
                                         meansigcritinv=meansigcritinv, $
                                         ntrial=ntrial, $
                                         fracatlens=fracatlens

  IF n_elements(scat) EQ 0 THEN BEGIN
      get_scat, 82, [2,3], scat
      photoz_dist_struct, scat, z, pofz, useind
  ENDIF 

  IF n_elements(ntrial) EQ 0 THEN ntrial=1
  nuse = n_elements(useind)

  FOR i=0L, ntrial-1 DO BEGIN 

      ;;  rand = randomu(seed, nuse)
      ;;  s = sort(rand)

      ;; randomize errors
      ;;  etanerr = sqrt(scat[useind[s]].e1e1err^2 + $
      ;;                 scat[useind[s]].e2e2err^2)
      
      etanerr = fltarr(nuse)
      
      rand = randomu(seed, nuse)
      s = sort(rand)
      
      ;;s=s[0:9999]
      ;;etanerr=etanerr[0:9999]
      
      IF n_elements(lcat) NE 0 THEN BEGIN 

          ;; generate a lens redshift for each source

          ;; get "main" galaxy sample
          w=where((lcat.primtarget AND 2L^6) NE 0, nlcat)

          ;; get nsource lenses randomly from this redshift distribution
          plothist,lcat[w].z,xhist,yhist,min=0.02,max=0.3,bin=0.02,/noplot
          genrand, yhist, xhist, nuse, zL

;          w2=where(lcat.m_e1[2] ne 0 AND lcat.m_e2[2] NE 0 AND $
;                  lcat.m_R[2] gt 0, nw2)

;          phase = arrscl( randomn(seed, nw2), 0.0, 2.*!pi, arrmin=0.,arrmax=1.)
;          lcat_etan = -(cos(2.*phase)*lcat[w2].m_e1[2] + $
;                        sin(2.*phase)*lcat[w2].m_e2[2])
;          e1e2err2 = sign(lcat[w2].m_e1e2err[2])*lcat[w2].m_e1e2err[2]^2
;          lcat_err = sqrt(lcat[w2].m_e1e1err[2]^2*cos(2.*phase)^2 + $
;                          lcat[w2].m_e2e2err[2]^2*sin(2.*phase)^2-$
;                          2.*e1e2err2*sin(2.*phase)*cos(2.*phase))

;          lrand = long( arrscl( randomu(seed, nuse), 0, nw2 ) )
;          input_evals = lcat_etan[lrand]
;          etanerr = lcat_err[lrand]

;          plothist, lcat_etan, exhist, eyhist, bin=0.03, min=-1.0, max=1.0,$
;                    /norm;,/noplot
;          genrand, eyhist, exhist, nuse, input_evals
;          plothist,input_evals,bin=0.03,/norm,/overplot,color=!green

      ENDIF ELSE BEGIN 
          zL = 0.15
      ENDELSE 

      photoz_z = scat[useind[s]].photoz_z
      photoz_zerr = scat[useind[s]].photoz_zerr

      phase = arrscl( randomn(seed, nuse), 0.0, 2.*!pi, arrmin=0.,arrmax=1.)
      erand1 = randomu(seed,nuse)
      es1 = sort(temporary(erand1))
      erand2 = randomu(seed,nuse)
      es2 = sort(temporary(erand2))

      input_evals = -(cos(2.*phase)*scat[useind[es1]].e1_recorr + $
                      sin(2.*phase)*scat[useind[es2]].e2_recorr)
      e1e2err2 = sign(scat[useind[es1]].e1e2err)*scat[useind[es1]].e1e2err^2
      etanerr = sqrt(scat[useind[es1]].e1e1err^2*cos(2.*phase)^2 + $
                     scat[useind[es2]].e2e2err^2*sin(2.*phase)^2-$
                     2.*e1e2err2*sin(2.*phase)*cos(2.*phase))

      ;; put some new sources at the lens redshift
      IF n_elements(fracatlens) NE 0 THEN BEGIN 

          IF n_elements(fracatlens) NE 0 AND n_elements(lcat) NE 0 THEN BEGIN 
              message,'Do not send both lcat and fracatlens. Not set up'+$
                      'to do that yet'
          ENDIF 

          N_atlens = long(nuse*fracatlens/(1.-fracatlens))
          print,'Adding '+ntostr(N_atlens)+' sources at zL = '+ntostr(zL)

          rand2 = randomu(seed, nuse)
          sf = sort(rand2)
          sf = sf[0:N_atlens-1]

          add_photoz_z = replicate(zL, N_atlens)
          add_photoz_zerr = photoz_zerr[sf]
          add_etanerr = etanerr[sf]

          photoz_z = [photoz_z, add_photoz_z]
          photoz_zerr = [photoz_zerr, add_photoz_zerr]
          etanerr  = [etanerr, add_etanerr]

          rand2 = randomu(seed, nuse+N_atlens)
          sf = sort(rand2)
          photoz_z = photoz_z[sf]
          photoz_zerr = photoz_zerr[sf]
          etanerr = etanerr[sf]

      ENDIF 

      test_denscont_likelihood_inputz, zL, $
                                       photoz_z, $
                                       photoz_zerr, $
                                       etanerr, doplot=doplot,$
                                       meansigcritinv=meansigcritinv,$
                                       fracatlens=fracatlens,input_evals=input_evals
  ENDFOR 

END 
