PRO make_pofsigcritinv, meanzs, sigzs, zlens, npts, sigcinv, sigcinvdist, WWi,$
                        median_siginv, mode_siginv, $
                        randzs_normal = randzs_normal

  !p.multi=[0,0,2]
  nrand = 500000L
  IF n_elements(randzs_normal) EQ 0 THEN BEGIN 
      print,'Creating normal'
      randzs_normal = randomu(seed, nrand, /normal)
      randzs = randzs_normal*sigzs + meanzs
  ENDIF ELSE BEGIN 
      randzs = randzs_normal*sigzs + meanzs
  ENDELSE 
  siginv = sigmacritinv(zlens, randzs)
  
  median_siginv = median(siginv)

  sigma_clip, siginv, meansc, sigsc, $
              nsig=4, niter=10, /silent
  ninsig = 20
  sbin = 2.*sigsc/ninsig
  plothist, siginv, sigcxhist, sigcyhist, bin=sbin, xcrange=xcrange
  
  ;; normalize
  minsc = min(sigcxhist)
  maxsc = max(sigcxhist)

  gauleg, minsc, maxsc, npts, sigcinv, WWi
  sigcinvdist = interpol(sigcyhist, sigcxhist, sigcinv)

  norm = total(sigcinvdist*WWi)
  sigcinvdist = sigcinvdist/norm

  plot,sigcinv,sigcinvdist,xrange=!x.crange,xstyle=1,psym=10

  w=where(sigcinvdist EQ max(sigcinvdist))
  mode_siginv = sigcinv[w[0]]
  !p.multi=0


END 

FUNCTION denscont_probability_sig, etan, ewidth, siginv, denscont

  maxarg = 15.0

  ;; prob of this ellipticity/2.: shear = e/2.
  efac = 1./sqrt(2.*!pi)/ewidth

  ;; uses h=1.0, pc^2/Msolar
  model = 2.*denscont*siginv

  earg = (etan - model)^2/2./ewidth^2 < maxarg
  pofe = efac*exp(-earg)

  return, pofe

END 

PRO denscont_likelihood_sig, zlens, meanzs, sigzs, $
                             etan, etanerr, shapenoise, npts, $
                             denscont, denscont_likelihood, $
                             silent=silent,$
                             randzs_normal = randzs_normal, $
                             med=med, mode=mode

  ndenscont = n_elements(denscont)
  denscont_likelihood = dblarr(ndenscont)
  ewidth = sqrt( etanerr^2 + shapenoise^2 )

  IF keyword_set(med) THEN BEGIN 

      sigcinv = sigmacritinv(zlens, meanzs)
      
      denscont_likelihood = denscont_probability_sig(etan, ewidth,$
                                                     sigcinv,$
                                                     denscont)

  ENDIF ELSE IF keyword_set(mode) THEN BEGIN 
      IF NOT keyword_set(silent) THEN BEGIN
          print,'Calculating weights, abscissa, and p(sigcinv)'
      ENDIF 
      make_pofsigcritinv, meanzs, sigzs, zlens, npts, $
                          sigcinv, sigcinvdist, WWi,$
                          median_siginv, mode_siginv, $
                          randzs_normal = randzs_normal

      denscont_likelihood = denscont_probability_sig(etan, ewidth,$
                                                     mode_siginv,$
                                                     denscont)

  ENDIF ELSE BEGIN 
      IF NOT keyword_set(silent) THEN BEGIN
          print,'Calculating weights, abscissa, and p(sigcinv)'
      ENDIF 
      make_pofsigcritinv, meanzs, sigzs, zlens, npts, $
                          sigcinv, sigcinvdist, WWi,median_siginv, $
                          randzs_normal = randzs_normal
      

            
      FOR i=0L, ndenscont-1 DO BEGIN 
          
          funcvals = sigcinvdist*denscont_probability_sig(etan, ewidth, $
                                                          sigcinv,$
                                                          denscont[i])
          
          denscont_likelihood[i] = total(funcvals*WWi)
                    
      ENDFOR 
  ENDELSE 

  denscont_likelihood = denscont_likelihood/max(denscont_likelihood)


END 
