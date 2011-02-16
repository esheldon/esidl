PRO fitgauss2im, image, cenx, ceny, sky, scale=scale, _extra=_extra

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: fitgauss2im, image, cenx, ceny [, sky, scale=scale, _extra=_extra]'
      return
  ENDIF 

  nterms = 4
  myusersym, 'fill_circle'
  psym=8

  IF n_elements(sky) EQ 0 THEN BEGIN 
      sky = median(image)
  ENDIF 

  IF n_elements(scale) EQ 0 THEN BEGIN 
      scale = 1.                
      addstr = ' pixels'
  ENDIF ELSE BEGIN 
      addstr = ' arcsec'
  ENDELSE 

  nsig = 5.

  sz=size(image)

  sx=sz[1]
  sy=sz[2]

  min = 0.
  max = min([ abs(sx-cenx), cenx])
  max = min([max, min([abs(sy-ceny), ceny]) ] )
  binsize = 1.                  ;pixels

  index=lindgen(sx*sy)
  x = index MOD sx
  y = index/sx

  r = sqrt( (x-float(cenx))^2 + (y-float(ceny))^2 )

  sigma_clip, image, mn, imsig, nsig=3., niter=4,/silent

  pmax = max(image); + 0.25*imsig
  w = where( image GT sky+nsig*imsig, nw)

  IF nw NE 0 THEN BEGIN 

      imhist = histogram(r[w], binsize=binsize, $
                         min = 0., max = max, reverse=rev_ind)

      nbin=n_elements(imhist)
      
      tmpr = fltarr(nbin)
      tmpvals = replicate(sky, nbin)
      tmpvdev = fltarr(nbin)

      FOR i=0, nbin-1 DO BEGIN 
          IF rev_ind(i) NE rev_ind(i+1) THEN BEGIN 
              
              wr = rev_ind( rev_ind[i]:rev_ind[i+1]-1 )
              
              ww = w[wr]
              nww = n_elements(ww)

              tmpvals[i] = mean_check(image[ww])
              tmpvdev[i] = sdev(image[ww])/sqrt(nww)
;              tmpvdev[i] = sqrt(tmpvals[i]-sky)
              tmpr[i] = mean_check(r[ww])

          ENDIF 
      ENDFOR 
  ENDIF ELSE BEGIN 
      print,'Cannot fit Gaussian'
      return
  ENDELSE 

  j=0
  go=1
  binvals = -1

  wgood = where(tmpvdev NE 0., ngood)

  IF ngood GT nterms THEN BEGIN 

      binvals = tmpvals[wgood]
      vdev = tmpvdev[wgood]
      rvals = tmpr[wgood]
  
      s=sort(rvals)
      binvals = binvals[s]
      rvals = rvals[s]
      vdev = vdev[s]

      newr = [-reverse(rvals), rvals]
      newbinvals = [reverse(binvals),binvals]
      newvdev = [reverse(vdev), vdev]

      gfit = gaussfit(newr, newbinvals, tAA, nterms=nterms)

      ;; get guesses from gaussfit, then get errors
      ;; from mpfitexpr

      AA = mpfitexpr('P(3) + P[0]*exp(- ( (X-P[1])/P[2] )^2/2. )', $
                     newr, newbinvals, newvdev, double(tAA), perror=AAerror)
      
      forprint,AA,AAerror
      sigma = AA[2]
      fwhm = 2.*sqrt(2.*alog(2))*sigma*scale
      fwhmerr = fwhm*(AAerror[2]/AA[2])

      maxr=max(rvals)
      xx = arrscl( findgen(100), -maxr, maxr )
      
      zz = ((xx-AA[1])/AA[2] )^2/2. < 10.
      gauss=AA[3] + AA[0]*exp(-zz )

      aploterror, 1, newr, newbinvals, newvdev,psym=psym, $
        xrange=[-maxr-1,maxr+1],xstyle=1, _extra=_extra
;      plot, newr, newbinvals, psym=psym,yrange=[sky,pmax]
      oplot, xx, gauss
      message = 'FWHM = '+ntostr( rnd(fwhm,2), 4)+$
        ' '+!tsym.plusminus+' '+ntostr( rnd(fwhmerr,2), 4)+addstr

      degfree = ngood-nterms+1  ;plus 1 because offset is not important

      ;;IF degfree GT 0. THEN BEGIN ;/2 because new has negative r-vals
      ;;    chisq = total( ( (gfit-newbinvals)/newvdev )^2 )/2.
      ;;    message = message + '  chisq/degfree = '+$
      ;;      ntostr(chisq)+'/'+ntostr(degfree)+' = '+ntostr(chisq/degfree)
      ;;ENDIF ELSE BEGIN
      ;;    message = message + '  no chisq'
      ;;ENDELSE 

      legend,message,box=0

  ENDIF ELSE BEGIN 
      print,'No good values found'
  ENDELSE 

  return
END 
