PRO stellar_locus_envelope, str

  ;; Plot the r vs mag, coloring the stars and galaxies and
  ;; drawing the envelope of the stellar locus

  IF n_elements(str) EQ 0 THEN read_tsobj,[3698,20,3],str,/corr,/all

  ;; magshift is just for comparison with deeper surveys, although
  ;; this overestimates the fraction of stars at the faint end

  magshift = 2.1
  petro = str.petrocounts[2] + magshift

  xrange = [14,22] + magshift
  
  minmag = 14.0 + magshift
  maxmag = 22.5 + magshift

  plot,petro,str.m_r_h[2],xrange=xrange,yrange=[0,2],psym=3, $
       xstyle=1, ystyle=1

  ws=where(str.objc_type eq 6 AND (str.m_r_h[2] GT 0) AND (str.m_r_h[2] LT 2))
  wg=where(str.objc_type EQ 3 AND (str.m_r_h[2] GT 0) AND (str.m_r_h[2] LT 2))
  oplot,petro[wg],str[wg].m_r_h[2],psym=8,symsize=0.5,color=!cyan
  oplot,petro[ws],str[ws].m_r_h[2],psym=8,symsize=0.5,color=!red

  
  oplot,xrange,[1,1]
  oplot,xrange,[0.8, 0.8], color=!green, thick=2


  magbin = 0.5
  h = histogram(petro[ws], bin=magbin, min=minmag, max=maxmag, $
                rev=rev)

  nh = n_elements(h)

  mnr = fltarr(nh)
  mnm = mnr
  rms = mnr
  envl = mnr
  envh = mnr
  cut = mnr

  FOR i=0L, nh-1 DO BEGIN 
      IF rev[i] NE rev[i+1] THEN BEGIN 

          w = rev[ rev[i]:rev[i+1] -1 ]
          w = ws[w]

          ;; calculate envelope 
          sigma_clip, str[w].m_r_h[2], tm, ts, nsig=3,  niter=4, /silent

          mnr[i] = tm
          rms[i] = ts

;          mom = moment(str[w].m_r_h[2])
;          mnr[i] = mom[0]
          mnm[i] = mean(petro[w])

;          rms[i] = sqrt(mom[1])

          envh[i] = 1 + 2.*rms[i]
          envl[i] = 1 - 2.*rms[i]

          cut[i] = min([0.8, envl[i]])

      ENDIF
      print,mnm[i],mnr[i],rms[i]
  ENDFOR 
  
  w=where(mnm NE 0)
  mnm = mnm[w]
  mnr = mnr[w]
  rms = rms[w]
  envl = envl[w]
  envh = envh[w]
  cut = cut[w]

  oplot, xrange, [mean(mnr), mean(mnr)], color=!blue
  oplot, mnm, envh, color=!blue, thick=2
  oplot, mnm, envl, color=!blue, thick=2
;  oplot, mnm, cut, color=!green, thick=2

return
END 
