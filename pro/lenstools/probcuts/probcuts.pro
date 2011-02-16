PRO probcuts, struct, rclr, $
              meanmag, $
              probcut, nobj, $
              s2n_probcut, s2n_nobj, s2n_purity, $
              hirata=hirata

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: probcuts, struct, rclr, $'
      print,'         meanmag, $'
      print,'         probcut, nobj, $'
      print,'         s2n_probcut, s2n_nobj, s2n_purity, $'
      print,'         hirata=hirata'
      return
  ENDIF 

  ;; test how to cut based on probability.  Use S/N: must be within
  ;; 0.99 of the maximum S/N based on the forulas:
  ;;
  ;; S = S_0 * (Ngal/Ntot)  Ngal = total(probgal)
  ;; N = sigma/sqrt(Ntot)

  ;; struct needs tags tags = ['m_r','m_r_h','objc_prob_psf','petrocounts']


  pcolor = !green
  pcolor2 = !cyan
  IF !d.name EQ 'PS' THEN BEGIN 
      pcolor=!blue
      pcolor2 = !magenta
  ENDIF 

  IF keyword_set(hirata) THEN rsmear = struct.m_r_h[rclr] $
  ELSE rsmear = struct.m_r[rclr]

  IF keyword_set(hirata) THEN rsmear = struct.m_r_avg_h $
  ELSE rsmear = struct.m_r_avg

  ;; make a hard cut on smear polarizability first
  hardrcut = 0.8

  ;; test a simple hard cut on probability
  hardprobcut = 0.8

  ;; test S/N > hards2ncut*max(S/N) and S > hardscut*S0
  hards2ncut = 0.99
  hardscut = 0.99

  w=where(rsmear GT 0.0 AND rsmear LT hardrcut AND $
          struct.objc_prob_psf GE 0.0)

  pclr = 2
  binsize = 0.25
  hist = histogram(struct[w].petrocounts[pclr], bin=binsize, $
                   min=18.0, max=22.0,$
                   rev=rev_ind)

  nhist = n_elements(hist)

  s2n_probcut = fltarr(nhist)
  s2npure_probcut = fltarr(nhist)

  s2n_nobj = lonarr(nhist)
  s2npure_Nobj = lonarr(nhist)

  s2n_purity = fltarr(nhist)
  s2npure_purity = fltarr(nhist)

  meanmag = fltarr(nhist)

  hardcut_purity = fltarr(nhist)
  hardcut_nobj = lonarr(nhist)

  !p.multi=[0,0,2]

  FOR i=0L, nhist-1 DO BEGIN 

      IF rev_ind[i] NE rev_ind[i+1] THEN BEGIN 

          w2=rev_ind[ rev_ind[i]:rev_ind[i+1]-1 ]
          w2 = w[w2]
          meanmag[i] = mean_check( struct[w2].petrocounts[pclr] )

          probgal = 1. - struct[w2].objc_prob_psf

          ;; reverse sorted
          s = reverse( sort( probgal ) )

          nprob = n_elements(probgal)

          cumprob = total(probgal[s], /cumulative)
          num = findgen(nprob) + 1.

          SoverS0 = cumprob/num

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; an estimate of the S/N propto Ngal/Ntot/sqrt(Ntot)
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
          s2n = cumprob/sqrt(num)
          s2n = s2n/max(s2n)

          nplot = 1000

          xx = arrscl(findgen(nplot), 0.0, 1.0)
          yy = interpol(s2n, probgal[s], xx)
          yy2 = interpol(SoverS0, probgal[s], xx)

          plot, xx, yy, $
                xtitle='Prob Gal', ytitle='S/N,  S', $
                yrange=[0,1.2],/ystyle,$
                title='petrocounts['+!colors[pclr]+'] = '+ntostr(meanmag[i],5)
          oplot, xx, yy2, color=!blue

          ws = where(s2n GT hards2ncut, nws)
          ws2 = where(s2n GT hards2ncut AND $
                      SoverS0 GT hardscut, nws2)

          IF nws NE 0 THEN BEGIN 
              ws = (max(ws))[0]
              s2n_probcut[i] = probgal[s[ws]]

              s2n_nobj[i] = num[ws]
              s2n_purity[i] = SoverS0[ws]
          ENDIF 

          IF nws2 NE 0 THEN BEGIN 
              ws2 = (max(ws2))[0]
              s2npure_probcut[i] = probgal[s[ws2]]

              s2npure_nobj[i] = num[ws2]
              s2npure_purity[i] = SoverS0[ws2]
          ENDIF 


          oplot, [s2n_probcut[i], s2n_probcut[i]], $
                 [0, 1000], color=!red

          oplot, [s2npure_probcut[i], s2npure_probcut[i]], $
                 [0, 1000], color=!green

          print,s2n_probcut[i],s2npure_probcut[i]

;          legend, ['petrocounts['+!colors[pclr]+'] = '+ntostr(meanmag[i]),$
;                   'S/N Cut: '+ntostr(s2n_probcut[i]),$
;                   'S/N+Purity Cut: '+ntostr(s2npure_probcut[i])],charsize=1,$
;                  /bottom,/right

          legend,['S/N','S'],line=[0,0],color=[!p.color,!blue],$
                 /bottom,/right,charsize=1,thick=[!p.thick,!p.thick],/clear

          ;;;;;;;;;;;;;;;;;;;;;
          ;; hard cut
          ;;;;;;;;;;;;;;;;;;;;;

          wh = where( probgal GT hardprobcut, nwh)
          hardcut_purity[i] = total(probgal[wh])/nwh
          hardcut_nobj[i] = nwh

          IF !d.name EQ 'X' THEN key=get_kbrd(1)

      ENDIF 

  ENDFOR 

  !p.multi=[0,0,2]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plot the probability cut
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  hsnstr = ntostr(hards2ncut, 4)
  hsstr = ntostr(hardscut, 4)
  pcutstr = ntostr(hardprobcut, 4)

  lsnstr = 'S/N > '+hsnstr
  lsstr = 'Purity > '+hsstr
  lsn_lsstr = lsnstr + ' '+lsstr
  lpstr = 'Prob > '+pcutstr

  xtit = 'petrocunts['+!colors[pclr]+']'
  ww=where(meanmag NE 0)

  plot, meanmag[ww], s2n_probcut[ww], psym=8, $
        xtitle=xtit,$
        ytitle='Prob Gal Cut', yrange=[-0.2, 1.0], /ysty

  oplot, meanmag[ww], s2n_probcut[ww]

  oplot, meanmag[ww], s2npure_probcut[ww], psym=8, color=pcolor
  legend,[lsnstr,lsn_lsstr], $
         psym=[8,8], color=[!p.color,pcolor],$
         charsize=1

  oplot, meanmag[ww], s2npure_probcut[ww], color=pcolor

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plot the purity
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  minp = min([s2n_purity[ww], hardcut_purity[ww]], max=maxp)
  yrange = [0.95*minp, 1.05*maxp]
  plot,meanmag[ww], hardcut_purity[ww], psym=8, $
       xtitle=xtit,$
       ytitle='Purity', yrange=yrange, /ysty

  oplot, meanmag[ww], hardcut_purity[ww]

  oplot, meanmag[ww], s2n_purity[ww], psym=8, color=pcolor
  oplot, meanmag[ww], s2n_purity[ww], color=pcolor

  oplot, meanmag[ww], s2npure_purity[ww], psym=8, color=pcolor2
  oplot, meanmag[ww], s2npure_purity[ww], color=pcolor2

  legend,[lpstr,lsnstr,lsn_lsstr], psym=[8,8,8], $
         color=[!p.color,pcolor,pcolor2], /right,$
         charsize=1

  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  fac = 1.0
  tmax = max( [hardcut_nobj[ww], s2n_nobj[ww]] )
  yrange = [0.0, 1.1*tmax]/fac
  plot,  meanmag[ww], hardcut_nobj[ww]/fac, psym=8, xtitle=xtit, $
         ytit='# of objects', yrange=yrange
  oplot, meanmag[ww], hardcut_nobj[ww]/fac

  oplot, meanmag[ww], s2n_nobj[ww]/fac, psym=8, color=pcolor
  oplot, meanmag[ww], s2n_nobj[ww]/fac, color=pcolor

  oplot, meanmag[ww], s2npure_nobj[ww]/fac, psym=8, color=pcolor2
  oplot, meanmag[ww], s2npure_nobj[ww]/fac, color=pcolor2

  legend,[lpstr,lsnstr,lsn_lsstr], psym=[8,8,8], $
         color=[!p.color,pcolor,pcolor2],$
         charsize=1

  !p.multi=0



END 
