PRO test_kcorr_meangr, lcat, usemean=usemean

  default = -1000.

  cindex = 2
  IF n_elements(lcat) EQ 0 THEN get_spectra_lcat, 10, lcat

  gr=lcat.counts_model[1]-lcat.counts_model[2]

  w=where(gr LT 2 AND gr GT 0 AND $
          lcat.petrocounts[2] GT 0.0 AND $
          lcat.petrocounts[2] LT 17.5, nw)

  minz = 0.05
  maxz = 0.2
  nz = 30
  medz = arrscl( findgen(nz), minz, maxz )

  FOR i=0L, nz-1 DO BEGIN 
    
      
      wtheta_absmag, medz[i], cindex, lcat[w].petrocounts[2], gr[w], $
        absmag, lum_solar, kcorr=kcorr

      mink = min(kcorr) + 0.01
      maxk = max(kcorr) - 0.01

      w2 = where( ( kcorr GT mink ) AND $
                  ( kcorr LT maxk ) AND $
                  ( kcorr NE default), nw2)
      plothist, kcorr[w2], xhist, yhist, bin=0.005, max=maxk, min=mink,/noplot
      
      print,'Z = ',medz[i]

      print,nw2
      
      ;; find mean kcorr
      IF keyword_set(usemean) THEN meank = mean(kcorr[w2]) $
      ELSE meank = median(kcorr[w2])
      errk = sdev(kcorr[w2])/sqrt(nw2)

      ;;meank = mean(kcorr[w2])

      print,'Mean k = ',meank,' +/- ',errk

      ;; find mean g-r
      IF keyword_set(usemean) THEN meangr = mean(gr[w[w2]]) $
      ELSE meangr = median(gr[w[w2]])
      ;meangr = mean(gr[w[w2]])
      ;;meangr = total( kcorr[w2]*gr[w[w2]] )/total(kcorr[w2])
      print,'Mean g-r = ',meangr
      
      wtheta_absmag, medz[i], cindex, lcat[w[0]].petrocounts[2], meangr, $
        absmag, lum_solar, kcorr=kcorr2
      
      print,'kcorr(<g-r>) = ',kcorr2

      ;key=get_kbrd(1)

      add_arrval, meank, meankcorr
      add_arrval, errk, errkcorr
      add_arrval, meangr, meangmr
      add_arrval, kcorr2, kcorrmeangmr

  ENDFOR 

;  !p.multi=[0,0,2]
  aploterror, !gratio, medz, meankcorr, errkcorr, psym=1, yrange=[-0.3, 0],$
    xrange = [min(medz)-0.02, max(medz)+0.02], $
    xtitle='Z',ytitle='Kcorr['+!colorsp[cindex]+']'
  oplot, medz, kcorrmeangmr, psym=4, color=!red
  legend, ['<kcorr>','kcorr(<g-r>)'],$
    psym=[1,4],colors=[!p.color,!red],thick=[!p.thick,!p.thick],/right

;  diff = abs(meankcorr-kcorrmeangmr)/meankcorr
;  differr = errkcorr/meankcorr
;  ploterror, medz, diff, differr, psym=1,$
;    xrange = [min(medz)-0.02, max(medz)+0.02], yrange=[-0.15,0.05],$
;    xtitle='Z',ytitle='% diff'
;  oplot,[-100,100],[0,0]

;  !p.multi=0

END   
