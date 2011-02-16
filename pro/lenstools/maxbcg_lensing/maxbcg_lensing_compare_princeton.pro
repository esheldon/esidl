
;; Make comparisions for samples using the two different shape corrections
;; schemes

PRO maxbcg_lensing_compare_princeton, $
  match=match, $
  corr=corr, $
  lold=lold, lnew=lnew, bcg=bcg

  ;; older style of shape correction
  mbold = obj_new('maxbcg_lensing',8)
  ;; regaussianizatioin method
  mbnew = obj_new('maxbcg_lensing',9)

  ngals_nbin = 12

  meanarr = fltarr(ngals_nbin)
  errarr  = fltarr(ngals_nbin)

  IF keyword_set(match) THEN BEGIN 
      IF n_elements(bcg) EQ 0 THEN bcg = mbold->get()
      IF n_elements(lold) EQ 0 THEN lold = mbold->lensoutput_get()
      IF n_elements(lnew) EQ 0 THEN lnew = mbnew->lensoutput_get()

      match, lold.zindex, lnew.zindex, mold, mnew

      mbold->ngals_bins, ngals_nbin, lowlim, highlim

  ENDIF 

;  !p.multi = [0,0,2]

  delvarx, allratio, allratioerr
  FOR ngals_bin=0L, ngals_nbin-1 DO BEGIN 


      IF keyword_set(match) THEN BEGIN 
          w = where(bcg[lold[mold].zindex].tngals GE lowlim[ngals_bin] AND $
                    bcg[lold[mold].zindex].tngals LE highlim[ngals_bin])

          cold = combine_lensum(lold[mold[w]],/silent)
          cnew = combine_lensum(lnew[mnew[w]],/silent)
      ENDIF ELSE BEGIN 
          cold = mbold->combined_get(ngals_nbin=ngals_nbin, ngals_bin=ngals_bin,$
                                 corr=corr)
          cnew = mbnew->combined_get(ngals_nbin=ngals_nbin, ngals_bin=ngals_bin,$
                                 corr=corr)
      ENDELSE 
 
      nrad = n_elements(cold.meanr)
      IF n_elements(allratio) EQ 0 THEN BEGIN 
          allratio = fltarr(ngals_nbin, nrad)
          allratioerr = fltarr(ngals_nbin, nrad)
      ENDIF 

;      wc = where(cold.meanr GT 100)
      wc = lindgen(nrad)
      calc_ratio_cov, cnew.sigma[wc], cnew.sigmaerr[wc], $
                      cold.sigma[wc], cold.sigmaerr[wc], $
                      ratio, ratioerr

      allratio[ngals_bin, *] = ratio
      allratioerr[ngals_bin, *] = ratioerr

      ;; mean ratio and errors
      wmom, ratio, ratioerr, wmean, wsig, werr
      print,'Mean = ',wmean,werr
      meanarr[ngals_bin] = wmean
      errarr[ngals_bin] = werr

      erase & multiplot, [1,2]
      plot_density_contrast, cold, /log
      oploterror, cnew.meanr, cnew.sigma, cnew.sigmaerr,psym=4, $
                  color=!green, errc=!green
      multiplot

      ploterror, cold.meanr[wc], ratio, ratioerr, psym=8, /xlog, $
                 yrange = wmean + [-1.0, 1.0]*5*wsig, ystyle=1+2, $
                 xrange = 10.0^!x.crange, xstyle=1
      oplot, [1,50000], [1,1]

      mean_error_legend, 'Mean', wmean, werr, /right, box=0

      multiplot,/reset

      key = prompt_kbrd()
      IF key EQ 'q' THEN return

  ENDFOR 


  ratio_sum = total(allratio/allratioerr^2, 1)
  ratio_wsum = total(1.0/allratioerr^2, 1)
  allratio = ratio_sum/ratio_wsum
  allratioerr = sqrt(1.0/ratio_wsum)

  ploterror, cold.meanr, allratio, allratioerr, psym=8,/xlog, xstyle=1+2
  oplot, [1, 50000], [1,1]

  obj_destroy, mbold, mbnew



  wmom, meanarr, errarr, wmean, wsig, werr

  oplot,[1,50000],[wmean,wmean], color=!green
  polyfill,10.0^!X.crange([0,1,1,0]),[wmean-werr,wmean-werr,wmean+werr,wmean+werr],$
           color=!green, /line_fill, orientation=45

  forprint,meanarr,errarr
  print,'------------------------'
  print,wmean,werr
END 
