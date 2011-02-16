PRO calc_xigg_xigm_bias, rgg, xigg, xigg_err, rgm, xigm, cov_xigm, $
                         rcommon, $
                         bias, bias_cov, bias_err, $
                         cbias, cbias_err, $
                         bias_inv, bias_inv_cov, bias_inv_err, $
                         cbias_inv, cbias_inv_err, $
                         _extra=_extra

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Take the ratio of xigg/xigm and calculate the mean bias, by fitting
  ;; to a constant.  Do with diagonal errors and full covariance matrix
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; We will interpolate the xigg measurements to the xigm positions (since
  ;; it should be higher S/N.  Assume we only have diagonal errors on wgg

  minrgg = min(rgg, max=maxrgg)

  wm = where(rgm LE maxrgg AND rgm GE minrgg, nm)
  rcommon = rgm[wm]

  xigg_interp = interpol(xigg, rgg, rcommon)
  xigg_err_interp = interpol(xigg_err, rgg, rcommon)

  ;; now calculate the covariance matrix for the bias
  cov_xigg = diagonal_array(xigg_err_interp^2)

  sendcov = cov_xigm[ wm[0]:wm[nm-1], wm[0]:wm[nm-1] ]

  help, xigg_interp, cov_xigg, xigm[wm], sendcov

  calc_ratio_cov, xigg_interp, cov_xigg, xigm[wm], sendcov, $
                  bias, bias_cov
  calc_ratio_cov, xigm[wm], sendcov, xigg_interp, cov_xigg, $
                  bias_inv, bias_inv_cov

;  xigm_err = fltarr(nm)
;  FOR i=0L, nm-1 DO xigm_err[i] = sqrt(cov_xigm[wm[i],wm[i]])
;  tbias_err = bias*sqrt( (xigg_err_interp/xigg_interp)^2 + $
;                         (xigm_err/xigm[wm])^2 )
;  tbias_inv_err = bias_inv*sqrt( (xigg_err_interp/xigg_interp)^2 + $
;                             (xigm_err/xigm[wm])^2 )

;  FOR i=0L, nm-1 DO print,tbias_err[i],sqrt(bias_cov[i,i])
;  print
;  FOR i=0L, nm-1 DO print,tbias_inv_err[i],sqrt(bias_inv_cov[i,i])
;stop

  bias_err = fltarr(nm)
  bias_inv_err = bias_err
  FOR i=0L, nm-1 DO BEGIN
      bias_err[i] = sqrt(bias_cov[i,i])
      bias_inv_err[i] = sqrt(bias_inv_cov[i,i])
  ENDFOR 

  ;; now straight mean
  wmom, bias, bias_err, wbias, wbias_sig, wbias_err
  wmom, rcommon, bias_err, rmean, rsig, rerr,/calc
  wmom, bias_inv, bias_inv_err, wbias_inv, wbias_inv_sig, wbias_inv_err
  wmom, rcommon, bias_inv_err, rmean_inv, rsig_inv, rerr_inv,/calc

  !p.multi=[0,0,2]

  setup_mystuff
  crange = [-10.0, 10.0]
  nc = 10000

  ;; Fit for bias
  fit_const, bias, bias_cov, crange, nc, chisq, $
             cbias, cerrhigh=cbias_errh,cerrlow=cbias_errl,/noplot
  cbias_err = cbias_errh[0]

  ;; Fit for 1/b
  fit_const, bias_inv, bias_inv_cov, crange, nc, chisq, $
             cbias_inv, cerrhigh=cbias_inv_errh,cerrlow=cbias_inv_errl,/noplot
  cbias_inv_err = cbias_inv_errh[0]

  tt = 1./cbias_inv
  tterr = tt*cbias_inv_err/cbias_inv

  print,'---------------------------------------------------------'
  print,'<b> = '+ntostr(cbias)+' '+!plusminus+' '+ntostr(cbias_err)+$
        ' Full Cov'
  print,'<b> = '+ntostr(wbias)+' '+!plusminus+' '+ntostr(wbias_err)+$
        ' Weighted Mean'
  print,'Effective Radius: '+ntostr(rmean)+' '+!plusminus+' '+ntostr(rerr)
  print
  print,'<1/b> = '+ntostr(cbias_inv)+' '+!plusminus+' '+ntostr(cbias_inv_err)+$
        ' Full Cov'
  print,'<1/b> = '+ntostr(wbias_inv)+' '+!plusminus+' '+ntostr(wbias_inv_err)+$
        ' Weighted Mean'
  print,'Effective Radius: '+ntostr(rmean_inv)+' '+!plusminus+' '+ntostr(rerr_inv)
  print,'<1/b>^-1 = '+ntostr(tt)+' '+!plusminus+' '+ntostr(tterr)+$
        ' Full Cov'
  print,'---------------------------------------------------------'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot b
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  setup_mystuff

  biasyt = 'b!Dgm!N = '+!csym.xi+'!Dgg!N/'+!csym.xi+'!Dgm!N'
  aploterror, 1.7, rcommon, bias, bias_err, /xlog, xstyle=2, $
              xtitle = !xigmxtitle, ytitle=biasyt, psym=8,/ynozero, $
              _extra=_extra

  minr = min(rcommon, max=maxr)

  xx = [minr, maxr, maxr, minr]
  yy = [cbias+cbias_err,cbias+cbias_err,$
        cbias-cbias_err,cbias-cbias_err]
  polyfill, xx, yy, /line_fill, orientation=45

  mess = '<b!Dgm!N> = '+$
    ntostr(cbias,4,/round)+!csym.plusminus+$
    ntostr(cbias_err,4,/round)
  legend, mess, /right, box=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now 1/b
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  biasyt = 'b!S!Dgm!R!U'+!csym.minus+'1!N = '+!csym.xi+'!Dgm!N/'+!csym.xi+'!Dgg!N'
  aploterror, 1.7, rcommon, bias_inv, bias_inv_err, /xlog, xstyle=2, $
              xtitle = !xigmxtitle, ytitle=biasyt, psym=8,/ynozero, $
              _extra=_extra

  xx = [minr, maxr, maxr, minr]
  yy = [cbias_inv+cbias_inv_err,cbias_inv+cbias_inv_err,$
        cbias_inv-cbias_inv_err,cbias_inv-cbias_inv_err]
  polyfill, xx, yy, /line_fill, orientation=45

  mess = '<b!S!Dgm!R!U'+!csym.minus+'1!N> = '+$
    ntostr(cbias_inv,4,/round)+!csym.plusminus+$
    ntostr(cbias_inv_err,4,/round)
  legend, mess, /right, box=0


  !p.multi=0

END 
