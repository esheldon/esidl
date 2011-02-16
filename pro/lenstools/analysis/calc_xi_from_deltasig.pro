PRO calc_xi_from_deltasig_rebin_data, x, y, yerr, newx, newy, newyerr

  nx = n_elements(x)
  new_nx = nx/2

  nxmod = nx MOD 2

  newx = fltarr(new_nx)
  newy = newx
  newyerr = newx

  FOR i=0L, new_nx-1 DO BEGIN 

      IF nxmod NE 0 AND i EQ new_nx-1 THEN BEGIN 
          bin1 = i*2
          bin2 = i*2+2
      ENDIF ELSE BEGIN 
          bin1 = i*2
          bin2 = i*2 + 1
      ENDELSE 

      tx = x[bin1:bin2]
      ty = y[bin1:bin2]
      tyerr = yerr[bin1:bin2]

      wmom, ty, tyerr, wmean, wsig, werr

      newy[i] = wmean
      newyerr[i] = werr

      wmom, tx, tyerr, wmean, wsig, werr
      
      newx[i] = wmean

  ENDFOR 

END 

PRO calc_xi_from_deltasig_addstruct, t, r3, xi, covxi, corrfac, xi_interior, $
                                     newt

  nr = n_elements(r3)
  xierr = fltarr(nr)
  FOR i=0L,  nr-1 DO xierr[i] = sqrt(covxi[i,i])


  ;; The lower bound on xi should be xi_interior - xierr 
  ;; This will redefine xierr
  xierr2 = ( xi - (xi_interior - xierr) ) > xierr
  print
  print,'       OldErr       NewErr        Ratio'
  forprint, xierr, xierr2, xierr2/xierr
  xierr = xierr2

  ;; now remake the covariance matrix using this
  ;; new rescaling of the errors
  corrxi = cov2corr(covxi)
  covxi = corr2cov(corrxi, xierr)

  nrstr = ntostr(nr)
  add_tags, t, ['r3', 'xi', 'xierr', $
                'rmax_corrfac','xi_interior',$
                'covxi','corrxi'], $
            ['fltarr('+nrstr+')', $
             'fltarr('+nrstr+')', $
             'fltarr('+nrstr+')', $
             'fltarr('+nrstr+')', $
             'fltarr('+nrstr+')', $
             'fltarr('+nrstr+','+nrstr+')',$
             'fltarr('+nrstr+','+nrstr+')'], newt




  newt.r3 = r3
  newt.xi = xi
  newt.xierr = xierr
  newt.rmax_corrfac = corrfac
  newt.xi_interior = xi_interior
  newt.covxi = covxi
  newt.corrxi = corrxi

END 

PRO calc_xi_from_deltasig, type, skip=skip, rebin=rebin, noerror=noerror

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; For creating xi from delta sigma.  For plots use plot_xi
; First you pick the range of radii for fitting the power law.  This
; Power law is for extrapolation, used to calculate xi.  Especially
; only use the outer values if there is possibility of different slope
; at smaller radii. 
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  max_radius = 7.0              ;Mpc

  dir = '/net/cheops1/data3/lensout/combstripe/comb/'

  IF !d.name EQ 'X' THEN BEGIN 
      clrs = [!green, !red, !cyan]
      ticklen = 0.04
  ENDIF ELSE BEGIN 
      clrs = [!p.color, !blue, !red]
      ticklen = 0.04
  ENDELSE 

  IF n_elements(type) EQ 0 THEN type = ''

  IF keyword_set(skip) THEN GOTO, jump

  ;; comoving?
  comove_str = 'comoving_'

  IF type EQ 'lum' OR type EQ 'redlum' THEN GOTO, lumjump

  ;; doing the single type
  


  ext = comove_str + 'N1.fit'
  IF type EQ 'eclass1_' THEN ext = comove_str + 'N2.fit'
  IF type EQ 'eclass2_' THEN ext = comove_str + 'N2.fit'

  file = dir + type+'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_'+ext
  outfile = dir + type+'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_'+ext

  print,'Reading file: ',file
  t = mrdfits(file, 1)

  all = mrdfits(dir + 'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_'+comove_str+'N1.fit',1)

  ;; We will use the rebinned stuff for blue galaxies
  IF type EQ 'eclass2_' OR type EQ 'gmr1_' THEN BEGIN 
     
      ;; need to redo the rebinning
;      fac = 4
;      t.sigmaerr[4] = t.sigmaerr[4]*fac
;      t.sigmaerr[5] = t.sigmaerr[5]*fac

      corr = cov2corr(all.covariance)
      cov = corr2cov(corr, t.sigmaerr)
      t.covariance = cov

      rebin_data, t.meanr, t.sigma, meanr, sigma, $
                  yerr = t.sigmaerr, newyerr=sigmaerr, $
                  cov = t.covariance, newcov = cov

      meanr = meanr/1000
      t.sigma_rebin = sigma
      t.sigmaerr_rebin = sigmaerr

;      nrebin = n_elements(t.meanr_rebin)
;      meanr = t.meanr_rebin/1000.
;      sigma = t.sigma_rebin
;      cov = fltarr(nrebin,nrebin)
;      FOR i=0L, nrebin-1 DO cov[i,i] = t.sigmaerr_rebin[i]^2

  ENDIF ELSE BEGIN 
      meanr = t.meanr/1000.
      sigma = t.sigma
      cov = t.covariance
  ENDELSE 

  ;; dave's code
  print,'Calling delta_sigma2der_sigma'
  delta_sigma2der_sigma, meanr, sigma, dersig, cov, $
                         covder, derr
  print,'Calling der_sigma2xi'
  der_sigma2xi, meanr, dersig, r3, xi, corrfac, covder, covxi, $
                xi_interior=xi_interior

  wplot = where(r3 LT max_radius)
  plot, r3, corrfac/xi_interior,xtitle='r',ytitle='C/xi_interior',/xlog
  oplot, r3, corrfac/xi_interior, psym=8
  oplot, r3[wplot], corrfac[wplot]/xi_interior[wplot], psym=8, color=!green
  key=get_kbrd(1)

  !p.multi=0

  nr = n_elements(r3)
  calc_xi_from_deltasig_addstruct, t, r3, xi, covxi, corrfac, xi_interior, newt

  t = newt
  xierr = t.xierr

  miny = min(xi-xierr) > 0.1
  maxy = max(xi+xierr)
  xrange = [0.015, 7.0]

  ;; rebin
  IF type EQ 'eclass2_' OR type EQ 'gmr1_' THEN BEGIN 
      ;; in this case, we used the already rebinned data
      newr3 = r3
      newxi = xi
      newxierr = xierr
  ENDIF ELSE BEGIN 
      calc_xi_from_deltasig_rebin_data, t.r3[wplot], t.xi[wplot], $
                                        t.xierr[wplot], $
                                        newr3, newxi, newxierr
  ENDELSE 

  nnew = ntostr(n_elements(newr3))
  sarr = 'fltarr('+nnew+')'
  add_tags, t, $
            ['r3_rebin', 'xi_rebin', 'xierr_rebin'], $
            [sarr,sarr,sarr], newt
  t = newt
  t.r3_rebin = newr3
  t.xi_rebin = newxi
  t.xierr_rebin = newxierr

  aploterror,1,r3[wplot],xi[wplot],xierr[wplot],/xlog,/ylog,$
             ystyle=1+2,yrange=[miny,maxy],xstyle=1+2,xrange=xrange,psym=8, $
             ytickf='loglabels',xtickf='loglabels'


;  oplot, r3, yfit, color=!green
;  oplot, r3[wplot], (r0/r3[wplot])^gamma, color=!green
;  if !d.name eq 'X' then key=get_kbrd(1)


  print
  print,'           r3           xi      corrfac'
  colprint,r3,xi,corrfac

  oploterror, newr3, newxi, newxierr, psym=4, color=!green, errc=!green
  print
  print,'Output file: ',outfile
  mwrfits, t, outfile, /create

return

lumjump:

;  clr = 2

  FOR clr = 0,4 DO BEGIN 

      cstr = !colors[clr]
      
      FOR lumnum = 1,3 DO BEGIN 
          lstr = ntostr(lumnum)
          
          IF type EQ 'redlum' THEN BEGIN 
              file = dir + 'sublum/'+cstr+'/redlum'+lstr+'threebin_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_'+comove_str+'N1.fit'
              outfile = dir + 'sublum/'+cstr+'/redlum'+lstr+'threebin_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_'+comove_str+'N1.fit'

          ENDIF ELSE BEGIN 
              file = dir + 'sublum/'+cstr+'/lum'+lstr+'threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_'+comove_str+'N1.fit'
              outfile = dir + 'sublum/'+cstr+'/lum'+lstr+'threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_'+comove_str+'N1.fit'
          ENDELSE 
          print,'Reading file: ',file
          t = mrdfits(file, 1)
          
;      delta_sigma2xi, t.meanr/1000., t.sigma, t.covariance, r3,xi,covxi,xierr
          
          delta_sigma2der_sigma, t.meanr/1000, t.sigma, dersig, t.covariance, $
                                 covder, derr
          der_sigma2xi, t.meanr/1000, dersig, r3, xi, corrfac, covder, covxi,$
                        xi_interior=xi_interior

          wplot = where(r3 LT max_radius)
          plot, r3, corrfac/xi_interior,xtitle='r',ytitle='C/xi_interior',/xlog
          oplot, r3, corrfac/xi_interior, psym=8
          oplot, r3[wplot], corrfac[wplot]/xi_interior[wplot], psym=8, $
                 color=!green
          key=get_kbrd(1)

          !p.multi=0


          nr = n_elements(r3)
          calc_xi_from_deltasig_addstruct, t, r3, xi, $
                                           covxi, corrfac, xi_interior, newt
          t = newt
          xierr = newt.xierr
          
          miny = min(xi-xierr) > 0.1
          maxy = max(xi+xierr)
          yrange = [miny,maxy]
          xrange = [0.02, 10.0]
          
          calc_xi_from_deltasig_rebin_data, t.r3[wplot], t.xi[wplot], $
                                            t.xierr[wplot], $
                                            newr3, newxi, newxierr
          nnew = ntostr(n_elements(newr3))
          sarr = 'fltarr('+nnew+')'
          add_tags, t, $
                    ['r3_rebin', 'xi_rebin', 'xierr_rebin'], $
                    [sarr,sarr,sarr], newt
          t = newt
          t.r3_rebin = newr3
          t.xi_rebin = newxi
          t.xierr_rebin = newxierr

          aploterror,1,r3[wplot],xi[wplot],xierr[wplot],/xlog,/ylog,$
                     ystyle=1+2,yrange=yrange,xstyle=1+2,$
                     xrange=xrange,psym=8, $
                     ytickf='loglabels',xtickf='loglabels'
          oploterror, newr3, newxi, newxierr, psym=4, color=!green, errc=!green

          if !d.name eq 'X' then key=get_kbrd(1)
          
          print
          print,'Output file: ',outfile
          mwrfits, t, outfile, /create
          
      ENDFOR 

  ENDFOR 

jump:

return

  erase & multiplot, [3,2], /square

  yrange = [0.5,1.e5]
  xrange = [0.015, 7.0]

  yt = !csym.xi+'!Dgm!N(r)'
  xt = 'r [ h!U'+!csym.minus+'1!N Mpc ]'

  FOR clr=0,4 DO BEGIN 

      cstr = !colors[clr]
      
      FOR lumnum = 1,3 DO BEGIN 
          lstr = ntostr(lumnum)
          
          file = dir + 'sublum/'+cstr+'/lum'+lstr+'threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb'+comove_str+'_xi_N1.fit'
          
          print,'Reading file: ',file
          t = mrdfits(file, 1)
          
          ;;miny = min(t.xi-t.xierr) > 0.1
          ;;maxy = max(t.xi+t.xierr)
          ;;yrange = [miny,maxy]
          
          IF keyword_set(rebin) THEN BEGIN 
              r3 = t.r3_rebin
              xi = t.xi_rebin
              xierr = t.xierr_rebin
          ENDIF ELSE BEGIN 
              wplot = where(t.r3 LT max_radius)
              r3 = t.r3[wplot]
              xi = t.xi[wplot]
              xierr = t.xierr[wplot]
          ENDELSE 
          
          delvarx, ytickf, xtickf, ytitle, xtitle
          CASE clr OF 
              0: BEGIN 
                  ytitle = yt
                  ytickf = 'loglabels'
              END 
              2: BEGIN 
                  xtitle = xt
                  xtickf = 'loglabels'
              END 
              3: BEGIN 
                  xtitle = xt
                  ytitle = yt
                  xtickf='loglabels' 
                  ytickf=xtickf 
              END 
              4: BEGIN 
                  xtitle = xt
                  xtickf = 'loglabels'
              END 
              ELSE: 
          ENDCASE 

          CASE lumnum OF
              1: psym = 8
              2: psym = 1
              3: psym = 5
              ELSE: message,'What?'
          ENDCASE 

          IF lumnum EQ 1 THEN BEGIN 
              IF NOT keyword_set(noerror) THEN BEGIN 
                  ploterror,r3,xi,xierr,/xlog,/ylog,$
                            ystyle=1+2,yrange=yrange,xstyle=1+2,$
                            xrange=xrange,psym=psym, $
                            ytickf=ytickf,xtickf=xtickf, $
                            xtitle=xtitle, ytitle=ytitle,$
                            xticklen=ticklen,yticklen=ticklen,charsize=1.25
              ENDIF ELSE BEGIN 
                  plot,r3,xi,/xlog,/ylog,$
                       ystyle=1+2,yrange=yrange,xstyle=1+2,$
                       xrange=xrange,psym=psym, $
                       ytickf=ytickf,xtickf=xtickf, $
                       xtitle=xtitle, ytitle=ytitle,$
                       xticklen=ticklen,yticklen=ticklen,charsize=1.25
              ENDELSE 

          ENDIF 
          oplot, r3, (xi > 0.01), color=clrs[lumnum-1]
          IF NOT keyword_set(noerror) THEN BEGIN 
              oploterror, r3, xi, xierr,psym=psym,color=clrs[lumnum-1],errc=clrs[lumnum-1]
          ENDIF ELSE BEGIN 
              oplot, r3, xi, psym=psym,color=clrs[lumnum-1]
          ENDELSE 
;          oplot, r3, (t.r0/r3)^t.gamma, color=clrs[lumnum-1],line=2
          

      ENDFOR

      legend, !colors[clr], /right, box=0

      CASE clr OF 
          0: multiplot
          1: multiplot
          2: multiplot
          3: multiplot
          ELSE: 
      ENDCASE 

 
  ENDFOR 

  multiplot,/reset

END 
