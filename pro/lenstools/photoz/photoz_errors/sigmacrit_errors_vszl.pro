PRO sigmacrit_errors_vszl, overwrite=overwrite

  ;; runs sigmacrit_errors for various zL and fits a line to the
  ;; resulting norm and power law indices 
  ;; (zerr/sigmacrit_err) = norm(zL)*zsource^( -pow(zL) )
  ;;
  ;; norm(zL) = a_norm + b_norm*zL
  ;; pow(zL) = a_pow + b_pow*zL

  ;; see run_sigmacrit_errors for the zLarr values
  ;;zLarr = [0.050,0.075,0.100,0.125,0.150,0.175,0.200]

  zLarr = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, $
           0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, $
           0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30]

  nzL = n_elements(zLarr)

  ;; these are the files that were created
  dir = '/net/cheops1/data0/esheldon/pzerr2scriterr/'
  files = dir + 'sigerr_zL'+ntostr(zLarr,5)+'.fit'
  outfile = dir + 'sigerr_normpow_vszl.fit'
  outps   = dir + 'sigerr_normpow_vszl.ps'

  ;; read in the output files and fit a line to
  ;; norm vs zL and pow vs. zL
  
  norm = fltarr(nzL)
  pow  = fltarr(nzL)
  FOR i=0L, nzL-1 DO BEGIN 

      IF fexist(files[i]) AND (NOT keyword_set(overwrite)) THEN BEGIN 
          tmp = mrdfits(files[i],1)
          tmp.pow = abs(tmp.pow)

          norm[i] = tmp.norm
          pow[i]  = tmp.pow
      ENDIF ELSE BEGIN 
          sigmacrit_errors, zLarr[i], tnorm, tpow, /dops
          norm[i] = tnorm
          pow[i] = tpow
      ENDELSE 

  ENDFOR 

  ;; fit power law to norm vs zL
  ;;result = linfit(zLarr, norm)
  ;;a_norm = result[0]
  ;;b_norm = result[1]

  aguess = [1.0,1.0]
  fitpower, zLarr, norm, replicate(1.0, n_elements(zLarr)), aguess, $
            yfit, result
  ap_norm = result[0]
  bp_norm = result[1]

  ;; fit line to pow 
  result = linfit(zLarr, pow)
  a_pow = result[0]
  b_pow = result[1]

  !p.multi=[0,0,2]
  begplot, name=outps

  xrange = [0.9*min(zLarr), 1.1*max(zLarr)]
  yrange = [0.9*min(norm), 1.1*max(norm)]
  title = 'Norm of power law [ '+$
    !csym.sigma+'('+!csym.sigma_cap+'!Dcrit!N) / '+$
    !csym.sigma_cap+'!Dcrit!N ] / '+$
    !csym.sigma+'(z!Ds!N) vs z!Ds!N'
  xtitle = 'z!DL!N'
  ytitle = 'Norm'
  plot, zLarr, norm, psym=4, /ynozero, $
        xrange=xrange, yrange=yrange, $
        title=title,xtitle=xtitle,ytitle=ytitle,charsize=1
  ;; show it goes through zero
  tx = arrscl( findgen(100), 0.0, max(zLarr) )
  ty = ap_norm*tx^bp_norm
  oplot, tx, ty
  message = ['norm = a'+!csym.times+'z!DL!Ub!N', $
             'a = '+ntostr(ap_norm,4,/round), $
             'b = '+ntostr(bp_norm,4,/round)]
  legend,message,/left,charsize=1

  yrange = [0.9*min(pow), 1.1*max(pow)]
  title = 'Index of power law [ '+$
    !csym.sigma+'('+!csym.sigma_cap+'!Dcrit!N) / '+$
    !csym.sigma_cap+'!Dcrit!N ] / '+$
    !csym.sigma+'(z!Ds!N) vs z!Ds!N'
  xtitle = 'z!DL!N'
  ytitle = 'Power Law Index'
  plot, zLarr, pow, psym=4, /ynozero, xrange=xrange,$
        title=title,xtitle=xtitle,ytitle=ytitle,charsize=1
  oplot, zLarr, a_pow + b_pow*zLarr
  message = ['pow = a + b'+!csym.times+'z!DL!N', $
             'a = '+ntostr(a_pow,4,/round), $
             'b = '+ntostr(b_pow,4,/round)]
  legend,message,/left,charsize=1

  endplot

  str = create_struct('zl', zLarr, $
                      'norm', norm, $
                      'pow', pow, $
                      'ap_norm', ap_norm, $
                      'bp_norm', bp_norm, $
                      'a_pow',  a_pow, $
                      'b_pow',  b_pow)

  mwrfits, str, outfile, /create

  !p.multi=0

END 
