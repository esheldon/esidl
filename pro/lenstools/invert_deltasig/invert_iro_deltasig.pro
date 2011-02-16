PRO invert_iro_deltasig, type, dops=dops, replace=replace, comoving=comoving

  dir = '/net/cheops1/data3/lensout/iro/'
  IF keyword_set(comoving) THEN addstr='_comoving' ELSE addstr=''

  IF type EQ 1 THEN BEGIN 
      fitfile = dir+'iro_invert_18.0_to_21.6_sdssw_z0.1'+addstr+'.fit'
      deltasig_file = dir+'gm2d_18.0_to_21.6_sdssw_z0.1.dat'
      xi_file = dir+'gm3d_18.0_to_21.6_sdssw_z0.1.dat'
      psfile = dir+'iro_invert_18.0_to_21.6_sdssw_z0.1'+addstr+'.ps'
      yrange = [0.1, 1.e5]
  ENDIF ELSE BEGIN 
      fitfile = dir + 'iro_invert_21.6_to_22.0_sdssw_z0.1'+addstr+'.fit'
      deltasig_file = dir+'gm2d_21.6_to_22.0_sdssw_z0.1.dat'
      xi_file = dir+'gm3d_21.6_to_22.0_sdssw_z0.1.dat'
      psfile = dir + 'iro_invert_21.6_to_22.0_sdssw_z0.1'+addstr+'.ps'
      yrange = [0.1, 2.e5]
  ENDELSE 

  IF fexist(fitfile) AND (NOT keyword_set(replace)) THEN BEGIN 
      print,'Reading file: ',fitfile
      struct = mrdfits(fitfile, 1)
  ENDIF ELSE BEGIN 
      readcol, deltasig_file, meanr, deltasig, junk, junk, $
               dsigerr_poiss, junk, dsigerr_jack, dsigerr_max

      readcol, xi_file, r3, xi, xierr_poisson, junk, xierr_jack, xierr_max
      
      ;; Will send this to der_sigma2xi so we get xi at the mean redshift.  
      ;; Then convert to comoving by shifting the scale
      zmean = 0.1

      ;; Their box has omega_m = 0.3 at redshift 0
      omega_m = 0.3
      
      ;; make covariance
      covariance = diagonal_array(dsigerr_max^2)

      ;; dave's code
      print,'Calling delta_sigma2der_sigma'
      delta_sigma2der_sigma, meanr, deltasig, dersig, covariance, $
                             covder, derr
      print,'Calling der_sigma2xi'
      der_sigma2xi, meanr, dersig, r3out, xiout, corrfac, covder, covxiout, $
                    xi_interior=xi_interior, zmean=zmean, omega_m=omega_m
      !p.multi=0

      nr3 = n_elements(r3out)
      xierrout = fltarr(nr3)
      FOR i=0L, nr3-1 DO xierrout[i] = sqrt(covxiout[i,i])

      IF keyword_set(comoving) THEN BEGIN 

          print,'Converting to comoving'
          meanr = meanr*(1+zmean)
          r3 = r3*(1+zmean)
          deltasig = deltasig/(1+zmean)^2
          dsigerr_max = dsigerr_max/(1+zmean)^2
          r3out = r3out*(1+zmean)

;          xiout = xiout/(1+zmean)
      ENDIF 

      struct=create_struct('meanr', meanr, $
                           'deltasig', deltasig, $
                           'deltasig_err', dsigerr_max, $
                           'r3', r3, $
                           'xi', xi, $
                           'xierr', xierr_max, $
                           'r3_invert', r3out, $
                           'xi_invert', xiout, $
                           'xierr_invert', xierrout, $
                           'xicov_invert', covxiout)

      print,'Writing file: ',fitfile
      mwrfits, struct, fitfile, /create

  ENDELSE 

  IF keyword_set(dops) THEN begplot,name=psfile, /color
  setup_mystuff

  xicolor = !p.color
  IF !d.name EQ 'PS' THEN invcolor = !blue ELSE invcolor = !green

  psym=8
  oplotsym=4

  ratio = struct.xi_invert/struct.xi
  ratioerr = ratio*sqrt( (struct.xierr_invert/struct.xi_invert)^2 + $
                         (struct.xierr/struct.xi)^2 )
  plottwoerr, struct.r3, struct.xi, struct.xierr, $
              struct.r3, ratio, ratioerr, $
              psym=psym, oplotsym=oplotsym, $
              frac1=0.75, topaspect=1, /topylog, /xlog, toplinestyle=0, $
              xoplot=struct.r3_invert, $
              yoplot=struct.xi_invert, oploterr=struct.xierr_invert, $
              oplotcolor=invcolor, $
              oplotlinestyle=2, $
              xrange=[0.01,10], xstyle=1+2, /ynozero, $
              xtitle=!xigmxtitle, topytitle=!xigmytitle, $
              botytitle=!csym.xi+'!DINV!N /'+!csym.xi, $
              xwindow1=xwindow1, ywindow1=ywindow1, $
              botyrange = [0.5,2.5]

  oplot, [0.001, 100], [1,1]

;  aplot, 1, [0], [0], /nodata, psym=8, /xlog, /ylog, $
;         yrange = yrange, xrange=[0.01, 20.0], $
;         ystyle=1+2, xstyle=1+2
;  oploterror, struct.r3, struct.xi, struct.xierr, psym=8, $
;              color=xicolor, errc=xicolor
;  oplot, struct.r3, struct.xi, color=xicolor

;  oploterror, struct.r3_invert, struct.xi_invert, struct.xierr_invert, $
;              psym=8, color=invcolor, errc=invcolor

  xwold = !x.window & !x.window=xwindow1
  ywold = !y.window & !y.window=ywindow1

  legend, ['3d', 'Inversion'], psym=[psym,oplotsym], $
          color=[xicolor, invcolor], $
          /right, box=0

  !x.window=xwold & !y.window=ywold
;  oplot, r3, (r0/r3)^gam, color=oclr

  colprint,ratio

  IF keyword_set(dops) THEN endplot

return
END 
