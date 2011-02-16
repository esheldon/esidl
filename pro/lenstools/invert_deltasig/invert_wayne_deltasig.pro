PRO invert_wayne_deltasig, dops=dops, replace=replace

  dir = '~/halo_model/'
  
  fitfile = dir+'wayne_invert.fit'
  deltasig_file = dir+'sig.dat'
  xi_file = dir+'xi.dat'
  psfile = dir+'wayne_invert.ps'
  yrange = [0.1, 1.e5]

  IF fexist(fitfile) AND (NOT keyword_set(replace)) THEN BEGIN 
      print,'Reading file: ',fitfile
      struct = mrdfits(fitfile, 1)
  ENDIF ELSE BEGIN 
      readcol, deltasig_file, meanr, deltasig, junk1, junk2, junk3

      readcol, xi_file, r3, ximm, xigg, xigm, bgm, bgg
      
      ;; make covariance
      ;; give all same S/N
      SN = 100
      deltasig_err = deltasig/SN
      xierr = xigm/SN

      covariance = diagonal_array(deltasig_err)

      ;; dave's code
      print,'Calling delta_sigma2der_sigma'
      delta_sigma2der_sigma, meanr, deltasig, dersig, covariance, $
                             covder, derr
      print,'Calling der_sigma2xi'
      der_sigma2xi, meanr, dersig, r3out, xiout, corrfac, covder, covxiout, $
                    xi_interior=xi_interior
      !p.multi=0

      nr3 = n_elements(r3out)
      xierrout = fltarr(nr3)
      FOR i=0L, nr3-1 DO xierrout[i] = sqrt(covxiout[i,i])

      ;; Wayne used 0.27
;      xiout = xiout*0.3/0.27

      struct=create_struct('meanr', meanr, $
                           'deltasig', deltasig, $
                           'deltasig_err', deltasig_err, $
                           'r3', r3, $
                           'xi', xigm, $
                           'xierr', xierr, $
                           'r3_invert', r3out, $
                           'xi_invert', xiout, $
                           'xierr_invert', xierrout, $
                           'xicov_invert', covxiout)

      print,'Writing file: ',fitfile
      mwrfits, struct, fitfile, /create

  ENDELSE 


  IF keyword_set(dops) THEN begplot,name=psfile, /color
      
  setup_mystuff

  xicolor = !red
  IF !d.name EQ 'PS' THEN invcolor = !blue ELSE invcolor = !green

  plottwo, struct.r3, struct.xi, struct.r3, struct.xi_invert/struct.xi, $
           frac1=0.75, topaspect=1, /topylog, /xlog, toplinestyle=0, $
           xoplot=struct.r3, yoplot=struct.xi_invert, oplotcolor=invcolor, $
           oplotlinestyle=2, $
           xrange=[0.01,10], xstyle=1+2, /ynozero, $
           xtitle=!xigmxtitle, topytitle=!xigmytitle, $
           botytitle=!csym.xi+'!DINV!N /'+!csym.xi, $
           xwindow1=xwindow1, ywindow1=ywindow1

;  aplot, 1, [0], [0], /nodata, psym=3, /xlog, /ylog, $
;         yrange = yrange, xrange=[0.01, 20.0], $
;         ystyle=1+2, xstyle=1+2
;  oplot, struct.r3, struct.xi, psym=3, $
;              color=xicolor
;  oplot, struct.r3, struct.xi, color=xicolor

;  oplot, struct.r3_invert, struct.xi_invert, $
;              psym=3, color=invcolor

  xwold = !x.window & !x.window=xwindow1
  ywold = !y.window & !y.window=ywindow1
  
  legend, ['3d', 'Inversion'], line=[0,2], color=[!p.color, invcolor], $
          /right, box=0, thick=[!p.thick,!p.thick]

  !x.window=xwold & !y.window=ywold
;  oplot, r3, (r0/r3)^gam, color=oclr

  IF keyword_set(dops) THEN endplot


return
END 
