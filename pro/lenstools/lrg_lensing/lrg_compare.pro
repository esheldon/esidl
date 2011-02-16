PRO lrg_compare, dops=dops

  plotdir = '~/plots/lrg/'

  file = plotdir + 'compare_samples.ps'
  IF keyword_set(dops) THEN begplot, name=file, /color

  l1 = obj_new('lrg_lensing', 1) ; normal, princeton regauss
  l2 = obj_new('lrg_lensing', 2) ; zbuffer of 0.1
  l3 = obj_new('lrg_lensing', 3) ; zbuffer of 0.2
  l4 = obj_new('lrg_lensing', 4) ; analytic corr
  l5 = obj_new('lrg_lensing', 5) ; more well resolved r > 0.585
  l6 = obj_new('lrg_lensing', 6) ; less well resolved r <= 0.585
  l7 = obj_new('lrg_lensing', 7) ; seeing < 1.5
  l8 = obj_new('lrg_lensing', 8) ; min angle 20 arcsec
  l9 = obj_new('lrg_lensing', 9) ; min angle 10 arcsec
  l10 = obj_new('lrg_lensing', 10) ; not blended
  l11 = obj_new('lrg_lensing', 11) ; r > 1/3
  l12 = obj_new('lrg_lensing', 12) ; r > 2/3

  t1 = l1->combined_get()
  t2 = l2->combined_get()
  t3 = l3->combined_get()
  t4 = l4->combined_get()
  t5 = l5->combined_get()
  t6 = l6->combined_get()
  t7 = l7->combined_get()
  t8 = l8->combined_get()
  t9 = l9->combined_get()
  t10 = l10->combined_get()
  t11 = l11->combined_get()
  t12 = l12->combined_get()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Compare redshift bins
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  z1_0 = l1->combined_get(z_nbin=2, z_bin=0)
  z1_1 = l1->combined_get(z_nbin=2, z_bin=1)

  oclr = !darkGreen
  plot_density_contrast, z1_0, /log, /mpc, aspect=1
  oploterror, z1_1.meanr/1000, z1_1.sigma, z1_1.sigmaerr, psym=4, color=oclr, errc=oclr


  mess = ['Zlens < 0.3', 'Zlens > 0.3']
  legend, mess, psym=[8,4], color=[!p.color, oclr], $
    /right, box=0, charsize=1


  key = prompt_kbrd()


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Compare normal with zbuffer of 0.1
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  oclr = !darkGreen
  plot_density_contrast, t1, /log, /mpc, aspect=1
  oploterror, t2.meanr/1000, t2.sigma, t2.sigmaerr, psym=4, color=oclr, errc=oclr


  mess = ['Standard','z-buffer = 0.1']
  legend, mess, psym=[8,4], color=[!p.color, oclr], $
    /right, box=0, charsize=1


  key = prompt_kbrd()


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Compare normal with zbuffer of 0.2
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  oclr = !darkGreen
  plot_density_contrast, t1, /log, /mpc, aspect=1
  oploterror, t3.meanr/1000, t3.sigma, t3.sigmaerr, psym=4, color=oclr, errc=oclr


  mess = ['Standard','z-buffer = 0.2']
  legend, mess, psym=[8,4], color=[!p.color, oclr], $
    /right, box=0, charsize=1


  key = prompt_kbrd()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Compare more versus less well resolved
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  oclr = !darkGreen
;  plot_density_contrast, t5, /log, /mpc, aspect=1
;  oploterror, t6.meanr/1000, t6.sigma, t6.sigmaerr, psym=4, color=oclr, errc=oclr

;  oplot, t1.meanr/1000, t1.sigma, color=!red
;  mess = ['Res > 0.585','Res < 0.585']
;  legend, mess, psym=[8,4], color=[!p.color, oclr], $
;    /right, box=0, charsize=1

;  key = prompt_kbrd()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Compare all to well resolved
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  oclr = !darkGreen
  plot_density_contrast, t1, /log, /mpc, aspect=1
  oploterror, t5.meanr/1000, t5.sigma, t5.sigmaerr, psym=4, color=oclr, errc=oclr

  mess = ['Standard','Res > 0.585']
  legend, mess, psym=[8,4], color=[!p.color, oclr], $
    /right, box=0, charsize=1

  key = prompt_kbrd()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Compare all to well resolved
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  oclr = !darkGreen
  plot_density_contrast, t1, /log, /mpc, aspect=1
  oploterror, t11.meanr/1000, t11.sigma, t11.sigmaerr, psym=4, color=oclr, errc=oclr

  mess = ['Standard','Res > 1/3']
  legend, mess, psym=[8,4], color=[!p.color, oclr], $
    /right, box=0, charsize=1

  key = prompt_kbrd()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Compare all to well resolved 2/3
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  oclr = !darkGreen
  plot_density_contrast, t1, /log, /mpc, aspect=1
  oploterror, t12.meanr/1000, t12.sigma, t12.sigmaerr, psym=4, color=oclr, errc=oclr

  mess = ['Standard','Res > 2/3']
  legend, mess, psym=[8,4], color=[!p.color, oclr], $
    /right, box=0, charsize=1

  key = prompt_kbrd()


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Compare more versus less well resolved with ratio
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  oclr = !darkGreen
  calc_ratio_cov, t5.sigma, t5.sigmaerr, t6.sigma, t6.sigmaerr, ratio, ratioerr
  plottwoerr, $
    t5.meanr/1000, t5.sigma, t5.sigmaerr, $
    t5.meanr/1000, ratio, ratioerr, $
    xoplot=t6.meanr/1000, yoplot=t6.sigma, oploterr=t6.sigmaerr, oplotcolor=oclr, $
    oplotsym=4, $
    psym=8, botpsym=8, $
    /xlog, /ylog, botyrange=[0.5, 50], botystyle=3, $
    topaspect=1, frac1=0.75, ytickformat1='loglabels', xtickformat='loglabels',$
    xwindow1=xwindow1, ywindow1=ywindow1, /center, $
    xtitle=!mpcxtitle2, topytitle=!deltaytitle, botytitle='ratio'

  oplot, [0.001, 100], [1,1]

  xwold = !x.window & !x.window=xwindow1
  ywold = !y.window & !y.window=ywindow1
  

  mess = ['Res > 0.585','Res < 0.585']
  legend, mess, psym=[8,4], color=[!p.color, oclr], $
    /right, box=0, charsize=1

  
  !x.window=xwold & !y.window=ywold



;  mess = ['Res > 0.585','Res < 0.585']
;  legend, mess, psym=[8,4], color=[!p.color, oclr], $
;    /right, box=0, charsize=1

  key = prompt_kbrd()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Compare seeing
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  oclr = !darkGreen
  plot_density_contrast, t1, /log, /mpc, aspect=1
  oploterror, t7.meanr/1000, t7.sigma, t7.sigmaerr, psym=4, color=oclr, errc=oclr


  mess = ['Standard','Seeing < 1.5']
  legend, mess, psym=[8,4], color=[!p.color, oclr], $
    /right, box=0, charsize=1

  key = prompt_kbrd()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Minimum angle of 10 arcsec
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  oclr = !darkGreen
  plot_density_contrast, t9, /log, /mpc, aspect=1, psym=4
  oploterror, t1.meanr/1000, t1.sigma, t1.sigmaerr, psym=8
  oploterror, t9.meanr/1000, t9.sigma, t9.sigmaerr, psym=4, color=oclr, errc=oclr


  mess = ['Standard','Minangle 10 arcsec']
  legend, mess, psym=[8,4], color=[!p.color, oclr], $
    /right, box=0, charsize=1

  key = prompt_kbrd()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Minimum angle of 20 arcsec
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  oclr = !darkGreen
  xrange = [0.5*min(t1.meanr), 1.5*max(t1.meanr)]/1000.
  plot_density_contrast, t8, /log, /mpc, aspect=1, psym=4, xrange=xrange
  oploterror, t1.meanr/1000, t1.sigma, t1.sigmaerr, psym=8
  oploterror, t8.meanr/1000, t8.sigma, t8.sigmaerr, psym=4, color=oclr, errc=oclr


  mess = ['Standard','Minangle 20 arcsec']
  legend, mess, psym=[8,4], color=[!p.color, oclr], $
    /right, box=0, charsize=1

  key = prompt_kbrd()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Not blended
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  oclr = !darkGreen
  plot_density_contrast, t10, /log, /mpc, aspect=1, psym=4
  oploterror, t1.meanr/1000, t1.sigma, t1.sigmaerr, psym=8
  oploterror, t10.meanr/1000, t10.sigma, t10.sigmaerr, psym=4, color=oclr, errc=oclr


  mess = ['Standard','Not blended']
  legend, mess, psym=[8,4], color=[!p.color, oclr], $
    /right, box=0, charsize=1




  obj_destroy, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12
  IF keyword_set(dops) THEN endplot

END 
