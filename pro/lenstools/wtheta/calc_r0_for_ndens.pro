PRO calc_r0_for_ndens

  ;; space density of galaxies # per Mpc^2 from
  ;; integrating Blanton et al. Luminosity functions
  ;; from 0.1*Lstar
  ;;
  ;; need to look at errors

  dens = [0.111, 0.051, 0.033, 0.031, 0.031]

  indir='/sdss5/data0/lensout/stripe10/'

  r=mrdfits(indir+'matchrN8_wthetalumweq_stripe10_rw_N5.fit',1)
  rr=mrdfits(indir+'matchrN8_wthetarandlumweq_stripe10_rw_N5.fit',1)

  i=mrdfits(indir+'matchrN8_wthetalumweq_stripe10_iw_N5.fit',1)
  ii=mrdfits(indir+'matchrN8_wthetarandlumweq_stripe10_iw_N5.fit',1)

  z=mrdfits(indir+'matchrN8_wthetalumweq_stripe10_zw_N5.fit',1)
  zz=mrdfits(indir+'matchrN8_wthetarandlumweq_stripe10_zw_N5.fit',1)

  !p.multi=[0,1,3]

  tdiff = r.meanlum - rr.meanlum
  rdiff = r.npair - rr.npair
  rerr = sqrt(r.meanlumerr^2 + rr.meanlumerr^2)*(rdiff/tdiff)

  tdiff = i.meanlum - ii.meanlum
  idiff = i.npair - ii.npair
  ierr = sqrt(i.meanlumerr^2 + ii.meanlumerr^2)*(idiff/tdiff)

  tdiff = z.meanlum - zz.meanlum
  zdiff = z.npair - zz.npair
  zerr = sqrt(z.meanlumerr^2 + zz.meanlumerr^2)*(zdiff/tdiff)

  aguess = [-0.8, 2.0]
  xrange=[100., 2300.]
  yrange=[1.0, 20.0]
  ytitle = 'Density  # Mpc!U'+!tsym.minus+'2!N'

  erase & multiplot, [1,3], /square
  print,'r-band measurements'
  w=where(r.meanr GT 110.)
  ploterror, r.meanr[w], rdiff[w], rerr[w], psym=1, /xlog, /ylog,$
             xrange=xrange, yrange=yrange, /ystyle,/xstyle,$
             ytitle=ytitle
  fitpower, r.meanr[w]/1000., rdiff[w], rerr[w], aguess, ryfit, rout
  oplot, r.meanr[w], ryfit

  print
  print

  multiplot
  print,'i-band measurements'
  w=where(i.meanr GT 110.)
  ploterror, i.meanr[w], idiff[w], ierr[w], psym=1, /xlog, /ylog,$
             xrange=xrange, yrange=yrange, /ystyle,/xstyle,$
             ytitle=ytitle
  fitpower, i.meanr[w]/1000., idiff[w], ierr[w], aguess, iyfit, iout
  oplot, i.meanr[w], iyfit

  print
  print

  multiplot
  print,'z-band measurements'
  w=where(z.meanr GT 110.)
  ploterror, z.meanr[w], zdiff[w], zerr[w], psym=1, /xlog, /ylog,$
              xrange=xrange, yrange=yrange, /ystyle,/xstyle,$
             ytitle=ytitle,xtitle=xtitle
  fitpower, z.meanr[w]/1000., zdiff[w], zerr[w], aguess, zyfit, zout
  oplot, z.meanr[w], zyfit

  multiplot,/reset

  ;; now get r0
  print
  print,'r-band r0'
  calc_r0, rout[0], 0.0, abs(rout[1]), 0.0, dens[2], 0.0, r_r0
  print,'i-band r0'
  calc_r0, iout[0], 0.0, abs(iout[1]), 0.0, dens[3], 0.0, i_r0
  print,'z-band r0'
  calc_r0, zout[0], 0.0, abs(zout[1]), 0.0, dens[4], 0.0, z_r0


END 
