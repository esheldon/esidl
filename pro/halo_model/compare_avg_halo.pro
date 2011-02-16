PRO compare_avg_halo, struct, wzstruct

  ;; compare the co-moving averaged deltasig to the
  ;; physical averaged deltasig, and then corrected to 
  ;; co-moving based on the mean redshift

  avg_halo_out, struct, wzstruct, r,  deltasig,  xigm, meanz=meanz
  avg_halo_out, struct, wzstruct, rp, deltasigp, xigmp, /phys

  !p.multi=[0,0,2]

  xrange = [0.01, 10.0]

  ;; The physical values, corrected based on mean z
  rp_fix = rp*(1+meanz)
  deltasigp_fix = deltasigp/(1+meanz)^2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plot first the co-moving stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  setup_mystuff
  comcolor = !p.color
  IF !d.name EQ 'PS' THEN fixcolor = !blue ELSE fixcolor=!green
  physcolor = !red

  line1 = 0
  line2 = 2
  line3 = 1

  plot, r,deltasig,/xlog,/ylog, xrange=xrange, $
        xtitle=!mpcxtitle2, ytitle=!deltaytitle, $
        xtickf='loglabels', ytickf='loglabels'
  ;; overplot the physical, corrected to comoving
  oplot, rp_fix, deltasigp_fix,color=fixcolor, line=line2
  ;; also overplot the uncorrected, physical measurements
  oplot, rp, deltasigp,color=physcolor, line=line3

  legend, ['Comoving', 'Phys. Fix', 'Phys'], $
          line=[line1,line2,line3], $
          color=[!p.color, fixcolor, physcolor], /right, box=0,$
          thick=replicate(!p.thick,3)

  ;; Compare in ratios
  deltasig_int = interpol(deltasig, r, rp)
  deltasig_int_fix = interpol(deltasig, r, rp_fix)

  ratio = deltasig_int/deltasigp
  ratio_fix = deltasig_int_fix/deltasigp_fix
  
  ;;yrange = prange(ratio, ratio_fix, /noerror)
  yrange = [0.95, 1.05]
  plot, [0], [0], /nodata, /xlog, $
        xrange=xrange, yrange=yrange, xtitle=!mpcxtitle2, ytitle='ratio', $
        xtickf='loglabels'
  
  oplot,rp_fix,deltasig_int_fix/deltasigp_fix, $
        color=fixcolor, line=line2
  ;;oplot, rp, ratio, color=physcolor

  oplot, [0.01, 1000], [1,1]

  !p.multi=0

END 
