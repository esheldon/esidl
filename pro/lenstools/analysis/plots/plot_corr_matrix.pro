PRO plot_corr_matrix, all, dops=dops, encap=encap

  ;; must have run inversion code first

  indir = '~/lensout/combstripe/comb/'
  all = mrdfits(indir+'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_xi_comoving_N1.fit',1)



;  begplot,name='~/tmp/test.ps',/color
;  loadct,13

  IF keyword_set(encap) THEN fend = '.eps' ELSE fend = '.ps'
  dir = '~/plots/'
  file_proj = dir + 'corr_matrix_deltasig'+fend
  file_3d   = dir + 'corr_matrix_xi'+fend

  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          begplot, name=file_proj, /encap, xsize=7, ysize=7
      ENDIF ELSE BEGIN 
          begplot, name=file_proj
      ENDELSE 
  ENDIF 

  xtitle = 'log(R [h!U'+!csym.minus+'1!N Mpc])'

  logr = alog10(all.meanr/1000.)
  xrange = [min(logr), max(logr)]

  tvim2, all.corr, /scale, xrange=xrange, yrange=xrange,$
         xtitle=xtitle,ytitle=xtitle

  tmp=prompt_kbrd()

  IF keyword_set(dops) THEN BEGIN 
      endplot
      IF keyword_set(encap) THEN BEGIN  
          set_bbox, file_proj, '%%BoundingBox: 10 30 470 445'
      ENDIF
  ENDIF 

  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          begplot, name=file_3d, /encap, xsize=7, ysize=7
      ENDIF ELSE BEGIN 
          begplot, name=file_3d
      ENDELSE 
  ENDIF 

  xtitle = 'log(r [h!U'+!csym.minus+'1!N Mpc])'
  logr = alog10(all.r3)
  xrange = [min(logr), max(logr)]
  tvim2, all.corrxi, /scale, xrange=xrange, yrange=xrange,$
         xtitle=xtitle,ytitle=xtitle

  tmp=prompt_kbrd()

  IF keyword_set(dops) THEN BEGIN 
      endplot
      IF keyword_set(encap) THEN BEGIN  
          set_bbox, file_3d, '%%BoundingBox: 10 30 483 445'
      ENDIF
  ENDIF 

  ;; Both together

  !p.multi=[0,0,2]

  logr = alog10(all.meanr/1000.)
  xrange = [min(logr), max(logr)]
  xtitle = 'log(R [h!U'+!csym.minus+'1!N Mpc])'

  tvim2, all.corr, /scale, range=[-1,1], xrange=xrange, yrange=xrange,$
         xtitle=xtitle,ytitle=xtitle

  xtitle = 'log(r [h!U'+!csym.minus+'1!N Mpc])'
  logr = alog10(all.r3)
  xrange = [min(logr), max(logr)]
  tvim2, all.corrxi, /scale, range=[-1,1], xrange=xrange, yrange=xrange,$
         xtitle=xtitle,ytitle=xtitle

;  print,max(all.corr),min(all.corr)
;  print,max(all.corrxi),min(all.corrxi)
;  endplot
  !p.multi=0

END 
