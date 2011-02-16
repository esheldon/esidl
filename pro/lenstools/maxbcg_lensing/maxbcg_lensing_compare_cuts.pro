PRO maxbcg_lensing_compare_cuts, dops=dops

  IF keyword_set(dops) THEN BEGIN 
      psfile = '~/plots/MaxBCG/compare_source_cuts.ps'
      begplot,name=psfile,/color
  ENDIF 

  ngals_nbin = 12
  ngals_bin = lindgen(ngals_nbin)
  n = n_elements(ngals_bin)


  m11 = obj_new('maxbcg_lensing',11) ; standard


  m14 = obj_new('maxbcg_lensing',14) ; R > 1/3
  m15 = obj_new('maxbcg_lensing',15) ; R > 0.4
  m13 = obj_new('maxbcg_lensing',13) ; R > 0.585
  m16 = obj_new('maxbcg_lensing',16) ; R > 2/3

  simpctable, colorlist=clist

  FOR i=0L, n-1 DO BEGIN 

      bin = ngals_bin[i]

      t11_11 = m11->combined_get(ngals_nbin=ngals_nbin, ngals_bin=bin)
      t13_11 = m13->combined_get(ngals_nbin=ngals_nbin, ngals_bin=bin)
      t14_11 = m14->combined_get(ngals_nbin=ngals_nbin, ngals_bin=bin)
      t15_11 = m15->combined_get(ngals_nbin=ngals_nbin, ngals_bin=bin)
      t16_11 = m16->combined_get(ngals_nbin=ngals_nbin, ngals_bin=bin)
      
      plot_density_contrast, t11_11, /log, /mpc, aspect=1, $
        title = 'ngals bin = '+ntostr(bin)
      oplot, t11_11.meanr/1000*1.05, t11_11.sigma
      
;  oploterror, t13_11.meanr/1000*1.05, t13_11.sigma, t13_11.sigmaerr, $
;    psym=4, color=clist[1], errc=clist[1]
      
      oplot, t14_11.meanr/1000*1.05, t14_11.sigma, $
        color=clist[1]
      oplot, t15_11.meanr/1000*1.05, t15_11.sigma, $
        color=clist[2]
      oplot, t13_11.meanr/1000*1.05, t13_11.sigma, $
        color=clist[3]
;      oplot, t16_11.meanr/1000*1.05, t16_11.sigma, $
;        color=clist[4]
      
;      mess = ['Standard', 'R > 1/3', 'R > 0.4', 'R > 0.585','R > 2/3']
      mess = ['R > 1/4', 'R > 1/3', 'R > 0.4', 'R > 0.585']
 
      nmess = n_elements(mess)
      colors = clist[0:nmess-1]
      lines = replicate(0, nmess)
      legend, mess, $
        colors=colors, lines=lines, $
        /right, box=0, charsize=1

      IF i NE n-1 AND !d.name EQ 'X' THEN BEGIN 
          key=prompt_kbrd('hit a key')
          IF key EQ 'q' THEN return
      ENDIF 

  ENDFOR 
  obj_destroy, m11, m13, m14, m15, m16

  IF keyword_set(dops) THEN endplot


END 
