PRO correct_clustering, cfiles, files


  combine_like, cfiles, mean=cmean, err=cerr, input_denscont=input_denscont
  combine_like, files,  mean=mean, err=err

;  fac = 10.0/mean
;  mean = mean*fac
;  err = err*fac

  ;; correct
  ncf = n_elements(cfiles)
  nf  = n_elements(files)

  corr = float(nf)/ncf*(err/cerr)^2

  mean_corr = cmean*corr
  err_corr = cerr*corr

  print
  print,'Corrected: '+ntostr(mean_corr)+' '+!plusminus+' '+ntostr(err_corr)+$
        '  (Nsig = '+ntostr( (mean_corr-input_denscont)/err_corr )+')'
  
  nsig = 4.5
  cmind = cmean - nsig*cerr
  cmaxd = cmean + nsig*cerr

  mind_corr = mean_corr - nsig*err_corr
  maxd_corr = mean_corr + nsig*err_corr

  mind = min([cmind, mind_corr])
  maxd = max([cmaxd, maxd_corr])

  d = arrscl(findgen(1000), mind, maxd)

  clike = gaussprob(d, cmean, cerr)
  clike = clike/max(clike)

  like_corr = gaussprob(d, mean_corr, err_corr)
  like_corr = like_corr/max(like_corr)

  setup_mystuff

  xtitle = !deltaytitle
  ytitle = 'Likelihood'
  aplot, !gratio, d, clike, yrange=[0,1.2], xtitle=xtitle, ytitle=ytitle
  oplot, d, clike, color=!red
  oplot, d, like_corr, color=!blue

  oplot,[input_denscont, input_denscont], [0, 1000], color=!green

  IF !d.name EQ 'X' THEN charsize = 1 ELSE charsize = 0.7

  legend,['Correction Factor = '+ntostr(corr,5,/round), $
          'Corr. Mean = '+ntostr(mean_corr,6,/round)+!csym.plusminus+$
          ntostr(err_corr,6,/round)], $
         charsize=charsize
  legend, ['Uncorrected', 'Corrected'], line=[0,0], colors=[!red, !blue],$
          charsize = charsize, /right

END 
