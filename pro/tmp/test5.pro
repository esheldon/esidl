PRO test5, norm, pow, yynorm, yypow, justplot=justplot, equal=equal

  yypow = [0.1, 0.2, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0,$
           2.7,3.0,3.5,4.0]
  npow = n_elements(yypow)
  yynorm = replicate(1.0, npow)

  ntrial = 1L

  IF keyword_set(justplot) THEN GOTO,jump

  pow=fltarr(npow, ntrial)
  norm=fltarr(npow, ntrial)

  FOR i=0L, npow-1 DO BEGIN 

      FOR j=0L, ntrial-1 DO BEGIN 

          test_binning, yynorm[i], yypow[i], tnorm, tpow, /nodisplay, equal=equal
          pow[i,j] = tpow & norm[i,j] = tnorm

      ENDFOR 

  ENDFOR 

jump:
  IF ntrial EQ 1 THEN BEGIN

      !p.multi=[0,0,2]
      IF keyword_set(equal) THEN BEGIN
          title='Equal Number Each Bin'
          name='test_binning_nonoise_equal.ps'
      ENDIF ELSE BEGIN
          title='Graded Binning'
          name='test_binning_nonoise.ps'
      ENDELSE 
      begplot,name=name

      xrange=[0,1.1*max(yypow)]
      yrange=[0,1.1*max(pow)]
      plot,yypow,pow,psym=4,xrange=xrange,yrange=yrange,$
        ytitle='Measured Power',xtitle='Input Power',title=title
      oplot,yypow,yypow

      yrange=[0,1.1*max(norm)]
      plot,yypow,norm,psym=4,xrange=xrange,yrange=yrange,$
        xtitle='Input Power',ytitle='Measured Norm',title=title
      legend,['Measured Norm','Input Norm'],psym=[4,0],/left,$
        thick=[!p.thick,!p.thick],box=0
      oplot,yypow,yynorm

      endplot 

  ENDIF ELSE BEGIN 

      !p.multi=[0,0,2]
      IF keyword_set(equal) THEN BEGIN
          title='Equal Number Each Bin'
          name='test_binning_equal.ps'
      ENDIF ELSE BEGIN
          title='Graded Binning'
          name='test_binning.ps'
      ENDELSE 
      begplot,name=name
      FOR i=0L, npow-1 DO BEGIN 
          
          plothist, pow[i, *], bin=0.01, xtitle='Power', xrange=[0,2],$
            title=title
          legend, 'Power = '+ntostr(yypow[i])
          oplot,[yypow[i],yypow[i]],[0,100000], line=2
          
          plothist, norm[i, *], bin=0.05, xtitle='Norm', xrange=[0,2],$
            title=title
          legend, 'Norm = '+ntostr(yynorm[i])
          oplot,[yynorm[i],yynorm[i]],[0,100000], line=2
          
      ENDFOR 

      endplot
  ENDELSE 

  !p.multi=0

END 
