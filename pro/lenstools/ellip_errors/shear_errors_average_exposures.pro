PRO shear_errors_average_exposures

  IF !d.name EQ 'PS' THEN BEGIN 
      colors = [!p.color, !darkGreen, !red, !blue, !magenta]
  ENDIF ELSE BEGIN 
      colors = [!p.color, !green, !red, !dodgerBlue, !yellow]
  ENDELSE 

  shapenoise = 0.32
  names = !csym.sigma+'!Umeas!N = '+ ['0.2','0.3','0.4','0.5','0.6']
  sigmeas = [0.2,0.3,0.4,0.5,0.6]
  Nmeas = lindgen(10)+1


  yrange = [0.2,1.1]
  sig_gamma = !csym.sigma+'!D'+!csym.gamma+'!N'
  ytitle =  sig_gamma + '(1)/'+sig_gamma+'(N!DMeas!N)'
  xtitle = 'N!DMeas!N'

  errwrong = 1.0/sqrt(Nmeas)

  plot, $
    Nmeas, errwrong, line=1, yrange=yrange, xtitle=xtitle, ytitle=ytitle

  nsigmeas = n_elements(sigmeas)
  FOR i=0L, nsigmeas-1 DO BEGIN 

      base = sqrt( shapenoise^2 + sigmeas[i]^2 )
      err = sqrt( shapenoise^2 + sigmeas[i]^2/Nmeas )

      oplot, Nmeas, err/base, color=colors[i]

  ENDFOR 

  legend,names,color=colors,line=replicate(0,nsigmeas),/bottom,box=0,charsize=1

END 
