PRO compare_combine, new_e1, new_e2, new_e1e1err, new_e1e2err, new_e2e2err, $
                     gnew_e1, gnew_e2, gnew_e1e1err, gnew_e1e2err, gnew_e2e2err

  ;; compare combined ellip with and without the g-band
  !p.multi=[0,0,2]

  n=n_elements(new_e)
  new_e = sqrt(new_e1^2 + new_e2^2)

  new_eerr_tot = sqrt( $
                       (new_e1/new_e)^2*new_e1e1err^2 + $
                       (new_e2/new_e)^2*new_e2e2err^2 + $
                       2.*new_e1*new_e2/new_e^2*sign(new_e1e2err)*new_e1e2err^2 $
                     )
  
  gnew_e = sqrt(gnew_e1^2 + gnew_e2^2)

  gnew_eerr_tot = sqrt( $
                       (gnew_e1/gnew_e)^2*gnew_e1e1err^2 + $
                       (gnew_e2/gnew_e)^2*gnew_e2e2err^2 + $
                       2.*gnew_e1*gnew_e2/gnew_e^2*sign(gnew_e1e2err)*gnew_e1e2err^2 $
                     )

  xrange=[0,2]
  bin=0.01
  ytitle='N'
  xtitle='Combined Ellipticity'
  title='Run: 1458'
  plothist, new_e, bin=bin, xtitle=xtitle, ytitle=ytitle, title=title,$
            xrange=xrange
  plothist, gnew_e, bin=bin, color=!green,/overplot

  legend,['Without '+!colorsp[1],'With '+!colorsp[1]],$
         color=[!p.color,!green],/right,line=[0,0],charsize=1

  ratio = new_eerr_tot/gnew_eerr_tot

  xrange=[0.5,2]
  xtitle = 'err/err with '+!colorsp[1]
  plothist, ratio, bin=0.01, min=0, max=2, $
            xtitle=xtitle,ytitle=ytitle,title=title, $
            xrange=xrange

  !p.multi=0

END 
