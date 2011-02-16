PRO rsmear_cuts_plotall, run, rerun, purity, clr, hirata=hirata

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: rsmear_cuts_plotall, run, rerun, purity, clr, hirata=hirata'
      return
  ENDIF 

  rsmear_cuts_files, run, rerun, purity, clr, fitfile, psfile, $
                     comb_psfile, comb_fitfile, hirata=hirata
  all = mrdfits(fitfile,1)

  rsmear_cuts_files, run, rerun, purity, clr, fitfile, psfile, camcol=1, hirata=hirata
  col1 = mrdfits(fitfile,1)

  rsmear_cuts_files, run, rerun, purity, clr, fitfile, psfile, camcol=2, hirata=hirata
  col2 = mrdfits(fitfile,1)

  rsmear_cuts_files, run, rerun, purity, clr, fitfile, psfile, camcol=3, hirata=hirata
  col3 = mrdfits(fitfile,1)

  rsmear_cuts_files, run, rerun, purity, clr, fitfile, psfile, camcol=4, hirata=hirata
  col4 = mrdfits(fitfile,1)

  rsmear_cuts_files, run, rerun, purity, clr, fitfile, psfile, camcol=5, hirata=hirata
  col5 = mrdfits(fitfile,1)

  ;; file name
  rsmear_cuts_files, run, rerun, purity, clr, fitfile, psfile, camcol=6, $
                     hirata=hirata
  col6 = mrdfits(fitfile,1)

  meanmag = all.meanmag

  begplot, name=comb_psfile,/color
  
  simpctable
  colors = [!p.color, !blue, !cyan, !green, !magenta, !red, !orange]

  xtit = 'r!Dpetrosian!N'
  ytit = 'R!Dsmear!N['+!colors[clr]+'] Cut'
;  aplot, 1, meanmag, all.rsmear_cuts, xtit=xtit,ytit=ytit,/ynozero, $
;         yrange=[0,1], xrange=[18,22]
  aplot, 1, [0],/nodata, xtit=xtit,ytit=ytit,/ynozero, $
         yrange=[0,1], xrange=[18,22]
  oplot, meanmag, col1.rsmear_cuts, color=colors[1]
  oplot, meanmag, col2.rsmear_cuts, color=colors[2]
  oplot, meanmag, col3.rsmear_cuts, color=colors[3]
  oplot, meanmag, col4.rsmear_cuts, color=colors[4]
  oplot, meanmag, col5.rsmear_cuts, color=colors[5]
  oplot, meanmag, col6.rsmear_cuts, color=colors[6]

  ;; get error bar from spread
  nmag = n_elements(all.rsmear_cuts)
  rcuts_err = fltarr(nmag)
  FOR i=0L, nmag-1 DO BEGIN 
      vals = [col1.rsmear_cuts[i], col2.rsmear_cuts[i], col3.rsmear_cuts[i], col4.rsmear_cuts[i], col5.rsmear_cuts[i], col6.rsmear_cuts[i]]
      rcuts_err[i] = sdev(vals);/sqrt(6.)
  ENDFOR 

  myusersym, 'fill_circle'
  oplot, meanmag, all.rsmear_cuts
  oploterror, meanmag, all.rsmear_cuts, rcuts_err, psym=8


  mess = ['All', 'col1', 'col2', 'col3', 'col4','col5', 'col6']
  mess2 = ['Bayesian',$
           'Purity Cut: '+ntostr(purity,4,/round),$
           'Run: '+ntostr(run)+' Rerun: '+ntostr(rerun)]

  legend, mess, /right, colors=colors, $
          line=replicate(0,7), charsize=0.7, $
          thick=replicate(!p.thick,7)
        
  legend, mess2, /bottom, charsize=0.7

  aplot, 1, all.meanmag, all.ngalcumul, $
         tit=tit, xtit=xtit, $
         ytit='Cumulative # of "galaxies"   R!Dsmear!N['+!colors[clr]+'] Cut', charsize=1.25
  legend, mess2, charsize=0.7

  endplot 

END 
