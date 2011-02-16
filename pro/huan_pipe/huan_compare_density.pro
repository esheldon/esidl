PRO huan_compare_density

  dir = '~/Huan/'
  cfhfile = dir + 'cfh/sensitivity/sensitivity.dat'
  goodsfile = dir + 'goods/sensitivity/cfh_sensitivity_30.0.dat'

  format = 'F,I,I,F,F,F,F,F'

  readcol, cfhfile, $
    cseeing, cnumgal, cwnumgal, $
    cdensity, cwdensity, ceffdensity, cshear_sens_sh, cshear_sens, $
    format=format

  readcol, goodsfile, $
    gseeing, gnumgal, gwnumgal, $
    gdensity, gwdensity, geffdensity, gshear_sens_sh, gshear_sens, $
    format=format
  
  ;; convolved data
  readcol, dir + 'goods/sensitivity/newseeing_and_density.dat', $
    seeing_conv, med_seeing, eff_dens

  ;; THe new file.  This doesn't seem to match the one above!
;  readcol, dir + 'goods/sensitivity/huan_goods_sensitivity_convolved.dat', $
;    seeing_conv, med_seeing, eff_dens


  ng = n_elements(gseeing)
  seeing = [gseeing[1:ng-1], med_seeing]
  effdensity = [geffdensity[1:ng-1], eff_dens]

;  seeing = gseeing
;  effdensity=geffdensity

  IF !d.name EQ 'PS' THEN BEGIN 
      oclr = !blue
      symsize = 1
  ENDIF ELSE BEGIN 
      oclr = !green
      symsize=2
  ENDELSE 

  xtitle = 'Seeing'
  ytitle = 'Effective Density'
  aplot, 1, seeing, effdensity, psym=8, xtitle=xtitle,ytitle=ytitle,$
    xstyle=2
  oplot, med_seeing, eff_dens, psym=8, color=!red
  oplot, [cseeing[0]], [ceffdensity[0]], psym=2, $
    color=oclr,symsize=symsize  

  ;; Also plot what we think is the density if we could combine
  ;; the three bandpasses with equal S/N

  eff_dens3 =

  ;; lsst numbers
  oplot,[0.632907, 0.732573],[16.6345, 14.3847],$
    psym=4,color=!magenta

;  oplot,cseeing, ceffdensity,psym=8,color=!magenta

  legend,['CFHT image',$
          'Convolved HST + noise',$
          'HST Analytic Approx',$
          'LSST 3 years (convolved)'],$
    psym=[2,8,8,4],$
    color=[oclr, !red, !p.color,!magenta],$
    box=0,/right, thick=replicate(!p.thick,4)

  ts = interpol(eff_dens, med_seeing, 0.9)
  print,0.9,ts

;  conv_seeing = [0.624853, 0.732204, 0.832780, 0.932661]
;  conv_dens   = [6.72843,  6.58228,  5.72743,  4.62274]

  

END 
