PRO huan_cfh_combine_plot

  dir = '~/Huan/cfh/sensitivity/'
  
  file = dir + 'sensitivity.dat'

  format = 'F,I,I,F,F,F,F,F'

  readcol, file, seeing, numgal, numgal_weighted, $
           raw_density, weighted_density, effective_density, $
           shear_sens_shonly, shear_sens, $
           format=format

  arcperpix = 0.206
  area = (2044.*4093)*arcperpix^2/60.0^2*12
  survey_area = 5000            ; square degrees
  survey_area_arcmin = survey_area*60.*60. ; square arcminutes

  begplot,name='~/plots/huan_sensitivity_effdens.ps',/color

  !p.multi=[0,0,2]

  line2 = 2
  color2 = !blue

  title = 'Shear Sensitivity'
  xt = 'seeing'
  yt = !csym.sigma+'('+!csym.gamma+') '+$
       !csym.times+' '+!csym.sqrt+'(5000 deg!U2!N/area) / 10!U'+!csym.minus+'5!N'
  survey_fac = sqrt(area/survey_area_arcmin)*1.e5

  shs = shear_sens*survey_fac
  plot,  seeing, shs, xtit=xt, ytit=yt, title=title
  oplot, seeing, shs, color=!blue

  w=where(seeing EQ 0.9)
  oplot, [0, seeing[w]], [shs[w], shs[w]], line=2
  oplot, [seeing[w], seeing[w]], [0, shs[w]], line=2

  title = 'Effective Density for '+!csym.sigma+'!Dsh!N = 0.32'
  xt = 'seeing'
  yt = '#/arcmin!U2!N'
  plot, seeing, effective_density, xtit=xt, ytit=yt, title=title,$
        yrange=[0,25]
  oplot,  seeing, effective_density, color=!blue

  oplot, [0, seeing[w]], [effective_density[w], effective_density[w]], line=2
  oplot, [seeing[w], seeing[w]], [0, effective_density[w]], line=2

  !p.multi=0

  endplot



END 
