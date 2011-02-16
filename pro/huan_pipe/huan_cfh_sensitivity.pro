PRO huan_cfh_sensitivity_shapenoise, struct, wgal

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now lets calculate stuff about the galaxies.  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Get shape noise from bright ones
  ;; look for low errors vs. magnitude
  magbin = 0.5
  binner, struct[wgal].mag_auto, struct[wgal].ellip_err, magbin, xo, yo, sig,$
          rev

  !p.multi=[0,0,2]

  ploterror, xo, yo, sig, psym=8

  nbin = n_elements(xo)
  shapenoise = fltarr(nbin)

  ntot = 0L
  FOR i=0L, nbin-1 DO BEGIN 
      IF rev[i] NE rev[i+1] THEN BEGIN 
          ww2 = rev[ rev[i]:rev[i+1]-1 ]

          ww2 = wgal[ww2]
          add_arrval, ww2, ww

          nw2 = n_elements(ww2)
          ntot = ntot + nw2

          angles = arrscl( randomu(seed, nw2), 0.0, 2*!pi )
          rotate_e1e2, angles, struct[ww2].e1, struct[ww2].e2, e1out, e2out

          add_arrval, e1out, e1

          shapenoise[i] = stdev( e1 )
      ENDIF 
  ENDFOR 

  plot, xo, shapenoise
  !p.multi=0

END 

PRO huan_cfh_sensitivity, struct, wgood, wstar, wgal, r

  plot_dir = '~/plots/'
  dir = '~/Huan/cfh/sensitivity/'

  ssh = 0.873

  ;; area of images (arcmin^2)
  arcperpix = 0.206
  area = (2044.*4093)*arcperpix^2/60.0^2*12
  survey_area = 5000            ; square degrees
  survey_area_arcmin = survey_area*60.*60. ; square arcminutes

;  IF n_elements(wgal) EQ 0 THEN BEGIN 
  rcut = 0.9

  minmag = 18.0
  maxmag = 26.0
  magstr = ntostr(maxmag,4)



;  begplot,name=plot_dir+'huan_gals.ps',/color
  huan_cfh_get_gals, $
    rcut, minmag, maxmag, struct, wgood, wstar, wgal, r, med_seeing
  key=prompt_kbrd()
;  endplot
      
;  ENDIF 

;  huan_cfh_sensitivity_shapenoise, struct, wgal
;return

  datfile = dir + 'sensitivity.dat'
;  begplot,name=plot_dir+'huan_sensitivity.ps',/color

  ;; get this for seeing in original images
  seeing_old = med_seeing


  ;; from above I got shapenoise = 0.25.  Interesting...

;  shapenoise = 0.25
  shapenoise=0.32
  shstr = ntostr(shapenoise, 4)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; New error
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ellip_err = struct.ellip_err/(1. - r)
  e1 = struct.e1/(1. - r)
  e2 = struct.e2/(1. - r)
  weights = 1./( shapenoise^2 + ellip_err[wgal]^2 ) ; [wgal]

  ;; save these.  We will add extra noise to these dilution
  ;; corrected ellipticities based on the new seeing

  e1old = e1
  e2old = e2
  ellip_err_old = ellip_err     ; all objects

  ;; raw galaxy density
  ngal = n_elements(wgal)
  raw_density = ngal/area

  ;; Now calculate the effective density
  weights_scaled = weights/max(weights)

  ngal_weighted = round( total(weights_scaled) )
  weighted_density = ngal_weighted/area

  print,'Raw density: ',raw_density
  print,'Weighted density: ',weighted_density

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Need to actually add error to the original ellipticities in
  ;; to account for this, otherwise the shear sensitivity will be
  ;; too good, despite our best attemptes at weighting.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wmom, e1[wgal], blah, we1mean, we1sig, we1err, inputweights = weights
  wmom, e2[wgal], blah, we2mean, we2sig, we2err, inputweights = weights

  ;; This is shear sensitivity for this field.
  shear_sens = 0.5*mean([we1err, we2err])/ssh

  shear_sens_shonly = 0.5*shapenoise/sqrt(ngal)

  print,'Shear sensitivity = '+ntostr(shear_sens)
  print,'Just using shapenoise: '+ntostr(shear_sens_shonly)
  ;; just using shapenoise, what is the effective source density?
  ;; shear sensitivity = 0.5*shapenoise/sqrt(density*area)
  effective_density = (shapenoise*0.5/shear_sens)^2/area
  print,'Effective density for sh='+shstr+': '+ntostr(effective_density)

  ;; Make a plot
  magbin = 0.25
  binner, struct[wgal].mag_auto, weights_scaled, magbin, xo, yo, sig,$
          rev, hist

  xt = 'IAB'
  yt = 'weight'
  xrange = [minmag,maxmag]
  aploterror, 1, xo, yo, sig, psym=8, yrange = [0, 1.2], $
              xrange=xrange, /xsty, xtit=xt, ytit=yt
  oplot, [xo, max(xo)+magbin/2.],  $
         [hist/float(max(hist))/2.,0.0],psym=10
;stop
  ;;legend,'IAB < '+magstr, box=0

  key=prompt_kbrd()

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now increase the seeing
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  seeing = $
    [0.650, 0.675, 0.700, 0.725, 0.750, 0.775, $
     0.800, 0.825, 0.850, 0.875, 0.900, 0.925, 0.950, 0.975, $
     1.000, 1.025, 1.050, 1.075, 1.100, 1.125, 1.150, 1.175, 1.20]

  seeing = [seeing_old, seeing]
  ns = n_elements(seeing)

  numgal            = replicate(ngal, ns)
  numgal_weighted   = replicate(ngal_weighted, ns)
  shear_sens        = replicate(shear_sens, ns)
  shear_sens_shonly = replicate(shear_sens_shonly, ns)
  weighted_density  = replicate(weighted_density, ns)
  effective_density = replicate(effective_density, ns)
  raw_density       = replicate(raw_density, ns)

  openw, lun, dir+'weights_vs_seeing.dat', /get_lun

  nrand=10
  FOR i=1L, ns-1 DO BEGIN 

      seeing_new = seeing[i]

      ;; ixx+iyy for this seeing
      size2_psf_new = 2.*(seeing_new/2.35/arcperpix)^2

      ;; The old ixx+iyy
      size2_psf_old = 2.*(seeing_old/2.35/arcperpix)^2

      ;; Same for objects. Correct image size for new seeing
      ;; assuming gaussians

      size2_old       = struct.ixx + struct.iyy
      size2_preseeing = size2_old - size2_psf_old

      size2_new = size2_preseeing + size2_psf_new
      
      ;; rsmear (not including 4th order) for each object in new seeing
      rr = size2_psf_new/size2_new

      ;; rsmear for old seeing
      rr_old = size2_psf_old/size2_old

      ;; How to scale the old r (again, not including 4th order)
      rscale = rr/rr_old

      r_rescaled = r*rscale

      print,seeing_old,seeing_new

      wgal2 = where( r_rescaled[wgal] LT rcut AND $
                     r_rescaled[wgal] GT 0, ngal2)
      wgal2 = wgal[wgal2]

      numgal[i] = ngal2

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; first try to account for more noise from seeing
      ;; This is due to increase in area of object due to seeing.
      ;; Only correct for sky-limited.  Also, the new dilution
      ;; correction will boost the noise as well
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      errscale2 = $
        ( size2_preseeing + size2_psf_new )/(size2_preseeing + size2_psf_old)
      errscale = sqrt(errscale2)

      ellip_err = $
        struct[wgal2].ellip_err*errscale[wgal2]/(1.-r_rescaled[wgal2])

      weights = 1./( shapenoise^2 + ellip_err^2 )
      weights_scaled = weights/max(weights)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Need to actually add error to the original ellipticities in
      ;; to account for this, otherwise the shear sensitivity will be
      ;; too good, despite our best attemptes at weighting.
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; These both include the dilution correct. The new one also
      ;; includes the effects of larger area due to seeing

      errdiff = sqrt( (ellip_err^2 - ellip_err_old[wgal2]^2) > 0.0 )

      huan_sens_addnoise, $
        nrand, $
        errdiff, $
        e1old[wgal2], e2old[wgal2], weights, $
        we1err, we2err

      ;; This is shear sensitivity for this field.
      shear_sens[i] = 0.5*mean([we1err, we2err])/ssh
      effective_density[i] = (shapenoise*0.5/shear_sens[i])^2/area

      ;; Some other types of estimators
      shear_sens_shonly[i] = 0.5*shapenoise/sqrt(ngal2)

      raw_density[i] = ngal2/area

      ngal_weighted = round( total(weights_scaled) )
      numgal_weighted[i] = ngal_weighted

      weighted_density[i] = ngal_weighted/area

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Make some plots
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      binner, struct[wgal2].mag_auto, weights_scaled, magbin, xo2, yo2, sig2,$
              rev, hist

      nbin = n_elements(xo2)
      colprint, replicate(seeing_new, nbin), xo2, yo2, sig2, hist, lun=lun

      xt = 'IAB'
      yt = 'weight'
      xrange = [minmag,maxmag]
      aploterror, 1, xo2, yo2, sig2, psym=8, $
                  yrange = [0.0, 1.2], xrange=xrange, /xsty, xtit=xt, ytit=yt
      oplot, [xo2, max(xo2)+magbin/2.],  $
             [hist/float(max(hist))/2.,0.0],psym=10

      legend,'IAB < '+magstr+' seeing = '+ntostr(seeing_new), box=0


  ENDFOR 

  free_lun, lun

  key=prompt_kbrd()

  !p.multi=[0,0,2]

  wcolor = !blue
  wline = 2

  xt = 'seeing'
  yt = '#/arcmin!U2!N'
  plot, seeing, raw_density, xtit=xt, ytit=yt
  oplot, seeing, effective_density, color=wcolor, line=wline
  legend,['density','effective density'],line=[0,wline],$
         color=[!p.color,wcolor],$
         /right, thick=[!p.thick,!p.thick], box=0
  ;;legend,'IAB < '+magstr, box=0, /bottom

  title = '5000 square degrees'
  yt = !csym.sigma+'('+!csym.gamma+') / 10!U'+!csym.minus+'5!N'
  survey_fac = sqrt(area/survey_area_arcmin)*1.e5
  plot, seeing, shear_sens*survey_fac, xtit=xt, ytit=yt, title=title, $
    line=wline
  oplot, seeing, shear_sens*survey_fac, color=wcolor, line=wline
  oplot, seeing, shear_sens_shonly*survey_fac
  legend,['shapenoise '+shstr,'true sensitivity'],line=[0,wline],$
         color=[!p.color,wcolor], thick=[!p.thick,!p.thick], box=0

  ;;legend,'IAB < '+magstr, box=0, /bottom,/right

  !p.multi=0

;  endplot

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Print to files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  file = dir + "weight.dat"
  print,'Writing to file: ',file
  openw, lun, file, /get_lun
  printf,lun,'# seeing=0.6"'
  printf,lun,'#     mag         weight'
  printf,lun,'#----------------------------'
  colprint,xo,yo,lun=lun
  free_lun,lun

  file = dir + "weight.tex"
  print,'Writing to file: ',file
  openw, lun, file, /get_lun
  printf,lun,"\begin{deluxetable}{ccccc}"
  printf,lun,"  \tabletypesize{\small}"
  printf,lun,"  \tablecaption{Relative Weight\label{tab:weight}}"
  printf,lun,"  \tablewidth{0pt}"
  printf,lun,"  \tablehead{"
  printf,lun,"    \colhead{$~~~~~~$IAB$~~~~~~$} &"
  printf,lun,"    \colhead{$~~~~~~$weight$~~~~~~$} "
  printf,lun,"  }"
  printf,lun,"  \"
  printf,lun,"  \startdata"

  nmag = n_elements(xo)
  FOR i=0L, nmag-1 DO BEGIN 

      strng = $
        "  "+$
        ntostr(xo[i], 5, /round) + ' & '+$
        ntostr(yo[i],5,/round)
      IF i NE nmag-1 THEN strng = strng + ' \\ '
      printf,lun,strng

  ENDFOR 
  printf,lun,"  \enddata"
  printf,lun,"  \tablecomments{Relative Weight as a function of IAB mag, accounting for measurement error and dilution.  The seeing for this data is 0.627\arcsec.}"
  printf,lun,"\end{deluxetable}"



  free_lun, lun


  print,'Writing to file: ',datfile
  openw, lun, datfile, /get_lun
  printf,lun,'     seeing         numgal      wnumgal    density      wdensity     effdensity  shear_sens_sh  shear_sens'
  printf,lun,'-----------------------------------------------------------------------------------------------------------------------------------------'

  FOR i=0L, ns-1 DO printf, lun, seeing[i], numgal[i], numgal_weighted[i], raw_density[i], weighted_density[i], effective_density[i], shear_sens_shonly[i], shear_sens[i], format='((F0,:,1X),2(I0,:,1X),5(F0,:,1X))'
  free_lun, lun

  ;; Just columns now

  file = dir + 'eff_dens.dat'
  print,'Writing to file: ',file
  openw, lun, file, /get_lun
  printf, lun,'#  seeing   eff_dens  shear_sens/(1.e-5 sqrt(5000 deg^2)'
  printf, lun,'#----------------------------------------------------------'
  colprint, seeing, effective_density, shear_sens*survey_fac, lun=lun
  free_lun, lun

  file = dir + "eff_dens.tex"
  print,'Writing to file: ',file
  openw, lun, file, /get_lun
  printf,lun,"\begin{deluxetable}{ccccc}"
  printf,lun,"  \tabletypesize{\small}"
  printf,lun,"  \tablecaption{Effective number density\label{tab:sens}}"
  printf,lun,"  \tablewidth{0pt}"
  printf,lun,"  \tablehead{"
  printf,lun,"    \colhead{seeing} &"
  printf,lun,"    \colhead{$n$} &"
  printf,lun,"    \colhead{\arcsec} & "
  printf,lun,"    \colhead{arcmin$^{-2}$} &"
  printf,lun,"    \colhead{arcmin$^{-2}$} "
  printf,lun,"  }"
  printf,lun,"  \"
  printf,lun,"  \startdata"

  FOR i=0L, ns-1 DO BEGIN 

      strng = $
        "  "+$
        ntostr(seeing[i], 5, /round) + ' & '+$
        ntostr(effective_density[i],6,/round)
      IF i NE ns-1 THEN strng = strng + ' \\ '
      printf,lun,strng

  ENDFOR 
  printf,lun,"  \enddata"
  printf,lun,"  \tablecomments{Effective number density $n$ as a function of seeing and magnitude cut. Defined such that shear sensitivity $\equiv 0.32/\sqrt{n*area}$}"
  printf,lun,"\end{deluxetable}"
  free_lun, lun



END 


