PRO huan_cfh_get_gals, rcut, minmag, maxmag, struct, wgood, wstar, wgal, r, $
                       med_seeing

  huan_dir = '~/Huan/cfh/'
  IF n_elements(struct) EQ 0 THEN BEGIN 
      files = findfile(huan_dir+'cat/*admom*')

      mrdfits_multi, files, struct
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; which aperture should we use
  ;; apertures in pixels: [5, 7, 10, 15, 20]
  ;; apertures in arcsec: 0.206 per pixel: [1.03, 1.442, 2.06, 3.09, 4.12]
  
  ;; crude s/g separation
  star_cut = 0.9
  
  ;; cuts
  badflags = 2L^0 + 2L^2 + 2L^3 + 2L^4 + 2L^5 + 2L^6 + 2L^7
  max_size2 = 100
  min_size2 = 2.8

  ;;mi = 1

  arcperpix = 0.206
  apertures = [5,7,10,15,20]
  apertures_arcsec = apertures*arcperpix
  ;;print,apertures
  ;;print,apertures_arcsec
;  print,'Using aperture: '+$
;        ntostr(apertures[mi])+' pix  ('+ntostr(apertures_arcsec[mi])+"'')"

  IF !d.name EQ 'X' THEN BEGIN 
      charsize=2
      star_color = !green
      bad_color = !red
  ENDIF ELSE BEGIN 
      charsize=1
      star_color = !blue
      bad_color = !red
  ENDELSE 

  xm = [20,6]
  ym = [8,4]

  pxm = [10,3]
  pym = [4,2]

  !x.margin = xm
  !y.margin = ym

;  plothist, struct.class_star, bin=0.01, $
;            xtitle='Star_Class', position=aspect(1.0)

;  key=prompt_kbrd()
  
  yt = 'ixx + iyy'
;  xt = 'Mag_Aper['+ntostr(mi)+']'
  xt = 'Mag_Auto'
  mrange = [15, 27]
  frange = [0, 20]

  size2 = struct.ixx + struct.iyy

  aplot, 1, struct.mag_auto, size2, psym=3, xmargin=xm, ymargin=ym, $
        xrange=mrange, yrange=frange, xtitle=xt, ytitle=yt

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Cuts:
  ;; flags we care about: skipping the deblend, 2L^1
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Checking cuts'
  help,where((struct.flags AND badflags ) NE 0) 
  help,where(size2 LT min_size2)
  help,where(size2 GT max_size2)
  help,where(struct.mag_auto LT minmag)
  help,where(struct.mag_auto GT maxmag)
  help,where(struct.whyflag NE 0)
  print

  OBJECT2_AMOMENT_FAINT =  '200000'X ; too faint for adaptive moments 
  OBJECT2_AMOMENT_UNWEIGHTED = '200000'X ; failed so tried unweighted mom
  OBJECT2_AMOMENT_SHIFT =  '400000'X ; centre moved too far while
                                ;determining adaptive moments 
  OBJECT2_AMOMENT_MAXITER = '800000'X ; Too many iterations while

  help,where((struct.whyflag AND OBJECT2_AMOMENT_FAINT) NE 0)
  help,where((struct.whyflag AND OBJECT2_AMOMENT_UNWEIGHTED) NE 0)
  help,where((struct.whyflag AND OBJECT2_AMOMENT_SHIFT) NE 0)
  help,where((struct.whyflag AND OBJECT2_AMOMENT_MAXITER) NE 0)

  wgood=where( (struct.flags AND badflags ) EQ 0 AND $
               size2 GT min_size2 AND $
               size2 LT max_size2 AND $
               struct.rho4 GT 0 AND $
               struct.mag_auto GT minmag AND $
               struct.mag_auto LT maxmag AND $
               struct.whyflag EQ 0, nw, $
               comp=bad, ncomp=nbad)

  help,struct,wgood

  IF nw NE 0 THEN BEGIN 
      oplot, struct[bad].mag_auto, size2[bad], $
             psym=8, symsize=0.25, color=bad_color

  ENDIF 

  key=prompt_kbrd()

  ;; just plot the good ones
  aplot, 1, struct.mag_auto, size2, psym=3, $
    xrange=mrange, yrange=frange, xtitle=xt, ytitle=yt, $
    xmargin=xm, ymargin=ym

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; stars
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; This is most liberal perhaps?  Really only use this
  ;; for "wrest"
  wstartmp = where(struct[wgood].class_star GT star_cut, nstar, $
                   comp=wrest, ncomp=nrest)
  wstartmp = wgood[wstartmp]
  wrest    = wgood[wrest]

  ;; Really cut for stars
  star_maxmag = 23
  star_minmag = 19.5
  star_max_size2 = 3.9
  star_min_size2 = min_size2
  wstar=where(struct[wgood].class_star GT star_cut AND $
              size2[wgood] LT star_max_size2 AND $
              size2[wgood] GT star_min_size2 AND $
              struct[wgood].mag_auto LT star_maxmag AND $
              struct[wgood].mag_auto GT star_minmag, nstar)

  wstar=wgood[wstar]
  
  help,wstar

;  oplot, struct[wstartmp].mag_auto, size2[wstartmp], $
;         psym=8, symsize=0.25,$
;         color=star_color
  plot_box,star_minmag,star_maxmag,star_min_size2,star_max_size2,color=!red

  key=prompt_kbrd()

  ;; now the rest
  aplot, 1, struct[wrest].mag_auto, size2[wrest], psym=3, $
        xrange=mrange, yrange=frange, xtitle=xt, ytitle=yt

  key=prompt_kbrd()

  ;; ellipticity
  ee = sqrt( struct.e1^2 + struct.e2^2 )
  e1 = struct.e1
  e2 = struct.e2
  momerr = struct.ellip_err

  ;; mean stellar size, to get rough value for rsmear
  psf_size2 = median(size2[wstar])
  med_seeing = sqrt(psf_size2/2.)*2.34*arcperpix
  print
  print,'Median seeing: '+ntostr(med_seeing)

  psf_r4val = median(4./struct[wstar].rho4 - 1.)

  r = (psf_size2/size2)*psf_r4val/(4./struct.rho4 - 1.)
  wgal = where(r[wgood] LT rcut AND r[wgood] GT 0, ngal)
  wgal = wgood[wgal]
  help,wgal

  !p.multi=[0,2,2]

  !x.margin = pxm
  !y.margin = pym

  bin = 0.02
  plothist,ee[wgood],bin=bin,xtitle='ellipticity', peak=1, charsize=charsize
  plothist,ee[wstar],bin=bin,/overplot,color=star_color, peak=1
  plothist,ee[wgal],bin=bin,/overplot,color=!red, peak=1

  plothist,e1[wgood],bin=bin,xtitle='e!D1!N', peak=1, charsize=charsize
  plothist,e1[wstar],bin=bin,/overplot,color=star_color, peak=1
  plothist,e1[wgal],bin=bin,/overplot,color=!red, peak=1

  plothist,e2[wgood],bin=bin,xtitle='e!D2!N', peak=1, charsize=charsize
  plothist,e2[wstar],bin=bin,/overplot,color=star_color, peak=1
  plothist,e2[wgal],bin=bin,/overplot,color=!red, peak=1

  !p.multi=0

  !x.margin = xm
  !y.margin = ym

  key=prompt_kbrd()

  aplot, 1, struct[wgood].mag_auto, r[wgood], psym=3, $
        xrange=mrange, yrange=[0,2], charsize=charsize, $
        xtitle=xt, ytitle='rsmear'
;  oplot, struct[wstartmp].mag_auto, r[wstartmp], psym=8, color=!OrangeRed, $
;         symsize=0.25
  oplot, struct[wstar].mag_auto, r[wstar], psym=8, color=star_color, $
         symsize=0.25
  oplot, [0, 50], [1,1]
  oplot, [0, 50], [rcut, rcut],color=!red

  print,mean(ee[wstar]),mean(ee[wgal])
  print,mean(e1[wstar]),mean(e1[wgal])
  print,mean(e2[wstar]),mean(e2[wgal])
  print,median(r[wstar]),median(r[wgal])
  
  !x.margin = pxm
  !y.margin = pym

END 
