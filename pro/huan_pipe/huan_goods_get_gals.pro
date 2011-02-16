PRO huan_goods_get_gals, struct, rcut, minmag, maxmag, $
                         wgood, wstar, wgal, r, med_seeing

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: huan_goods_get_gals, struct, rcut, minmag, maxmag, wgood, wstar, wgal, r, med_seeing'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; which aperture should we use
  ;; apertures in pixels: [5, 7, 10, 15, 20]
  ;; apertures in arcsec: 0.206 per pixel: [1.03, 1.442, 2.06, 3.09, 4.12]
    
  ;; cuts
  badflags = 2L^0 + 2L^2 + 2L^3 + 2L^4 + 2L^5 + 2L^6 + 2L^7
  max_fwhm = 2.0
  min_fwhm = 0.09

  arcperpix = 0.03              ; !!!

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

  mag = struct.mag_auto

  yt = 'FWHM (arcsec)'
  xt = 'i_AB (mag_auto)'
  mrange = [15, 30]
  frange = [0, 1]

  size2 = struct.ixx + struct.iyy
  afwhm = sqrt(size2/2.)*2.35*arcperpix

  aplot, 1, mag, afwhm, psym=3, /center, xmargin=xm, ymargin=ym, $
        xrange=mrange, yrange=frange, xtitle=xt, ytitle=yt

  w=where(struct.ixx GT 0)
  w2 = where(afwhm[w] LT 2)
  help,w,w2



  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Cuts:
  ;; flags we care about: skipping the deblend, 2L^1
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Checking cuts'
  help,where((struct.flags AND badflags) NE 0) 
  help,where(afwhm LT min_fwhm)
  help,where(afwhm GT max_fwhm)
  help,where(mag LT minmag)
  help,where(mag GT maxmag)
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

  wgood=where( (struct.flags AND badflags) EQ 0 AND $
               afwhm GT min_fwhm AND $
               afwhm LT max_fwhm AND $
               struct.rho4 GT 0 AND $
               mag GT minmag AND $
               mag LT maxmag AND $
               struct.whyflag EQ 0, nw, $
               comp=bad, ncomp=nbad)

  help,struct,wgood

  IF nw NE 0 THEN BEGIN 
      oplot, mag[bad], afwhm[bad], $
             psym=8, color=bad_color

  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; stars
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Really cut for stars
  star_maxmag = 23.0
  star_minmag = 20.5
  star_max_fwhm = 0.11
  star_min_fwhm = min_fwhm
  wstar=where(afwhm[wgood] LT star_max_fwhm AND $
              afwhm[wgood] GT star_min_fwhm AND $
              mag[wgood] LT star_maxmag AND $
              mag[wgood] GT star_minmag, nstar)


  wstar=wgood[wstar]
  help,wstar

  plot_box,star_minmag,star_maxmag,star_min_fwhm,star_max_fwhm,color=!red

  key=prompt_kbrd()

  ;; just plot the good ones
  aplot, 1, mag[wgood], afwhm[wgood], psym=3, /center, $
    xmargin=xm, ymargin=ym, $
    xrange=mrange, yrange=frange, xtitle=xt, ytitle=yt

  plot_box,star_minmag,star_maxmag,star_min_fwhm,star_max_fwhm,color=!red

  key=prompt_kbrd()

  ;; ellipticity
  ee = sqrt( struct.e1^2 + struct.e2^2 )
  e1 = struct.e1
  e2 = struct.e2
  momerr = struct.ellip_err

  ;; mean stellar size, to get rough value for rsmear
  psf_size2 = median(size2[wstar])
  med_seeing = median(afwhm[wstar])
  print
  print,'Median seeing: '+ntostr(med_seeing)
  print,'Median fwhm from cat: '+$
    ntostr(median(struct[wstar].fwhm_image)*arcperpix)

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
  legend,['All','Star','Gal'],lin=[0,0,0], $
    color=[!p.color,star_color,!red],/right,box=0

  plothist,e1[wgood],bin=bin,xtitle='e!D1!N', peak=1, charsize=charsize
  plothist,e1[wstar],bin=bin,/overplot,color=star_color, peak=1
  plothist,e1[wgal],bin=bin,/overplot,color=!red, peak=1

  plothist,e2[wgood],bin=bin,xtitle='e!D2!N', peak=1, charsize=charsize
  plothist,e2[wstar],bin=bin,/overplot,color=star_color, peak=1
  plothist,e2[wgal],bin=bin,/overplot,color=!red, peak=1

  fmax = 3.0
  plothist,afwhm[wgood],bin=bin,xtitle='adaptive FWHM [arcsec]', $
    peak=1, charsize=charsize, xrange=[0.0, 2.0], /xstyle
  plothist,afwhm[wstar],bin=bin, $
    /overplot,color=star_color, peak=1
  plothist,afwhm[wgal],bin=bin, $
    /overplot,color=!red, peak=1

  !p.multi=0

  !x.margin = xm
  !y.margin = ym

  key=prompt_kbrd()

  aplot, 1, mag[wgood], r[wgood], psym=3, xmargin=xm, ymargin=ym, $
        xrange=mrange, yrange=[0,2], charsize=charsize, $
        xtitle=xt, ytitle='rsmear', /center

  oplot, mag[wstar], r[wstar], psym=8, color=star_color, $
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
