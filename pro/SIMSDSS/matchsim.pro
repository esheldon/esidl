PRO matchsim, sexcat, gals, stars, plt=plt, wgals=wg,wstars=ws,psfile=psfile

  IF n_params() EQ 0 THEN BEGIN
    print,'-Syntax: matchsim, sexcat, gals, stars, plt=plt, wgals=wg,wstars=ws'
    return
  ENDIF 

  sexcat.psf_fwhm = stars[0].fwhm
  ;; close_match parameters
  allow = 1  ;; keep closest match
  tol = 5.0  ;; pixels

  nc = n_elements(sexcat)
  ng = n_elements(gals)
  ns = n_elements(stars)
  print,'Number in sexcat: ',nc
  print,'Number of galaxies: ',ng
  print,'Number of stars: ',ns

  xvals = [gals.x, stars.x]
  yvals = [gals.y, stars.y]

  stype = [intarr(ng)+1, intarr(ns)+2]  ;; gal=1, star=2
  smag = [gals.mag, stars.mag]
  saratio = [gals.aratio, stars.aratio]
  sposangle = [gals.posangle, stars.posangle ]
  se1 = [gals.e1, stars.e1]
  se2 = [gals.e2, stars.e2]

;;;;; Use close_match to find stars and galaxies
  close_match, sexcat.x_image, sexcat.y_image, xvals, yvals,$
               msex, m, tol, allow

  sexcat[ msex ].stype = stype[m]
  sexcat[ msex ].smag = smag[m]
  sexcat[ msex ].saratio = saratio[m]
  sexcat[ msex ].sposangle = sposangle[m]
  sexcat[ msex ].se1 = se1[m]
  sexcat[ msex ].se2 = se2[m]

  print,'------------------------------------------'
  wg=where(sexcat.stype EQ 1,nwg)
  ws=where(sexcat.stype EQ 2,nws)
  print,strtrim(string(nwg),2),' Galaxies Found'
  print,strtrim(string(nws),2),' Stars Found'

  wws=where(sexcat[ws].mag_best LE 20.0,nwws)
  print,strtrim(string(nwws),2),' Stars with mag_best < 20.0'
  print,'------------------------------------------'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Calculate the polarization since we know the real stars
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; pick objects that are well measured and were matched
  good = where(sexcat.e1_ad NE -10.0 AND $
               sexcat.rho4 GT .1 AND $
               sexcat.stype NE -10 AND $
               sexcat.uncert_ad LT 1.0 $
               ,ngood)
  IF (ngood NE 0) THEN BEGIN
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;; set goodflag=1 for objects that can be corrected
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    sexcat[good].goodflag = 1

    ;; find good bright stars
    gs = where(sexcat[good].stype EQ 2 AND $
               sexcat[good].mag_best LT 20.0, nstars)
    IF nstars NE 0 THEN BEGIN
      gs = good[gs]
      psfrho4 = median( sexcat[gs].rho4 )
      psfsize = median( sexcat[gs].x2_ad + sexcat[gs].y2_ad )
      ;; Use gary's correction factor
      psfsize = psfsize*(4.0/psfrho4 -1)
      allsize = sexcat[good].x2_ad + sexcat[good].y2_ad
      allsize = allsize*(4.0/sexcat[good].rho4 -1)
      sexcat[good].sr = psfsize/allsize
    ENDIF 
  ENDIF 

  IF keyword_set(plt) THEN BEGIN
    
    galsym = 4
    starsym = 1
    gsymsize= .4
    ssymsize = .4

    xtitle='Mag_best'
    ytitle='sqrt( ixx + iyy )'
    plot, sexcat[wg].mag_best, sqrt(sexcat[wg].x2_ad+sexcat[wg].y2_ad), $
      xrange=[14,26],yrange=[0,10],$
      xtitle=xtitle,ytitle=ytitle,psym=galsym,symsize=gsymsize
  
    oplot, sexcat[ws].mag_best, sqrt(sexcat[ws].x2_ad+sexcat[ws].y2_ad), $
      psym=starsym,symsize=ssymsize
    legend,['Galaxies','Stars'],psym=[galsym,starsym]

    IF keyword_set(psfile) THEN BEGIN
      makeps,psfile
      plot, sexcat[wg].mag_best, sqrt(sexcat[wg].x2_ad+sexcat[wg].y2_ad), $
        xrange=[14,26],yrange=[0,10],$
        xtitle=xtitle,ytitle=ytitle,psym=galsym,symsize=gsymsize
  
      oplot, sexcat[ws].mag_best, sqrt(sexcat[ws].x2_ad+sexcat[ws].y2_ad), $
        psym=starsym,symsize=ssymsize
      legend,['Galaxies','Stars'],psym=[galsym,starsym]
      ep
    ENDIF 

  ENDIF

  return
END

















