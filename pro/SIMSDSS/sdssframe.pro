pro  sdssframe, ngal, nstar, fwhm, aratio, posangle, image, gals, stars, $
                psf ,plt=plt

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME: sdssframe
;       
; PURPOSE: simulate a SDSS field
;	
;
; CALLING SEQUENCE:
;      sdssframe, ngal, nstar, fwhm, aratio, posangle, image, gals, stars, 
;                 psf ,plt=plt
;
; INPUTS: ngal:   The number of galaxies
;         nstar:  The number of stars
;         fwhm:   The fwhm of the psf gaussian in arcseconds
;         aratio: The axis ratio of the psf gaussian
;         posagnle: The position angle of the psf gaussian
;       
; OUTPUTS: image:  The simulated image = fltarr(2048,1489)
;          gals:   Structure holding information about the galaxies
;          stars:  Structure holding info about stars
;          psf:    Structure holding info about the psf
;
; INPUT KEYWORD PARAMETERS:
;          plt:  if /plt, then rdis is used to display the image
; 
; PROCEDURE: 
;	Uses empirically determined distributions for sizes and magnitudes.
;       Galaxies are exponential disks.
;       Begins with image on a grid twice as fine.  Places stars and galaxies.
;         Then rebins the image to (2048,1489)
;
; RECOMMENDATIONS:
;      ngal: ~
;      nstars: ~	
;
; REVISION HISTORY:
;	Erin Scott Sheldon   Umich  6/15/99
;       
;                                      
;-                                        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() EQ 0 THEN BEGIN
	print,'-Syntax: sdssframe, ngal, nstar, fwhm, aratio, posangle,  [image, gals, stars, psf, plt=plt]'
	return
  ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; some parameters
;;;; NOTE:  numbers followed by ;--- should be played with
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  s=systime(1)
  PI = !dpi
  COMMON seed, seed

; SDSS frames are 2048 X 1489, 0.4 arcseconds per pixel  
; We begin with a grid twice as fine

  sloanx = 2048
  sloany = 1489
  binfactor=2   ;---

; Must use binfactor^2 times the number of counts required for a 
; given magnitude because the image is rebinned.

  magfactor = binfactor^2
  gridx = binfactor*sloanx
  gridy = binfactor*sloany
  image = fltarr(gridx, gridy)

; make sky about same as run259 

  sky = 186.0  ;---
  gain = 5.25  ;--- this is a good average for red images.
  mzero = 30.0 ;--- Magnitude zero point.

; PSF parameters

  cfac = 2.*sqrt(2.*alog(2))  ;; fwhm = cfac*sigma for gaussian
  psfactor = 2*3.0  ;---      ;; make KERNEL 3 sigma in RADIUS => 6 sigma wide

  psfsigma = fwhm/cfac ;; psf scale length in arcseconds.
  psfsigma = psfsigma/.2       ;; in pixels.  .2 arcsec/pixel in fine grid
  psfsize = psfsigma*psfactor   ;; psf KERNEL size
  psfe1 = (1. - aratio^2)/(1. + aratio^2)*cos(2*posangle*PI/180.0)
  psfe2 = (1. - aratio^2)/(1. + aratio^2)*sin(2*posangle*PI/180.0)

  print,'------------------------------------------'
  print,'Seeing fwhm: ',strtrim(string(fwhm),2),' arcseconds'
  print,'Psf axis ratio: ',strtrim(string(aratio),2)
  print,'Psf position angle: ',strtrim(string(posangle),2), ' degrees'
  print,'Psf e1: ',strtrim(string(psfe1),2)
  print,'Psf e2: ',strtrim(string(psfe2),2)
  print,'------------------------------------------'

; Galaxy images will be [ galfactor*r0, galfactor*r0 ] in size
  galfactor = 20.0  ;---


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; make a structure to hold info about each galaxy, star and the psf
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  galstr = create_struct(name='ggalstr', $
                          'x', 0, $
                          'y', 0, $
                          'r0', 0.0, $
                          'mag', 0.0, $
                          'aratio', 0.0, $
                          'posangle', 0.0, $
                          'e1',0.0,$
                          'e2',0.0)
  starstr = create_struct(name='sstarstr', $
                          'x', 0, $
                          'y', 0, $
                          'mag', 0.0, $
                          'fwhm',0.0, $
                          'aratio',0.0,$
                          'posangle',0.0,$
                          'e1',0.0,$
                          'e2',0.0)
  psf = create_struct(name='psfstr', $
                          'fwhm',0.0, $
                          'sigma',0.0, $
                          'posangle',0.0, $
                          'aratio', 0.0, $
                          'e1', 0.0, $
                          'e2', 0.0)

  gals = replicate(galstr, ngal)
  stars = replicate(starstr, nstar)

;                          ---- GALAXIES ----
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;  Generate the magnitudes for each galaxy using a predetermined 
;;;;  emperical distribution based on actual sloan data (in maghist)
;;;;
;;;;  Then sample the distribution using genrand.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  galmaghist, rhist, rmag
  genrand, rhist, rmag, ngal, mag

  gals.mag = mag

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;  Generate sizes in a way similar to the magnitudes
;;;;
;;;;  NOTE: This outputs sqrt( ixx + iyy) from adaptive moments.  It must
;;;;        be scaled appropriately.  It is in pixels.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  galsizehist, yhist, xhist

  genrand, yhist, xhist, ngal, sizes

  ;;; first try at the scaling:
  gals.r0 = (binfactor-1)*sizes

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;  generate positions (uniformly), scale factors (uniformly), 
;;;;  axis ratios (sin(pi/2 * randomu ) ), position angles(uniformly)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR i=0,ngal-1 DO BEGIN 
    
    ;; generate positions on sloan sized image.  Multiply by 2 when placing.
    gals[i].x = fix( sloanx * randomu(seed) )
    gals[i].y = fix( sloany * randomu(seed) )

    ;;; Give it a chance to be really big and bright
;   if (gals[i].mag le 17. and gals[i].r0 ge .85) then begin
;     print,'Making a big one'
;     gals[i].r0 = 3.0 + (10.0 - 3.0)*randomu(seed)
;     gals[i].mag = 14 + randomn(seed)
;   endif
    gals[i].aratio = sin( randomu(seed)*PI/2.0 )
    gals[i].posangle = PI*randomu(seed)
    e0 = (1 - gals[i].aratio^2)/(1 + gals[i].aratio^2)
    gals[i].e1 = e0*cos(2*gals[i].posangle)
    gals[i].e2 = e0*sin(2*gals[i].posangle)
  ENDFOR

;                          ---- STARS ----
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;  Generate magnitudes for stars in similar way to galaxies
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  starmaghist, rhist, rmag
  genrand,rhist, rmag, nstar,mag

  stars.mag = mag

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;  Generate positions for each stars uniformly
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR i=0, nstar-1 DO BEGIN 
    stars[i].x = fix( sloanx * randomu(seed) )
    stars[i].y = fix( sloany * randomu(seed) )
  ENDFOR 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Put galaxies down on the image
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'------------------------------------------'
  print,'Adding Galaxies'
  numbad=0
  nwbad=0
  FOR i=0, ngal-1 DO BEGIN 

    counts = magfactor*10.0^( (mzero - gals[i].mag)/2.5 )
    r0 = gals[i].r0
    size = [ galfactor*r0, galfactor*r0 ] 

    ;; center on random point in center pixel to avoid wierd effects
    ;; factor of two is for finer grid
    cen = [(size[0]-1.)/2.+ (randomu(seed) -.5),$
           (size[1]-1.)/2. + (randomu(seed) -.5)]

    pa = gals[i].posangle*180.0/PI
    make_exp, disk, size, r0, counts=counts, theta=pa, $
                   aratio = gals[i].aratio, cen=cen
    
    j = indgen(size[0])
    k = indgen(size[1])

    ;;; positions on finer grid
    posx = binfactor*gals[i].x + (j - cen[0])
    posy = binfactor*gals[i].y + (k - cen[1] )
    wx = where(posx GE 0 AND posx LE gridx - 1, nwx)
    wy = where(posy GE 0 AND posy LE gridy - 1, nwy)
    IF (nwx NE 0 AND nwy NE 0) THEN BEGIN
      posx = posx[wx]
      posy = posy[wy]
      minposx = fix(min(posx))
      maxposx = fix(max(posx))
      minposy = fix(min(posy))
      maxposy = fix(max(posy))
      
      xnum=maxposx-minposx+1
      ynum=maxposy-minposy+1
      IF ( (nwx EQ xnum) AND (nwy EQ ynum) )THEN BEGIN
        image[minposx:maxposx, minposy:maxposy] = $
          image[minposx:maxposx, minposy:maxposy] + $
          disk[min(wx):max(wx), min(wy):max(wy)]
      ENDIF ELSE numbad = numbad+1
    ENDIF ELSE nwbad = nwbad+1
  ENDFOR

  IF (numbad NE 0) THEN print,'Number of incorrect sizes:  ',numbad
  IF (nwbad NE 0) THEN print,'Number of bad where statements: ',nwbad
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; put down stars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Adding Stars'
  print,'------------------------------------------'
  FOR i=0, nstar-1 DO BEGIN
    counts = magfactor*10.0^( (mzero - stars[i].mag)/2.5 )
    ;; use ;;; positions on finer grid
    image[binfactor*stars[i].x,binfactor*stars[i].y] = $
      image[binfactor*stars[i].x,binfactor*stars[i].y] + counts
  ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Define the parameters for psf (and hence stars)
;;;; Convolve the image with the psf
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'------------------------------------------'
  print,'Convolving with psf'
  print,'Kernel size is ',strtrim(string(fix(psfsize)),2),$
                          ' square on finer grid'

  psf.fwhm = fwhm  ;; arcseconds
  psf.sigma = fwhm/cfac  ;; arcseconds
  psf.aratio = aratio
  psf.posangle = posangle
  psf.e1 = psfe1
  psf.e2 = psfe2

  stars.fwhm = fwhm
  stars.aratio = aratio
  stars.posangle = posangle
  stars.e1 = psfe1
  stars.e2 = psfe2

  ;; center on random point in center pixel to avoid wierd effects
  ;; factor of two is for finer grid
  cen = [(psfsize-1.)/2.0 + (randomu(seed)-.5), $
         (psfsize-1.)/2.0 + (randomu(seed)-.5)]
  makegauss, psfim, [psfsize,psfsize], psfsigma, aratio=aratio, $
                 theta=posangle, cen=cen
  
  image = convol(image,psfim,/edge_truncate )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Add the sky, bin to sdss frame size, and add noise
;;;; Make sure noise is down by a factor of sqrt(gain)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
;;;  Right now this can make the wrong sky value.  However, since no 
;;;  negative values are added to the image up to this point, this shouldn't 
;;;  be a problem.

  print,'Adding sky value: ',strtrim(string(sky),2)
  min=min(image)
  image = image + ( min LE 0)*abs(min ) + sky
  
  print,'Rebinning Image'
  image = rebin(image,sloanx, sloany) 

  print,'Adding noise'
  ;;; make sure noise is down by factor of gain.
  add_noise, image, gain=gain

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; display if requested
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(plt) THEN BEGIN 
	rdis_setup,image, pls
	rdis,image,pls
  ENDIF 

  print
  print,'Execution was ',strtrim(string( (systime(1)-s)/60.0 ),2), ' min'
  print,'------------------------------------------'
return
END














