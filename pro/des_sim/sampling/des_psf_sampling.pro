;+
; NAME:
;  DES_PSF_SAMPLING
;
;
; PURPOSE:
;  Test the effects of pixel size on PSF shape measurements.  Create images
;  and place PSF, run shape measurement code on them, and return a structure
;  containing statistics on the resulting measurements.
;
;  PSF's are placed on a grid and the centers are perturbed by a uniform random
;  number within the pixel to test noise.
;
;
; CATEGORY:
;  Generic Shape Measurements
;
;
; CALLING SEQUENCE:
;  des_psf_sampling, seeing_fwhm, arcperpix, etot, posAngle, nstarx, $
;                    struct, $
;                    doplot=doplot
;
;
; INPUTS:
;  seeing_fwhm: The seeing FWHM in arcseconds
;  arcperpix: Pixel diameter in arcseconds
;  etot: The total ellipticity sqrt(e1^2 + e2^2) of the PSF. Can be an array
;  posAngle: Position angle of the PSF.  Can be an array.
;  nstarsx: The number of stars in the grid in the x direction.  The same
;           number will be placed in the y direction.
;    
;
; OPTIONAL INPUTS:
;  gridfac: Intiial PSF is sampled at gridfac times the resolution of final.
;
;
; KEYWORD PARAMETERS:
;  /doplot:  Make a histogram of the measured e1/e2 for each etot and posAngle
;
;
; OUTPUTS:
;  struct:  Structure containing statistics of recovered shapes for each etot
;          and posAngle
;
;
; OPTIONAL OUTPUTS:
;
;
; MODIFICATION HISTORY:
;  Created: Nov-2004 Erin Sheldon, UofChicago
;
;-

PRO des_psf_sampling, seeing_fwhm, arcperpix, etot, posAngle, nstarx, $
                      struct, $
                      gridfac=gridfac, $
                      doplot=doplot

  IF n_params() LT 6 THEN BEGIN 
      print,'-Syntax: des_psf_sampling, seeing_fwhm, arcperpix, etot, posAngle, nstarx, $'
      print,'               struct, $'
      print,'               gridfac=gridfac'
      print,'               /doplot, $'
      return
  ENDIF 


  IF keyword_set(doplot) THEN BEGIN 
      !p.multi=[0,0,2]
      !p.charsize=1
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  save_st = create_struct('arcperpix', arcperpix, $
                          'seeing_fwhm', seeing_fwhm, $
                          'psf_sigma_pixels', 0.0, $
                          'e1in', 0.0, $
                          'e2in', 0.0, $
                          'e1meas', 0.0, $
                          'e2meas', 0.0, $
                          'e1sdev', 0.0, $
                          'e2sdev', 0.0, $
                          'ixxmeas', 0.0, $
                          'iyymeas', 0.0, $
                          'ixymeas', 0.0, $
                          'npsf', 0L, $
                          'npsf_use', 0L)

  ;; Image will have no sky or sky noise
  sky = 0.0
  skysig = 0.0

  ;; Number of PSF stars
  nstary = nstarx
  nstar = long(nstarx)*nstary

  cfac = 2.*sqrt(2.*alog(2))    ; fwhm = cfac*sigma for gaussian
  psf_sigma = seeing_fwhm/cfac       ; arcsec

  seeing_pixels = seeing_fwhm/arcperpix ; fwhm in pixels
  psf_sigma_pixels = psf_sigma/arcperpix ; sigma in pixels
  pixscale_over_seeing = arcperpix/seeing_fwhm

  save_st.psf_sigma_pixels = psf_sigma_pixels

  counts = 10000.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The PSF image
  ;; want 4 sigma around each star => 8 sigma sized box  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  min_psf_size = 10 ;; pixels
  psf_size = round( [8*psf_sigma_pixels, 8*psf_sigma_pixels] ) > min_psf_size

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The overall image.  There is no space between the PSF images
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nx = psf_size[0]*nstarx
  ny = psf_size[1]*nstary

  image = fltarr(nx, ny)

  ;; print out some info
  help,seeing_fwhm,psf_sigma
  help,arcperpix
  help,pixscale_over_seeing
  help,seeing_pixels,psf_sigma_pixels
  help,nstarx,nstary,nx,ny
  print,'PSF SIZE = ',psf_size


  ;; We will perturb from this default center
  default_center = [(psf_size[0]-1.)/2.,(psf_size[1]-1.)/2.]

  print,'Default Center = ',default_center
  psfCenterX = fltarr(nstar)
  psfCenterY = fltarr(nstar)

  Netot = n_elements(etot)
  NposAngle = n_elements(posAngle)
  Ntot = Netot*NposAngle
  struct = replicate(save_st, Ntot)


  print,'Total number of ellipticities: ',Ntot

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over the ellipticities and position angles
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ii = 0L
  FOR ie=0L, Netot-1 DO BEGIN

      FOR ipos=0L, NposAngle-1 DO BEGIN 

          ;; Convert to e1 and e2
          e1in = etot[ie]*cos(2.0*posAngle[ipos]*!pi/180.)
          e2in = etot[ie]*sin(2.0*posAngle[ipos]*!pi/180.)

          print,'-------------------------------------------------'
          print,'test # '+ntostr(ii+1)+'/'+ntostr(Ntot)
          print,'e1 = '+ntostr(e1in)
          print,'e2 = '+ntostr(e2in)

          ;; The PSF code takes in axis ratio and position angle
          findABTheta, e1in, e2in, aratio, theta
          theta = theta*180.0/!pi

          ;; Place the PSF images into the overall image
          index = 0L
          begx = 0L
          image[*] = 0
          FOR ix=0L, nstarx-1 DO BEGIN 
              
              begy = 0L
              FOR iy=0L, nstary-1 DO BEGIN 
                  

                  ;; move the center by random number
                  rand_offset = randomu(seed,2)-0.5
                  cen = default_center[*] + rand_offset[*]

                  makegauss, $
                    psf, psf_size, psf_sigma_pixels, $
                    aratio=aratio, theta=theta, $
                    counts=counts, cen=cen, gridfac=gridfac


                  image[begx:begx+psf_size[0]-1, $
                        begy:begy+psf_size[1]-1 ] = psf[*,*]
                  
                  psfCenterX[index] = begx + cen[0]
                  psfCenterY[index] = begy + cen[1]
                  
                  index = index + 1
                  
                  begy = begy+psf_size[1]
              ENDFOR 
              begx = begx + psf_size[0]
              
          ENDFOR 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Run the PSF measurement code: Adaptive Moments
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          ;; We are giving it a good guess
          wguess = replicate( psf_sigma_pixels^2, nstar )
          
          run_admom, $
            image, psfCenterX, psfCenterY, sky, skysig, wguess, $
            ixx, iyy, ixy, rho4, err, whyflag
          
          e1 = replicate(-1000., nstar)
          e2 = e1

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Measure some statistics for the results
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          w=where(whyflag EQ 0, ngood)
          help,nstar,ngood,float(ngood)/nstar
          e1[w] = (ixx[w]-iyy[w])/(ixx[w]+iyy[w])
          e2[w] = 2.*ixy[w]/(ixx[w]+iyy[w])

          struct[ii].e1in = e1in
          struct[ii].e2in = e2in
          struct[ii].e1meas = median(e1[w])
          struct[ii].e2meas = median(e2[w])
          struct[ii].e1sdev = sdev(e1[w])
          struct[ii].e2sdev = sdev(e2[w])

          struct[ii].ixxmeas = median(ixx[w])
          struct[ii].iyymeas = median(iyy[w])
          struct[ii].ixymeas = median(ixy[w])

          struct[ii].npsf = nstar
          struct[ii].npsf_use = ngood

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Plot the results if requested
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF keyword_set(doplot) THEN BEGIN 

              title = $
                'pixscale/seeing = '+ntostr(pixscale_over_seeing,5,/round)
              e1meas = struct[ii].e1meas
              e2meas = struct[ii].e2meas
              
              e1sdev = struct[ii].e1sdev
              e2sdev = struct[ii].e2sdev

              nsig = 4

              miny = 0.1

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; plot e1 histogram
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              plotrange = e1meas + [-nsig*e1sdev, nsig*e1sdev]
              plotlines = e1meas + [-e1sdev, e1sdev]
              binsize = 4.0*e1sdev/sqrt(nstar)

              plothist, e1[w], xhist, yhist, $
                bin=binsize, xrange=plotrange, /noplot
              yrange = [miny, max(yhist)]
              plothist, e1[w], xhist, yhist, $
                bin=binsize, xrange=plotrange, $
                xtitle='e1',/ylog, yrange=yrange, $
                title=title,ytickf='loglabels'

;              plothist, e1[w], xhist, yhist, $
;                bin=binsize, xrange=plotrange, $
;                xtitle='e1'

              oplot, [e1meas, e1meas], [miny, 10*max(yhist)], color=!blue
              oplot, [plotlines[0],plotlines[0]], [miny, 10*max(yhist)], $
                color=!red, line=2
              oplot, [plotlines[1],plotlines[1]], [miny, 10*max(yhist)], $
                color=!red, line=2
              legend,'sdev = '+ntostr(e1sdev),/clear

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; plot e2 histogram
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              plotrange = e2meas + [-nsig*e2sdev, nsig*e2sdev]
              plotlines = e2meas + [-e2sdev, e2sdev]
              binsize = 4.0*e2sdev/sqrt(nstar)

              plothist, e2[w], xhist, yhist, $
                bin=binsize, xrange=plotrange, /noplot
              yrange = [miny, max(yhist)]
              plothist, e2[w], xhist, yhist, $
                bin=binsize, xrange=plotrange, $
                xtitle='e2',/ylog,yrange=yrange,ytickf='loglabels'

              oplot, [e2meas, e2meas], [miny, 10*max(yhist)], color=!blue

              oplot, [plotlines[0],plotlines[0]], [miny, 10*max(yhist)], $
                color=!red, line=2
              oplot, [plotlines[1],plotlines[1]], [miny, 10*max(yhist)], $
                color=!red, line=2

              legend,'sdev = '+ntostr(e2sdev),/clear

          ENDIF 

          ii = ii+1
      ENDFOR ;; over position angles

  ENDFOR ;; over total ellipticities

  IF keyword_set(doplot) THEN !p.multi=0


END 

