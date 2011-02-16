PRO wtheta_conv_old, type, clr=clr

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: wtheta_conv, clr=clr'
      print,' type=0 (high dens) type=1 (low dens) type=2 (tot)'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  Npts = 1000L                  ;Num points to use in convolution integral

  colors = ['u','g','r','i','z']
  dir = '/sdss4/data1/esheldon/GAL_GAL/DENSITY/'
  outdir = '/sdss4/data1/esheldon/TMP/conv/'
  IF n_elements(clr) EQ 0 THEN clr=2

  typestr = ['dense','low', 'tot']
  IF type EQ 0 THEN BEGIN
      print,'Doing High Density Stuff'
      file = dir + 'wtheta_dense_'+colors[clr]+'.dat'
      filtername = 'wtheta_dense_function'
  ENDIF ELSE IF type EQ 1 THEN BEGIN
      print,'Doing Low Density Stuff'
      file = dir + 'wtheta_low_'+colors[clr]+'.dat'
      filtername = 'wtheta_low_gaussian'
  ENDIF ELSE IF type EQ 2 THEN BEGIN
      print,'Doing tot'
      filtername = 'wtheta_tot_function'
  ENDIF ELSE BEGIN 
      print,'type = 0, 1, or 2'
      return
  ENDELSE 
  print


  Ssh = 1.13                    ;mean shear polarizeability
  fac = 1/3600.                 ;convert from #/(arcmin)^2 to arcsec

  t=systime(1)
  defsysv,'!conv_X0',exists=exists1
  IF (NOT exists1) THEN BEGIN
      defsysv,'!conv_X0',0.
      defsysv,'!conv_Y0',0.
      defsysv,'!conv_rmax',0.
      defsysv,'!conv_sigma',0.
      defsysv,'!conv_cutoff',0.
      defsysv,'!zsource',0.
      defsysv,'!zlens',0.
  ENDIF 

  IF type EQ 2 THEN BEGIN
  fff = '/sdss4/data1/esheldon/GAL_GAL/gal_gal_752_756_'+colors[clr]+'_N1.dat'
      rdobjshear,fff,str
  ENDIF ELSE BEGIN 
      rdwtheta, file, str
  ENDELSE 
  w=where(str.shear NE 0., nw)
  shear = str[w].shear
  shearerr = str[w].shearerr
  meanr = str[w].meanr

  w=where(shear GT 0.)
  shear[w] = shear[w]*Ssh

  IF clr NE 2. THEN GOTO,jump

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; inputs for functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; functions called by conv2d
  func = 'call'+typestr[type]
  pqlim = 'pq_limits'

  !conv_sigma = 1.              ;km/s

  frmax = max(meanr)
  !conv_rmax = 1200.

  !zsource = .4
  !zlens = .15
;  !zsource = .285                 ;makes sigcrit right...
;  !zlens = .172
  ab_limits = [-!conv_rmax, !conv_rmax]
  sigma_crit = sigmacrit(!zsource, !zlens, h=1.)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; where to evaluate
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  minx = 0.
  stepx = 10.
  nx = 61L
  xx = dblarr(nx)
  FOR i=0, nx-1 DO xx[i] = minx + i*stepx

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Values of cutoff radius
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  cutmin = 50.
  cutstep = 25.
  ncut = 35L

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot stuff and file stuff 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  message = '#     radius        kappa      kappa(<r)'

  pfile=outdir+'conv_'+typestr[type]+'_'+colors[clr]+'_plots.ps'
  title=typestr[type]+' regions'
  xtitle='radius (arcsec)'
  ytitle='shear(neighbors)'
  begplot,name=pfile

;; don't forget to change!
 

  FOR icut = 0, ncut-1 DO BEGIN

      !conv_cutoff = cutmin + icut*cutstep

      print
      print,'Using: sigma = ',ntostr(!conv_sigma)
      print,'Using: cutoff['+ntostr(icut)+'] = '+ntostr(!conv_cutoff)

      conv2d, func, pqlim, ab_limits, xx, Npts, tmpfunc

      tmpfunc = tmpfunc*fac
      intfunc, tmpfunc*2.*!pi*xx, xx, tmpfunc_int, /double
      outfunc_int = tmpfunc_int[1:nx-1]/!pi/xx[1:nx-1]^2

      outfunc = tmpfunc[1:nx-1]

      outfile = outdir + 'conv_'+typestr[type]+'_'+colors[clr]+'_cut'+$
        ntostr(long(!conv_cutoff),5)+'.txt'
      
      linestyle=icut MOD 2
      IF icut EQ 0 THEN BEGIN 
          IF (type EQ 0) OR (type EQ 2) THEN plot, xx[1:nx-1],$
            (outfunc_int - outfunc)*150.^2,$
            yrange=[0.,.0015],ystyle=1, line=linestyle, $
            title=title, ytitle=ytitle,xtitle=xtitle

          IF type EQ 1 THEN plot, xx[1:nx-1],(outfunc_int - outfunc)*150.^2,$
            yrange=[-.0004,0.],ystyle=1, line=linestyle,$
            title=title, ytitle=ytitle,xtitle=xtitle
      ENDIF ELSE BEGIN 
          IF (type EQ 0) OR (type EQ 2) THEN oplot, xx[1:nx-1],$
            (outfunc_int - outfunc)*150.^2,line=linestyle

          IF type EQ 1 THEN oplot, xx[1:nx-1],(outfunc_int - outfunc)*150.^2,$
            line=linestyle
      ENDELSE 

      openw, lun, outfile, /get_lun
      !textunit = lun
      printf, lun, '#cutoff_radius= ',!conv_cutoff
      printf, lun
      printf,lun,message
      forprint, TEXT=5, xx[1:nx-1], outfunc, outfunc_int, /silent
      free_lun, lun


  ENDFOR 
  endplot,/noprint

  jump:
  message='#     meanr     shear       shearerr'
  outfile1 = outdir + 'shear_'+typestr[type]+'_'+colors[clr]+'.txt'

  openw, lun, outfile1, /get_lun
  !textunit = lun

  printf,lun, message
  forprint, TEXT=5, meanr, shear, shearerr, /silent
  free_lun, lun

  ptime,systime(1)-t
return
END 
