PRO wtheta_conv, type, cutmin=cutmin, cutstep=cutstep, ncut=ncut, npts=npts

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: wtheta_conv, type, cutmin=cutmin, cutstep=cutstep, ncut=ncut, npts=npts'
      print,' type=0 (high dens) type=1 (low dens) type=2 (tot)'
      print,' cutmin=50d,cutstep=25d,ncut=47 fills [50,1200]'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;Num points to use in convolution integral
  IF n_elements(npts) EQ 0 THEN Npts = 1500L 


  colors = ['u','g','r','i','z']

  outdir = '/sdss5/data0/lensout/wtheta_conv/'
  typestr = ['high','low', 'tot']
  filters=['wtheta_high_function','wtheta_low_gaussian','wtheta_tot_function']

  defsysv,'!conv_X0',exists=exists1
  IF (NOT exists1) THEN BEGIN
      defsysv,'!conv_X0',0d
      defsysv,'!conv_Y0',0d
      defsysv,'!conv_rmax',0d
      defsysv,'!conv_sigma',0d
      defsysv,'!conv_cutoff',0d
      defsysv,'!atot', dblarr(2)
      defsysv,'!alow', dblarr(4)
      defsysv,'!ahigh', dblarr(2)
  ENDIF 

  ;; use sigma_v=1.0 then rescale later
  !conv_sigma = 1d              ;km/s

  ;; maximum radius to integrate to
  !conv_rmax = 4000d            ;kpc
  ab_limits = [-!conv_rmax, !conv_rmax]

  ;; This is with new tsgals approach
  ;; all galaxies (power law)
  !atot[0] = 0.922449           ;#/Mpc^2, radius in Mpc
  !atot[1] = -0.740816
  !atot[0] = !atot[0]/(1000d)^2 ;#/kpc^2 (for integration)
  !atot[0] = !atot[0]*(1000d)^(-!atot[1]) ;Now radius is in kpc


  ;; low-density regions (gaussian)
  !alow[0] = -5.53803           ;#/Mpc^2, radius in kpc
  !alow[0] = !alow[0]/(1000d)^2 ;#/kpc^2 (for integration)
  !alow[1] = -439.465
  !alow[2] = 397.702
  !alow[3] = -0.278953
  !alow[3] = !alow[3]/(1000d)^2 ;#/kpc^2 (for integration)

  ;; high-density regions (power law)
  !ahigh[0] = 1.95714           ;#/Mpc^2, radius in Mpc
  !ahigh[1] = -0.822449
  !ahigh[0] = !ahigh[0]/(1000d)^2 ;#/kpc^2 (for integration)
  !ahigh[0] = !ahigh[0]*(1000d)^(-!ahigh[1]) ;Now radius is in kpc


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; inputs for functions
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; functions called by conv2d
  print,'Doing '+typestr[type]
  filtername = filters[type]
  print
  func = 'call'+typestr[type]
  pqlim = 'pq_limits'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; where to evaluate the function
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  minx = 0d
  stepx = 20d
  nx = 125L
  xx = dblarr(nx)
  FOR i=0L, nx-1 DO xx[i] = minx + i*stepx

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Values of cutoff radius
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; this does cut in [50,1200] kpc in steps of 10
  IF n_elements(cutmin) EQ 0 THEN cutmin = 50d                  ;kpc
  IF n_elements(cutstep) EQ 0 THEN cutstep = 10d
  IF n_elements(ncut) EQ 0 THEN ncut = 116L

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  arrval = dblarr(nx-1)
  ss = create_struct('cutoff', 0d, $
                     'r', arrval, $
                     'density', arrval, $
                     'density_int', arrval)

  outstruct = replicate(ss, ncut)

  time=systime(1)
  FOR icut = 0L, ncut-1 DO BEGIN

      !conv_cutoff = cutmin + icut*cutstep

      print
      print,'Using: sigma = ',ntostr(!conv_sigma)
      print,'Using: cutoff['+ntostr(icut)+'] = '+ntostr(!conv_cutoff)

      conv2d, func, pqlim, ab_limits, xx, Npts, tmpfunc

      intfunc, tmpfunc*2.*!pi*xx, xx, tmpfunc_int, /double
      outfunc_int = tmpfunc_int[1:nx-1]/!pi/xx[1:nx-1]^2

      outfunc = tmpfunc[1:nx-1]

      outstruct[icut].cutoff = !conv_cutoff
      outstruct[icut].r = xx[1:nx-1]
      outstruct[icut].density = outfunc
      outstruct[icut].density_int = outfunc_int
      
  ENDFOR 
  dtime=systime(1)-time
  ptime,dtime
  print,dtime,' seconds'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot stuff and file stuff 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; plot the values for sigmav=150. km/s
  cutmax = max(outstruct.cutoff)
  cutstr = ntostr(long(cutmin))+'-'+ntostr(long(cutmax))
  pfile=outdir+'conv_'+typestr[type]+'_cut'+cutstr+'_plots.ps'
  print
  print,'Pfile: ',pfile
  title=typestr[type]+' regions   sigmav=150 km/s'
  xtitle='radius (kpc)'
  ytitle='Density Contrast (neighbors)'
  begplot,name=pfile,/landscape
  sigmav2=150d^2
  setupplot, 'PS'
  plot, outstruct[ncut-1].r, $
    (outstruct[ncut-1].density_int - outstruct[ncut-1].density)*sigmav2, $
    xtitle=xtitle, ytitle=ytitle, title=title, linestyle=linestyle
  FOR icut=0L, ncut-1 DO BEGIN 
      linestyle = (icut MOD 2)
      oplot, outstruct[icut].r, $
             (outstruct[icut].density_int-outstruct[icut].density)*sigmav2, $
             linestyle=linestyle
      
  ENDFOR 
  endplot,/noprint
  pslandfix,pfile
  setupplot,'X'
  
  ;; output structure
  fitsoutfile = outdir + 'conv_'+typestr[type]+'_cut'+cutstr+'.fit'

  print
  print,'Writing fits file: ',fitsoutfile
  mwrfits, outstruct, fitsoutfile, /create

return
END 
