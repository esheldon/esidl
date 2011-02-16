PRO make_wthetaconv_func, type, concat=concat

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: make_wthetaconv_func, type, concat=concat'
      print
      print,' type=0 (high dens) type=1 (low dens) type=2 (tot)'
      return
  ENDIF 

  typestr = ['high','low', 'tot']

  ;; read in two structs
  indir = '/sdss5/data0/lensout/wtheta_conv/'
  IF keyword_set(concat) THEN BEGIN 
      file1 = indir+'conv_'+typestr[type]+'_cut50-600.fit'
      file2 = indir+'conv_'+typestr[type]+'_cut625-1200.fit'
      
      st1 = mrdfits(file1,1)
      st2 = mrdfits(file2,1)
      
      concat_structs, st1, st2, dstruct
  ENDIF ELSE BEGIN
      file = indir+'conv_'+typestr[type]+'_cut50-1200.fit'
      dstruct = mrdfits(file,1)
  ENDELSE 
  cutoff = dstruct.cutoff
  ncut = n_elements(dstruct)

  ;; Make sure no NAN's
  fixwthetaconvnan, dstruct

  ;; now fill in different values of sigmav and save (can be redone later
  ;; with the output fits structure)
  nsig = 105L
  sigmin = 40d
;  sigstep = double(2.5)
;  sigmav = dblarr(nsig)
;  FOR isig=0L, nsig-1 DO sigmav[isig]=sigmin + isig*sigstep

  sigmin=50d
  sigmax=200d
  sigmav=arrscl(dindgen(nsig),sigmin,sigmax)

  radiuskpc = dstruct[0].r
  nx = n_elements(radiuskpc)
  model_density_contrast = dblarr(ncut, nsig, nx)
  neighbor_density_contrast = dblarr(ncut, nsig, nx)
  central_density_contrast = dblarr(ncut, nsig, nx)
  central_density = dblarr(ncut, nsig, nx)
  model_ratio = dblarr(ncut, nx)
  model_ratio_cum = model_ratio

  FOR icut=0L, ncut-1 DO BEGIN 
      ;; density contrast of central galaxy
      sig1 = sigmasis_trunc(1.0, cutoff[icut], radiuskpc, /core,/contrast)
      cendens = sigmasis_trunc(1.0, cutoff[icut], radiuskpc, /core)
      ;; density contrast of neighbors
      density_cont = dstruct[icut].density_int-dstruct[icut].density

      FOR isig=0L, nsig-1 DO BEGIN 
          neighbor_density_contrast[icut, isig, *] = density_cont*sigmav[isig]^2
          central_density_contrast[icut, isig, *] = sig1*sigmav[isig]^2
          central_density[icut, isig, *] = cendens*sigmav[isig]^2
      ENDFOR 

      ;; ratio
      model_ratio[icut,*] = density_cont/sig1

      intfunc, density_cont*2.*!pi*radiuskpc, radiuskpc, density_cont_int,$
        /double
      intfunc, sig1*2.*!pi*radiuskpc, radiuskpc, sig1_int,$
        /double
      model_ratio_cum[icut,*] = density_cont_int/sig1_int

      density_cont = density_cont + sig1


      FOR isig=0L, nsig-1 DO BEGIN 
          model_density_contrast[icut,isig,*] = (density_cont)*sigmav[isig]^2
      ENDFOR 

  ENDFOR 
;  return
  savefile = indir + 'conv_'+typestr[type]+'_cut50-1200'+'.sav'
  print,savefile
  save, cutoff, sigmav, radiuskpc, model_density_contrast, $
    neighbor_density_contrast, $
    central_density_contrast, $
    central_density, $
    model_ratio, $
    model_ratio_cum, filename=savefile, /verbose


END 
