PRO make_wthetaconv_func_lumw, clr, wclr, lumsolar, lensw, sigcritinv, retlum=retlum, sigonly=sigonly, noweight=noweight

  ;; model will now be a function of beta and sigma_star
  ;; we will fix the cutoff to be 230 kpc

  ;; Also!! We get power law same as the tot, so all we
  ;; need to do is change the normalization

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: make_wthetaconv_func_lumw, clr, wclr, lumsolar, lensw, sigcritinv, retlum=retlum'
      print
      print,'clr=lensing images   wclr=luminosity weight color'
      print,'/retlum just returs the luminosities and weights'
      return
  ENDIF 

  setup_mystuff

  Mstar = [-18.34,-20.04,-20.83,-21.26,-21.55]
  Msun = [6.39,5.07,4.62,4.52,4.48]
  Lstar = 10.0^((Mstar-Msun)/(-2.5))
  Lstar = Lstar/1.e10           ;in units of 10^10 solar lum

  maxlum = [10.0e10, 10.e10, 15.0e10, 30.0e10, 45.0e10]
  maxz = 0.6
  minz = 0.02

  type=0
  typestr = ['lum','lum1', 'lum2','lum3','lum4']

  ;; read in two structs
  indir = '/sdss5/data0/lensout/wtheta_conv/'

  file = indir+'conv_tot_cut50-1200.fit'
  dstruct = mrdfits(file,1)

  CASE 1 OF
      keyword_set(sigonly): betafile = indir+'fit_wtheta_tot_sigonly_'+!colors[wclr]+'w_lumw.fit'
      keyword_set(noweight): betafile = indir+'fit_wtheta_tot_noweight_'+!colors[wclr]+'w_lumw.fit'
      ELSE: betafile = indir+'fit_wtheta_tot_'+!colors[wclr]+'w_lumw.fit'
  ENDCASE 
  betafile = indir+'fit_wtheta_tot_'+!colors[wclr]+'w_lumw.fit'
  betastruct = mrdfits(betafile,1)
  beta = betastruct.beta
  nbeta = n_elements(betastruct.beta)

  stripes = [10,82,36,37,42,43]
;  stripes = 10
  nstripe = n_elements(stripes)
  stripestr = ['stripe10','stripe82','stripe36','stripe37','stripe42','stripe43','all']
;  stripestr = 'stripe10'

  delvarx,lumsolar,lensw,sigcritinv

  FOR i=0L, nstripe-1 DO BEGIN 

      tstripe = stripes[i]

      lensumfile = '/sdss5/data0/lensout/'+stripestr[i]+'/sublum/'+!colors[wclr]+$
        '/'+typestr[type]+'_zgal_gal_'+stripestr[i]+'_'+!colors[clr]+'_lensum_N1.fit'
      print,lensumfile
      lensum = mrdfits(lensumfile,1)
      print,min(lensum.lum[wclr]),max(lensum.lum[wclr])
      tsigcritinv = lensum.scritinv

      lumsolar = lensum.lum[wclr]/1.e10
      add_arrval, lumsolar/Lstar[wclr], lumw

      ttotpairs = lensum.totpairs
      add_arrval, ttotpairs, totpairs

      ;tDL = angdist_lambda(lensum.z1d)*1000.
      ;add_arrval, tDL, DL
      CASE 1 OF 
          keyword_set(sigonly): tlensw=tsigcritinv^2*1.e9
          keyword_set(noweight): tlensw=replicate(1.0, n_elements(lensum))
          ELSE: tlensw = tsigcritinv^2*ttotpairs*1.e6 ;1.e6 for order 1 numbers 
      ENDCASE 

      add_arrval, tlensw, lensw 
      add_arrval, tsigcritinv, sigcritinv

      lensum=0
  ENDFOR 
  lumsolar = lumw*Lstar[wclr]

  IF keyword_set(retlum) THEN return

  lenswtot = total(lensw)

;(
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; got this from tot
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  cutoff = 230.
  wcut = where(dstruct.cutoff EQ cutoff,nwcut)
  IF nwcut EQ 0 THEN message,'Cutoff '+ntostr(cutoff)+' not found!'

  ;; Make sure no NAN's
  fixwthetaconvnan, dstruct

  ;; normalizations: smoothly varying, so 
  ;; can interpolate
  nbeta = 200L
  maxb=max(betastruct.beta,min=minb)
  beta = arrscl( dindgen(nbeta), minb, maxb)
  norm=interpol(betastruct.norm,betastruct.beta,beta)

  ;; now fill in different values of sigmav and save (can be redone later
  ;; with the output fits structure)
  nsig = 105L
  sigmin = 40d

  sigmin=50d
  sigmax=200d
  sigmav=arrscl(dindgen(nsig),sigmin,sigmax)

  radiuskpc = dstruct[0].r
  nx = n_elements(radiuskpc)
  model_density_contrast = dblarr(nbeta, nsig, nx)
  neighbor_density_contrast = dblarr(nbeta, nsig, nx)
  central_density_contrast = dblarr(nbeta, nsig, nx)
  central_density = dblarr(nbeta, nsig, nx)
  model_ratio = dblarr(nbeta, nx)
  model_ratio_cum = model_ratio

  ;; density contrast of central galaxy
  sig1 = sigmasis_trunc(1.0, cutoff, radiuskpc, /core,/contrast)
  cendens = sigmasis_trunc(1.0, cutoff, radiuskpc, /core)
  ;; density contrast of neighbors
  density_cont = dstruct[wcut].density_int-dstruct[wcut].density
  ;; normalization used to get this density contrast (see wtheta_conv)
  normold = 0.922449           ;#/Mpc^2, radius in Mpc
  powuse = -0.74
  normold = normold/(1000d)^2   ;#/kpc^2 (for integration)
  normold = normold*(1000d)^(-powuse) ;Now radius is in kpc

  ;; divide out the old normalization
  density_cont = density_cont/normold

  FOR ibeta=0L, nbeta-1 DO BEGIN 

      betause = beta[ibeta]
      normuse = norm[ibeta]
      normuse = normuse/(1000d)^2 ;#/kpc^2 (for integration)
      normuse = normuse*(1000d)^(-powuse) ;Now radius is in kpc

      ;; now create model for foreground galaxies with luminosities
      fnorm = total( lensw*(lumw)^betause )/lenswtot
      
      FOR isig=0L, nsig-1 DO BEGIN 
          neighbor_density_contrast[ibeta, isig, *] = $
            normuse*density_cont*sigmav[isig]^2
          central_density_contrast[ibeta, isig, *] = fnorm*sig1*sigmav[isig]^2

          model_density_contrast[ibeta, isig, *] = $
            central_density_contrast[ibeta, isig, *] + $
            neighbor_density_contrast[ibeta, isig, *]

          central_density[ibeta, isig, *] = fnorm*cendens*sigmav[isig]^2
      ENDFOR 

      ;; ratio
      model_ratio[ibeta,*] = density_cont/sig1

      intfunc, density_cont*2.*!pi*radiuskpc, radiuskpc, density_cont_int,$
        /double
      intfunc, sig1*2.*!pi*radiuskpc, radiuskpc, sig1_int,$
        /double
      model_ratio_cum[ibeta,*] = density_cont_int/sig1_int

  ENDFOR 
;  return

  CASE 1 OF
      keyword_set(sigonly): savefile = indir + 'conv_'+!colors[wclr]+'w_'+$
        !colors[clr]+'_'+typestr[type]+'_sigonly_beta0-2.sav'
      keyword_set(noweight): savefile = indir + 'conv_'+!colors[wclr]+'w_'+$
        !colors[clr]+'_'+typestr[type]+'_noweight_beta0-2.sav'
      ELSE: savefile = indir + 'conv_'+!colors[wclr]+'w_'+$
        !colors[clr]+'_'+typestr[type]+'_beta0-2.sav'
  ENDCASE 

  print,savefile
  save, beta, sigmav, radiuskpc, model_density_contrast, $
    neighbor_density_contrast, $
    central_density_contrast, $
    central_density, $
    model_ratio, $
    model_ratio_cum, filename=savefile, /verbose


END 
