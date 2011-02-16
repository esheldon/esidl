PRO spectra_plot, stripes, clr, dops=dops, type=type_in, start=start, nf=nf, $
                  dojack=dojack, indir=indir, $
                  combined_stripes=combined_stripes, $
                  cstripe_strings=cstripe_strings, $
                  recorr=recorr, $
                  lumclr=lumclr, hirata=hirata, comoving=comoving, $
                  lrg_sources=lrg_sources, rlrg_sources=rlrg_sources, $
                  logbin=logbin

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: spectra_plot, stripes, clr, dops=dops, type=type, start=start, nf=nf, dojack=dojack, indir=indir, combined_stripes=combined_stripes, cstripe_strings=cstripe_strings, lumclr=lumclr, comoving=comoving, /lrg_sources, /rlrg_sources, /logbin'
      return
  ENDIF 

  IF n_elements(start) EQ 0 THEN strt=0 ELSE strt=start-1
  IF n_elements(nf) EQ 0 THEN nf=1
  IF n_elements(indir) EQ 0 THEN indir = esheldon_config("lensout_dir")
  IF n_elements(type) EQ 0 THEN type=''
  IF keyword_set(recorr) THEN recorrstr='_recorr' ELSE recorrstr=''
  IF keyword_set(hirata) THEN hirstr = '_h' ELSE hirstr = ''

  IF keyword_set(lrg_sources) THEN BEGIN 
      lrgstr = '_lrg' 
  ENDIF ELSE BEGIN 
      IF keyword_set(rlrg_sources) THEN lrgstr = '_rlrg' ELSE lrgstr=''
  ENDELSE 

  IF n_elements(type_in) EQ 0 THEN type='' ELSE type=type_in+'_'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; stripe strings
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(combined_stripes) THEN BEGIN 
      ;; input is strings
      stripestr = 'stripe'+cstripe_strings
  ENDIF ELSE BEGIN 
      stripestr = 'stripe'+ntostr(stripes)
  ENDELSE 
  nrstr=n_elements(stripestr)

  totstripestr = 'stripe'+stripearr2string(stripes)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; color string
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  clr_string = clrarr2string(clr)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; jackknifing?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(dojack) THEN jackstr = '_jack' ELSE jackstr=''

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; file names
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(lumclr) THEN BEGIN 
      faddstr = +'sublum/'+!colors[lumclr]+'/'
  ENDIF ELSE BEGIN 
      faddstr = ''
  ENDELSE 

  filedir = indir + stripestr + '/' + faddstr

  files =        filedir + type+ 'zgal_gal_'+stripestr+'_'+clr_string+lrgstr+recorrstr+hirstr
  lensumfiles =  filedir + type+ 'zgal_gal_'+stripestr+'_'+clr_string+lrgstr+recorrstr+hirstr+'_lensum'
  rfiles =       filedir + type+ 'zrand_'+stripestr+'_'+clr_string+lrgstr+recorrstr+hirstr
  rlensumfiles = filedir + type+ 'zrand_'+stripestr+'_'+clr_string+lrgstr+recorrstr+hirstr+'_lensum'
  zfiles =       filedir + type+ 'zgal_gal_'+stripestr+'_'+clr_string+lrgstr+recorrstr+hirstr+'_z'

  IF keyword_set(dojack) THEN BEGIN 
      tjackknife_file = filedir + 'jackknife_samples/jackknife_'+$
        type+ 'zgal_gal_'+stripestr+'_'+clr_string+lrgstr+recorrstr+hirstr+$
        '_lensum'
      trjackknife_file = filedir + 'jackknife_samples/jackknife_'+$
        type+ 'zrand_'+stripestr+'_'+clr_string+lrgstr+recorrstr+hirstr+$
        '_lensum'
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output file front end
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  outfront = type+ 'zgal_gal_'+totstripestr+'_'+clr_string+lrgstr+recorrstr+hirstr+jackstr
  routfront = type+ 'zrand_'+totstripestr+'_'+clr_string+lrgstr+recorrstr+hirstr+jackstr

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; ranges for fitting
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  cutrange=[100., 1300.]
  powrange=[0.3, 1.4]

  sigvrange = [50.0, 500.]

  IF keyword_set(comoving) THEN comove_str = '_comoving' ELSE comove_str = ''

  FOR i=strt, strt+nf-1 DO BEGIN 

      istr = ntostr(i+1)
      ext = '_N'+istr+'.fit'
      ps_ext = '_N'+istr+'.ps'

      out_ext = comove_str + ext
      out_ps_ext = comove_str + ps_ext 

      infiles = files + ext
      rinfiles = rfiles + ext
      zinfiles = zfiles + ext

      IF keyword_set(combined_stripes) THEN BEGIN 
          paramname = indir + 'combstripe/fits/'+faddstr+outfront + '_fitparam'+out_ext
          combname = indir + 'combstripe/comb/'+faddstr+outfront + '_comb'+out_ext
          rcombname = indir + 'combstripe/comb/'+faddstr+routfront + '_comb'+out_ext
          psfile    = indir + 'combstripe/figures/'+faddstr+outfront + out_ps_ext
      ENDIF ELSE BEGIN 
          paramname = indir + 'fits/'+faddstr+outfront + '_fitparam'+out_ext
          combname = indir + 'comb/'+faddstr+outfront + '_comb'+out_ext
          rcombname = indir + 'comb/'+faddstr+routfront + '_comb'+out_ext
          psfile    = indir + 'figures/'+faddstr+outfront + out_ps_ext
      ENDELSE 

      IF keyword_set(dojack) THEN BEGIN 
          jackknife_file = tjackknife_file + ext

          rjackknife_file = trjackknife_file + ext
      ENDIF 

;      IF keyword_set(dojack) THEN BEGIN
          lensuminfiles = lensumfiles + ext
;      ENDIF 
      
      IF keyword_set(dops) THEN begplot, name=psfile,/color
      setup_mystuff

      plot_zshear, infiles, rinfiles, paramname, combname, rcombname, $
        zfiles=zinfiles, $
        sigvrange=sigvrange, cutrange=cutrange, $
        normrange=normrange, powrange=powrange, $
        lensumfiles=lensuminfiles, $
        munit=1.e12, lumclr=lumclr, comoving=comoving,$
        logbin=logbin, $
        jackknife_file=jackknife_file,rjackknife_file=rjackknife_file

      IF  keyword_set(dops) THEN endplot
      setup_mystuff
      IF (NOT keyword_set(dops)) AND (i NE strt+nf-1) THEN BEGIN 
          key=get_kbrd(1)
          IF key EQ 'q' THEN return
      ENDIF 

  ENDFOR 


  return
END 
