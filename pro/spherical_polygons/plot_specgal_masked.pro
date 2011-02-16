PRO plot_specgal_masked_plotbounds, stripes

  primary_bound_multi, stripes, bound_arr, overall_bound

  thick = 2.*!p.thick
  nst = n_elements(stripes)
  FOR i=0, nst-1 DO BEGIN 
      plot_box, bound_arr[i].lammin, bound_arr[i].lammax, $
                bound_arr[i].etamin, bound_arr[i].etamax, color=!blue,thick=thick
  ENDFOR 
  
  plot_box, overall_bound.lammin, overall_bound.lammax, $
            overall_bound.etamin, overall_bound.etamax, color=!green,thick=thick

END 

PRO plot_specgal_masked, stripes, spec, phot, comp

  nst = n_elements(stripes)

  dir = sdssidl_config('SHAPECORR_DIR')+'combined/'

  IF n_elements(spec) EQ 0 THEN BEGIN 

      specfiles = strarr(nst)
      FOR i=0L,nst-1 DO specfiles[i] = $
        dir + 'stripe'+stripe2string(stripes[i])+'_spec.fit'

      mrdfits_multi, specfiles, spec, count=count

  ENDIF 

  IF n_elements(phot) EQ 0 THEN BEGIN 

      photfiles = strarr(nst)
      FOR i=0L,nst-1 DO photfiles[i] = $
        dir + 'stripe'+stripe2string(stripes[i])+'_srcgal_gri.fit'

      mrdfits_multi, photfiles, phot, count=count

  ENDIF 

  IF n_elements(comp) EQ 0 THEN BEGIN 
      sphpoly_completeness, spec.ra, spec.dec, comp
  ENDIF 

  eq2csurvey, spec.ra, spec.dec, lam, eta
  primary_bound_multi, stripes, bound_arr, overall_bound
  
  maxx = max(abs([overall_bound.lammin, overall_bound.lammax]))
  xrange = [-maxx,maxx]*1.1

;  pxmin = min(phot.clambda, max=pxmax)

;  sxmin = min(lam, max=sxmax)

;  xmin = min([pxmin, sxmin])
;  xmax = max([pxmax, sxmax])

  ;; symmetric
;  maxx = max(abs([xmin,xmax]))
;  xmin = -maxx
;  xmax =  maxx

  allgal_clr = !p.color
  specgal_clr = !cyan
  specgal_unmasked_clr = !red

  stripestring = stripearr2string(stripes)
  plot, phot.clambda, phot.ceta,psym=3,xtitle=$
        !csym.lambda+'!Dc!N', $
        ytitle=!csym.eta+'!Dc!N', $
        title = 'Stripes '+stripestring,/ynozero, $
        charsize = 1.3, xrange=xrange,/xstyle


  oplot, lam, eta, psym=3, color=specgal_clr

  nphot = n_elements(phot)
  nspec = n_elements(spec)
  w=where(comp GT 0,ngood)
  oplot, lam[w], eta[w], psym=3, color=specgal_unmasked_clr

  plot_specgal_masked_plotbounds, stripes

  mess = ['Gal: '+ntostr(nphot),'SpecGal: '+ntostr(nspec),'Unmasked: '+ntostr(ngood)]
  legend,mess,psym=8,color=[allgal_clr,specgal_clr,specgal_unmasked_clr],/right,/clear

  outdir = '~/tmp/'
  jpegfile = outdir + 'stripe'+stripestring+'.jpg'
  pngfile = outdir + 'stripe'+stripestring+'.png'

  print,'Writing jpeg file: ',jpegfile
  write_jpeg,jpegfile,tvrd(true=1),quality=75,true=1
  print,'Writing pngfile file: ',pngfile
  write_png,pngfile,tvrd(true=1)

END 
