PRO specgal_against_edge_mask, stripes, struct, maskflags, masked, unmasked, wbad_monte,wbad, completeness, maskfile=maskfile, redo=redo, poly=poly

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: specgal_against_edge_mask, stripes, struct, maskflags, masked, unmasked, wbad_monte,wbad, completeness, maskfile=maskfile, redo=redo, poly=poly'
      return
  ENDIF 

  FLAGS_MASKED = '1'X
  FLAGS_QUAD1_MASKED = '2'X
  FLAGS_QUAD2_MASKED = '4'X
  FLAGS_QUAD3_MASKED = '8'X
  FLAGS_QUAD4_MASKED = '10'X

  FLAGS_QUAD1_MASKED_MONTE = '20'X
  FLAGS_QUAD2_MASKED_MONTE = '40'X
  FLAGS_QUAD3_MASKED_MONTE = '80'X
  FLAGS_QUAD4_MASKED_MONTE = '100'X

  nstripe = n_elements(stripes)
  nst = n_elements(struct)
  IF n_elements(struct) EQ 0 THEN BEGIN 
      redo = 1
      indir = '/net/cheops1/data3/corrected.local/combined/'

      spfiles = strarr(nstripe)


      FOR i=0,nstripe-1 DO BEGIN
          sstr = stripe2string(stripes[i])
          spfiles[i] = indir+'stripe'+sstr+'_spec.fit'
      ENDFOR 

      mrdfits_multi, spfiles, struct

      s=sort(struct.clambda)
      rmclose_lameta, struct[s], keep, out
      struct = temporary(struct[s[keep]])
      nst = n_elements(struct)
  ENDIF 

  stripestring = ''
  stripestring2 = ''

  FOR i=0,nstripe-1 DO BEGIN
      sstr = stripe2string(stripes[i])
      
      stripestring = stripestring + sstr
      stripestring2 = stripestring2 + sstr
      IF i NE nstripe-1 THEN BEGIN 
          stripestring = stripestring + '_'
          stripestring2 = stripestring2 + ' '
      ENDIF 
  ENDFOR 


  ind = where(struct.clambda NE -9999 AND struct.z GT 0.04,nind)

  radius = 10.                  ;Mpc
  maxangle = radius/angdist_lambda(struct[ind].z)*180./!pi

  IF keyword_set(redo) THEN BEGIN 
      tt=systime(1)
      apply_pixel_mask, struct[ind].clambda, struct[ind].ceta, mm, umm, $
                        maskflags, maskfile=maskfile, maxangle=maxangle
      ptime,systime(1)-tt

      apply_sphpoly_mask, struct[ind].clambda, struct[ind].ceta, sphmasked, sphunmasked, completeness=completeness, /lameta

  ENDIF 

  compcut = 0.0
  sphunmasked = where(completeness GT compcut, comp=sphmasked )
  sphunmasked = ind[sphunmasked]
  sphmasked = ind[sphmasked]

  plot,struct[ind].clambda,struct[ind].ceta,psym=3,xrange=[-70,70],xstyle=1,/ynozero

  ;; look for objects which have two adjacent quadrants that are unmasked

  masked = where( ((maskflags AND FLAGS_MASKED) NE 0) , comp=unmasked)
  masked = ind[masked]
  unmasked = ind[unmasked]

  ;; this will subscript ind
  wgood = where( ((maskflags AND FLAGS_MASKED) EQ 0)  AND $
                 ( $
                   ((maskflags AND (FLAGS_QUAD1_MASKED+FLAGS_QUAD2_MASKED)) EQ 0) OR $
                   ((maskflags AND (FLAGS_QUAD2_MASKED+FLAGS_QUAD3_MASKED)) EQ 0) OR $
                   ((maskflags AND (FLAGS_QUAD3_MASKED+FLAGS_QUAD4_MASKED)) EQ 0) OR $
                   ((maskflags AND (FLAGS_QUAD4_MASKED+FLAGS_QUAD1_MASKED)) EQ 0) $
                 ), $
                 ngood, $
                 comp = wbad, ncomp = nbad)

  ;; will subscript wbad
  wbad_monte  = where( ((maskflags[wbad] AND FLAGS_QUAD1_MASKED_MONTE) NE 0) OR $
                       ((maskflags[wbad] AND FLAGS_QUAD2_MASKED_MONTE) NE 0) OR $
                       ((maskflags[wbad] AND FLAGS_QUAD3_MASKED_MONTE) NE 0) OR $
                       ((maskflags[wbad] AND FLAGS_QUAD4_MASKED_MONTE) NE 0), nwbad_monte)
  wbad_monte = wbad[wbad_monte]



  wgood = ind[wgood]
  wbad  = ind[wbad]
  wbad_monte = ind[wbad_monte]

  wgood_poly = where( ((maskflags AND FLAGS_MASKED) EQ 0)  AND $
                      (completeness GT compcut) AND $
                      ( $
                        ((maskflags AND (FLAGS_QUAD1_MASKED+FLAGS_QUAD2_MASKED)) EQ 0) OR $
                        ((maskflags AND (FLAGS_QUAD2_MASKED+FLAGS_QUAD3_MASKED)) EQ 0) OR $
                        ((maskflags AND (FLAGS_QUAD3_MASKED+FLAGS_QUAD4_MASKED)) EQ 0) OR $
                        ((maskflags AND (FLAGS_QUAD4_MASKED+FLAGS_QUAD1_MASKED)) EQ 0) $
                      ), $
                      ngood_poly, $
                      comp = wbad_poly, ncomp = nbad_poly)
  wgood_poly = ind[wgood_poly]
  wbad_poly = ind[wbad_poly]

  help,masked,wbad_monte,wbad,wgood,wgood_poly

  ;; plot "bad" ones
  oplot, struct[masked].clambda, struct[masked].ceta, psym=8,color=!red, symsize=0.5
  oplot, struct[wbad].clambda, struct[wbad].ceta, psym=8, color=!red, symsize=0.5

  ;; plot ones that were caught by the monte-carlo routine
  oplot, struct[wbad_monte].clambda, struct[wbad_monte].ceta, $
         psym=8, color=!green, symsize=0.5

  

  FOR i=0L, nstripe-1 DO BEGIN 
      ist = stripes[i]
      primary_bound, ist, bound
      plot_box, bound.lammin, bound.lammax, bound.etamin, bound.etamax,color=!blue
  ENDFOR 
  
  IF keyword_set(poly) THEN BEGIN 
      
      oplot, struct[sphmasked].clambda, struct[sphmasked].ceta,$
             psym=8,color=!yellow,symsize=0.5
      legend,['Edge Pixel', 'Pixel Monte Carlo', 'Masked Poly'], $
             psym=[8,8,8], color=[!red,!green,!yellow],/clear

      legend,['Radius = '+ntostr(long(radius))+' Mpc',$
              'Ngood = '+ntostr(ngood_poly)+'/'+ntostr(nind)],/right,/clear

      file = '~/tmp/stripe'+stripestring+'_edge_poly'

  ENDIF ELSE BEGIN 
      legend,['Edge Pixel', 'Pixel Monte Carlo'], psym=[8,8], color=[!red,!green],/clear
      legend,['Radius = '+ntostr(long(radius))+' Mpc',$
              'Ngood = '+ntostr(ngood)+'/'+ntostr(nind)],/right,/clear
      file = '~/tmp/stripe'+stripestring+'_edge'
  ENDELSE 

  legend,'Stripes: '+stripestring2,/right,/bottom,/clear
return
  print,'Writing files: '+file+'*'
  write_jpeg,file+'.jpg',tvrd(/true),/true
  write_png, file+'.png',tvrd(/true)

  return

END 
