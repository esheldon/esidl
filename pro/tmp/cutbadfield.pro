PRO cutbadfield, shstruct, nsig1, nsig2, wgood, wbad

  IF n_params() LT 2 THEN BEGIN 
      print,'test3, shstruct, nsig1, nsig2,wgood'
      return
  ENDIF 

                                ;Find out device
  IF (!d.flags AND 1) EQ 0 THEN doX=1 ELSE doX=0

  simpctable
  !p.background=!white
  !p.color=!black

  nf = n_elements(shstruct)
  
  maxe=0.05
  val = replicate(maxe, nf)

  step = .01
  add = 0.
  title = 'Run '+ntostr(shstruct[0].run)+' Rerun '+ntostr(shstruct[0].rerun)+$
          ' Bandpass '+ntostr(shstruct[0].bandpass)
  plot, [0], /nodata, xrange=[0,800], yrange=[0,.05], xtitle='field', ytitle='<e> PSF', title=title
  mean = fltarr(6)
  sig  = fltarr(6)
  meangal = mean
  siggal = mean
  meanstar = mean
  sigstar = mean

  pcol = 3

  FOR col=0, 5 DO BEGIN 

      ;; Stats on psf
      sigma_clip, shstruct.psfe[col], mn, sg,nsig=3.5,niter=4,/silent
      sigma_clip, shstruct.nstar[col], mnstar, sgstar, nsig=3.5, niter=4, /silent

      ;; Stats on gals
      sigma_clip, shstruct.ngal[col], mngal, sggal,nsig=3.5,niter=4,/silent

      

;      mean[col] = median(shstruct.psfe[col])
;      sig[col] = sqrt( total( (shstruct.psfe[col] - mean[col])^2 )/(nf-1) )

      mean[col] = mn
      sig[col] = sg
      meangal[col] = mngal
      siggal[col] = sggal
      meanstar[col] = mnstar
      sigstar[col] = sgstar

      maxe = mean[col] + nsig1*sig[col]

      w=where(shstruct.psfe[col] GT maxe, nw)
;      w = where(shstruct.psfe[col] GT maxe AND shstruct.ngal[col] GT meangal[col]-siggal[col],nw)
      IF col EQ pcol THEN wcol = w

;      oplot, shstruct.field, smooth(shstruct.psfe[col]+add,3)
      oplot, shstruct.field, shstruct.psfe[col]+add
      IF nw NE 0 THEN BEGIN 
          IF n_elements(wb) EQ 0 THEN wb = w ELSE wb = [wb,w]
;          oplot, shstruct[w].field, (smooth(shstruct.psfe[col]+add,3))[w], psym=4, color=!red
          oplot, shstruct[w].field, shstruct[w].psfe[col]+add, psym=4, color=!red
      ENDIF 

      add = add + step

  ENDFOR 

  wb = wb[ rem_dup(wb) ]
  nb =  n_elements(wb)
  corr = fltarr(nb)


  IF doX THEN key=get_kbrd(1)

  corr = fltarr(nf)
  corr2 = corr
  FOR i=0L, nf-1 DO BEGIN 
      FOR jj=0L, 5 DO BEGIN 
          corr[i] = corr[i] + (shstruct[i].psfe[jj]-mean[jj])
      ENDFOR 
      FOR jj= 0L, 4 DO BEGIN 
          FOR kk=jj+1, 5 DO BEGIN 
              corr2[i] = corr2[i] + (shstruct[i].psfe[jj]-mean[jj])*(shstruct[i].psfe[kk]-mean[kk])
          ENDFOR 
      ENDFOR 
  ENDFOR 

  ;; corr2
  sigma_clip, corr2, mn, sg,nsig=3.5,niter=4,/silent

  mean = replicate(mn, nf)
  sig = replicate(sg,nf)

  wbad = where( corr2[wb] GT mean[0]+nsig2*sig[0] , nbad)
  wbad = wb[wbad]
  print,nbad

  erase & multiplot, [1,2]

  title = 'Camcol = '+ntostr(pcol)
  plot,shstruct.field,shstruct.ngal[pcol],ytitle='Ngal', title=title
  oplot,shstruct.field,replicate(meangal[pcol], nf)
  oplot,shstruct.field,replicate(meangal[pcol], nf)-siggal[pcol]
  oplot,shstruct.field,replicate(meangal[pcol], nf)+siggal[pcol]
  oplot,shstruct[wcol].field, shstruct[wcol].ngal[pcol],psym=4,color=!red

  multiplot
  plot,shstruct.field,shstruct.nstar[pcol],ytitle='Nstar'
  oplot,shstruct.field,replicate(meanstar[pcol], nf)
  oplot,shstruct.field,replicate(meanstar[pcol], nf)-sigstar[pcol]
  oplot,shstruct.field,replicate(meanstar[pcol], nf)+sigstar[pcol]
  oplot,shstruct[wcol].field, shstruct[wcol].nstar[pcol],psym=4,color=!red
  multiplot,/reset

  IF doX THEN key=get_kbrd(1)

  fac=1.e-4
  plot, shstruct.field, corr2/fac, ytitle='Corr2',xtitle='Field'
  oplot,shstruct.field, mean/fac, color=!blue
  oplot,shstruct.field, (mean-nsig2*sig)/fac, color=!blue, line=2
  oplot,shstruct.field, (mean+nsig2*sig)/fac, color=!blue, line=2

  oplot, shstruct[wb].field, corr2[wb]/fac, psym=4, color=!red

  legend, 'out = '+ntostr(nbad),/right

  shstruct[wbad].fstatus = -1
  wgood = lindgen(nf)
  remove, wbad, wgood

END 
