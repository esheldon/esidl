PRO rdregress, file, meanr, sh1, sh1err, sh2, sh2err, shave, shaverr, old=old, _extra=extra, Sshinv=Sshinv

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: rdregress, file, meanr, sh1, sh1err, sh2, sh2err, shave, shaverr,old=old, _extra=extra, Sshinv=Sshinv'
      return
  ENDIF 

  IF n_elements(Sshinv) EQ 0 THEN Sshinv=1.0
  
  IF keyword_set(old) THEN BEGIN 
      readcol, file, meanr, sh1, sh1err, sh2, sh2err
  ENDIF ELSE BEGIN 
      readcol, file, meanr, sh1, sh1err, sh2, sh2err, ortho1, ortho2
  ENDELSE 

  IF keyword_set(old) THEN BEGIN 
      sh1 = -sh1*Sshinv
      sh2 = -sh2*Sshinv
  ENDIF ELSE BEGIN 
      sh1 = sh1*Sshinv
      sh2 = sh2*Sshinv
      ortho1=ortho1*Sshinv
      ortho2=ortho2*Sshinv
  ENDELSE 

  w1 = 1./sh1err^2
  w2 = 1./sh2err^2

  nbin = n_elements(sh1)
  shave = fltarr(nbin)
  shaverr = fltarr(nbin)
  orthoave = fltarr(nbin)

  FOR i=0, nbin-1 DO BEGIN 
      
      wsum = w1[i] + w2[i]
      shave[i] = ( sh1[i]*w1[i] + sh2[i]*w2[i])/wsum
      
      shaverr[i] = sqrt(1./wsum)
      
      IF NOT keyword_set(old) THEN BEGIN 
          orthoave[i] = ( ortho1[i]*w2[i] + ortho2[i]*w1[i] )/wsum
      ENDIF 

  ENDFOR 

  erase & multiplot, [1,3]

  xtitle = 'radius(arcsec)'
  ytitle = 'shear1'
  yrange=[-.001,.007]

  ploterr, meanr, sh1, sh1err, psym=1, ytitle=ytitle,yrange=yrange,ystyle=1, $
    _extra=extra
  oplot,[0,10000],[0,0]
  multiplot

  ytitle='shear2'
  ploterr,meanr, sh2, sh2err, psym=1, ytitle=ytitle,yrange=yrange,ystyle=1
  oplot,[0,10000],[0,0]
  multiplot

  ytitle='Shear Avg'
  ploterr, meanr, shave, shaverr, psym=1, ytitle=ytitle,xtitle=xtitle,$
    yrange=yrange,ystyle=1
  oplot,[0,10000],[0,0]
  multiplot, /reset
  
  IF NOT keyword_set(old) THEN BEGIN 
      
      key = get_kbrd(1)
      erase & multiplot, [1,3]

      xtitle = 'radius(arcsec)'
      ytitle = 'Ortho-tangential Shear1'
      yrange=[-.004,.004]

      ploterr, meanr, ortho1, sh2err, psym=1, $
        ytitle=ytitle,yrange=yrange,ystyle=1, _extra=extra
      oplot,[0,10000],[0,0]
      multiplot

      ytitle = 'Ortho-tangential Shear2'
      ploterr, meanr, ortho2, sh1err, psym=1, $
        ytitle=ytitle,yrange=yrange,ystyle=1
      oplot,[0,10000],[0,0]
      multiplot

      ytitle='Ortho Avg'
      ploterr, meanr, orthoave, shaverr, psym=1, ytitle=ytitle,xtitle=xtitle,$
        yrange=yrange,ystyle=1
      oplot,[0,10000],[0,0]
      multiplot, /reset
      
  ENDIF 

  colprint,meanr,sh1,sh2,shave


return
END
