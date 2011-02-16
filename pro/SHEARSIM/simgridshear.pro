PRO simgridshear, imtan, imrad, imtanerr, imraderr, imtans2n, gridx, gridy, $
                  sdec=sdec,sra=sra,$
                  density=density, $
                  sigma=sigma, zlens=zlens, zsource=zsource, $
                  rmin=rmin, rmax=rmax, binsize=binsize, $
                  gridsize=gridsize, nx=nx, ny=ny, $
                  write=write, noise=noise, rotate=rotate,noprompt=noprompt,$
                  _extra=e

  IF n_params() EQ 0 THEN BEGIN
      print,'-Syntax: simgridshear, imtan, imrad, imtanerr, imraderr, imtans2n, gridx, gridy, density=density, sigma=sigma, zlens=zlens, zsource=zsource, rmin=rmin, rmax=rmax, binsize=binsize, gridsize=gridsize, nx=nx, ny=ny, write=write, noise=noise, rotate=rotate,_extra=e'
      return
  ENDIF 

  time=systime(1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some Parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  COMMON seed,seed
  IF n_elements(sdec) EQ 0 THEN sdec = 1. ;degrees
  IF n_elements(sra) EQ 0 THEN sra = 1.
  IF n_elements(density) EQ 0 THEN density = 7500. ; per square degree
  IF n_elements(sigma) EQ 0 THEN sigma = 1000.  ;km/s
  IF n_elements(zlens) EQ 0 THEN zlens = .15
  IF n_elements(zsource) EQ 0 THEN zsource = .4
  IF n_elements(rmin) EQ 0 THEN rmin = 10.  ;arcsec
  IF n_elements(rmax) EQ 0 THEN rmax = 611. ;arcsec
  IF n_elements(binsize) EQ 0 THEN binsize = 100. ;arcsec
  IF n_elements(gridsize) EQ 0 THEN gridsize = 401. ;arcsec

  tx = n_elements(nx)
  ty = n_elements(ny)
  IF tx EQ 0 AND ty EQ 0 THEN BEGIN
      nx = 61L && ny = 61L
  ENDIF ELSE BEGIN
      IF tx EQ 0 AND ty NE 0 THEN nx=ny
      IF tx NE 0 AND ty EQ 0 THEN ny=nx
  ENDELSE 

  IF NOT keyword_set(noise) THEN noise=0 ELSE noise=1
  IF NOT keyword_set(write) THEN write=0
  IF NOT keyword_set(noprompt) THEN noprompt=0

  radmin = rmin/3600.
  radmax = rmax/3600.
  bsize = binsize/3600.

  gsize = gridsize/3600.
  gridmax_x  = (sdec + gsize)/2.
  gridmin_x  = (sdec - gsize)/2.
  gridmax_y  = (sra  + gsize)/2.
  gridmin_y  = (sra  - gsize)/2.

  IF ( gridmax_x + radmax GT sdec OR $
       gridmin_x - radmax LT 0. OR $
       gridmax_y + radmax GT sra OR $
       gridmin_y - radmax LT 0. ) THEN BEGIN 
      print,'Cannot accomodate this grid'
      return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up the output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF write THEN BEGIN 
      outdir = '/sdss4/data1/esheldon/GRIDSHEAR/'
      prename = outdir+'sim_grid_s'+ntostr(long(sigma))

      psfile = prename+'_N1.ps'
      rtfile = prename+'_rot_N1.ps'
      tanfit = prename+'_tan_N1.fit'
      radfit = prename+'_rad_N1.fit'
      s2nfit = prename+'_s2n_N1.fit'
      WHILE exist(psfile) DO BEGIN
          psfile = newname(psfile)
          rtfiel = newname(rtfile)
          tanfit = newname(tanfit)
          radfit = newname(radfit)
          s2nfit = newname(s2nfit)
      ENDWHILE 
      print,'Postscript file: ',psfile
  ENDIF 
  IF NOT keyword_set(rotate) THEN rotate=0
  IF noprompt THEN rotate=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Define the grid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  gridx = arrscl( findgen(nx), gridmin_x, gridmax_x )
  gridy = arrscl( findgen(ny), gridmin_y, gridmax_y )

  print,' sigma: ',ntostr(long(sigma)),'  zlens: ',ntostr(zlens,4),$
        '  zsource: ',ntostr(zsource,4)
  print,' rmin: ',rmin,' rmax: ',rmax
  print,' binsize: ',binsize
  
  print,' grid center: ',median(gridx), median(gridy)
  print,' grid size: ',gridsize
  print,' grid spacing: ',abs( gridx[1]-gridx[0] )*3600.
  
  addstring = ' < '+ntostr(long(rmax))
  IF noise THEN BEGIN
      print
      print,'Adding Noise'
      addstring = addstring+' (Noise)'
  ENDIF
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Create the simulation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sis_shear, zlens, zsource, sigma, sdec, sra, cat, cen, $
    density=density, noise=noise

  gridshear, gridx, gridy, cat, radmin, radmax, bsize, $
             imtan, imrad, imtanerr, imraderr, imtans2n

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; View the images, make outputs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ptime,systime(1)-time

  IF write THEN BEGIN 
      writefits, tanfit, imtan
      writefits, radfit, imrad
      writefits, s2nfit, imtans2n
  ENDIF 

  usex = (gridx - cen[0])*3600.
  usey = (gridy - cen[1])*3600.

  ntick=5
  
  dummy=strarr(ntick)
  FOR i=0, ntick-1 DO dummy[i]=' '

  xtickn=strarr(ntick)
  minx=min(gridx)
  maxx=max(gridx)
  stepx = (maxx-minx)/(ntick-1)
  FOR i=0, ntick-1 DO xtickn[i]=ntostr( (minx + i*stepx-cen[0])*3600., 6)

  ytickn=strarr(ntick)
  miny=min(gridy)
  maxy=max(gridy)
  stepy=(maxy-miny)/(ntick-1)
  FOR i=0, ntick-1 DO ytickn[i]=ntostr( ( miny + i*stepy -cen[1])*3600., 6)

  IF write THEN begplot, name=psfile
  tvim2,/silent,imtan,title='Tangential Shear'+addstring, /noframe, _extra=e
  axis, xaxis=0, xticks=ntick-1, xtickn=xtickn
  axis, xaxis=1, xticks=ntick-1, xtickn=dummy
  axis, yaxis=0, yticks=ntick-1, ytickn=ytickn
  axis, yaxis=1, yticks=ntick-1, ytickn=dummy

  IF NOT write AND NOT noprompt THEN key=get_kbrd(1)
  surface, imtan, usex, usey, $
    title='Tangential Shear'+addstring,charsize=2., _extra=e



  IF NOT write AND NOT noprompt THEN key=get_kbrd(1)
  tvim2,/silent,imrad,title='Radial Shear'+addstring, /noframe, _extra=e
  axis, xaxis=0, xticks=ntick-1, xtickn=xtickn
  axis, xaxis=1, xticks=ntick-1, xtickn=dummy
  axis, yaxis=0, yticks=ntick-1, ytickn=ytickn
  axis, yaxis=1, yticks=ntick-1, ytickn=dummy

  IF NOT write AND NOT noprompt THEN key=get_kbrd(1)
  surface,imrad, usex, usey, title='Radial Shear'+addstring, charsize=2., $
    _extra=e



  IF NOT write AND NOT noprompt THEN key=get_kbrd(1)
  tvim2,/silent,imtans2n,title='Tan S/N'+addstring, /noframe, _extra=e
  axis, xaxis=0, xticks=ntick-1, xtickn=xtickn
  axis, xaxis=1, xticks=ntick-1, xtickn=dummy
  axis, yaxis=0, yticks=ntick-1, ytickn=ytickn
  axis, yaxis=1, yticks=ntick-1, ytickn=dummy  

  IF NOT write AND NOT noprompt THEN key=get_kbrd(1)
  surface,imtans2n, usex, usey, title='Tan S/N'+addstring, charsize=2., $
    _extra=e


  IF write THEN ep

  IF rotate THEN BEGIN 
      IF write THEN BEGIN 
          rotate_plot, imtan, usex, usey, step=5., $
            title='Tangential Shear'+addstring,psfile=rtfile,charsize=2., $
            _extra=e
      ENDIF ELSE BEGIN 
          key=get_kbrd(1)
          rotate_plot, imtan, usex, usey, step=5., $
            title='Tangential Shear'+addstring,charsize=2., _extra=e
      ENDELSE 
  ENDIF 

  return
END 
