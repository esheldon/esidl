PRO clustergrid, run1, run2, clr, clust, $
                 imtan, imrad, imtanerr, imraderr, imtans2n, $
                 gridx, gridy, $
                 gridsize=gridsize, rmin=rmin, rmax=rmax, binsize=binsize, $
                 nx=nx, ny=ny, $
                 scat=scat, $
                 write=write, rotate=rotate,noprompt=noprompt,_extra=e

  IF n_params() EQ 0 THEN BEGIN
      print,'-Syntax: clustergrid, run1, run2, clr, clust, imtan, imrad, gridx, gridy, gridsize=gridsize, rmin=rmin, rmax=rmax, binsize=binsize, nx=nx, ny=ny, scat=scat, write=write, rotate=rotate,_extra=e'
      return
  ENDIF 

  time=systime(1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some Parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  COMMON seed,seed
  IF NOT keyword_set(write) THEN write=0
  IF NOT keyword_set(noprompt) THEN noprompt=0
  IF n_elements(gridsize) EQ 0 THEN gridsize=200. ;arcsec
  IF n_elements(rmin) EQ 0 THEN rmin = 10.        ;arcsec
  IF n_elements(rmax) EQ 0 THEN rmax = 600.       ;arcsec
  IF n_elements(binsize) EQ 0 THEN binsize = 100. ;arcsec

  tx = n_elements(nx)
  ty = n_elements(ny)
  IF tx EQ 0 AND ty EQ 0 THEN BEGIN
      nx = 60 && ny = 60
  ENDIF ELSE BEGIN
      IF tx EQ 0 AND ty NE 0 THEN nx=ny
      IF tx NE 0 AND ty EQ 0 THEN ny=nx
  ENDELSE 


  r1str = ntostr(run1)
  r2str = ntostr(run2)
  colors = ['u','g','r','i','z']

  ;; Define x-position of cluster
  cenx = clust.dec
  ceny = clust.ra

  ;; convert to degrees
  radmin = rmin/3600.
  radmax = rmax/3600.
  bsize = binsize/3600.
  gsize = gridsize/3600.

  xmax  = cenx + gsize/2.
  xmin  = cenx - gsize/2.
  ymax  = ceny + gsize/2.
  ymin  = ceny - gsize/2.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Get the background catalog
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(scat) EQ 0 THEN BEGIN
      dir='/sdss4/data1/esheldon/CORRECTED/'
      fend='_srcgal_'+colors[clr]+'_overlap.fit'
      file=dir+'run'+r1str+'_'+r2str+fend
      IF NOT exist(file) THEN BEGIN
          file=dir+'run'+r2str+'_'+r1str+fend
          IF NOT exist(file) THEN BEGIN
              print,'No '+colors[clr]+' band overlap file exists for ', $
                    'runs '+r1str+' and '+r2str
              return
          ENDIF 
      ENDIF 

      scat=mrdfits(file, 1, hdr)
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up the output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  zstr = ntostr(clust.z,4)

  IF write THEN BEGIN 
      outdir = '/sdss4/data1/esheldon/GRIDSHEAR/SHEAR/'
      prename=outdir+clust.name+'_R'+ntostr(long(rmax))+'_Z'+zstr

      add1 = '_N1.ps'
      add2 = '_N1.fit'
      psfile = prename+add1
      rtfile = prename+'_rot'+add1
      tanfit = prename+'_tan'+add2
      radfit = prename+'_rad'+add2
      s2nfit = prename+'_s2n'+add2
      WHILE exist(psfile) DO BEGIN
          add1=newname(add1)
          add2=newname(add2)
          psfile = prename+add1
          rtfile = prename+'_rot'+add1
          tanfit = prename+'_tan'+add2
          radfit = prename+'_rad'+add2
          s2nfit = prename+'_s2n'+add2
      ENDWHILE 
      print,'Postscript file: ',psfile
  ENDIF 
  IF NOT keyword_set(rotate) THEN rotate=0
  IF noprompt THEN rotate=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Define the grid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  gridx = arrscl( findgen(nx), xmin, xmax )
  gridy = arrscl( findgen(ny), ymin, ymax )

  print,' Cluster location   DEC: ',ntostr(cenx),'  RA: ',ntostr(ceny)
  print,' grid center: ',median(gridx), median(gridy)
  print,' grid xrange: ',min(gridx),max(gridx)
  print,' grid yrange:    ',min(gridy),max(gridy)
  print,' grid size: ', gridsize
  print,' grid spacing: ',abs( gridx[1]-gridx[0] )*3600.

  print,' rmin: ',rmin
  print,' rmax: ',rmax
  print,' binsize: ',binsize
    
  dstr=ntostr(clust.dec,6)
  rstr=ntostr(clust.ra,6)
  addstring = ' < '+ntostr(long(rmax))+'  Z='+zstr+'  '+clust.name
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Select nearby background galaxies and look at their distribution.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF max(scat.ra)  LT ymax+radmax OR min(scat.ra)  GT ymin-radmax OR $
     max(scat.dec) LT xmax+radmax OR min(scat.dec) GT xmin-radmax  THEN BEGIN
      
      print,'Cannot accomodate this grid'
      return
  ENDIF 
  w1 = where(scat.ra GT ymin-radmax AND scat.ra LT ymax+radmax, nw1)
  w2 = where(scat[w1].dec GT xmin-radmax AND scat[w1].dec LT xmax+radmax, nw2)
  w=w1[w2]

  ;; NOTE: plotting ra as x.
  IF write THEN begplot,name=psfile

  print
  print,ntostr(nw2),' Background galaxies found in ', $
        ntostr( 2*rmax+gridsize ),' arcsecond box'
  plot, scat[w].ra, scat[w].dec, psym=3, ystyle=1, xstyle=1, $
        yrange=[xmin-radmax,xmax+radmax],xrange=[ymin-radmax,ymax+radmax]
  plot_box, ymin,ymax,xmin,xmax
  oplot, [ceny], [cenx], psym=4, symsize=2.

  rdannis, /silent
  oplot, [clust.ra], [clust.dec], psym=6, symsize=2 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate the tangential shear around each of the grid points
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  gridshear, gridx, gridy, scat[w], radmin, radmax, bsize, $
             imtan, imrad, imtanerr, imraderr, imtans2n

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make some more outputs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ptime,systime(1)-time

  IF write THEN BEGIN 
      writefits, tanfit, imtan
      writefits, radfit, imrad
      writefits, s2nfit, imtans2n
  ENDIF 

  usex = (gridx - cenx)*3600.
  usey = (gridy - ceny)*3600.

  ntick=5
  
  dummy=strarr(ntick)
  FOR i=0, ntick-1 DO dummy[i]=' '

  xtickn=strarr(ntick)
  minx=min(gridx)
  maxx=max(gridx)
  stepx = (maxx-minx)/(ntick-1)
  FOR i=0, ntick-1 DO xtickn[i]=ntostr( (minx + i*stepx-cenx)*3600., 6)

  ytickn=strarr(ntick)
  miny=min(gridy)
  maxy=max(gridy)
  stepy=(maxy-miny)/(ntick-1)
  FOR i=0, ntick-1 DO ytickn[i]=ntostr( ( miny + i*stepy -ceny)*3600., 6)

  tvim2,/silent,imtan,title='Tangential Shear'+addstring, /noframe, _extra=e
  axis, xaxis=0, xticks=ntick-1, xtickn=xtickn
  axis, xaxis=1, xticks=ntick-1, xtickn=dummy
  axis, yaxis=0, yticks=ntick-1, ytickn=ytickn
  axis, yaxis=1, yticks=ntick-1, ytickn=dummy

  IF NOT write AND NOT noprompt THEN key=get_kbrd(1)
  surface, imtan, usex, usey, $
    title='Tangential Shear'+addstring,charsize=2., _extra=e

  IF NOT write  AND NOT noprompt THEN key=get_kbrd(1)
  tvim2,/silent,imrad,title='Radial Shear'+addstring, /noframe, _extra=e
  axis, xaxis=0, xticks=ntick-1, xtickn=xtickn
  axis, xaxis=1, xticks=ntick-1, xtickn=dummy
  axis, yaxis=0, yticks=ntick-1, ytickn=ytickn
  axis, yaxis=1, yticks=ntick-1, ytickn=dummy

  IF NOT write  AND NOT noprompt THEN key=get_kbrd(1)
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


  IF write THEN endplot,/noprint

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
