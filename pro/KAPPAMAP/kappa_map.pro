PRO kappa_map, run1, run2, clr, lenses, $
               slength=slength, $
               rfac=rfac, $
               gridsize=gridsize, $
               stepfac=stepfac, $
               yoverx=yoverx, $
               scat=scat, $
               fgal=fgal, $
               donoise=donoise, $
               check = check, $
               write=write, $
               allign=allign, $
               abs=abs, $
               rotate=rotate, surface=surface, $
               noprompt=noprompt, verbose=verbose, ngood=ngood, $
               outdir=outdir, $
               no_ps=no_ps, $
               kapfit=kapfit, kerrfit=kerrfit, densfit=densfit, $
               wsource=wsource, $
               _extra=extra

  IF n_params() EQ 0 THEN BEGIN
      print,'-Syntax:  kappa_map, run1, run2, clr, lenses, '
      print,'           slength=slength, '
      print,'           rfac=rfac, '
      print,'           gridsize=gridsize, '
      print,'           stepfac=stepfac, '
      print,'           yoverx=yoverx, '
      print,'           scat=scat, '
      print,'           donoise=donoise, '
      print,'           check = check, '
      print,'           write=write, '
      print,'           allign=allign, '
      print,'           rotate=rotate, surface=surface, '
      print,'           noprompt=noprompt, verbose=verbose, ngood=ngood, '
      print,'           outdir=outdir, '
      print,'           no_ps=no_ps, '
      print,'           wsource=wsource, '
      print,'           _extra=extra'
      return
  ENDIF 

  time=systime(1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some Parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  COMMON seed,seed
  check = 0
  IF NOT keyword_set(write) THEN write=0
  IF NOT keyword_set(noprompt) THEN noprompt=0
  IF NOT keyword_set(surface) THEN surface = 0
  IF NOT keyword_set(donoise) THEN donoise=0
  IF NOT keyword_set(allign) THEN allign=0
  IF NOT keyword_set(abs) THEN abs=0
  IF NOT keyword_set(no_ps) THEN no_ps=0

  IF n_elements(fgal) NE 0 THEN dodens = 1 ELSE dodens = 0
  IF n_elements(verbose) EQ 0 THEN verbose = 2 ;Maximum verbosity
  IF verbose EQ 2 THEN silent = 0 ELSE silent = 1
  IF n_elements(outdir) EQ 0 THEN outdir = $
                                  '/sdss4/data1/esheldon/TMP/'
  IF n_elements(stepfac) EQ 0 THEN stepfac = 2.
  IF n_elements(slength) EQ 0 THEN slength = 120. ;arcsec
  IF n_elements(rfac) EQ 0 THEN rfac = 10. ; Gets >90% of S/N
  rmax = rfac*slength
  step = slength/stepfac

  IF n_elements(gridsize) EQ 0 THEN gridsize=1000. ;size of long direction
  IF n_elements(yoverx) EQ 0 THEN BEGIN
      xgridsize = gridsize & ygridsize = xgridsize
  ENDIF ELSE BEGIN
      xgridsize = gridsize & ygridsize = double(yoverx)*xgridsize
  ENDELSE 
  nx = long(xgridsize/step)
  ny = long(ygridsize/step)
  IF nx MOD 2 EQ 0 THEN nx = nx+1
  IF ny MOD 2 EQ 0 THEN ny = ny+1

  r1str = ntostr(run1)
  r2str = ntostr(run2)
  colors = ['u','g','r','i','z']

  step = 300L

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Get the catalogs
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
  nscat =  n_elements(scat)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check if certain tags are defined for structure 'lenses'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  a=tag_names(lenses)

  wz = where(a EQ 'Z', nz)      ;Check for redshift tag
  IF nz EQ 0 THEN BEGIN
      wz = where(a EQ 'PHOTOZ', nz)
  ENDIF
  wn = where(a EQ 'NAME', nwn)  ;Check for a tag called "NAME"
  IF nwn NE 0 THEN name = ntostr(lenses[0].name) ELSE name = 'galaxies'
  we = where(a EQ 'E1', nwe)
  IF nwe NE 0 THEN lense = 1 ELSE lense = 0
  IF lense AND allign THEN donoise = 0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Define positions of lens centers and ranges around them
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ntot = n_elements(lenses)
  cenx = lenses.dec
  ceny = lenses.ra

  ;; convert to degrees
  s = slength/3600.
  radmax = rmax/3600.
  xsize = xgridsize/3600.
  ysize = ygridsize/3600.

  gfac = 1.
  IF (lense AND allign) THEN BEGIN
      IF verbose EQ 2 THEN BEGIN 
          print
          print,'Alligning grids with lens ellipticities'
      ENDIF 
      gfac = sqrt(2.)
      theta = .5*atan( lenses.e2, lenses.e1 )
      cos = cos(theta)
      sin = sin(theta)
  ENDIF 
  ;; sqrt(2.) is overestimate in most cases.
  xmax  = cenx + xsize/2.*gfac  
  xmin  = cenx - xsize/2.*gfac
  ymax  = ceny + ysize/2.*gfac
  ymin  = ceny - ysize/2.*gfac

  gridxmax = cenx + xsize/2.
  gridxmin = cenx - xsize/2.
  gridymax = ceny + ysize/2.
  gridymin = ceny - ysize/2.

  maxscatra = max(scat.ra)
  minscatra = min(scat.ra)
  maxscatdec = max(scat.dec)
  minscatdec = min(scat.dec)

  ly1 = maxscatra - (ymax+radmax) ; > 0 to accomodate grid.
  ly2 = (ymin-radmax) - minscatra
  lx1 = maxscatdec - (xmax+radmax)
  lx2 = (xmin-radmax) - minscatdec  

  tol = -20./3600.

  goodlenses = where( (ly1 GE tol) AND $
                      (ly2 GE tol) AND $
                      (lx1 GE tol) AND $
                      (lx2 GE tol) , nlens)


  IF nlens EQ 0 THEN BEGIN      
      print,'Cannot accomodate grid'
      check = 1
      return
  ENDIF ELSE BEGIN
      ndiff = ntot - nlens
      IF (ndiff NE 0) AND (verbose NE 0) THEN BEGIN
          print,'Using ',ntostr(nlens),' lenses'
          print,'Threw out ',ntostr(ndiff),' Lenses for Accomodation'
      ENDIF 
  ENDELSE 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up the output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sstr = ntostr(long(slength))
  rstr = ntostr(long(rmax))
  mid = '_S'+sstr+'_R'+rstr
  IF write  THEN BEGIN 
      IF nlens EQ 1 THEN BEGIN
          IF nz NE 0 THEN zstr='_Z'+ntostr(lenses.(wz[0]),4) ELSE zstr = ''
          prename=outdir+name[0]+mid+zstr
      ENDIF ELSE BEGIN
          IF name EQ 'galaxies' THEN BEGIN
              prename=outdir+'galaxies'+mid
          ENDIF ELSE BEGIN 
              prename=outdir+'lensesgrp'+mid
          ENDELSE 
      ENDELSE 

      add1 = '_'+colors[clr]+'_N.ps'
      add2 = '_'+colors[clr]+'_N.fit'
      add3 = '_'+colors[clr]+'_N.dat'
      continue = 1
      WHILE continue DO BEGIN 
          add1=newname(add1) & add2=newname(add2) & add3=newname(add3)
          psfile  = prename+add1          & rtfile  = prename+'_rot'+add1
          kapfit  = prename+'_kappa'+add2 & kerrfit = prename+'_kerr'+add2
          radfit  = prename+'_rad'+add2   & rerrfit = prename+'_rerr'+add2
          densfit = prename+'_dens'+add2  & datfile = prename+add3
          continue = exist(kapfit)
      ENDWHILE 
      print,'Kappa File: ',kapfit
  ENDIF 
  IF NOT keyword_set(rotate) THEN rotate=0
  IF noprompt THEN rotate=0


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Group the lenses for speed.  This is faster if lenses are sorted by ra.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF verbose EQ 2 THEN BEGIN 
      print
      IF dodens THEN print,' Finding Smoothed Galaxy Density'
      print,' x grid size: ',xgridsize,'  y grid size: ',ygridsize
      print,' x grid spacing: ',xgridsize/(nx-1),$
            ' y grid spacing: ',ygridsize/(ny-1)
      print,' nx: ',ntostr(nx),'  ny: ',ntostr(ny)
      print,' Smoothing length:  ',ntostr(slength),$
        ' rmax: ',rmax
      print,' Weight Max:  ', ntostr(slength*1.8938227)
      IF donoise THEN print,' Randomizing'
      print
  ENDIF 

  step = 300L                   ;Number of lenses in group
  nstepOld = nlens/step
  nstep = nstepOld
  left = nlens MOD step

  ;; To account for leftover stuff
  IF nstepOld EQ 0 THEN BEGIN
      nstepOld = -10
      step = left
      nstep = 1
  ENDIF ELSE BEGIN
      IF left NE 0 THEN nstep = nstepOld + 1
  ENDELSE 

  wlens = goodlenses
  skip = 0
  FOR group = 0L, nstep-1 DO BEGIN ; Loop over groups

      IF group EQ nstepOld THEN BEGIN ; Last group
          ii = wlens( group*step: group*step+left-1  )
          step = left
      ENDIF ELSE ii = wlens[ group*step : (group+1)*step -1 ]

      tmpmaxy = max( ymax[ii] ) ; Select region of interest for this group
      tmpminy = min( ymin[ii] )
      tmpmaxx = max( xmax[ii] )
      tmpminx = min( xmin[ii] )

      wsrc = where( scat.ra LE (tmpmaxy + radmax) AND $
                    scat.ra GE (tmpminy - radmax), nwsrc)
      IF dodens THEN wdens1 = where(fgal.ra LE (tmpmaxy + radmax) AND $
                                    fgal.ra GE (tmpminy - radmax), ndens1)
      IF (nwsrc NE 0) THEN BEGIN
          wsrc2 = where( scat[wsrc].dec LE (tmpmaxx + radmax) AND $
                         scat[wsrc].dec GE (tmpminx - radmax), nwsrc2)

          IF nwsrc2 NE 0 THEN BEGIN 
              wsrc = wsrc[ wsrc2 ]

              IF dodens THEN BEGIN
                  IF ndens1 NE 0 THEN BEGIN
                      wdens = where($
                               fgal[wdens1].dec LE (tmpmaxx + radmax) AND $
                               fgal[wdens1].dec GE (tmpminx - radmax), ndens2)
                      IF ndens2 NE 0 THEN sendf = fgal[wdens1[wdens]]
                  ENDIF 
              ENDIF 

              IF verbose NE 0 THEN print,'Lens Group = ',$
                ntostr(group+1)+'/'+ntostr(nstep)
              FOR gi = 0, step-1 DO BEGIN 

                  jj = ii[gi]
                  
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; Define the grid
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  gridx = arrscl( dindgen(nx), gridxmin[jj], gridxmax[jj] )
                  gridy = arrscl( dindgen(ny), gridymin[jj], gridymax[jj] )

                  IF nz NE 0 THEN BEGIN
                      zstr = ntostr(lenses[jj].(wz[0]) )
                      dstr = ntostr(angdist_lambda(lenses[jj].z))
                  ENDIF ELSE BEGIN
                      zstr = ''
                      dstr = ''
                  ENDELSE 

                  IF verbose EQ 2 THEN BEGIN 
                      ;print
                      print,' Lens Name: ',name
                      print,' Lens Redshift~ ',zstr,'  Distance: ',dstr,' Mpc'
                      print,' BCG   DEC: ',ntostr(cenx[jj]),'  RA: ',$
                        ntostr(ceny[jj])
                      ;print
                  ENDIF 
                  ymx = ymax[jj] + radmax
                  ymn = ymin[jj] - radmax
                  xmx = xmax[jj] + radmax
                  xmn = xmin[jj] - radmax

                  w1 = where(scat[wsrc].ra GE ymn AND $
                             scat[wsrc].ra LE ymx, nw1)
                  IF nw1 NE 0 THEN BEGIN
                      w1 = wsrc[w1]
                      w2 = where(scat[w1].dec GE xmn AND $
                                 scat[w1].dec LE xmx, nw)
                      IF nw NE 0 THEN BEGIN
                          w = w1[w2]
                         
                          ;; Calculate the signal.
                          IF lense AND allign THEN BEGIN 
                              rotgridkappa, scat[w].e1, scat[w].e2, $
                                scat[w].uncert, $
                                scat[w].dec, scat[w].ra, $
                                gridx, gridy, $
                                cenx[jj], ceny[jj], $
                                cos[jj], sin[jj], $
                                radmax, s, $
                                kappa, radial, wsum, npair, $
                                kerr1, kerr2, rerr1, rerr2, err3, $
                                silent=silent
                          ENDIF ELSE BEGIN
                              IF donoise THEN BEGIN
                                  
                                  ;; Scramble ellipticities
                                  a=randomu(seed, nw)
                                  new_w=w[ sort(a) ]
                                  
                                  ;; Randomly rotate source ellipticities.
                                  rth = arrscl(randomu(seed,nw), $
                                               0, !pi, arrmin=0., arrmax=1.)
                                  randcos = cos(2.*rth)
                                  randsin = sin(2.*rth)
                                  
                                  se1 =  scat[new_w].e1*randcos + $
                                         scat[new_w].e2*randsin
                                  se2 = -scat[new_w].e1*randsin + $
                                         scat[new_w].e2*randcos
                                  gridkappa, se1, se2, $
                                    scat[new_w].uncert, $
                                    scat[w].dec, scat[w].ra, $
                                    gridx, gridy, radmax, s, $
                                    kappa, radial, wsum, npair, $
                                    kerr1, kerr2, rerr1, rerr2, err3, $
                                    fgal=sendf, dens=dens, $
                                    silent=silent
                              ENDIF ELSE BEGIN 
                                  gridkappa, scat[w].e1, scat[w].e2, $
                                    scat[w].uncert, $
                                    scat[w].dec, scat[w].ra, $
                                    gridx, gridy, radmax, s, $
                                    kappa, radial, wsum, npair, $
                                    kerr1, kerr2, rerr1, rerr2, err3, $
                                    fgal=sendf, dens=dens, $
                                    silent=silent
                              ENDELSE 
                          ENDELSE 
                          
                      ENDIF ELSE skip=1
                  ENDIF ELSE skip=1
                  IF skip EQ 1 THEN BEGIN
                      skip = 0
                      IF nlens EQ 1 THEN BEGIN
                          print,'No Data!'
                          check=1
                          return
                      ENDIF ELSE BEGIN
                          IF verbose EQ 2 THEN print,$
                                            'No data, skipping this lens'
                          IF n_elements(badlenses) EQ 0 THEN badlenses=[jj] $
                          ELSE badlenses = [badlenses, jj]
                      ENDELSE 
                  ENDIF 
              ENDFOR 
          ENDIF ELSE skip=2
      ENDIF ELSE skip=2
      IF skip EQ 2 THEN BEGIN
          skip = 0
          IF n_elements(badlenses) EQ 0 THEN badlenses = ii $
          ELSE badlenses = [badlenses, ii]
          IF verbose NE 0 THEN print,'No data for this lens group'
      ENDIF 
  ENDFOR 
  nbad = n_elements(badlenses)
  IF nbad NE 0 THEN remove, badlenses, goodlenses
  ngood = n_elements(goodlenses)
  IF (nlens NE 1) AND (verbose NE 0) THEN BEGIN 
      print
      print,'Used ',ntostr(ngood),' lenses'
      print,'Threw out ',ntostr(nbad),' lenses for lack of data'
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Calculate final values
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  IF n_elements(kappa) EQ 0 THEN BEGIN
      print,'No lenses could fit in grid'
      return
  ENDIF 

  area = !pi*radmax^2

  ;; revaluate
  wcheck=where(dens NE 0, ncheck)
  print,'Ncheck = ',ncheck
  IF ncheck NE 0 THEN dodens = 1 ELSE dodens = 0
;  IF dodens THEN BEGIN 
;      dlength = s
;      darea = !pi*(3.*dlength)^2
;      dens = npair/darea
;  ENDIF 

  kappa  = area*kappa/wsum      ; Signal
  radial = area*radial/wsum

  kerr1 = area^2*kerr1          ; Calculate the Uncertainty
  kerr2 = -2.*area*kappa*kerr2

  rerr1 = area^2*rerr1
  rerr2 = -2.*area*radial*rerr2

  kerr3  = kappa^2*err3
  rerr3  = radial^2*err3

  kerror = sqrt( kerr1 + kerr2 + kerr3 )/wsum
  rerror = sqrt( rerr1 + rerr2 + rerr3 )/wsum

  ksig = kappa/kerror
  rsig = radial/rerror

  IF verbose NE 0 THEN BEGIN 
      print
      print,'Mean kerror: ',mean(kerror)
      print,'Mean rerror: ',mean(rerror)
      print,'Max kappa S/N = ',max(ksig)
      print,'Max radial S/N = ',max(rsig)
      print
  ENDIF 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
;
; Make some outputs
;
;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF write THEN BEGIN 
      writefits, kapfit, kappa
      writefits, radfit, radial
      writefits, kerrfit, kerror
      writefits, rerrfit, rerror
      IF dodens THEN writefits, densfit, dens
  ENDIF 

  IF verbose NE 0 THEN ptime,systime(1)-time

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Print some general info.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF write THEN BEGIN           
      openw, lun1, datfile, /get_lun
      !textunit = lun1

      printf, lun1, 'Lenses Used: ',ngood
      printf, lun1, 'Lenses Thrown out: ',nbad
      printf, lun1

      IF nz NE 0 AND nwn NE 0 THEN BEGIN 
          goodnames = ntostr(lenses[goodlenses].name)
          IF nz NE 0 THEN BEGIN
              goodz  = ntostr(lenses[goodlenses].z, 5)
          ENDIF ELSE BEGIN
              goodz  = strarr(ngood)
          ENDELSE 

          message = '   name      z     '
          printf, lun1, message
          forprint, goodnames, goodz, /silent, TEXT=5
      ENDIF 
      close, lun1
      free_lun,  lun1
  ENDIF 


  IF write OR (NOT noprompt) THEN doplots = 1 ELSE doplots=0
  IF write AND no_ps THEN doplots = 0
  
  IF doplots THEN BEGIN ;Check if we should even make the plot

      IF nz NE 0 THEN BEGIN 
          mz = mean_check( lenses[goodlenses].(wz[0]) )
          zstr = ntostr(mz, 5)
          zdist = angdist_lambda(mz) ;in Mpc
          zdist = zdist*1000.   ;kpc
      ENDIF ELSE BEGIN
          zdist = 0. & zstr = '0'
      ENDELSE 

      IF write THEN begplot,name=psfile,/invbw

      IF nlens EQ 1 THEN BEGIN 
          addstring='  Z='+zstr+'  '+name+' s='+ntostr(long(slength))
      ENDIF ELSE BEGIN 
          addstring='  Z='+zstr+'  '+'lensesgrp'+' s='+ntostr(long(slength))
      ENDELSE 
      addstring = addstring

      ;; Define grid and axis labels.
      IF (NOT allign) AND abs THEN BEGIN 
          cx = gridx
          cy = gridy
      ENDIF 
      usex = (gridx - cenx[jj])*3600.
      usey = (gridy - ceny[jj])*3600.
      ntick=5
  
      xangle=strarr(ntick)
      xdist = xangle

      minx=min(gridx)
      maxx=max(gridx)
      stepx = (maxx-minx)/(ntick-1)
      FOR i=0, ntick-1 DO BEGIN
          xangle[i]=ntostr( (minx + i*stepx-cenx[jj])*3600., 6)
          xdist[i] =ntostr( (minx + i*stepx-cenx[jj])*!pi/180.*zdist, 6)
      ENDFOR 

      yangle=strarr(ntick)
      ydist = yangle

      miny=min(gridy)
      maxy=max(gridy)
      stepy=(maxy-miny)/(ntick-1)
      FOR i=0, ntick-1 DO BEGIN
          yangle[i]=ntostr( ( miny + i*stepy -ceny[jj])*3600., 6)
          ydist[i] =ntostr( (miny + i*stepy-ceny[jj])*!pi/180.*zdist, 6)
      ENDFOR 

      ;; Contour levels (sigma levels)
      nnl = 20
      IF max(ksig) GT 5 THEN BEGIN
          ltmp = indgen(nnl)+1  ; 1 sig levels
          levels = [-reverse(ltmp), ltmp]
      ENDIF ELSE BEGIN
          ltmp = indgen(nnl)/2. + .5 ;half sig levels
          levels = [-reverse(ltmp), ltmp]
      ENDELSE 

      c_labels=levels*0. + 1.
      c_linestyle = (levels LT 0.)
      c_colors = c_linestyle
      c_colors = c_colors*(!d.n_colors-1)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Make the plots
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
      style=4+1                 ; This uses exact axis sizes and also
                                ; suppresses the axes
      
      xt1='kpc'
      xt2='Kappa  '+addstring
      yt1='kpc'
      yt2='arcsec'
      image_contour, ksig, $
                 ntick=ntick, xtickn1=xdist, xtickn2=xangle, $
                 ytickn1=ydist, ytickn2=yangle, $
                 xtitle1 = xt1, xtitle2=xt2, $
                 ytitle1 = yt1, ytitle2=yt2, $
                 levels=levels, c_linestyle=c_linestyle, c_labels=c_labels, $
                 c_colors=c_colors, position=position
      IF abs AND (NOT allign) THEN BEGIN
          IF (NOT write) AND (NOT noprompt) THEN key=get_kbrd(1)
          contour, ksig, cx, cy, ystyle=1, xstyle=1, $
                 title = xt2, $
                 levels=levels, c_linestyle=c_linestyle, $
                 c_labels=c_labels, $
                 position=position
      ENDIF 
      IF surface THEN BEGIN 
          IF (NOT write) AND (NOT noprompt) THEN key=get_kbrd(1)
          surface, kappa, usex, usey, $
            title='Kappa'+addstring,charsize=2., _extra=extra
      ENDIF 
      
      IF (NOT write) AND (NOT noprompt) THEN key=get_kbrd(1)
      xt2 = 'Radial  '+addstring
      image_contour, rsig, $
                 ntick=ntick, xtickn1=xdist, xtickn2=xangle, $
                 ytickn1=ydist, ytickn2=yangle, $
                 xtitle1 = xt1, xtitle2=xt2, $
                 ytitle1 = yt1, ytitle2=yt2, $
                 levels=levels, c_linestyle=c_linestyle, c_labels=c_labels, $
                 c_colors=c_colors,position=position

      IF abs AND (NOT allign) THEN BEGIN
          IF (NOT write) AND (NOT noprompt) THEN key=get_kbrd(1)
          contour, rsig, cx, cy, ystyle=1, xstyle=1, $
                 title = xt2, $
                 levels=levels, c_linestyle=c_linestyle, $
                 c_labels=c_labels, $
                 position=position
      ENDIF 

      IF surface THEN BEGIN 
          IF (NOT write)  AND (NOT noprompt) THEN key=get_kbrd(1)
          surface, radial, usex, usey, title='Radial'+addstring, charsize=2., $
            _extra=extra
      ENDIF 

      IF dodens THEN BEGIN 
          IF (NOT write) AND (NOT noprompt) THEN key=get_kbrd(1)
          xt2 = 'Smoothed Density  '+addstring
          ;; First attempt at this number
          dvar = sdev(dens)
          md = median(dens)
          image_contour, (dens-md)/dvar, $
                 ntick=ntick,xtickn1=xdist,xtickn2=xangle, $
                 ytickn1=ydist, ytickn2=yangle, $
                 xtitle1 = xt1, xtitle2=xt2, $
                 ytitle1 = yt1, ytitle2=yt2, $
                 levels=levels, c_linestyle=c_linestyle, c_labels=c_labels, $
                 c_colors=c_colors
  
          IF surface THEN BEGIN
              IF (NOT write)  AND (NOT noprompt) THEN key=get_kbrd(1)
              surface, dens, usex, usey, title='Density'+addstring, $
                charsize=2., $
                _extra=extra
          ENDIF 
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Plot source gal distribution, annis' clusters and position of 
      ;; our cluster
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF nlens EQ 1 THEN BEGIN 
          IF (NOT write) AND (NOT noprompt) THEN key=get_kbrd(1)
          print
          print,ntostr(nw),' Background galaxies found in ', $
            ntostr( 2*rmax+gridsize ),' arcsecond box'
          plot, scat[w].dec, scat[w].ra, psym=3, ystyle=1, xstyle=1, $
            xrange=[xmin[goodlenses]-radmax,xmax[goodlenses]+radmax],$
            yrange=[ymin[goodlenses]-radmax,ymax[goodlenses]+radmax],$
            position=position
          plot_box, xmin[goodlenses], xmax[goodlenses], $
                    ymin[goodlenses], ymax[goodlenses]
          oplot, [cenx[jj]], [ceny[jj]], psym=4, symsize=2.
      ENDIF 

;      IF nlens LE 100 THEN BEGIN 
;          IF (NOT write) AND (NOT noprompt) THEN key=get_kbrd(1)
;          rdannis, /silent          ;obsolote?
;          oplot, [lenses[goodlenses].dec],[lenses[goodlenses].ra],$
;            psym=6,symsize=2 
;          IF nlens EQ 1 THEN plot_box, xmin[goodlenses], xmax[goodlenses], $
;                                       ymin[goodlenses], ymax[goodlenses]
;      ENDIF 
      
      IF write THEN endplot,/noprint

      IF (NOT write) AND (NOT noprompt) THEN key=get_kbrd(1)
      IF rotate THEN BEGIN 
          IF write THEN BEGIN 
              rotate_plot, kappa, usex, usey, step=5., $
                title='Kappa'+addstring,psfile=rtfile,charsize=2., $
                _extra=extra
          ENDIF ELSE BEGIN 
              key=get_kbrd(1)
              rotate_plot, kappa, usex, usey, step=5., $
                title='Kappa'+addstring,charsize=2., _extra=extra
          ENDELSE 
      ENDIF 
  ENDIF                         ;Checked if should plot

  wsource = wsrc

  return
END 
