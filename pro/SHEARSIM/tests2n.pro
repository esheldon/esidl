PRO tests2n, diffpos, s2n, maxx=maxx, maxy=maxy

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time = systime(1)

  sigma = 1200.                 ; km/s

  ;; number of times for each bin.
  nsample = 10L

  ;; Number of binsizes and Number of rmax
  nbsize = 10L
;  nbsize=2L
  nrmax = nbsize
  nx = 61L
  ny = nx

  gridsize = 401.
  fac = gridsize/nx

  biggestr = 1200.
  smallestr = 600.
  addmin = 10.5
  rstep = (biggestr - smallestr)/(nrmax-1)
  rmax = findgen(nrmax)*rstep+smallestr+addmin
  print,rmax

  bmax = 400.
  bmin = 50.
  bstep = (bmax-bmin)/(nbsize-1)

  binsize = findgen(nbsize)*bstep+bmin
  print,binsize

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  dir = '/sdss4/data1/esheldon/GRIDSHEAR/'
  psname = dir+'tests2n_N1.ps'
  s2nfitfile = dir+'tests2n_N1.fit'
  posfitfile = dir+'testposerr_N1.fit'
  WHILE exist(psname) OR exist(s2nfitfile) OR exist(posfitfile) DO BEGIN
      psname         = newname(psname)
      s2nfitfile = newname(s2nfitfile)
      posfitfile     = newname(posfitfile)
  ENDWHILE 

  print,'Ps file: ',psname
  print,'S/N map   : ',s2nfitfile
  print,'Poserr    : ',posfitfile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Define some arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  maxx = fltarr(nbsize, nrmax)
  maxy = maxx
  s2n = maxx
  diffpos = maxx

  tmpx = fltarr(nsample)
  tmpy = tmpx
  tmpmt = tmpx

  index = lindgen(nx*ny)
  makegauss, filter, [3,3], 3.0

  FOR nb=0, nbsize-1 DO BEGIN
      FOR nr=0, nrmax-1 DO BEGIN 
          FOR j=0, nsample-1 DO BEGIN 

              ;; First make a measure of the noise

              IF j EQ 0 THEN BEGIN 

                  print,'Measuring The Noise'
                  simgridshear, imtan2, imrad2,sigma=100., $
                    binsize=binsize[nb], rmax=rmax[nr], nx=nx, ny=ny, $
                    /noise, /noprompt
                  
                  noise1 = sqrt(variance(imtan2))
                  simgridshear, imtan2, imrad2, sigma=100., $
                    binsize=binsize[nb], rmax=rmax[nr], nx=nx, ny=ny, $
                    /noise, /noprompt
                  
                  noise2 = sqrt(variance(imtan2))
                  noise = (noise1 + noise2)/2.

              ENDIF 
                  
              simgridshear, imtan, imrad, sigma=sigma, $
                binsize=binsize[nb], rmax=rmax[nr], nx=nx, ny=ny, $
                /noise, /noprompt
          
              t=convol( imtan, filter, /edge_truncate)
              ii1 = (nx-1)/4.
              ii2 = (ny-1)/4.
              tmpmt[j] = max( t[ii1:3*ii1, ii2:3*ii2] )
              w=where(t EQ tmpmt[j], nw)

                      
              tmpx[j] = index[w[0]] MOD nx
              tmpy[j] = index[w[0]]/nx

              print
              print,'At Sample: ',j
              print,'Found middle at: ',tmpx[j]*fac,tmpy[j]*fac
              print

          ENDFOR
          IF nsample EQ 1 THEN BEGIN 
              maxx[nb, nr] = tmpx[0]
              maxy[nb, nr] = tmpy[0]
          ENDIF ELSE BEGIN 
              maxx[nb, nr] = mean( tmpx )
              maxy[nb, nr] = mean( tmpy )
          ENDELSE 

          print,'Average Center:  ',fac*maxx[nb, nr],fac*maxy[nb, nr]

          ;; Some measure of signal to noise in map.
          s2n[nb, nr] = median( tmpmt )/noise
      ENDFOR 
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; gridsize = 400. arcseconds across, using 60x60 pixel grid
  ;; maxx, maxy are in pixels so each pixel is 400./60. = 6.67 arcsec
  ;;
  ;; center is at .5 degrees, but our box doesn't cover the whole area
  ;; so center for us is at gridsize/2.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  diffx = fac*(maxx - (nx-1)/2.)
  diffy = fac*(maxy - (ny-1)/2.)

  diffpos = sqrt( diffx^2 + diffy^2 )

  ;; Open postscript file
  begplot,name=psname

  ntick = 2L*nbsize+1

  xdummy = strarr(ntick)
  ydummy = xdummy
  xtickn = xdummy
  ytickn = xdummy
  minbin=min(binsize)
  maxbin=max(binsize)
  minrmax = min(rmax)
  maxrmax = max(rmax)
  stepx = ( maxbin - minbin )/(nbsize-1)
  stepy = ( maxrmax - minrmax )/(nbsize-1)

  kk = 0
  FOR j=0, ntick-1 DO BEGIN
      IF j MOD 2 EQ 0 THEN BEGIN 
          xtickn[j] = ' '
          ytickn[j] = ' '
      ENDIF ELSE BEGIN 
          xtickn[j]=ntostr( long( minbin + kk*stepx  ) )
          ytickn[j]=ntostr( long( minrmax + kk*stepy ) )
          kk=kk+1
      ENDELSE 
      xdummy[j] = ' '
      ydummy[j] = ' '
  ENDFOR 


  tit = 'S/N map'
  tvim2,s2n, title=tit,/noframe
  axis, xaxis=0, xticks=ntick-1, xtickn=xtickn, xtitle=xtit
  axis, xaxis=1, xticks=ntick-1, xtickn=xdummy
  axis, yaxis=0, yticks=ntick-1, ytickn=ytickn, ytitle=ytit
  axis, yaxis=1, yticks=ntick-1, ytickn=ydummy
  

  tit = 'Pos Diff'
  tvim2,diffpos, title=tit,/noframe
  axis, xaxis=0, xticks=ntick-1, xtickn=xtickn, xtitle=xtit
  axis, xaxis=1, xticks=ntick-1, xtickn=xdummy
  axis, yaxis=0, yticks=ntick-1, ytickn=ytickn, ytitle=ytit
  axis, yaxis=1, yticks=ntick-1, ytickn=ydummy
  


  endplot,/noprint
  
  mwrfits, s2n, s2nfitfile, /create
  mwrfits, diffpos, posfitfile, /create


  ptime,systime(1)-time
  return
END 
