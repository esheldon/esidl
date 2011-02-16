pro binner2, xi,yi,ysig,xo,yo,sig,bin
  
  histo=histogram(xi,bin=0.01,reverse_indices=r)
  nn=n_elements(histo)
  xo=fltarr(nn)
  yo=fltarr(nn)
  sig=fltarr(nn)
  kk=0
  for i=0,nn-1 do begin
      if(r(i+1)-r(i) gt 2) then BEGIN
          
          wmom, yi[ r(r(i):r(i+1)-1) ], ysig[ r(r(i):r(i+1)-1) ], $
            ymean, ysdev, yerr
          yo[kk] = ymean
          sig[kk] = yerr

          result=moment(xi(r(r(i):r(i+1)-1)))
          xo[kk]=result[0]
          kk=kk+1

      endif
  endfor
  yo=yo[0:kk-1]
  xo=xo[0:kk-1]
  sig=sig[0:kk-1]
  
  return
end

pro binner, xi,yi,xo,yo,sig,bin
  
  histo=histogram(xi,bin=0.01,reverse_indices=r)
  nn=n_elements(histo)
  xo=fltarr(nn)
  yo=fltarr(nn)
  sig=fltarr(nn)
  kk=0
  for i=0,nn-1 do begin
    if(r(i+1)-r(i) gt 2) then begin
      result=moment(yi(r(r(i):r(i+1)-1)))
      yo[kk]=result[0]
      sig[kk]=sqrt(result[1]/(r(i+1)-r(i)))
      result=moment(xi(r(r(i):r(i+1)-1)))
      xo[kk]=result[0]
;    print,xo[kk],yo[kk],sig[kk],r(i+1)-r(i)
      kk=kk+1
    endif
  endfor
  yo=yo[0:kk-1]
  xo=xo[0:kk-1]
  sig=sig[0:kk-1]
  
  return
  
end

pro test_rho4_corshape, run, rerun, camcol, clr, iorder, start=start, nframes=nframes,$
                        minmag=minmag, maxmag=maxmag, $
                        rcut=rcut, ermscut=ermscut, errmax=errmax, printfile=printfile

;   ermscut = cut off frames with psfe rms gt ermscut
  
  IF n_params() EQ 0 THEN BEGIN 
      print,'-syntax corshape, run, rerun, camcol, clr, iorder, start=start, nframes=nframes, minmag=minmag, maxmag=maxmag, rcut=rcut, ermscut=ermscut, printfile=printfile'
      return
  ENDIF 
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Written by Phil Fischer.  Adapted to new file formats E.S.S.
; Also, since stars are identified in the adat files, can save 
; factor of ~two in time by not reading in stars separately
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time=systime(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  newfront = 'adatc'            ;The front of files with adaptive moments

  colors=['u','g','r','i','z']
  bands = [1,2,3]               ; we will process g,r,i
  nband = n_elements(bands)
  
  rstr = run2string(run)
  rrstr = ntostr(rerun)
  cstr = ntostr(camcol)

  IF n_elements(maxmag) EQ 0 THEN maxmag=22.
  IF n_elements(minmag) EQ 0 THEN minmag=0.
  IF n_elements(rcut) EQ 0 THEN rcut=0.8
  IF n_elements(ermscut) EQ 0 THEN ermscut = 1000.
  IF n_elements(errmax) EQ 0 THEN errmax = 1000.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Setup the column
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  fetch_dir, run, camcol, rerun, dir, atldir, corrdir=corrdir,$
    corratldir = fitdir
  fetch_file_list, dir, files, fnums, start=start, nframes=nframes, $
                   fieldmin=fieldmin, fieldmax=fieldmax
  ntot = fieldmax-fieldmin+1

  IF n_elements(printfile) EQ 0 THEN BEGIN
      printfile = fitdir+'testcorshape_'+rstr+'_'+cstr+'_'+colors[clr]+'_N1.ps'
      WHILE exist(printfile) DO printfile=newname(printfile)
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up new filenames. 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  rename_tsobj, files, corrdir, newfront, adatcfiles, nchar

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; setup moment arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  momname, run, camcol, clr, iorder, psffile1
  testfile = fitdir+'test'+repstr(psffile1,'.txt','.fit')
  psffile1 = fitdir + psffile1
  print,psffile1

  print
  print,'Reading in moment arrays'
  readcol,psffile1,sjunk,junk,a,b,c,d,e,f,format='A,I,D,D,D,D,D,D,D'
  testfit = mrdfits(testfile, 1)
  
  galflag   = 0
  starflag  = 0
  ebv       = fltarr(nframes)
  ff        = fltarr(nframes)
  psfmean   = fltarr(nframes)
  testpsfmean = fltarr(nframes)
  psfsig    = fltarr(nframes)
  psfe1mean = fltarr(nframes)
  psfe1sig  = fltarr(nframes)
  psfe2mean = fltarr(nframes)
  psfe2sig  = fltarr(nframes)
  foredens  = fltarr(nframes)
  gl        = fltarr(nframes)
  gb        = fltarr(nframes)
  
  print,'-----------------------------------------------------'
  print,' Run: ',rstr,' Rerun: ',rrstr,' Camcol: ',cstr,' Band: ',colors[clr]
  print,'-----------------------------------------------------'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over fields
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR ic = 0L, nframes-1 DO BEGIN 
      infile = adatcfiles[ic]
      field = fnums[ic]
      fstr = ntostr(field)
      ff[ic]=field

      ind  = ic*4
      ind1 = ind+1
      ind2 = ind+2
      ind3 = ind+3

      ;; Fields have already been trimmed of overlap.  Don't need to
      ;; use read_tsobjmin
      openr, lun, infile, /get_lun
      IF ic EQ 0 THEN BEGIN 
          pstruct = mrdfits3(lun, 1, 0, /silent)
      ENDIF ELSE BEGIN
          pstruct = mrdfits3(lun, 1, 0, /silent, /deja_vu)
      ENDELSE 
      free_lun, lun

                                ;Initialize arrays for field
      nps = n_elements(pstruct)
      gal = lindgen(nps)

      x = pstruct.colc[clr]
      y = pstruct.rowc[clr]
      objsize = pstruct.ixx[clr]+pstruct.iyy[clr]
      testobjsize = pstruct.petrorad[clr]^2

      psfsize = a[ind] + b[ind]*x   $
                       + c[ind]*y   $
                       + d[ind]*x*y $
                       + e[ind]*x*x $
                       + f[ind]*y*y

      testpsfsize = testfit[ic].szcoeff[0] + testfit[ic].szcoeff[1]*x $
                                           + testfit[ic].szcoeff[2]*y $
                                           + testfit[ic].szcoeff[3]*x*y $
                                           + testfit[ic].szcoeff[4]*x*x $
                                           + testfit[ic].szcoeff[5]*y*y

      psfrho4 = a[ind3] + b[ind3]*x   $
                        + c[ind3]*y   $
                        + d[ind3]*x*y $
                        + e[ind3]*x*x $
                        + f[ind3]*y*y

      stars = where(pstruct.starflag[clr] EQ 1., nstars)
      IF nstars NE 0 THEN BEGIN
          remove, stars, gal
      ENDIF 
      
      zmag = pstruct.petrocounts[clr]-pstruct.reddening[clr]
                                ;Throw out stuff we don't want for gals
      ww=where( (objsize[gal] gt 0) and (zmag[gal] le maxmag) $ 
                and (zmag[gal] ge minmag) AND $
                (pstruct[gal].petroraderr[clr] GT 0.) AND $ ;!ADDED FOR TESTING
                (pstruct[gal].petrorad[clr] GT 0.) AND $
                (pstruct[gal].momerr[clr] LE errmax), nww)

      IF nww NE 0 THEN BEGIN
          gal = gal[ww]
      ENDIF ELSE BEGIN
          GOTO, JUMP1
      ENDELSE 
    
      cor1 = psfsize[gal]/objsize[gal]*( 4./psfrho4[gal]-1. )/ $
                                       ( 4./pstruct[gal].rho4[clr]-1. )
      cor2 = psfsize[gal]/objsize[gal]
      cor3 = testpsfsize[gal]/testobjsize[gal]

      ;; We find galaxies from their smear polarizeability
      ww = 0
      ww=where(cor1 le rcut AND cor2 LE rcut AND cor3 LE rcut, nww)
      IF nww NE 0 THEN BEGIN 
          gal = gal[ww]
          cor1 = cor1[ww]
          cor2 = cor2[ww]
          cor3 = cor3[ww]
          ww = 0
          ww=where(psfsize[gal] GT 0., nww)
          IF nww NE 0 THEN BEGIN 
              gal = gal[ww]        
              cor1 = cor1[ww]
              cor2 = cor2[ww]
              cor3 = cor3[ww]
              psfe1 = (  a[ind1] + b[ind1]*x[gal] + c[ind1]*y[gal] $
                       + d[ind1]*x[gal]*y[gal] + e[ind1]*x[gal]*x[gal] $
                       + f[ind1]*y[gal]*y[gal] )/psfsize[gal]

              psfe2 = 2.*(  a[ind2] + b[ind2]*x[gal] + c[ind2]*y[gal] $
                          + d[ind2]*x[gal]*y[gal] + e[ind2]*x[gal]*x[gal] $
                          + f[ind2]*y[gal]*y[gal] )/psfsize[gal]
        
              IF nww GT 1 THEN BEGIN 
                  ramean = mean( pstruct[gal].ra)
                  decmean = mean( pstruct[gal].dec)

                  testpsfmean[ic] = mean(testpsfsize[gal])

                  tt            = moment(psfsize[gal])
                  psfmean[ic]   = tt[0]
                  psfsig[ic]    = sqrt(tt[1])
                  tt            = moment(pstruct[gal].reddening[clr])
                  ebv[ic]       = tt[0]
                  tt            = moment(psfe2)
                  psfe2mean[ic] = tt[0]
                  psfe2sig[ic]  = sqrt(tt[1])
                  tt            = moment(psfe1)
                  psfe1mean[ic] = tt[0]
                  psfe1sig[ic]  = sqrt(tt[1])
                  foredens[ic]  = nww
              ENDIF ELSE BEGIN 
                  ramean        = pstruct[gal[0]].ra
                  decmean       = pstruct[gal[0]].dec

                  testpsfmean[ic] = testpsfsize[0]
                  psfmean[ic]   = psfsize[0]
                  psfsig[ic]    = 0.
                  ebv[ic]       = pstruct[gal[0]].reddening[clr]
                  psfe2mean[ic] = psfe2[0]
                  psfe2sig[ic]  = 0.
                  psfe1mean[ic] = psfe1[0]
                  psfe1sig[ic]  = 0.
                  foredens[ic]  = 1
              ENDELSE 

              glactc, ramean/15., decmean, 2000., gll , gbb, 1
              gl[ic]=gll
              gb[ic]=gbb

              IF (ic MOD 20 EQ 0) OR (ic EQ 0) OR (ic EQ nframes-1) THEN BEGIN
                  print,'Field: ',fstr,' Objects: ',ntostr(nps)
                  print,'Number = ',ntostr(nww),sqrt(psfmean[ic]/2.)*2.35*0.4,sqrt(testpsfmean[ic])*2.35*0.4
              ENDIF 

              ;; Original e1, e2 measurements
              obje1 = (pstruct[gal].ixx[clr] $
                                 -pstruct[gal].iyy[clr])/objsize[gal]
              obje2 = 2.*pstruct[gal].ixy[clr]/objsize[gal]

              
              IF (psfe1sig[ic] le ermscut and psfe2sig[ic] le ermscut) then begin
                  IF (galflag EQ 0) THEN BEGIN 
                      galflag   = 1
                      galpsfe1  = psfe1
                      galpsfe2  = psfe2
                      galobjerr = pstruct[gal].momerr[clr]
                      galobje1o = pstruct[gal].e1[clr]
                      galobje2o = pstruct[gal].e2[clr]
                      galobje1test2 = obje1 - cor2*psfe1
                      galobje2test2 = obje2 - cor2*psfe2
                      galobje1test3 = obje1 - cor3*psfe1
                      galobje2test3 = obje2 - cor3*psfe2
                      galobje1  = obje1
                      galobje2  = obje2
                      galmag    = zmag[gal]
                      galcor    = cor1
                  ENDIF ELSE BEGIN 
                      galpsfe1  = [ temporary(galpsfe1),  psfe1 ]
                      galpsfe2  = [ temporary(galpsfe2),  psfe2 ]
                      galobjerr = [ temporary(galobjerr), $
                                    pstruct[gal].momerr[clr]]
                      galobje1o = [ temporary(galobje1o), $
                                    pstruct[gal].e1[clr] ]
                      galobje2o = [ temporary(galobje2o), $
                                    pstruct[gal].e2[clr] ]
                      galobje1test2 = [ temporary(galobje1test2), $
                                       obje1 - cor2*psfe1 ]
                      galobje2test2 = [ temporary(galobje2test2), $
                                       obje2 - cor2*psfe2 ]
                      galobje1test3 = [ temporary(galobje1test3), $
                                       obje1 - cor3*psfe1 ]
                      galobje2test3 = [ temporary(galobje2test3), $
                                       obje2 - cor3*psfe2 ]
                      galobje1  = [ temporary(galobje1),  obje1 ]
                      galobje2  = [ temporary(galobje2),  obje2 ]
                      galmag    = [ temporary(galmag), zmag[gal] ]
                      galcor    = [ temporary(galcor), cor1]
                  ENDELSE 
              ENDIF 
          ENDIF 
      ENDIF                     ;End of galaxy stuff
                                ;Begin star stuff

      JUMP1:
      IF nstars EQ 0 THEN GOTO,JUMP2
                                ;Free up some memory
      cor1 = 0
      obje1 = 0
      obje2 = 0
      psfe1 = 0
      psfe2 = 0
      zmag = 0
      ww = 0

      ww = where(objsize[stars] GT 0, nww)
      IF nww GT 1 THEN BEGIN
          stars = stars[ww]
      ENDIF ELSE BEGIN
          GOTO,JUMP2
      ENDELSE 

      cor1 = psfsize[stars]/objsize[stars]*( 4./psfrho4[stars]-1. )/ $
                                       ( 4./pstruct[stars].rho4[clr]-1. )
      
      psfe1 = (  a[ind1] + b[ind1]*x[stars] + c[ind1]*y[stars] $
                 + d[ind1]*x[stars]*y[stars] + e[ind1]*x[stars]*x[stars] $
                 + f[ind1]*y[stars]*y[stars] )/psfsize[stars]
      
      psfe2 = 2.*(  a[ind2] + b[ind2]*x[stars] + c[ind2]*y[stars] $
                    + d[ind2]*x[stars]*y[stars] + e[ind2]*x[stars]*x[stars] $
                    + f[ind2]*y[stars]*y[stars] )/psfsize[stars]

      obje1 = (pstruct[stars].ixx[clr] $
               -pstruct[stars].iyy[clr])/objsize[stars]
      obje2 = 2.*pstruct[stars].ixy[clr]/objsize[stars]
   
      if(starflag eq 0) then begin
          starflag = 1
          tpsfe1   = psfe1
          tpsfe2   = psfe2
          tobjerr  = pstruct[stars].momerr[clr]
          tobje1o  = pstruct[stars].e1[clr]
          tobje2o  = pstruct[stars].e2[clr]
          tobje1   = obje1
          tobje2   = obje2
          mag      = pstruct[stars].petrocounts[clr]
          tcor     = cor1
          size     = sqrt(objsize[stars]/2.)*2.35*0.40
      endif else begin
          tpsfe1  = [ temporary(tpsfe1),  psfe1 ]
          tpsfe2  = [ temporary(tpsfe2),  psfe2 ]
          tobjerr = [ temporary(tobjerr), pstruct[stars].momerr[clr] ]
          tobje1o = [ temporary(tobje1o), pstruct[stars].e1[clr] ]
          tobje2o = [ temporary(tobje2o), pstruct[stars].e2[clr] ]
          tobje1  = [ temporary(tobje1),  obje1 ]
          tobje2  = [ temporary(tobje2),  obje2 ]
          mag     = [ temporary(mag),     pstruct[stars].petrocounts[clr] ]
          tcor    = [ temporary(tcor),    cor1 ]
          size    = [ temporary(size),    sqrt(objsize[stars]/2.)*2.35*0.40 ]
      endelse
    
      JUMP2:

      pstruct = 0               ;Free up memory
      x = 0
      y = 0
      psfsize = 0
      testpsfsize = 0
      psfrho4 = 0
      objsize = 0
      obje1 = 0
      obje2 = 0
      psfe1 = 0
      psfe2 = 0
      cor1 = 0
      ww = 0

  ENDFOR 
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Almost no modification beyond this point E.S.S.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  set_plot,'X'
  if keyword_set(printfile) then begin
;  BEGPLOT,NAME=PRINTFILE
    set_plot,'ps'
    device,file=printfile,ysize=20,yoffset=3
  endif
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Plots of galaxy stuff
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                ;Extinction vs field
  !p.multi=0
  plot, ff, ebv, xtitle='Field', ytitle='Extinction (mag)',$
    title=psffile1
  psfmean=sqrt(psfmean/2.)*2.35*0.4
  psfsig=sqrt(psfsig/2.)*2.35*0.4
  !p.multi=[0,1,3]
                                ;Mean psf vs field
  plot,ff, psfmean, xtitle='Field', ytitle='FWHM (arcsec)',$
    title=psffile1
                                ;Mean psf_e1, psf_e2 vs field 
  plot, ff, psfe1mean, xtitle='Field', ytitle='e_1'
  plot, ff, psfe2mean, xtitle='Field',ytitle='e_2'


                                ;psf_e1, psf_e2 variance vs field
  !p.multi=[0,1,2]
  plot,ff,psfe1sig,xtitle='Field',ytitle='RMS(e_1)',title=psffile1
  xx=[ff[0],ff[nframes-1]]
  yy=[ermscut,ermscut]
  oplot,xx,yy
  
  plot,ff, psfe2sig, xtitle='Field',ytitle='RMS(e_2)'
  oplot,xx,yy
                                ; Number of gals vs galactic lattitude
                                ; shows contamination
  !p.multi=0
  plot, gb, foredens, xtitle='Latitude',ytitle='Number'
                                ; Fit to the number above
  fsig=sqrt(foredens)
  fitlin, gb, foredens, fsig, aa, siga, bb, sigb
  print, aa, siga, bb, sigb

                                ;Plots of residual gal shape vs. the 
                                ;psf shape near them after correction
  yrange=[-0.1, 0.1]

  !p.multi=[0,2,2]
  xmin = -0.4
  xmax = 0.4
  e1err = moment(galobje1o)
  print, 'Mean e1', e1err[0], sqrt(e1err[1])
  e2err = moment(galobje2o)
  print, 'Mean e2', e2err[0], sqrt(e2err[1])
  erre1 = sqrt(galobjerr^2+0.32^2)
  xx=[-1,1]
  binner2, galpsfe1, galobje1o, erre1, xo, yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e1 '+colors[clr],[-1,1],yrange
  print,'e1PSF vs e1'
;  forprint,xo,yo,sig

  binner2,galpsfe1,galobje2o,erre1, xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e2 '+colors[clr],[-1,1],yrange

  binner2,galpsfe2,galobje1o,erre1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e1 '+colors[clr],[-1,1],yrange

  binner2,galpsfe2,galobje2o,erre1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e2 '+colors[clr],[-1,1],yrange

                                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                ;; Same but with only size used for correction
                                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  !p.multi=[0,2,2]
  xmin = -0.4
  xmax = 0.4
  e1err = moment(galobje1test2)
  print, 'Mean e1', e1err[0], sqrt(e1err[1])
  e2err = moment(galobje2test2)
  print, 'Mean e2', e2err[0], sqrt(e2err[1])
  xx=[-1,1]
  binner2, galpsfe1, galobje1test2, erre1, xo, yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e1 test'+colors[clr],[-1,1],yrange
  print,'e1PSF vs e1 test'
;  forprint,xo,yo,sig

  binner2,galpsfe1,galobje2test2,erre1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e2 test'+colors[clr],[-1,1],yrange

  binner2,galpsfe2,galobje1test2,erre1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e1 test'+colors[clr],[-1,1],yrange

  binner2,galpsfe2,galobje2test2,erre1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e2 test'+colors[clr],[-1,1],yrange

                                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                ;; Same but with only size used for correction
                                ;; and size is an alternate size
                                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  !p.multi=[0,2,2]
  xmin = -0.4
  xmax = 0.4
  e1err = moment(galobje1test3)
  print, 'Mean e1', e1err[0], sqrt(e1err[1])
  e2err = moment(galobje2test3)
  print, 'Mean e2', e2err[0], sqrt(e2err[1])
  xx=[-1,1]
  binner2, galpsfe1, galobje1test3, erre1,xo, yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e1 test'+colors[clr],[-1,1],yrange
  print,'e1PSF vs e1 test'
;  forprint,xo,yo,sig

  binner2,galpsfe1,galobje2test3,erre1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e2 test'+colors[clr],[-1,1],yrange

  binner2,galpsfe2,galobje1test3,erre1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e1 test'+colors[clr],[-1,1],yrange

  binner2,galpsfe2,galobje2test3,erre1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e2 test'+colors[clr],[-1,1],yrange

                                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                ;; same as above, except before correction
                                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  !p.multi=[0,2,2]
  xmin=-0.4
  xmax=0.4
  e1err=moment(galobje1)
  print,'Mean e1',e1err[0],sqrt(e1err[1])
  e2err=moment(galobje2)
  print,'Mean e2',e2err[0],sqrt(e2err[1])
  xx=[-1,1]
  tt=moment(galobje1)
  print,tt[0],tt[1]
  binner2,galpsfe1,galobje1,erre1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e1 '+colors[clr],[-1,1],yrange

  binner2,galpsfe1,galobje2,erre1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e2 '+colors[clr],[-1,1],yrange

  binner2,galpsfe2,galobje1,erre1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e1 '+colors[clr],[-1,1],yrange

  binner2,galpsfe2,galobje2,erre1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e2 '+colors[clr],[-1,1],yrange

                                ; magnitude histogram
  !p.multi=0
  histo=histogram(galmag,binsize=0.5,min=15.,max=25.)
  bin=fltarr(1)
  bin[0]=15.25
  for i=15.75,24.75,0.5 do begin
    bin=[bin,i]
  endfor
  plot,bin,histo,xtitle=colors[clr],ytitle='N'
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Star plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                ;Size-mag diagram for stars
  plot,mag,size,psym=3,xtitle=colors[clr],ytitle="FWHM (arcsec)"
  !p.multi=[0,2,2]

                                ;Residual star shape vs. psf shape
                                ;after correction
  xmin = -0.4
  xmax = 0.4
  ymin = -0.1
  ymax = 0.1
  e1err = moment(tobje1o)
  print,'Mean star e1',e1err[0],sqrt(e1err[1])
  e2err = moment(tobje2o)
  print,'Mean star e2',e2err[0],sqrt(e2err[1])
  erre1 = sqrt(tobjerr^2+0.32^2)
  xx=[-1,1]
  binner,tpsfe1,tobje1o,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e1 '+colors[clr],[-1,1],[ymin,ymax]

  binner,tpsfe1,tobje2o,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e2 '+colors[clr],[-1,1],[ymin,ymax]

  binner,tpsfe2,tobje1o,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e1 '+colors[clr],[-1,1],[ymin,ymax]

  binner,tpsfe2,tobje2o,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e2 '+colors[clr],[-1,1],[ymin,ymax]


  if keyword_set(printfile) then begin
    device,/close
  endif
  set_plot,'X'

  
  ptime, systime(1)-time

end

