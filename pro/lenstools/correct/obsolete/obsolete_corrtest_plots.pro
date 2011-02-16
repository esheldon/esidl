PRO corrtest_plots_getrange, x1, y1, x2, y2, nperbin, xrange, yrange

  binner_bynum, x1, y1, nperbin, txo1, tyo1, tsig1, tnum1
  binner_bynum, x2, y2, nperbin, txo2, tyo2, tsig2, tnum2

  w1=where(tnum1 EQ nperbin)
  w2=where(tnum2 EQ nperbin)
  xrange=prange(txo1[w1],txo2[w2],/noerror,slack=0.3,/symmetric)
  yrange=prange(tyo1[w1],tyo2[w2],/noerror,slack=0.4,/symmetric)

END 

pro corrtest_plots, run, rerun, camcol, clr, start=start, nframes=nframes,$
                    minmag=minmag, maxmag=maxmag, $
                    rcut=rcut, ermscut=ermscut, errmax=errmax, psfile=psfile, $
                    maxseeing=maxseeing, offset=offset, onefile=onefile, status=status,$
                    overwrite=overwrite

  ;; status is 1 unless we reach end
  status=1

;   ermscut = cut off frames with psfe rms gt ermscut
  
  IF n_params() EQ 0 THEN BEGIN 
      print,'-syntax corshape, run, rerun, camcol, clr, start=start, nframes=nframes, minmag=minmag, maxmag=maxmag, rcut=rcut, ermscut=ermscut, psfile=psfile, maxseeing=maxseeing, offset=offset, onefile=onefile, status=status,overwrite=overwrite'
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
  nchar=25
;  newfront = 'adat'
;  nchar=24
  
  selectclr = 2                 ;r-band

  colors=['u','g','r','i','z']
  bands = [1,2,3]               ; we will process g,r,i
  nband = n_elements(bands)
  
  rstr = run2string(run)
  rrstr = ntostr(rerun)

  ;; cuts used in make_corrected_files for stars
  stmaxmag = [0.0,21.,20.,20.,0.0]
  IF n_elements(maxmag) EQ 0 THEN maxmag=22.
  IF n_elements(minmag) EQ 0 THEN minmag=0.
  IF n_elements(rcut) EQ 0 THEN rcut=0.8
  IF n_elements(ermscut) EQ 0 THEN ermscut = 1000.
  IF n_elements(errmax) EQ 0 THEN errmax = 1000.
  IF n_elements(maxseeing) EQ 0 THEN maxseeing=3.0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Setup the column
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  cstr = ntostr(camcol)
  
  fetch_dir, run, camcol, rerun, dir, atldir, corrdir=corrdir,$
    corratldir = fitdir
  fetch_file_list, corrdir, adatcfiles, fnums, start=start, nframes=nframes, $
    fieldmin=fieldmin, fieldmax=fieldmax, front=newfront, nchar=nchar

  nfields = n_elements(adatcfiles)

  IF keyword_set(offset) THEN addstr='offset' ELSE addstr=''
  IF n_elements(psfile) EQ 0 THEN BEGIN
      psfile = fitdir+addstr+'corshape_'+rstr+'_'+cstr+'_'+colors[clr]+'_N1.ps'
      
      IF NOT keyword_set(overwrite) THEN BEGIN 
          WHILE exist(psfile) DO psfile=newname(psfile)
      ENDIF 
  ENDIF 
  fitsfile = repstr(psfile, '.ps','.fit')

  print
  print,'psfile: ',psfile
  print,'fitsfile: ',fitsfile
  print

  print
  print,'corrected files are '+newfront+'-*'
  print
  ;rename_tsobj, files, corrdir, newfront, adatcfiles, nchar
  
  galflag   = 0
  starflag  = 0
  ebv       = fltarr(nfields)
  ff        = fltarr(nfields)
  psfmean   = fltarr(nfields)
  psfsig    = fltarr(nfields)
  psfe1mean = fltarr(nfields)
  psfe1sig  = fltarr(nfields)
  psfe2mean = fltarr(nfields)
  psfe2sig  = fltarr(nfields)
  foredens  = fltarr(nfields)
  gl        = fltarr(nfields)
  gb        = fltarr(nfields)
  
  IF keyword_set(onefile) THEN BEGIN 
      psfname, run, camcol, fname
      print
      print,'Reading psf fit file ',fname
      psfstr = mrdfits(fitdir+fname,1)
      print
      error=0
  ENDIF 

  print,'-----------------------------------------------------'
  print,' Run: ',rstr,' Rerun: ',rrstr,' Camcol: ',cstr,' Band: ',colors[clr]
  print,'-----------------------------------------------------'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over fields
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR ic = 0L, nfields-1 DO BEGIN 
      infile = adatcfiles[ic]
      field = fnums[ic]
      fstr = ntostr(field)
      ff[ic]=field

      ;; Fields have already been trimmed of overlap.  Don't need to
      ;; use read_tsobjmin
      read_tsobj, corrdir, pstruct, start=field, nframes=1,$
                  /everytag, tsobjstr=tsobjstr, verbose=0, front='adatc'
;      openr, lun, infile, /get_lun, ERROR=error1
;      IF error1 EQ 0 THEN BEGIN
      IF n_elements(pstruct) NE 0 THEN BEGIN 

;          IF ic EQ 0 THEN BEGIN 
;              pstruct = mrdfits3(lun, 1, 0, /silent)
;          ENDIF ELSE BEGIN
;              pstruct = mrdfits3(lun, 1, 0, /silent, /deja_vu)
;          ENDELSE 
;          free_lun, lun
          
          ;; we want this for everything, will use for stars later
          objsize = pstruct.ixx[clr]+pstruct.iyy[clr]
          psfsize = pstruct.psfixx[clr] + pstruct.psfiyy[clr]
          zmag = pstruct.petrocounts[selectclr]-pstruct.reddening[selectclr]
                                ;Get good objects
          mobj=where( (objsize GT 0)                  AND $
                      (psfsize GT 0.)                 AND $
                      (zmag LE maxmag)                AND $ 
                      (zmag GE minmag)                AND $
                      (pstruct.e1[clr] NE 1.e10)      AND $
                      (pstruct.momerr[clr] LE errmax) AND $
                      (pstruct.r[clr] GT 0.0 )        AND $
                      (abs(pstruct.ixy[clr]) LT 2. )  AND $
                      (pstruct.seeing[clr] LT maxseeing), nobj)
          
          IF nobj NE 0 THEN BEGIN 

              ;; uncomment this to get all objects with r > rcut
;              stars = where( pstruct[mobj].r[clr] GT rcut, nstars)

              stars = where( ( (pstruct[mobj].starflag AND 2b^clr) NE 0) AND $
                             ( pstruct[mobj].r[clr] GT rcut ) AND $
                             ( (pstruct[mobj].type[2] EQ 6) AND $
                               (pstruct[mobj].type[1] EQ 6) ), nstars)
              

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; galaxies
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              
              gal = where( (pstruct[mobj].r[clr] LT rcut), ngal)
              
              IF ngal NE 0 THEN BEGIN 
              
                  gal = mobj[gal]
                  ngal = n_elements(gal)
                  cor1 = psfsize[gal]/objsize[gal]*(4./pstruct[gal].psfrho4[clr]-1.)/(4./pstruct[gal].rho4[clr]-1)
                  cor2 = psfsize[gal]/objsize[gal]
                  
                  psfe1 = (pstruct[gal].psfixx[clr] - $
                           pstruct[gal].psfiyy[clr])/psfsize[gal]
                  psfe2 = 2.*pstruct[gal].psfixy[clr]/psfsize[gal]
                  
                  IF ngal GT 1 THEN BEGIN 
                      ramean = mean( pstruct[gal].ra)
                      decmean = mean( pstruct[gal].dec)
                      
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
                      foredens[ic]  = ngal
                  ENDIF ELSE BEGIN 
                      ramean        = pstruct[gal[0]].ra
                      decmean       = pstruct[gal[0]].dec
                      
                      psfmean[ic]   = psfsize[gal[0]]
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
                  
                  IF (ic MOD 20 EQ 0) OR (ic EQ 0) OR (ic EQ nfields-1) THEN BEGIN
                      print,'Field: ',fstr,' Objects: ',ntostr(nobj)
                      print,'Number = ',ntostr(ngal),sqrt(psfmean[ic]/2.)*2.35*0.4
                  ENDIF 
                  
                  ;; Original e1, e2 measurements
                  obje1 = (pstruct[gal].ixx[clr] $
                           -pstruct[gal].iyy[clr])/objsize[gal]
                  obje2 = 2.*pstruct[gal].ixy[clr]/objsize[gal]
                  obje1_2 = obje1 - cor2*psfe1
                  obje2_2 = obje2 - cor2*psfe2
                  
                  IF (psfe1sig[ic] le ermscut and psfe2sig[ic] le ermscut) then BEGIN
                      add_arrval, psfe1, galpsfe1
                      add_arrval, psfe2, galpsfe2
                      add_arrval, pstruct[gal].momerr[clr], galobjerr
                      add_arrval, pstruct[gal].e1[clr], galobje1o
                      add_arrval, pstruct[gal].e2[clr], galobje2o
                      add_arrval, obje1, galobje1
                      add_arrval, obje2, galobje2
                      add_arrval, obje1_2, galobje1o2
                      add_arrval, obje2_2, galobje2o2
                      add_arrval, pstruct[gal].petrocounts[clr]-$
                                  pstruct[gal].reddening[clr], galmag 
                      add_arrval, cor1, galcor
                  ENDIF 
              ENDIF ELSE BEGIN ;; ngal ne 0
                  print,'-- No Galaxies'
              ENDELSE 

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; begin stars
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          
              IF nstars NE 0 THEN BEGIN 
                  stars = mobj[stars]
                                ;Free up some memory
                  setzero, cor1, obje1, obje2, psfe1, psfe2, zmag, ww, tmpmag, $
                           galpsf, gal
                  
                  
                  psfe1 = (pstruct[stars].psfixx[clr] - pstruct[stars].psfiyy[clr])/psfsize[stars]
                  psfe2 = 2.*pstruct[stars].psfixy[clr]/psfsize[stars]
                  
                  cor1 = pstruct[stars].R[clr]
                  
                  obje1 = (pstruct[stars].ixx[clr] - pstruct[stars].iyy[clr])/objsize[stars]
                  obje2 = 2.*pstruct[stars].ixy[clr]/objsize[stars]
                  
                  ;; try correcting with corr=1.0!
                  obje1_cor2 = obje1 - psfe1
                  obje2_cor2 = obje2 - psfe2
                  
                  add_arrval, psfe1, tpsfe1
                  add_arrval, psfe2, tpsfe2
                  add_arrval, pstruct[stars].momerr[clr], tobjerr
                  add_arrval, pstruct[stars].e1[clr], tobje1o
                  add_arrval, pstruct[stars].e2[clr], tobje2o
                  add_arrval, obje1_cor2, tobje1o_cor2
                  add_arrval, obje2_cor2, tobje2o_cor2
                  add_arrval, obje1, tobje1
                  add_arrval, obje2, tobje2
                  add_arrval, pstruct[stars].petrocounts[clr], mag
                  add_arrval, cor1, tcor
                  add_arrval, cor1-pstruct[stars].r[clr], diff
                  add_arrval, sqrt(objsize[stars]/2.)*2.35*0.40, size
                  
              ENDIF ELSE BEGIN ;; nstar ne 0
                  print,'-- No Stars'
              ENDELSE 
         
          ENDIF ELSE BEGIN ;; nobj ne 0
              print,'Field: ',fstr,'  No objects passed mom cuts and seeing cut of ',maxseeing
          ENDELSE  
      ENDIF ELSE BEGIN ;; field not empty
          ;print,!ERR_STRING
          ;free_lun, lun
      ENDELSE 
      delvarx, pstruct
      setzero, psfsize, objsize, obje1, obje2, $
               psfe1, psfe2, cor1, ww, stpsf, stars, obje1_cor2, obje2_cor2
  ENDFOR 


  IF n_elements(psfile) NE 0 THEN begplot,name=psfile,/color
  !p.charsize=1.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Plots of galaxy stuff
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                ;Extinction vs field
  !p.multi=0

  myusersym, 'fill_circle'
  fpsym=8
  aplot, !gratio, ff, ebv, xtitle='Field', ytitle='Extinction (mag)',$
    title=psffile1,psym=fpsym,charsize=1.75
  psfmean=sqrt(psfmean/2.)*2.35*0.4
  psfsig=sqrt(psfsig/2.)*2.35*0.4
  !p.multi=[0,1,3]
                                ;Mean psf vs field
  tsize=2.0
  plot,ff, psfmean, xtitle='Field', ytitle='FWHM (arcsec)',$
    title=psffile1,charsize=tsize,psym=fpsym
                                ;Mean psf_e1, psf_e2 vs field 
  plot, ff, psfe1mean, xtitle='Field', ytitle='e!D1!N',charsize=tsize,psym=fpsym
  plot, ff, psfe2mean, xtitle='Field',ytitle='e!D2!N',charsize=tsize,psym=fpsym


                                ;psf_e1, psf_e2 variance vs field
  !p.multi=[0,1,2]
  tsize=1.75
  plot,ff,psfe1sig,xtitle='Field',ytitle='RMS(e!D1!N)',title=psffile1,psym=fpsym,$
       charsize=tsize
  xx=[ff[0],ff[nfields-1]]
  yy=[ermscut,ermscut]
  oplot,xx,yy
  
  plot,ff, psfe2sig, xtitle='Field',ytitle='RMS(e!D2!N)',psym=fpsym,charsize=tsize
  oplot,xx,yy
                                ; Number of gals vs galactic lattitude
                                ; shows contamination
  !p.multi=0
  aplot, !gratio,gb, foredens, xtitle='Latitude',ytitle='Number',psym=fpsym,$
         charsize=tsize
                                ; Fit to the number above
  fsig=sqrt(foredens)

  ;;;;;;;;;;;;;
  ;; binsizes
  ;;;;;;;;;;;;;

  ebin = 0.015
  erbin = ebin/2.
  rbin = 0.05

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plots of egal vs epsf
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,"-------------------------------------------"

  !p.multi=[0,2,2]

  e1err = moment(galobje1o)
  print, 'Mean corrected gal e1', e1err[0], sqrt(e1err[1])
  e2err = moment(galobje2o)
  print, 'Mean corrected gal e2', e2err[0], sqrt(e2err[1])
  erre1 = sqrt(galobjerr^2+0.32^2)
  erre2 = sqrt(galobjerr^2+0.32^2)
  xx=[-1,1]

  ;;; about 8 bins
  frac = 1./8.
  nperbin = long( round(frac*n_elements(galpsfe1)) )

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find out xrange for gals
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  corrtest_plots_getrange, galpsfe1, galobje1, galpsfe2, galobje2, nperbin,$
                           xrange, yrange

  erase & multiplot, [2,2], /square

  plot_fitlin, galpsfe1, galobje1o, nperbin, $
               intercept_ge1pe1, intercept_ge1pe1_error, $
               slope_ge1pe1, slope_ge1pe1_error,$
               ytitle="e!D1!N Gal",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band corrected";,/plotnum

  cyrange=!y.crange 
  plothist,galpsfe1,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe2, galobje1o, nperbin, $
               intercept_ge1pe2, intercept_ge1pe2_error, $
               slope_ge1pe2, slope_ge1pe2_error,$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Gal"

  cyrange=!y.crange 
  plothist,galpsfe2,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe1, galobje2o, nperbin, $
               intercept_ge2pe1, intercept_ge2pe1_error, $
               slope_ge2pe1, slope_ge2pe1_error,$
               ytitle="e!D2!N Gal",xtitle="e!D1!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, galpsfe2, galobje2o, nperbin, $
               intercept_ge2pe2, intercept_ge2pe2_error, $
               slope_ge2pe2, slope_ge2pe2_error,$
               xtitle="e!D2!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D2!N Gal"

  multiplot,/reset

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Same but with only size used for correction
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,"-------------------------------------------"
  !p.multi=[0,2,2]
  e1err = moment(galobje1o2)
  print, 'Mean corrected gal (size only) e1', e1err[0], sqrt(e1err[1])
  e2err = moment(galobje2o2)
  print, 'Mean corrected gal (size only) e2', e2err[0], sqrt(e2err[1])
  erre1 = sqrt(galobjerr^2+0.32^2)
  erre2 = sqrt(galobjerr^2+0.32^2)
  xx=[-1,1]

  erase & multiplot, [2,2], /square

  plot_fitlin, galpsfe1, galobje1o2, nperbin, $
               ytitle="e!D1!N Gal",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band corrected    Size Only";,/plotnum

  cyrange=!y.crange 
  plothist,galpsfe1,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe2, galobje1o2, nperbin, $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Gal"

  cyrange=!y.crange 
  plothist,galpsfe2,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe1, galobje2o2, nperbin, $
               ytitle="e!D2!N Gal",xtitle="e!D1!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, galpsfe2, galobje2o2, nperbin, $
               xtitle="e!D2!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D2!N Gal"

  multiplot,/reset

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; same as above, except before correction
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,"-------------------------------------------"
  !p.multi=[0,2,2]
  e1err=moment(galobje1)
  print,'Mean uncorrected gal e1',e1err[0],sqrt(e1err[1])
  e2err=moment(galobje2)
  print,'Mean uncorrected gal e2',e2err[0],sqrt(e2err[1])
  erre1=sqrt(galobjerr^2+0.32^2)
  erre2=sqrt(galobjerr^2+0.32^2)
  xx=[-1,1]


  erase & multiplot, [2,2], /square

  plot_fitlin, galpsfe1, galobje1, nperbin, $
               ytitle="e!D1!N Gal",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Uncorrected";,/plotnum

  cyrange=!y.crange 
  plothist,galpsfe1,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe2, galobje1, nperbin, $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Gal"

  cyrange=!y.crange 
  plothist,galpsfe2,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe1, galobje2, nperbin, $
               ytitle="e!D2!N Gal",xtitle="e!D1!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, galpsfe2, galobje2, nperbin, $
               xtitle="e!D2!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D2!N Gal"

  multiplot,/reset

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; magnitude histogram
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  !p.multi=0
  histo=histogram(galmag,binsize=0.5,min=15.,max=25.)
  bin=fltarr(1)
  bin[0]=15.25
  for i=15.75,24.75,0.5 do begin
    bin=[bin,i]
  endfor
  aplot,!gratio,bin,histo,xtitle=colors[clr],ytitle='N',charsize=tsize
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; e vs smear polarizability
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  !p.multi=[0,0,2]
  
  pold=!p.charsize
  !p.charsize=1.25
  binner, galcor, galobje1o, rbin, xo, yo, sig
  plot1,xo,yo,sig,intercept_ge1R,slope_ge1R,$
    'R smear','e1',[0,0.8],[-0.05,0.05],!colors[clr]+"-band corrected",$
    siga=intercept_ge1R_error,sigb=slope_ge1R_error

  cyrange=!y.crange 
  plothist,galcor,xhist,yhist,bin=rbin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  binner, galcor, galobje2o, rbin, xo, yo, sig
  plot1,xo,yo,sig,intercept_ge2R,slope_ge2R,$
    'R smear','e2',[0,0.8],[-0.05,0.05],$
    siga=intercept_ge2R_error,sigb=slope_ge2R_error

  !p.charsize=pold
  !p.multi=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; e vs smear polarizability*epsf
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  erase & multiplot, [2,2], /square

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find out xrange for gals including galcor
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  corrtest_plots_getrange, galpsfe1*galcor, galobje1, galpsfe2*galcor, $
                           galobje2, nperbin,$
                           xrange, tyrange ;tyrange not used

  plot_fitlin, galpsfe1*galcor, galobje1o, nperbin, $
               intercept_ge1pe1R, intercept_ge1pe1R_error, $
               slope_ge1pe1R, slope_ge1pe1R_error,$
               ytitle="e!D1!N Gal",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band corrected";,/plotnum

  cyrange=!y.crange 
  plothist,galpsfe1*galcor,xhist,yhist,bin=erbin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe2*galcor, galobje1o, nperbin, $
               intercept_ge1pe2R, intercept_ge1pe2R_error, $
               slope_ge1pe2R, slope_ge1pe2R_error,$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Gal"

  cyrange=!y.crange 
  plothist,galpsfe2*galcor,xhist,yhist,bin=erbin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe1*galcor, galobje2o, nperbin, $
               intercept_ge2pe1R, intercept_ge2pe1R_error, $
               slope_ge2pe1R, slope_ge2pe1R_error,$
               ytitle="e!D2!N Gal",xtitle="e!D1!N PSF*R",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, galpsfe2*galcor, galobje2o, nperbin, $
               intercept_ge2pe2R, intercept_ge2pe2R_error, $
               slope_ge2pe2R, slope_ge2pe2R_error,$
               xtitle="e!D2!N PSF*R",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D2!N Gal"

  multiplot,/reset

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Star plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  yrange=[-0.2,0.2]
  pold=!p.charsize
  !p.charsize=1.75
                                ;Size-mag diagram for stars
  ;;plot,mag,size,psym=3,xtitle=colors[clr],ytitle="FWHM (arcsec)"
  ploth,mag,size,xtitle=colors[clr],ytitle="FWHM (arcsec)",$
    nxbins=50,nybins=50,/silent
                                ;smear polarizability for all objects
  yrr=[0.,2.0]
  xrr=[12.,maxmag]
  ;;plot,galmag,galcor,psym=3,xtitle=colors[clr],ytitle="Rsm",yrange=yrr,xrange=xrr
  ;;oplot,mag,tcor,psym=4,symsize=.5
  ploth,[galmag,mag],[galcor,tcor],$
    xtitle=colors[clr],ytitle="Rsm",yrange=yrr,xrange=xrr,$
    nxbins=100,nybins=100,/silent,/sqrt
  oplot, [0.0, 1000.], [rcut, rcut], color=!grey50
  ;oplot, [stmaxmag[clr], stmaxmag[clr]], [rcut, 1000.], color=!grey50
  !p.charsize=pold

  sigma_clip, tcor, meantcor, sigmatcor, nsig=3.5, niter=3, /silent
  momtmp = moment(tcor)
  meantcor2 = momtmp[0]
  sigmatcor2 = sqrt(momtmp[1])
  print
  print,"-------------------------------------------"
  print,'Mean star smear polarizability: ',meantcor, sigmatcor
  print,'straight mean: ',meantcor2,sigmatcor2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;star shape vs. psf shape before correction
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,"-------------------------------------------"
  e1err = moment(tobje1)
  print,'Mean star e1',e1err[0],sqrt(e1err[1])
  e2err = moment(tobje2)
  print,'Mean star e2',e2err[0],sqrt(e2err[1])
  erre1 = sqrt(tobjerr^2+0.32^2)
  erre2 = sqrt(tobjerr^2+0.32^2)
  xx=[-1,1]

  sfrac = 1./8.
  snperbin = long( round(frac*n_elements(tpsfe1)) )

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find out xrange for stars
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  corrtest_plots_getrange, tpsfe1, tobje1, tpsfe2, tobje2, nperbin,$
                           xrange, yrange

  erase & multiplot, [2,2], /square

  plot_fitlin, tpsfe1, tobje1, snperbin, $
               ytitle="e!D1!N Star",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Uncorrected";,/plotnum

  cyrange=!y.crange 
  plothist,tpsfe1,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, tpsfe2, tobje1, snperbin, $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Star"

  cyrange=!y.crange 
  plothist,tpsfe2,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, tpsfe1, tobje2, snperbin, $
               ytitle="e!D2!N Star",xtitle="e!D1!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, tpsfe2, tobje2, snperbin, $
               xtitle="e!D2!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D2!N Star"

  multiplot,/reset

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Residual star shape vs. psf shape after correction
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,"-------------------------------------------"
  e1err = moment(tobje1o)
  print,'Mean corrected star e1',e1err[0],sqrt(e1err[1])
  e2err = moment(tobje2o)
  print,'Mean corrected star e2',e2err[0],sqrt(e2err[1])
  erre1 = sqrt(tobjerr^2+0.32^2)
  erre2 = sqrt(tobjerr^2+0.32^2)
  xx=[-1,1]

  erase & multiplot, [2,2], /square

  plot_fitlin, tpsfe1, tobje1o, snperbin, $
               ytitle="e!D1!N Star",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Corrected";,/plotnum

  cyrange=!y.crange 
  plothist,tpsfe1,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, tpsfe2, tobje1o, snperbin, $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Star"

  cyrange=!y.crange 
  plothist,tpsfe2,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, tpsfe1, tobje2o, snperbin, $
               ytitle="e!D2!N Star",xtitle="e!D1!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, tpsfe2, tobje2o, snperbin, $
               xtitle="e!D2!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D2!N Star"

  multiplot,/reset

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Try correcting with cor=1.0
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,"-------------------------------------------"
  e1err = moment(tobje1o_cor2)
  print,'Mean corrected star e1 (R=1)',e1err[0],sqrt(e1err[1])
  e2err = moment(tobje2o_cor2)
  print,'Mean corrected star e2 (R=1)',e2err[0],sqrt(e2err[1])
  erre1 = sqrt(tobjerr^2+0.32^2)
  erre2 = sqrt(tobjerr^2+0.32^2)
  xx=[-1,1]

  erase & multiplot, [2,2], /square

  plot_fitlin, tpsfe1, tobje1o_cor2, snperbin, $
               ytitle="e!D1!N Star",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Corrected (R=1)";,/plotnum

  cyrange=!y.crange 
  plothist,tpsfe1,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, tpsfe2, tobje1o_cor2, snperbin, $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Star"

  cyrange=!y.crange 
  plothist,tpsfe2,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, tpsfe1, tobje2o_cor2, snperbin, $
               ytitle="e!D2!N Star",xtitle="e!D1!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, tpsfe2, tobje2o_cor2, snperbin, $
               xtitle="e!D2!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D2!N Star"

  multiplot,/reset

  IF n_elements(psfile) NE 0 THEN endplot

  print
  print,"-------------------------------------------"
;  if keyword_set(psfile) then begin
;    device,/close
;  endif
;  set_plot,'X'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; For output
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  os=create_struct('galpsfe1', galpsfe1,$
                   'galpsfe2', galpsfe2,$
                   'gale1', galobje1,$
                   'gale2', galobje2,$
                   'corr',galcor,$
                   'gale1_corrected', galobje1o,$
                   'gale2_corrected', galobje2o,$
                   $ ;; linear fits to galaxy ellip vs psf ellip
                   'intercept_ge1pe1', intercept_ge1pe1, $
                   'slope_ge1pe1',slope_ge1pe1,$
                   'intercept_ge1pe1_error', intercept_ge1pe1_error, $
                   'slope_ge1pe1_error',slope_ge1pe1_error,$
                   $
                   'intercept_ge2pe1', intercept_ge2pe1, $
                   'slope_ge2pe1',slope_ge2pe1,$
                   'intercept_ge2pe1_error', intercept_ge2pe1_error, $
                   'slope_ge2pe1_error',slope_ge2pe1_error,$
                   $
                   'intercept_ge1pe2', intercept_ge1pe2, $
                   'slope_ge1pe2',slope_ge1pe2,$
                   'intercept_ge1pe2_error', intercept_ge1pe2_error, $
                   'slope_ge1pe2_error',slope_ge1pe2_error,$
                   $
                   'intercept_ge2pe2', intercept_ge2pe2, $
                   'slope_ge2pe2',slope_ge2pe2,$
                   'intercept_ge2pe2_error', intercept_ge2pe2_error, $
                   'slope_ge2pe2_error',slope_ge2pe2_error,$
                   $ ;; linear fits to gal ellip vs smear polarizability
                   'intercept_ge1R',intercept_ge1R,$
                   'slope_ge1R',slope_ge1R,$
                   'intercept_ge1R_error',intercept_ge1R_error,$
                   'slope_ge1R_error',slope_ge1R_error,$
                   $
                   'intercept_ge2R',intercept_ge2R,$
                   'slope_ge2R',slope_ge2R,$
                   'intercept_ge2R_error',intercept_ge2R_error,$
                   'slope_ge2R_error',slope_ge2R_error)
  os2=create_struct($
                   $ ;; linear fits to gal ellip vs smear*psfellip
                   'intercept_ge1pe1R',intercept_ge1pe1R,$
                   'slope_ge1pe1R',slope_ge1pe1R,$
                   'intercept_ge1pe1R_error',intercept_ge1pe1R_error,$
                   'slope_ge1pe1R_error',slope_ge1pe1R_error,$
                   $
                   'intercept_ge2pe1R',intercept_ge2pe1R,$
                   'slope_ge2pe1R',slope_ge2pe1R,$
                   'intercept_ge2pe1R_error',intercept_ge2pe1R_error,$
                   'slope_ge2pe1R_error',slope_ge2pe1R_error,$
                   $
                   'intercept_ge1pe2R',intercept_ge1pe2R,$
                   'slope_ge1pe2R',slope_ge1pe2R,$
                   'intercept_ge1pe2R_error',intercept_ge1pe2R_error,$
                   'slope_ge1pe2R_error',slope_ge1pe2R_error,$
                   $
                   'intercept_ge2pe2R',intercept_ge2pe2R,$
                   'slope_ge2pe2R',slope_ge2pe2R,$
                   'intercept_ge2pe2R_error',intercept_ge2pe2R_error,$
                   'slope_ge2pe2R_error',slope_ge2pe2R_error )

  os=create_struct(os,os2)

  mwrfits2, os, fitsfile, /create, /destroy

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; try re-correcting
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  corshape_again, fitsfile

  !p.multi=0
  ptime, systime(1)-time

  status=0

end

