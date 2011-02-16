pro binner, xi,yi,xo,yo,sig,bin
  
  IF n_elements(bin) EQ 0 THEN bin=0.01
  histo=histogram(xi,bin=bin,reverse_indices=r)
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

pro corshape, run, rerun, camcol, clr, start=start, nframes=nframes,$
              minmag=minmag, maxmag=maxmag, $
              rcut=rcut, ermscut=ermscut, errmax=errmax, psfile=psfile, $
              maxseeing=maxseeing, offset=offset, onefile=onefile, status=status

  ;; status is 1 unless we reach end
  status=1

;   ermscut = cut off frames with psfe rms gt ermscut
  
  IF n_params() EQ 0 THEN BEGIN 
      print,'-syntax corshape, run, rerun, camcol, clr, start=start, nframes=nframes, minmag=minmag, maxmag=maxmag, rcut=rcut, ermscut=ermscut, psfile=psfile, maxseeing=maxseeing, offset=offset, onefile=onefile, status=status'
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
      WHILE exist(psfile) DO psfile=newname(psfile)
  ENDIF 
  fitsfile = repstr(psfile, '.ps','.fit')

  print
  print,'psfile: ',psfile
  print,'fitsfile; ',fitsfile
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

      ;; psf fit name
      psfname2, run, rerun, camcol, field, psffitname
      psffitname = fitdir + psffitname

      ;; Fields have already been trimmed of overlap.  Don't need to
      ;; use read_tsobjmin
      openr, lun, infile, /get_lun, ERROR=error1
      IF error1 NE 0 THEN BEGIN
          print,!ERR_STRING
          free_lun, lun
          GOTO,jump2
      ENDIF 



      IF ic EQ 0 THEN BEGIN 
          pstruct = mrdfits3(lun, 1, 0, /silent)
      ENDIF ELSE BEGIN
          pstruct = mrdfits3(lun, 1, 0, /silent, /deja_vu)
      ENDELSE 
      free_lun, lun

      IF NOT keyword_set(onefile) THEN BEGIN 
          openr, lun2, psffitname, /get_lun, ERROR=error2
          IF error2 EQ 0 THEN BEGIN 
              IF ic EQ 0 THEN BEGIN 
                  IF error2 EQ 0 THEN psfstr = mrdfits4(lun2, 1, 0, /silent)
              ENDIF ELSE BEGIN 
                  IF error2 EQ 0 THEN psfstr = mrdfits4(lun2, 1, 0, /silent, /deja_vu )
              ENDELSE 
          ENDIF ELSE BEGIN 
              free_lun, lun2
              print,!ERR_STRING
              GOTO,jump2
          ENDELSE 
          free_lun, lun2
      ENDIF 

      IF keyword_set(onefile) THEN BEGIN 
          wf=where(psfstr.field EQ field)
          photo_match, $
            pstruct.run,pstruct.rerun, pstruct.camcol, pstruct.field, $
            pstruct.id,$
            psfstr[wf].run,psfstr[wf].rerun,psfstr[wf].camcol,psfstr[wf].field,  $
            psfstr[wf].id, $
            mobj, mpsf
          mpsf = wf[mpsf]
      ENDIF ELSE BEGIN 
          photo_match, $
            pstruct.run,pstruct.rerun, pstruct.camcol, pstruct.field, $
            pstruct.id,$
            psfstr.run,psfstr.rerun,psfstr.camcol,psfstr.field,  $
            psfstr.id, $
            mobj, mpsf
      ENDELSE 
;help,mobj,mpsf

;      photo_match, $
;        pstruct[mobj].run, pstruct[mobj].rerun, $
;        pstruct[mobj].camcol, pstruct[mobj].field, pstruct[mobj].id,$
;        psfstr[mpsf].run,  psfstr[mpsf].rerun,  $
;        psfstr[mpsf].camcol,  psfstr[mpsf].field,  psfstr[mpsf].id, $
;        blahobj,blahpsf
;help,blahobj,blahpsf

                                ;Initialize arrays for field
      nps = n_elements(mobj)

      ;; we want this for everything, will use for stars later
      objsize = pstruct.ixx[clr]+pstruct.iyy[clr]
      psfsize = psfstr[mpsf].psfixx[clr] + psfstr[mpsf].psfiyy[clr]
      zmag = pstruct.petrocounts[selectclr]-pstruct.reddening[selectclr]
                                ;Get good objects
      ww=where( (objsize[mobj] gt 0)   AND $
                (psfsize GT 0.)        AND $
                (zmag[mobj] le maxmag) AND $ 
                (zmag[mobj] ge minmag) AND $
                (pstruct[mobj].e1[clr] NE 1.e10) AND $
                (pstruct[mobj].momerr[clr] LE errmax) AND $
                (pstruct[mobj].r[clr] GT 0.0 ) AND $
                (abs(pstruct[mobj].ixy[clr]) LT 2. ), nww)
      IF nww NE 0 THEN BEGIN
          mobj=mobj[ww]
          mpsf=mpsf[ww]
      ENDIF ELSE BEGIN
          GOTO, JUMP1
      ENDELSE 
      ww=where(sqrt(psfsize[mpsf]/2.)*2.35*0.4  LT maxseeing, nww)
      IF nww NE 0 THEN BEGIN
          mobj=mobj[ww]
          mpsf=mpsf[ww]
      ENDIF ELSE BEGIN
          print,'Field: ',fstr,'  No objects passed seeing cut of ',maxseeing
          GOTO,jump2
      ENDELSE 
;help,mobj,mpsf      
      ;; maybe should do this with r[clr] > 0.8 or whatever
;      stars = where(pstruct[mobj].starflag[clr] EQ 1 AND $
;                    pstruct[mobj].r[clr] GT rcut AND $
;                    ( (pstruct[mobj].type[2] EQ 6) AND $
;                      (pstruct[mobj].type[1] EQ 6) ), nstars)
      stars = where( ( (pstruct[mobj].starflag AND 2b^clr) NE 0) AND $
                     ( pstruct[mobj].r[clr] GT rcut ) AND $
                     ( (pstruct[mobj].type[2] EQ 6) AND $
                       (pstruct[mobj].type[1] EQ 6) ), nstars)
      gal=mobj
      galpsf = mpsf
      IF nstars NE 0 THEN BEGIN
          remove, stars, gal, galpsf
          ;; stars and their psf's
          starsold = stars
          stars = mobj[starsold]
          stpsf = mpsf[starsold]
      ENDIF
;help,gal,galpsf
;help,stars,stpsf,nstars
      psfsize = psfstr[galpsf].psfixx[clr] + psfstr[galpsf].psfiyy[clr]
      psfrho4 = psfstr[galpsf].psfrho4[clr]
      cor1 = psfsize/objsize[gal]*(4./psfrho4-1.)/(4./pstruct[gal].rho4[clr]-1)
      cor2 = psfsize/objsize[gal]
;      cor1 = pstruct[gal].R[clr]
      ;; We find resolved galaxies from their smear polarizeability

      ww = 0
      ww=where(cor1 le rcut, nww)
      IF nww NE 0 THEN BEGIN 
          gal = gal[ww]
          galpsf = galpsf[ww]
;help,gal,galpsf
          cor1 = cor1[ww] & cor2=cor2[ww]
          psfsize = psfsize[ww] & psfrho4 = psfrho4[ww]

          ww = 0
          ww=where(psfsize GT 0., nww)
          IF nww NE 0 THEN BEGIN 
              gal = gal[ww]     
              galpsf = galpsf[ww]
;help,gal,galpsf

              cor1 = cor1[ww] & cor2 = cor2[ww] 
              psfsize = psfsize[ww] & psfrho4 = psfrho4[ww]

              psfe1 = (psfstr[galpsf].psfixx[clr] - $
                       psfstr[galpsf].psfiyy[clr])/psfsize
              psfe2 = 2.*psfstr[galpsf].psfixy[clr]/psfsize

              IF nww GT 1 THEN BEGIN 
                  ramean = mean( pstruct[gal].ra)
                  decmean = mean( pstruct[gal].dec)

                  tt            = moment(psfsize)
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

              IF (ic MOD 20 EQ 0) OR (ic EQ 0) OR (ic EQ nfields-1) THEN BEGIN
                  print,'Field: ',fstr,' Objects: ',ntostr(nps)
                  print,'Number = ',ntostr(nww),sqrt(psfmean[ic]/2.)*2.35*0.4
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
          ENDIF                 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ENDIF                     ;End of galaxy stuff
                                ;Begin star stuff
                                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      JUMP1:
      IF nstars EQ 0 THEN GOTO,JUMP2
                                ;Free up some memory
      setzero, cor1, obje1, obje2, psfe1, psfe2, zmag, ww, tmpmag, $
        psfsize, psfrho4, galpsf, gal
     
;      photo_match, $
;        pstruct[stars].run, pstruct[stars].rerun, pstruct[stars].camcol, $
;        pstruct[stars].field, pstruct[stars].id,$
;        psfstr[stpsf].run,  psfstr[stpsf].rerun,  psfstr[stpsf].camcol,  $
;        psfstr[stpsf].field,  psfstr[stpsf].id, $
;        objtmp, psftmp
;help,stars,stpsf
;help,objtmp,psftmp
;return
      psfsize = psfstr[stpsf].psfixx[clr] + psfstr[stpsf].psfiyy[clr]
      psfrho4 = psfstr[stpsf].psfrho4[clr]
      psfe1 = (psfstr[stpsf].psfixx[clr] - psfstr[stpsf].psfiyy[clr])/psfsize
      psfe2 = 2.*psfstr[stpsf].psfixy[clr]/psfsize
;      cor1=psfsize/objsize[stars]*(4./psfrho4-1.)/(4./pstruct[stars].rho4[clr]-1.0)
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
    
      JUMP2:

      setzero, pstruct, psfsize, psfrho4, objsize, obje1, obje2, $
        psfe1, psfe2, cor1, ww, stpsf, stars, obje1_cor2, obje2_cor2
      IF NOT keyword_set(onefile) THEN psfstr=0
      
  ENDFOR 
print,'Max/min: ',max(tobje1o),min(tobje1o)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Almost no modification beyond this point E.S.S.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  set_plot,'X'
;  if keyword_set(psfile) then begin
;  BEGPLOT,NAME=PSFILE
;    set_plot,'ps'
;    device,file=psfile,ysize=20,yoffset=3
;  endif
  
  IF n_elements(psfile) NE 0 THEN begplot,name=psfile
  !p.charsize=1

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
  xx=[ff[0],ff[nfields-1]]
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
  ;fitlin, gb, foredens, fsig, aa, siga, bb, sigb, /silent
  ;print, aa, siga, bb, sigb

                                ;Plots of residual gal shape vs. the 
                                ;psf shape near them after correction
  yrange=[-0.1,0.1]
  xrange=[-.2,.2]
  !p.multi=[0,2,2]

  e1err = moment(galobje1o)
  print, 'Mean corrected gal e1', e1err[0], sqrt(e1err[1])
  e2err = moment(galobje2o)
  print, 'Mean corrected gal e2', e2err[0], sqrt(e2err[1])
  erre1 = sqrt(galobjerr^2+0.32^2)
  erre2 = sqrt(galobjerr^2+0.32^2)
  xx=[-1,1]

  print,'e1PSF vs e1'
;  forprint,xo,yo,sig

  binner, galpsfe1, galobje1o, xo, yo,sig  
  plot1,xo,yo,sig,intercept_ge1pe1,slope_ge1pe1,$
    'e1 PSF '+colors[clr],'e1 '+colors[clr],xrange,yrange,$
    'corrected',$
    siga=intercept_ge1pe1_error,sigb=slope_ge1pe1_error

  binner,galpsfe1,galobje2o,xo,yo,sig  
  plot1,xo,yo,sig,intercept_ge2pe1,slope_ge2pe1,$
    'e1 PSF '+colors[clr],'e2 '+colors[clr],xrange,yrange,$
    siga=intercept_ge2pe1_error,sigb=slope_ge2pe1_error

  binner,galpsfe2,galobje1o,xo,yo,sig  
  plot1,xo,yo,sig,intercept_ge1pe2,slope_ge1pe2,$
    'e2 PSF '+colors[clr],'e1 '+colors[clr],xrange,yrange,$
    siga=intercept_ge1pe2_error,sigb=slope_ge1pe2_error

  binner,galpsfe2,galobje2o,xo,yo,sig  
  plot1,xo,yo,sig,intercept_ge2pe2,slope_ge2pe2,$
    'e2 PSF '+colors[clr],'e2 '+colors[clr],xrange,yrange,$
    siga=intercept_ge2pe2_error,sigb=slope_ge2pe2_error



                                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                ;; Same but with only size used for correction
                                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  !p.multi=[0,2,2]
  e1err = moment(galobje1o2)
  print, 'Mean corrected gal (size only) e1', e1err[0], sqrt(e1err[1])
  e2err = moment(galobje2o2)
  print, 'Mean corrected gal (size only) e2', e2err[0], sqrt(e2err[1])
  erre1 = sqrt(galobjerr^2+0.32^2)
  erre2 = sqrt(galobjerr^2+0.32^2)
  xx=[-1,1]
  binner, galpsfe1, galobje1o2, xo, yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e1 '+colors[clr],xrange,yrange,$
    'corrected size only'
  print,'e1PSF vs e1'

  binner,galpsfe1,galobje2o2,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e2 '+colors[clr],xrange,yrange

  binner,galpsfe2,galobje1o2,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e1 '+colors[clr],xrange,yrange

  binner,galpsfe2,galobje2o2,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e2 '+colors[clr],xrange,yrange

                                ;same as above, except before correction
  !p.multi=[0,2,2]
  e1err=moment(galobje1)
  print,'Mean uncorrected gal e1',e1err[0],sqrt(e1err[1])
  e2err=moment(galobje2)
  print,'Mean uncorrected gal e2',e2err[0],sqrt(e2err[1])
  erre1=sqrt(galobjerr^2+0.32^2)
  erre2=sqrt(galobjerr^2+0.32^2)
  xx=[-1,1]
  tt=moment(galobje1)
  print,tt[0],tt[1]
  binner,galpsfe1,galobje1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e1 '+colors[clr],xrange,yrange

  binner,galpsfe1,galobje2,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e2 '+colors[clr],xrange,yrange

  binner,galpsfe2,galobje1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e1 '+colors[clr],xrange,yrange

  binner,galpsfe2,galobje2,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e2 '+colors[clr],xrange,yrange

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
  ;; e vs smear polarizability
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  !p.multi=[0,0,2]
  
  binner, galcor, galobje1o, xo, yo, sig,0.025
  plot1,xo,yo,sig,intercept_ge1R,slope_ge1R,$
    'R smear','e1',[0,0.8],[-0.05,0.05],$
    siga=intercept_ge1R_error,sigb=slope_ge1R_error

  binner, galcor, galobje2o, xo, yo, sig,0.025
  plot1,xo,yo,sig,intercept_ge2R,slope_ge2R,$
    'R smear','e2',[0,0.8],[-0.05,0.05],$
    siga=intercept_ge2R_error,sigb=slope_ge2R_error

  !p.multi=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; e vs smear polarizability*epsf
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  !p.multi=[0,2,2]

  binner, galpsfe1*galcor, galobje1o, xo, yo,sig  
  plot1,xo,yo,sig,intercept_ge1pe1R,slope_ge1pe1R,$
    'e1 PSF*R '+colors[clr],'e1 '+colors[clr],xrange,yrange,$
    'corrected',$
    siga=intercept_ge1pe1R_error,sigb=slope_ge1pe1R_error

  binner,galpsfe1*galcor,galobje2o,xo,yo,sig  
  plot1,xo,yo,sig,intercept_ge2pe1R,slope_ge2pe1R,$
    'e1 PSF*R '+colors[clr],'e2 '+colors[clr],xrange,yrange,$
    siga=intercept_ge2pe1R_error,sigb=slope_ge2pe1R_error

  binner,galpsfe2*galcor,galobje1o,xo,yo,sig  
  plot1,xo,yo,sig,intercept_ge1pe2R,slope_ge1pe2R,$
    'e2 PSF*R '+colors[clr],'e1 '+colors[clr],xrange,yrange,$
    siga=intercept_ge1pe2R_error,sigb=slope_ge1pe2R_error

  binner,galpsfe2*galcor,galobje2o,xo,yo,sig  
  plot1,xo,yo,sig,intercept_ge2pe2R,slope_ge2pe2R,$
    'e2 PSF*R '+colors[clr],'e2 '+colors[clr],xrange,yrange,$
    siga=intercept_ge2pe2R_error,sigb=slope_ge2pe2R_error

  !p.multi=0


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Star plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                ;Size-mag diagram for stars
  ;;plot,mag,size,psym=3,xtitle=colors[clr],ytitle="FWHM (arcsec)"
  ploth,mag,size,xtitle=colors[clr],ytitle="FWHM (arcsec)",$
    nxbins=50,nybins=50
                                ;smear polarizability for all objects
  yrr=[0.,2.0]
  xrr=[12.,maxmag]
  ;;plot,galmag,galcor,psym=3,xtitle=colors[clr],ytitle="Rsm",yrange=yrr,xrange=xrr
  ;;oplot,mag,tcor,psym=4,symsize=.5
  ploth,[galmag,mag],[galcor,tcor],$
    xtitle=colors[clr],ytitle="Rsm",yrange=yrr,xrange=xrr,$
    nxbins=100,nybins=100


  sigma_clip, tcor, meantcor, sigmatcor, nsig=3.5, niter=3
  momtmp = moment(tcor)
  meantcor2 = momtmp[0]
  sigmatcor2 = sqrt(momtmp[1])
  print
  print,'Mean star smear polarizability: ',meantcor, sigmatcor
  print,'straight mean: ',meantcor2,sigmatcor2
  print

;  plothist,diff,bin=.01,xtitle='cor1-r[clr]'
;  plothist,galdiff,bin=.01,xtitle='galcor1-gal.r[clr]'

  !p.multi=[0,2,2]

                                ;star shape vs. psf shape
                                ;before correction

  e1err = moment(tobje1)
  print,'Mean star e1',e1err[0],sqrt(e1err[1])
  e2err = moment(tobje2)
  print,'Mean star e2',e2err[0],sqrt(e2err[1])
  erre1 = sqrt(tobjerr^2+0.32^2)
  erre2 = sqrt(tobjerr^2+0.32^2)
  xx=[-1,1]
  binner,tpsfe1,tobje1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e1 '+colors[clr],xrange,yrange

  binner,tpsfe1,tobje2,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e2 '+colors[clr],xrange,yrange

  binner,tpsfe2,tobje1,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e1 '+colors[clr],xrange,yrange

  binner,tpsfe2,tobje2,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e2 '+colors[clr],xrange,yrange


  !p.multi=[0,2,2]

                                ;Residual star shape vs. psf shape
                                ;after correction

  e1err = moment(tobje1o)
  print,'Mean corrected star e1',e1err[0],sqrt(e1err[1])
  e2err = moment(tobje2o)
  print,'Mean corrected star e2',e2err[0],sqrt(e2err[1])
  erre1 = sqrt(tobjerr^2+0.32^2)
  erre2 = sqrt(tobjerr^2+0.32^2)
  xx=[-1,1]
  binner,tpsfe1,tobje1o,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e1 '+colors[clr],xrange,yrange

  binner,tpsfe1,tobje2o,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e2 '+colors[clr],xrange,yrange

  binner,tpsfe2,tobje1o,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e1 '+colors[clr],xrange,yrange

  binner,tpsfe2,tobje2o,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e2 '+colors[clr],xrange,yrange


                                ;Try correcting with cor=1.0
  e1err = moment(tobje1o_cor2)
  print,'Mean corrected star e1 (R=1)',e1err[0],sqrt(e1err[1])
  e2err = moment(tobje2o_cor2)
  print,'Mean corrected star e2 (R=1)',e2err[0],sqrt(e2err[1])
  erre1 = sqrt(tobjerr^2+0.32^2)
  erre2 = sqrt(tobjerr^2+0.32^2)
  xx=[-1,1]
  binner,tpsfe1,tobje1o_cor2,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e1 '+colors[clr],xrange,yrange

  binner,tpsfe1,tobje2o_cor2,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e1 PSF '+colors[clr],'e2 '+colors[clr],xrange,yrange

  binner,tpsfe2,tobje1o_cor2,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e1 '+colors[clr],xrange,yrange

  binner,tpsfe2,tobje2o_cor2,xo,yo,sig  
  plot1,xo,yo,sig,aa,bb,'e2 PSF '+colors[clr],'e2 '+colors[clr],xrange,yrange

  IF n_elements(psfile) NE 0 THEN endplot


;  if keyword_set(psfile) then begin
;    device,/close
;  endif
;  set_plot,'X'

  ;; For output
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

  corshape_again, fitsfile

  !p.multi=0
  ptime, systime(1)-time

  status=0

end

