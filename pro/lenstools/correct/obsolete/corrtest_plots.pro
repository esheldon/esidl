PRO corrtest_plots_plot1,x,y,sig,a,b,xtitle,ytitle,xrange,yrange,title,$
          siga=siga,sigb=sigb

  myusersym, 'fill_circle'
  psym=8

  xx=[-1,1]

  srt=sort(x)

  fitlin,x,y,sig,a,siga,b,sigb,/silent
  aploterror,1,x[srt],y[srt],sig[srt],xtitle=xtitle,ytitle=ytitle,yrange=yrange,xrange=xrange,title=title,psym=psym
  yy=a+b*xx
  oplot,xx,yy
  sa=strmid(strtrim(string(a,format='(e8.1)'),2),0,8)
  ssiga=strmid(strtrim(string(abs(siga/a)),2),0,5)
  sb=strmid(strtrim(string(b,format='(e8.1)'),2),0,8)
  ssigb=strmid(strtrim(string(abs(sigb/b)),2),0,5)
  stta='a='+sa+' '+ssiga
  sttb='b='+sb+' '+ssigb
;  xyouts,xrange[0]+0.01,yrange[0]+0.8*(yrange[1]-yrange[0]),stt,font=-1
  legend, [stta,sttb], /right, /clear
;  print,'yaya',a,b,strtrim(string(a),2),strtrim(string(b,format='(e8.1)'),2)
  
return

end


PRO corrtest_plots_getrange, x1, y1, x2, y2, nperbin, xrange, yrange

  binner_bynum, x1, y1, nperbin, txo1, tyo1, tsig1, tnum1
  binner_bynum, x2, y2, nperbin, txo2, tyo2, tsig2, tnum2

  w1=where(tnum1 EQ nperbin)
  w2=where(tnum2 EQ nperbin)
  xrange=prange(txo1[w1],txo2[w2],/noerror,slack=0.3,/symmetric)
  yrange=prange(tyo1[w1],tyo2[w2],/noerror,slack=0.4,/symmetric)

END 

pro corrtest_plots, run, rerun, camcol, clr, $
                    start=start, nframes=nframes,$
                    minmag=minmag, maxmag=maxmag, $
                    hardrcut=hardrcut, hardprobcut=hardprobcut, $
                    star_hardrcut=star_hardrcut, $
                    ermscut=ermscut, errmax=errmax, $
                    maxseeing=maxseeing, $
                    offset=offset, onefile=onefile, status=status,$
                    no_overwrite=no_overwrite, hirata=hirata, varcut=varcut,$
                    outdir=outdir

  ;; status is 1 unless we reach end
  status=1

;   ermscut = cut off frames with psfe rms gt ermscut
  
  IF n_params() EQ 0 THEN BEGIN 
      print,'-syntax corrtest_plots, run, rerun, camcol, clr, $'
      print,'            start=start, nframes=nframes,$'
      print,'            minmag=minmag, maxmag=maxmag, $'
      print,'            hardrcut=hardrcut, hardprobcut=hardprobcut, $'
      print,'            star_hardrcut=star_hardrcut, $'
      print,'            ermscut=ermscut, errmax=errmax, $'
      print,'            maxseeing=maxseeing, $'
      print,'            offset=offset, onefile=onefile, status=status,$'
      print,'            no_overwrite=no_overwrite, hirata=hirata, varcut=varcut,$'
      print,'            outdir=outdir'

      return
  ENDIF 
  
  time=systime(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  newfront = 'adatc'            ;The front of files with adaptive moments

  selectclr = 2                 ;r-band

  colors=['u','g','r','i','z']
  bands = [1,2,3]               ; we will process g,r,i
  nband = n_elements(bands)
  
  rstr = run2string(run)
  rrstr = ntostr(rerun)

  ;; cuts used in make_corrected_files for stars
  stmaxmag = [0.0,21.,20.,20.,0.0]
  IF n_elements(maxmag) EQ 0 THEN maxmag=22.
  IF n_elements(minmag) EQ 0 THEN minmag=18.
  IF n_elements(hardrcut) EQ 0 THEN hardrcut = !hardrcut
  IF n_elements(hardprobcut) EQ 0 THEN hardprobcut = !hardprobcut
  IF n_elements(star_hardrcut) EQ 0 THEN star_hardrcut = 0.9
  IF n_elements(ermscut) EQ 0 THEN ermscut = 1000.
  IF n_elements(errmax) EQ 0 THEN errmax = 0.4
  IF n_elements(maxseeing) EQ 0 THEN maxseeing=3.0

  print
  print,'Hard Rcut: '+ntostr(hardrcut)
  print,'Hard Probcut: '+ntostr(hardprobcut)
  print,'Errmax: '+ntostr(errmax)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; For selecting good objects
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  make_flag_struct, fs
  fs.AMOMENT_FAINT = 'N'
  fs.AMOMENT_UNWEIGHTED = 'N'
  fs.AMOMENT_SHIFT = 'N'
  fs.AMOMENT_MAXITER = 'N'
  fs.DEBLENDED_AS_PSF = 'N'

  edef = -9999.                 ;default value for m_e[12]_corr[*]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; For star-galaxy separation
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  IF keyword_set(varcut) THEN BEGIN 
      purity = 0.99
      rsmear_cuts_files, run, rerun, purity, clr, rsmear_file, $
                         camcol=camcol, hirata=hirata
      IF NOT file_test(rsmear_file) THEN BEGIN 
          print
          print,'Attempting to make rsmear files for camcol: '+ntostr(camcol)
          rsmear_cuts, run, rerun, purity, 1, tstruct, $
                       camcol=camcol, /overwrite, hirata=hirata
          rsmear_cuts, run, rerun, purity, 2, tstruct, $
                       camcol=camcol, /overwrite, hirata=hirata
          rsmear_cuts, run, rerun, purity, 3, tstruct, $
                       camcol=camcol, /overwrite, hirata=hirata
          IF NOT file_test(rsmear_file) THEN BEGIN 
              print
              print,'Failed to make rsmear files for camcol: '+ntostr(camcol)
              status = 1
              return
          ENDIF 
      ENDIF 
      read_rsmear_cuts, run, rerun, purity, clr, rcutstr, $
                        camcol=camcol, status=status, hirata=hirata
      IF status NE 0 THEN return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Setup the column
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  cstr = ntostr(camcol)
  
  fetch_dir, run, camcol, rerun, dir, atldir, corrdir=corrdir,$
    corratldir = fitdir
  fetch_file_list, corrdir, adatcfiles, fnums, start=start, nframes=nframes, $
    fieldmin=fieldmin, fieldmax=fieldmax, front=newfront

  IF n_elements(outdir) NE 0 THEN fitdir = outdir

  nfields = n_elements(adatcfiles)

  IF keyword_set(offset) THEN addstr='offset' ELSE addstr=''

  IF keyword_set(hirata) THEN rsmearstr = '_h' ELSE rsmearstr=''
  fitsfile = fitdir+addstr+'corshape_'+rstr+'_'+cstr+'_'+$
    !colors[clr]+rsmearstr+'_N1.fit'

  IF keyword_set(no_overwrite) THEN BEGIN 
      WHILE fexist(fitsfile) DO fitsfile=newname(fitsfile)
  ENDIF 
  psfile = repstr(fitsfile, '.fit','.ps')

  print
  print,'psfile: ',psfile
  print,'fitsfile: ',fitsfile
  print

  print
  print,'corrected files are '+newfront+'-*'
  print
  
  galflag   = 0
  starflag  = 0
  ebv       = fltarr(nfields)
  ff        = fltarr(nfields)
  seeing    = fltarr(nfields)
  seeingsig = fltarr(nfields)
  psfe1mean = fltarr(nfields)
  psfe1sig  = fltarr(nfields)
  psfe2mean = fltarr(nfields)
  psfe2sig  = fltarr(nfields)
  foredens  = fltarr(nfields)
  gl        = fltarr(nfields)
  gb        = fltarr(nfields)
  
  print,'-----------------------------------------------------'
  print,' Run: ',rstr,' Rerun: ',rrstr,' Camcol: ',cstr,' Band: ',!colors[clr]
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
                  tsobjstr=tsobjstr, verbose=0, front='adatc'

      IF n_elements(pstruct) NE 0 THEN BEGIN 
          
;          IF keyword_set(hirata) THEN BEGIN 
;              compea4_struct, pstruct, clr, $
;                              e1_out, e2_out, R_out, compflags
;              pstruct.m_e1_corr_h[clr] = e1_out
;              pstruct.m_e2_corr_h[clr] = e2_out
;              pstruct.m_R_h[clr] = R_out
;              pstruct.CompEA4flags[clr] = compflags
;          ENDIF 

          ;; we want this for everything, will use for stars later
          objsize = pstruct.m_rr_cc[clr]
          psfsize = pstruct.m_rr_cc_psf[clr]

          zmag = pstruct.petrocounts[selectclr]-pstruct.reddening[selectclr]

          flag_select, pstruct, fs, clr, good

          IF good[0] NE -1 THEN BEGIN 
             
              ;; Further cuts
              IF keyword_set(hirata) THEN BEGIN 
                  ;; remember hirata's stuff is already 
                  ;; dilution corrected
                  e1err = pstruct[good].m_e1e1err[clr]
                  e1err = e1err/(1. - pstruct[good].m_R_h[clr])
                  e2err = pstruct[good].m_e2e2err[clr]
                  e2err = e2err/(1. - pstruct[good].m_R_h[clr])

                  good2=where( (zmag[good] LE maxmag)            AND $ 
                               (zmag[good] GE minmag)            AND $
                               (pstruct[good].m_R_h[clr] GT 0.0) AND $
                               (pstruct[good].CompEA4flags[clr] EQ 0) AND $
                               abs(pstruct[good].m_e1_corr_h[clr]) LT 2 AND $
                               abs(pstruct[good].m_e2_corr_h[clr]) LT 2 AND $
                               (e1err LE errmax)                  AND $
                               (e2err LE errmax)                  AND $
                               (pstruct[good].seeing[clr] LT maxseeing), nobj)
              ENDIF ELSE BEGIN 
                  e1err = pstruct[good].m_e1e1err[clr]
                  e1err = e1err/(1. - pstruct[good].m_R[clr])
                  e2err = pstruct[good].m_e2e2err[clr]
                  e2err = e2err/(1. - pstruct[good].m_R[clr])

                  good2=where( (zmag[good] LE maxmag)          AND $ 
                               (zmag[good] GE minmag)          AND $
                               (pstruct[good].m_R[clr] GT 0.0) AND $
                               abs(pstruct[good].m_e1_corr[clr]) LT 2 AND $
                               abs(pstruct[good].m_e2_corr[clr]) LT 2 AND $
                               (e1err LE errmax)                AND $
                               (e2err LE errmax)                AND $
                               (pstruct[good].seeing[clr] LT maxseeing), nobj)
              ENDELSE 
              
              IF nobj NE 0 THEN BEGIN 
                  good = good[good2]

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; galaxies
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              
                  IF keyword_set(varcut) THEN BEGIN 
                      gal_rcut = interpol(rcutstr.rsmear_cuts, $
                                          rcutstr.meanmag, $
                                          zmag[good]) < hardrcut
                  ENDIF ELSE BEGIN 
                      gal_rcut = hardrcut
                  ENDELSE 

                  IF keyword_set(hirata) THEN BEGIN 
                      gal = where( (pstruct[good].objc_prob_psf GE 0.0) AND $
                                   (pstruct[good].objc_prob_psf LT (1.-hardprobcut)) AND $
                                   (pstruct[good].m_r_h[clr] LT gal_rcut), $
                                   ngal, complement=compgal, ncomp=ncompgal)
                  ENDIF ELSE BEGIN 
                      gal = where( (pstruct[good].objc_prob_psf GE 0.0) AND $
                                   (pstruct[good].objc_prob_psf LT (1.-hardprobcut)) AND $
                                   (pstruct[good].m_r[clr] LT gal_rcut), $
                                   ngal, complement=compgal, ncomp=ncompgal)
                  ENDELSE 

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; stars
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  IF ncompgal EQ 0 THEN BEGIN 
                      nstars = 0
                  ENDIF ELSE BEGIN 
                      compgal = good[compgal]
                      IF keyword_set(varcut) THEN BEGIN 
                          star_rcut = interpol(rcutstr.rsmear_cuts, $
                                               rcutstr.meanmag, $
                                               zmag[compgal]) > star_hardrcut
                      ENDIF ELSE BEGIN 
                          star_rcut = star_hardrcut
                      ENDELSE 

                      IF keyword_set(hirata) THEN BEGIN 
                          stars = where( ( pstruct[compgal].m_R_h[clr] GT star_rcut ) AND $
                                         ( zmag[compgal] LT stmaxmag[clr] )           AND $
                                         ( (pstruct[compgal].type[2] EQ 6)            AND $
                                           (pstruct[compgal].type[1] EQ 6) )          AND $
                                         (abs(pstruct[compgal].m_e1[clr] LT 2))       AND $
                                         (abs(pstruct[compgal].m_e2[clr] LT 2)), nstars)
                      ENDIF ELSE BEGIN 
                          stars = where( ( pstruct[compgal].m_R[clr] GT star_rcut ) AND $
                                         ( zmag[compgal] LT stmaxmag[clr] )         AND $
                                         ( (pstruct[compgal].type[2] EQ 6)          AND $
                                           (pstruct[compgal].type[1] EQ 6) )        AND $
                                         (abs(pstruct[compgal].m_e1[clr]) LT 2)     AND $
                                         (abs(pstruct[compgal].m_e2[clr]) LT 2), nstars)
                      ENDELSE 

                  ENDELSE 
              
                  IF ngal NE 0 THEN BEGIN 
                  
                      gal = good[gal]
                      ngal = n_elements(gal)
                      IF keyword_set(hirata) THEN BEGIN 
                          corr = pstruct[gal].m_R_h[clr]
                      ENDIF ELSE BEGIN 
                          corr = pstruct[gal].m_R[clr]
                      ENDELSE 
                      corr2 = psfsize[gal]/objsize[gal]
                      
                      psfe1 = pstruct[gal].m_e1_psf[clr]
                      psfe2 = pstruct[gal].m_e2_psf[clr]
                      
                      IF ngal GT 1 THEN BEGIN 
                          ramean = mean( pstruct[gal].ra)
                          decmean = mean( pstruct[gal].dec)
                          
                          tt            = moment(pstruct[gal].seeing[clr])
                          seeing[ic]    = tt[0]
                          seeingsig[ic] = sqrt(tt[1])
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
                          
                          seeing[ic]    = pstruct[gal[0]].seeing[clr]
                          seeingsig[ic] = 0.
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
                      
                      IF ((ic MOD 20 EQ 0) OR (ic EQ 0) OR $
                          (ic EQ nfields-1) )THEN BEGIN
                          print,'Field: ',fstr,' Objects: ',ntostr(nobj)
                          print,'Number = ',ntostr(ngal),seeing[ic]
                      ENDIF 

                      ;; either dilution correct or not
                      ;;fac = 1.0
                      fac = 1./(1. - pstruct[gal].m_R[clr])

                      ;; Original e1, e2 measurements
                      e1 = pstruct[gal].m_e1[clr]
                      e2 = pstruct[gal].m_e2[clr]

                      ;; corrected (not hirata)
                      e1_corr = pstruct[gal].m_e1_corr[clr]*fac
                      e2_corr = pstruct[gal].m_e2_corr[clr]*fac

                      ;; corrected with just size info
                      e1_corr2 = (e1 - corr2*psfe1)*fac
                      e2_corr2 = (e2 - corr2*psfe2)*fac
                      

                      IF (psfe1sig[ic] LE ermscut AND $
                          psfe2sig[ic] LE ermscut) THEN BEGIN

                          add_arrval, psfe1, galpsfe1
                          add_arrval, psfe2, galpsfe2
                          add_arrval, pstruct[gal].m_e1e1err[clr], gale1err
                          add_arrval, pstruct[gal].m_e2e2err[clr], gale2err
                          add_arrval, e1, gale1
                          add_arrval, e2, gale2

                          IF keyword_set(hirata) THEN BEGIN 
                              add_arrval, pstruct[gal].m_e1_corr_h[clr], gale1_corr
                              add_arrval, pstruct[gal].m_e2_corr_h[clr], gale2_corr
                          ENDIF ELSE BEGIN 
                              add_arrval, e1_corr, gale1_corr
                              add_arrval, e2_corr, gale2_corr
                          ENDELSE 

                          add_arrval, e1_corr2, gale1_corr2
                          add_arrval, e2_corr2, gale2_corr2
                          add_arrval, pstruct[gal].petrocounts[clr]-$
                                      pstruct[gal].reddening[clr], galmag 
                          add_arrval, corr, galcorr
                      ENDIF 
                  ENDIF ELSE BEGIN ;; ngal ne 0
                      print,'-- No Galaxies'
                  ENDELSE 

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; begin stars
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          
                  IF nstars NE 0 THEN BEGIN 
                      stars = compgal[stars]
                                ;Free up some memory
                      setzero, corr, e1, e2, e1_corr2, e2_corr2,$
                               psfe1, psfe2, $
                               ww, tmpmag, galpsf, gal
                      
                      
                      psfe1 = pstruct[stars].m_e1_psf[clr]
                      psfe2 = pstruct[stars].m_e2_psf[clr]
                      
                      IF keyword_set(hirata) THEN BEGIN 
                          corr = pstruct[stars].m_R_h[clr]
                      ENDIF ELSE BEGIN 
                          corr = pstruct[stars].m_R[clr]
                      ENDELSE 

                      e1 = pstruct[stars].m_e1[clr]
                      e2 = pstruct[stars].m_e2[clr]

                      ;; try correcting with corr=1.0!
                      e1_corr2 = e1 - psfe1
                      e2_corr2 = e2 - psfe2
                  
                      add_arrval, psfe1, starpsfe1
                      add_arrval, psfe2, starpsfe2
                      ;; add_arrval, pstruct[stars].momerr[clr], tobjerr
                      add_arrval, pstruct[stars].m_e1e1err[clr], stare1err
                      add_arrval, pstruct[stars].m_e2e2err[clr], stare2err
                      IF keyword_set(hirata) THEN BEGIN 
                          add_arrval, pstruct[stars].m_e1_corr_h[clr], $
                                      stare1_corr
                          add_arrval, pstruct[stars].m_e2_corr_h[clr], $
                                      stare2_corr
                      ENDIF ELSE BEGIN 
                          add_arrval, pstruct[stars].m_e1_corr[clr], $
                                      stare1_corr
                          add_arrval, pstruct[stars].m_e2_corr[clr], $
                                      stare2_corr
                      ENDELSE 
                      add_arrval, e1_corr2, stare1_corr2
                      add_arrval, e2_corr2, stare2_corr2
                      add_arrval, e1, stare1
                      add_arrval, e2, stare2
                      add_arrval, zmag[stars], starmag
                      add_arrval, corr, starcorr
                      add_arrval, corr-pstruct[stars].m_r[clr], diff
                      add_arrval, sqrt(objsize[stars]/2.)*2.35*0.40, starsize
                  
                  ENDIF ELSE BEGIN ;; nstar eq 0
                      print,'-- No Stars'
                  ENDELSE 

                  add_arrval, zmag[good], allmag
                  IF keyword_set(hirata) THEN BEGIN 
                      add_arrval, pstruct[good].m_R_h[clr], allcorr
                  ENDIF ELSE BEGIN 
                      add_arrval, pstruct[good].m_R[clr], allcorr
                  ENDELSE 

              ENDIF ELSE BEGIN ;; nobj eq 0
                  print,'Field: ',fstr,$
                        '  No objects passed mom cuts and seeing cut of ',$
                        maxseeing
              ENDELSE  
          ENDIF ELSE BEGIN ;; AMOMENT checks passed no objects
              print,'Field: ',fstr,'  No objects passed AMOMENT flag cuts'
          ENDELSE 
      ENDIF ELSE BEGIN ;; field empty
          ;print,!ERR_STRING
          ;free_lun, lun
      ENDELSE 
      delvarx, pstruct
      setzero, psfsize, objsize, e1, e2, $
               psfe1, psfe2, corr, ww, stars, e1_corr2, e2_corr2, zmag
  ENDFOR 


  IF n_elements(psfile) NE 0 THEN begplot,name=psfile,/color
  setup_mystuff
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
  !p.multi=[0,1,3]
                                ;Mean psf vs field
  tsize=2.0
  plot,ff, seeing, xtitle='Field', ytitle='FWHM (arcsec)',$
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
  aplot, !gratio,gb, foredens, xtitle='Latitude',ytitle='# of Galaxies',$
         psym=fpsym,$
         charsize=tsize
                                ; Fit to the number above
  fsig=sqrt(foredens)

  ;;;;;;;;;;;;;
  ;; binsizes
  ;;;;;;;;;;;;;

  ebin = 0.015
  erbin = ebin/2.
  rbin = 0.05

  ngal = n_elements(gale1_corr)
  nstar = n_elements(stare1_corr)

  print
  print,"-------------------------------------------"
  print,'# of galaxies: '+ntostr(ngal)
  print,'# of stars:    '+ntostr(nstar)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plots of egal vs epsf
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,"-------------------------------------------"

  !p.multi=[0,2,2]

  e1mom = moment(gale1_corr)
  print, 'Mean corrected gal e1', e1mom[0], sqrt(e1mom[1]/ngal)
  e2mom = moment(gale2_corr)
  print, 'Mean corrected gal e2', e2mom[0], sqrt(e2mom[1]/ngal)
  erre1 = sqrt(gale1err^2+0.32^2)
  erre2 = sqrt(gale2err^2+0.32^2)
  xx=[-1,1]

  ;;; about 8 bins
  frac = 1./9.
  nperbin = long( round(frac*n_elements(galpsfe1)) )

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find out xrange for gals
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  corrtest_plots_getrange, galpsfe1, gale1, galpsfe2, gale2, nperbin,$
                           xrange, yrange

  erase & multiplot, [2,2], /square

  plot_fitlin, galpsfe1, gale1_corr, nperbin, $
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
  plot_fitlin, galpsfe2, gale1_corr, nperbin, $
               intercept_ge1pe2, intercept_ge1pe2_error, $
               slope_ge1pe2, slope_ge1pe2_error,$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Gal"

  cyrange=!y.crange 
  plothist,galpsfe2,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe1, gale2_corr, nperbin, $
               intercept_ge2pe1, intercept_ge2pe1_error, $
               slope_ge2pe1, slope_ge2pe1_error,$
               ytitle="e!D2!N Gal",xtitle="e!D1!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, galpsfe2, gale2_corr, nperbin, $
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
  e1mom = moment(gale1_corr2)
  print, 'Mean corrected gal (size only) e1', e1mom[0], sqrt(e1mom[1]/ngal)
  e2mom = moment(gale2_corr2)
  print, 'Mean corrected gal (size only) e2', e2mom[0], sqrt(e2mom[1]/ngal)
  erre1 = sqrt(gale1err^2+0.32^2)
  erre2 = sqrt(gale2err^2+0.32^2)
  xx=[-1,1]

  erase & multiplot, [2,2], /square

  plot_fitlin, galpsfe1, gale1_corr2, nperbin, $
               ytitle="e!D1!N Gal",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band corrected    Size Only";,/plotnum

  cyrange=!y.crange 
  plothist,galpsfe1,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe2, gale1_corr2, nperbin, $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Gal"

  cyrange=!y.crange 
  plothist,galpsfe2,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe1, gale2_corr2, nperbin, $
               ytitle="e!D2!N Gal",xtitle="e!D1!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, galpsfe2, gale2_corr2, nperbin, $
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
  e1mom=moment(gale1)
  print,'Mean uncorrected gal e1',e1mom[0],sqrt(e1mom[1]/ngal)
  e2mom=moment(gale2)
  print,'Mean uncorrected gal e2',e2mom[0],sqrt(e2mom[1]/ngal)
  erre1=sqrt(gale1err^2+0.32^2)
  erre2=sqrt(gale2err^2+0.32^2)
  xx=[-1,1]


  erase & multiplot, [2,2], /square

  plot_fitlin, galpsfe1, gale1, nperbin, $
               ytitle="e!D1!N Gal",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Uncorrected";,/plotnum

  cyrange=!y.crange 
  plothist,galpsfe1,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe2, gale1, nperbin, $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Gal"

  cyrange=!y.crange 
  plothist,galpsfe2,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe1, gale2, nperbin, $
               ytitle="e!D2!N Gal",xtitle="e!D1!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, galpsfe2, gale2, nperbin, $
               xtitle="e!D2!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D2!N Gal"

  multiplot,/reset

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; magnitude histogram
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  !p.multi=[0,0,2]
  plothist, galmag, galmagxhist, galmagyhist, bin=0.1, min=minmag, max=maxmag, $
            ytitle='N', xtitle=!colors[clr]+'-mag'
  legend, 'N!Dtotal!N = '+ntostr(long(total(galmagyhist))), /right, box=0

;  histo=histogram(galmag,binsize=0.5,min=15.,max=25.)
;  bin=fltarr(1)
;  bin[0]=15.25
;  for i=15.75,24.75,0.5 do begin
;    bin=[bin,i]
;  endfor
;  aplot,!gratio,bin,histo,xtitle=!colorsp[clr],ytitle='N',$
;        charsize=tsize, psym=10

  !p.multi=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; e vs smear polarizability
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  !p.multi=[0,0,2]
  
  pold=!p.charsize
  !p.charsize=1.25
  binner, galcorr, gale1_corr, rbin, xo, yo, sig
  corrtest_plots_plot1,xo,yo,sig,intercept_ge1R,slope_ge1R,$
    'R smear','e1',[0,0.8],[-0.05,0.05],!colors[clr]+"-band corrected",$
    siga=intercept_ge1R_error,sigb=slope_ge1R_error

  cyrange=!y.crange 
  plothist,galcorr,xhist,yhist,bin=rbin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  binner, galcorr, gale2_corr, rbin, xo, yo, sig
  corrtest_plots_plot1,xo,yo,sig,intercept_ge2R,slope_ge2R,$
    'R smear','e2',[0,0.8],[-0.05,0.05],$
    siga=intercept_ge2R_error,sigb=slope_ge2R_error

  !p.charsize=pold
  !p.multi=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; e vs smear polarizability*epsf
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  erase & multiplot, [2,2], /square

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find out xrange for gals including galcorr
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  corrtest_plots_getrange, galpsfe1*galcorr, gale1, galpsfe2*galcorr, $
                           gale2, nperbin,$
                           xrange, tyrange ;tyrange not used

  plot_fitlin, galpsfe1*galcorr, gale1_corr, nperbin, $
               intercept_ge1pe1R, intercept_ge1pe1R_error, $
               slope_ge1pe1R, slope_ge1pe1R_error,$
               ytitle="e!D1!N Gal",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band corrected";,/plotnum

  cyrange=!y.crange 
  plothist,galpsfe1*galcorr,xhist,yhist,bin=erbin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe2*galcorr, gale1_corr, nperbin, $
               intercept_ge1pe2R, intercept_ge1pe2R_error, $
               slope_ge1pe2R, slope_ge1pe2R_error,$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Gal"

  cyrange=!y.crange 
  plothist,galpsfe2*galcorr,xhist,yhist,bin=erbin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, galpsfe1*galcorr, gale2_corr, nperbin, $
               intercept_ge2pe1R, intercept_ge2pe1R_error, $
               slope_ge2pe1R, slope_ge2pe1R_error,$
               ytitle="e!D2!N Gal",xtitle="e!D1!N PSF*R",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, galpsfe2*galcorr, gale2_corr, nperbin, $
               intercept_ge2pe2R, intercept_ge2pe2R_error, $
               slope_ge2pe2R, slope_ge2pe2R_error,$
               xtitle="e!D2!N PSF*R",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D2!N Gal"

  multiplot,/reset

  pold=!p.charsize
  !p.charsize=1.75

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Size-mag diagram for stars
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  yrange=[-0.2,0.2]

  loadct,0
  ploth,starmag,starsize,xtitle=!colors[clr]+'-band mag',ytitle="FWHM (arcsec)",$
    nxbins=50,nybins=50,/silent,title='Size of stars'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; smear polarizability vs. mag for all objects
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  yrr=[0.,2.0]
  xrr=[12.,maxmag]

  ploth,allmag,allcorr,$
    xtitle=!colors[clr]+'-band mag',ytitle="R!Dsm!N",yrange=yrr,$
    nxbins=100,nybins=100,/silent,/sqrt

  oplot, [0.0, 1000.], [hardrcut, hardrcut];, color=!grey50
;  oplot, [0.0, 1000.], [star_hardrcut, star_hardrcut];, color=!grey50
;  oplot, [stmaxmag[clr],stmaxmag[clr]], [star_hardrcut, 1000.];, color=!grey50
  simpctable

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Star plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  !p.charsize=pold

  sigma_clip, starcorr, meanstarcorr, sigmastarcorr, nsig=3.5, niter=3, /silent
  momtmp = moment(starcorr)
  meanstarcorr2 = momtmp[0]
  sigmastarcorr2 = sqrt(momtmp[1])
  print
  print,"-------------------------------------------"
  print,'Mean star smear polarizability: ',meanstarcorr, sigmastarcorr
  print,'straight mean: ',meanstarcorr2,sigmastarcorr2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;star shape vs. psf shape before correction
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,"-------------------------------------------"
  e1mom = moment(stare1)
  print,'Mean star e1',e1mom[0],sqrt(e1mom[1]/nstar)
  e2mom = moment(stare2)
  print,'Mean star e2',e2mom[0],sqrt(e2mom[1]/nstar)
  erre1 = sqrt(stare1err^2+0.32^2)
  erre2 = sqrt(stare2err^2+0.32^2)
  xx=[-1,1]

  sfrac = 1./9.
  snperbin = long( round(frac*n_elements(starpsfe1)) )

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find out xrange for stars
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  corrtest_plots_getrange, starpsfe1, stare1, starpsfe2, stare2, snperbin,$
                           xrange, yrange

  erase & multiplot, [2,2], /square

  plot_fitlin, starpsfe1, stare1, snperbin, $
               ytitle="e!D1!N Star",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Uncorrected";,/plotnum

  cyrange=!y.crange 
  plothist,starpsfe1,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, starpsfe2, stare1, snperbin, $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Star"

  cyrange=!y.crange 
  plothist,starpsfe2,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, starpsfe1, stare2, snperbin, $
               ytitle="e!D2!N Star",xtitle="e!D1!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, starpsfe2, stare2, snperbin, $
               xtitle="e!D2!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D2!N Star"

  multiplot,/reset

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Residual star shape vs. psf shape after correction
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,"-------------------------------------------"
  e1mom = moment(stare1_corr)
  print,'Mean corrected star e1',e1mom[0],sqrt(e1mom[1]/nstar)
  e2mom = moment(stare2_corr)
  print,'Mean corrected star e2',e2mom[0],sqrt(e2mom[1]/nstar)
  erre1 = sqrt(stare1err^2+0.32^2)
  erre2 = sqrt(stare2err^2+0.32^2)
  xx=[-1,1]

  erase & multiplot, [2,2], /square

  plot_fitlin, starpsfe1, stare1_corr, snperbin, $
               ytitle="e!D1!N Star",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Corrected";,/plotnum

  cyrange=!y.crange 
  plothist,starpsfe1,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, starpsfe2, stare1_corr, snperbin, $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Star"

  cyrange=!y.crange 
  plothist,starpsfe2,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, starpsfe1, stare2_corr, snperbin, $
               ytitle="e!D2!N Star",xtitle="e!D1!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, starpsfe2, stare2_corr, snperbin, $
               xtitle="e!D2!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D2!N Star"

  multiplot,/reset

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Try correcting with cor=1.0
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,"-------------------------------------------"
  e1mom = moment(stare1_corr2)
  print,'Mean corrected star e1 (R=1)',e1mom[0],sqrt(e1mom[1]/nstar)
  e2mom = moment(stare2_corr2)
  print,'Mean corrected star e2 (R=1)',e2mom[0],sqrt(e2mom[1]/nstar)
  erre1 = sqrt(stare1err^2+0.32^2)
  erre2 = sqrt(stare2err^2+0.32^2)
  xx=[-1,1]

  erase & multiplot, [2,2], /square

  plot_fitlin, starpsfe1, stare1_corr2, snperbin, $
               ytitle="e!D1!N Star",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1,$
               title=!colors[clr]+"-band  Corrected (R=1)";,/plotnum

  cyrange=!y.crange 
  plothist,starpsfe1,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, starpsfe2, stare1_corr2, snperbin, $
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1
  axis, yaxis=1, ystyle=1, ytitle="e!D1!N Star"

  cyrange=!y.crange 
  plothist,starpsfe2,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, starpsfe1, stare2_corr2, snperbin, $
               ytitle="e!D2!N Star",xtitle="e!D1!N PSF",$
               xrange=xrange,yrange=yrange,ystyle=1,xstyle=1

  multiplot
  plot_fitlin, starpsfe2, stare2_corr2, snperbin, $
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
                   'gale1', gale1,$
                   'gale2', gale2,$
                   'corr',galcorr,$
                   'gale1_corrected', gale1_corr,$
                   'gale2_corrected', gale2_corr,$
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

