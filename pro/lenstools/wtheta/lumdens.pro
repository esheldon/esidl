PRO lumdens, files, rfiles, clrs, lumnum=lumnum, galtype=galtype, create=create, kaa=kaa, hkaaerr=kaaerrh, lkaaerr=kaaerrl, doplot=doplot, nodofits=nodofits, fourbin=fourbin, fixlum=fixlum, basedir=basedir, outdir=outdir, nrad=nrad, corrfunc=corrfunc
  
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: lumdens, files, rfiles, [, clrs, lumnum=lumnum, galtype=galtype, create=create, kaa=kaa, hkaaerr=kaaerrh, lkaaerr=kaaerrl, doplot=doplot, nodofits=nodofits, fourbin=fourbin, fixlum=fixlum, basedir=basedir, outdir=outdir, nrad=nrad]'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; just input rw file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; basedir is base directory to read from

  ;; Note fixlum here means multiply the luminosity 
  ;; and output
  ;; in calc_m2l_omega_lumdens it means multiply the normalizations
  
  ;; nrad: what is index of first inner bin to use?
  IF n_elements(nrad) EQ 0 THEN nrad=0

  ;; Blanton's luminosity densities
  blanton_lumdensity, [0,1,2,3,4], bld, blderr

  ;; to account for the missing luminosity below 0.1*Lstar
  lumfix = [1.314430862, 1.235023618, 1.194514446, 1.227677096, 1.220581657]

;  IF n_params() LT 1 THEN BEGIN 
;      print,'-Syntax: lumdens, wclr, lum, rlum, lumlensum, rlumlensum'
;      return
;  ENDIF 

  nfile = n_elements(files)
  mlumfile = strarr(nfile)
  mrlumfile = strarr(nfile)

  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss5/data0/lensout/mass2light/'

  IF n_elements(clrs) EQ 0 THEN clrs=[1,2,3,4]
  nclrs=n_elements(clrs)

  IF (n_elements(lumnum) NE 0) AND (n_elements(galtype) NE 0) THEN BEGIN
      message,'Do not send lumnum and galtype at the same time'
  ENDIF 

  ;; N3 is 0.1 L* stuff
  ;; N2 is magcut stuff
  ;; N1 is 0.1 L* stuff to 2 Mpc
  ;; should not do fixlum on magcut stuff since we don't know
  ;; selection 

  ;;IF keyword_set(fixlum) THEN BEGIN
  ;;    fstr = '_fixlum_'
  ;;ENDIF ELSE BEGIN
  ;;    fstr = '_'
  ;;ENDELSE 
  ;;psend = fstr+extstr+'.ps'
  ;;fitend = fstr+extstr+'.fit'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output files, same extension as input
  ;; files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;outps = outdir + 'lumdense_allplot'+psend
  ;;outfits = outdir + 'lumdense_allplot'+fitend
  ;;IF n_elements(lumnum) NE 0 THEN BEGIN
  ;;   IF NOT keyword_set(fourbin) THEN ttstr='twobin' ELSE ttstr=''
  ;;    type = 'lum'+ntostr(long(lumnum))+ttstr+'_'
  ;;    outps = outdir+'lumdense_lumnum'+ntostr(long(lumnum))+ttstr+psend
  ;;    outfits = outdir+'lumdense_lumnum'+ntostr(long(lumnum))+ttstr+fitend
  ;;ENDIF ELSE type=''
  type=''
  ;;IF n_elements(galtype) NE 0 THEN BEGIN
  ;;    IF (galtype NE 'spiral') AND (galtype NE 'ellip') AND $
  ;;      (galtype NE 'low') AND (galtype NE 'high') THEN $
  ;;      message,'galtype '+ntostr(galtype)+' unknown'
  ;;    type = galtype+'_'
  ;;    outps = outdir+'lumdense_'+galtype+psend
  ;;    outfits = outdir+'lumdense_'+galtype+fitend
  ;;ENDIF 

  dirsep, files[0], tmpdir, tmlf

  outfits = outdir+'lumdens_'+tmlf
  outps = repstr(outfits,'.fit','.ps')

  print,'Output ps file: ',outps
  print,'Output fits file: ',outfits
  print

  IF keyword_set(doplot) THEN begplot,name=outps
  setup_mystuff
  dtype = display_type()

  CASE nclrs OF
      5: BEGIN & erase & multiplot,[2,3] & END 
      4: BEGIN & erase & multiplot,[2,2], /square & END 
      3: BEGIN & erase & multiplot,[2,2] & END 
      2: BEGIN & erase & multiplot,[1,2] & END 
      ELSE:
  ENDCASE 

  kaa = fltarr(5, 2)
  kaaerrh = fltarr(5,2)
  kaaerrl = fltarr(5,2)

  IF n_elements(basedir) EQ 0 THEN BEGIN
      dirsep, files[0], basedir, throwaway
  ENDIF 

  FOR i=0,nclrs-1 DO BEGIN 
      wclr=clrs[i]
      
      
      IF ( (n_elements(lumnum) NE 0) OR $
           (n_elements(galtype) NE 0) ) THEN substr='sublum/'+!colors[wclr]+'/'$
      ELSE substr=''

      front = 'wthetalumw'
      rfront = 'wthetarandlumw'

;      front = 'wthetalumwlg'
;      rfront = 'wthetarandlumwlg'
      
      indir=basedir+substr
;indir='/sdss6/data0/wtheta/'
;IF wclr EQ 2 THEN indir = '/sdss4/data1/esheldon/tmp/'

;      indir='/sdss4/data1/esheldon/tmp/'+substr
      
      clrstr = !colors[wclr]
      
      ;; create file names for this bandpass
      FOR fi=0L, nfile-1 DO BEGIN 
          mlumfile[fi] = repstr(files[fi],'rw',clrstr+'w')
          mrlumfile[fi] = repstr(rfiles[fi],'rw',clrstr+'w')
      ENDFOR 
      combine_wthetalumw, mlumfile, lum
      combine_wthetalumw, mrlumfile, rlum

;      mlumfile=indir+type+front+'_stripe'+stripestr+'_sum_'+!colors[wclr]+'w_'+ext
;      mrlumfile=indir+type+rfront+'_stripe'+stripestr+'_sum_'+!colors[wclr]+'w_'+ext
;      lum=mrdfits(mlumfile,1,status=ok)
;      IF ok ne 0 THEN message,'file not found'
;      rlum=mrdfits(mrlumfile,1,status=ok)
;      IF ok ne 0 THEN message,'file not found'


      nbin=n_elements(lum.meanr)

      meanlum=lum.meanlum
      error=lum.meanlumerr
      tmeanlum=lum.tmeanlum
      terror=lum.tmeanlumerr

      rmeanlum=rlum.meanlum
      rerror=rlum.meanlumerr
      rtmeanlum=rlum.tmeanlum
      rterror=rlum.tmeanlumerr

      ;; last one might be zero
      w=where(lum.meanr NE 0,nw)
      wfit=w[nrad:nw-1]

      ;; add in extra luminosity below 0.1 lstar
      IF keyword_set(fixlum) THEN BEGIN 
          meanlum = meanlum*lumfix[wclr] & error=error*lumfix[wclr]
          rmeanlum = rmeanlum*lumfix[wclr] & rerror=rerror*lumfix[wclr]

          tmeanlum = tmeanlum*lumfix[wclr] & terror=terror*lumfix[wclr]
          rtmeanlum = rtmeanlum*lumfix[wclr] & rterror=rterror*lumfix[wclr]
      ENDIF 
      ;;forprint,meanlum,lum.meanlum

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; plots of mean lum around lenses, random points. Then plot
      ;; of difference
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      yrange=prange(meanlum[wfit],rmeanlum[wfit],error[wfit],rerror[wfit])
      yrange[0]=0
      xtitle=!mpcxtitle2
      ytitle = 'Mean '+!colorsp[wclr]+' lum [10!U10!N L'+sunsymbol()+']'

      lumdiff = meanlum - rmeanlum
      differror=sqrt(error^2 + rerror^2)
      tlumdiff = tmeanlum - rtmeanlum
      tdifferror=sqrt(terror^2 + rterror^2)

      command = !colors[wclr]+'lumdiff = lumdiff & ' + $
                !colors[wclr]+'differror = differror & ' + $
                !colors[wclr]+'tlumdiff = tlumdiff & ' + $
                !colors[wclr]+'tdifferror = tdifferror & '
      IF NOT execute(command) THEN message,' Doh!'

      IF nclrs EQ 1 THEN BEGIN 
          erase & multiplot,[1,3]
          ploterror,lum.meanr[wfit],meanlum[wfit],error[wfit],psym=1,yrange=yrange,$
            ytitle='around lenses'
          multiplot
          ploterror,rlum.meanr[wfit],rmeanlum[wfit],rerror[wfit],psym=1,yrange=yrange,$
            ytitle='around rand'

          multiplot
          ploterror,lum.meanr[wfit],lumdiff[wfit],differror[wfit],psym=1,$
            ytitle=ytitle,xtitle=xtitle
          multiplot,/reset
      ENDIF 
      diffdense = lumdiff
      diffdenseerror = differror

      if (dtype eq 'X') AND (nclrs EQ 1) then key=get_kbrd(1)
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; luminosity density and best fitting power law
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      yt=''
      xt=''
      ytitle = !colorsp[wclr]+!lumytitle
      IF nclrs NE 1 THEN BEGIN
          ;;!p.ticklen=0.05
          ytitle = !lumytitle
      ENDIF 
      IF nclrs EQ 1 THEN BEGIN 
          xt=xtitle 
          yt=ytitle 
      END 

      IF (nclrs NE 1) AND (i EQ 2) THEN yt=ytitle
      IF (nclrs NE 1) AND ( (i EQ 2) OR (i EQ 3) ) THEN xt=xtitle
                                ;IF wclr GT 0 THEN ycrange=!y.crange

      ymin=0.5
      minfac = 0.7
      maxfac = 1.3
      xrange=[minfac*min(lum.meanr[wfit]/1000.), $
              maxfac*max(lum.meanr[wfit])/1000.]
      xrange[0] = xrange[0] < 0.1
;      xrange[0] = xrange[0] < 0.01

      CASE type OF
          'spiral_': ycrange=[ymin,50]
          'ellip_': ycrange=[ymin,50]
          'low': ycrange=[ymin,20]
          'high': ycrange=[ymin,20]
          ELSE: BEGIN 
              IF n_elements(lumnum) EQ 0 THEN ycrange=[ymin,20] $
              ELSE BEGIN 
                  CASE lumnum OF
                      1: ycrange=[ymin,100]
                      2: ycrange=[ymin,100]
                      3: ycrange=[ymin,100]
                      4: ycrange=[ymin,110]
                  ENDCASE 
              ENDELSE 
          END 
      ENDCASE 

      myusersym, 'fill_circle', symsize=0.8
      psym=8
      IF (nclrs NE 1) AND (i GT 0) THEN multiplot, doxaxis=doxaxis
      ploterror,lum.meanr[wfit]/1000.,diffdense[wfit],diffdenseerror[wfit],$
        psym=psym,$
        ytitle=yt,xtitle=xt,/xlog,/ylog,$
        xrange=xrange,xstyle=1,yrange=ycrange, ystyle=1,$
        xticklen=0.05,yticklen=0.05, hat=0
;      ploterror,lum.meanr/1000.,diffdense,diffdenseerror,psym=psym,$
;        ytitle=yt,xtitle=xt,/xlog,/ylog,$
;        xrange=xrange,xstyle=1,yrange=ycrange, ystyle=1,$
;        xticklen=0.05,yticklen=0.05, hat=0
      IF nclrs GT 1 THEN BEGIN
          legend,!colorsp[wclr],/right,box=0
          ;;!p.ticklen=0.02
      ENDIF 
      ;IF (nclrs EQ 4) AND ( (i EQ 2) OR (i EQ 3) ) THEN add_labels,xtickv=[.2,2]
      ;IF (nclrs EQ 4) AND ( (i EQ 0) OR (i EQ 2) ) THEN add_labels,ytickv=[2,3,4,5]
      
;      ENDIF ELSE BEGIN 
;          aploterror,!gratio,lum.meanr[wfit],diffdense[wfit],diffdenseerror[wfit],psym=1,$
;            ytitle=yt,xtitle=xt
;      ENDELSE 
      
;      wfit=w


;      fitpower, lum.meanr[wfit]/1000., diffdense[wfit], diffdenseerror[wfit],$
;        [1.0, -1.0], yfit, aa
;      norm=aa[0]
;      pow=aa[1]

      npow = 400L
      nnorm = 400L
      IF n_elements(lumnum) EQ 0 THEN BEGIN 
          IF n_elements(galtype) EQ 0 THEN BEGIN 
                  powmin=-1.0 & powmax = -0.4
                  normmin=0.0 & normmax = 4.0
          ENDIF ELSE BEGIN 
              CASE galtype OF 
                  'ellip': BEGIN 
                      powmin=-1.0 & powmax = -0.4
                      normmin=1.0 & normmax = 5.
                  END 
                  'spiral': BEGIN 
                      powmin=-1.0 & powmax = -0.4
                      normmin=1.0 & normmax = 3.5
                  END 
                  ELSE : BEGIN
                      powmin=-1.0 & powmax = -0.4
                      normmin=1.0 & normmax = 5.
                  END 
              ENDCASE 
          ENDELSE 
      ENDIF ELSE BEGIN 
          CASE lumnum OF 
              1: BEGIN 
                  powmin=-1.0 & powmax = -0.4
                  normmin=1.0 & normmax = 4.0
              END 
              2: BEGIN 
                  powmin=-1.0 & powmax = -0.4
                  normmin=1.0 & normmax = 5.5
              END 
              3: BEGIN 
                  powmin=-1.5 & powmax = -0.4
                  normmin=0.0 & normmax = 6.0
              END 
              4: BEGIN 
                  IF clrs[i] EQ 0 THEN BEGIN 
                      powmin=-1.0 & powmax = 0.0
                  ENDIF ELSE BEGIN 
                      powmin=-1.5 & powmax = -0.4
                  ENDELSE 
                  normmin=1.0 & normmax = 8.0
              END 
              ELSE: message,'unknown lumnum'+ntostr(lumnum)
          ENDCASE 
      ENDELSE 
      IF clrs[i] EQ 0 THEN powmax = 0.0

      IF nclrs EQ 1 THEN BEGIN
          nodisplay=0 
          IF (dtype EQ 'X') THEN key=get_kbrd(1)
      ENDIF ELSE nodisplay=1
      pow_chisq_conf_gen, $
        lum.meanr[wfit]/1000., diffdense[wfit], diffdenseerror[wfit],$
        [powmin,powmax], [normmin,normmax], npow, nnorm, $
        chisq_surf, $
        pow, norm, powlow, powhigh, normlow, normhigh,$
        nodisplay=nodisplay, yfit=lumdens

      range2error, powlow[0], pow, powhigh[0], peh, pel
      range2error, normlow[0], norm, normhigh[0], neh, nel
      
      IF nclrs EQ 1 THEN BEGIN 
          IF (dtype EQ 'X') THEN key=get_kbrd(1)
          ploterror,lum.meanr[wfit],diffdense[wfit],diffdenseerror[wfit],psym=1,$
            ytitle=yt,xtitle=xt,/xlog,/ylog,$
            xrange=xrange,xstyle=1,yrange=ycrange, ystyle=1
      ENDIF 
      oplot, lum.meanr[wfit]/1000., lumdens

      print
      print,'Norm: ',ntostr(norm),' + ',ntostr(neh),' - ',ntostr(nel)
      print,'Pow:  ',ntostr(pow),' + ',ntostr(peh),' - ',ntostr(pel)
      print
;      mess1=['A = '+ntostr( rnd(norm,2), 5 ),$
;             !tsym.alpha +' = '+ntostr( rnd(pow,2), 5 )]
      mess1=['b = '+ntostr( rnd(norm,2), 4 )+$
             '!S!U+'+ntostr(rnd(neh,2),4)+'!R!D'+!tsym.minus+ntostr(rnd(nel,2),4), $
             '',$
             !tsym.beta +' = '+ntostr( rnd(abs(pow),3), 5 )+$
             '!S!U+'+ntostr(rnd(peh,3),5)+'!R!D'+!tsym.minus+ntostr(rnd(pel,3),5) ]

      IF wclr EQ 0 THEN legend,mess1,/top,/left, box=0, charsize=1 $
      ELSE legend,mess1,/bottom,/left, box=0, charsize=1

      IF (dtype EQ 'X') AND (i LT nclrs-1) AND (nclrs EQ 1) THEN key=get_kbrd(1)
      
      aa = [norm, pow]
;      kaa[i,*] = aa
;      kaaerrh[i,0] = normhigh[0]
;      kaaerrl[i,0] = normlow[0]
;      kaaerrh[i,1] = powhigh[0]
;      kaaerrl[i,1] = powlow[0]

      kaa[wclr,*] = aa
      kaaerrh[wclr,0] = normhigh[0]
      kaaerrl[wclr,0] = normlow[0]
      kaaerrh[wclr,1] = powhigh[0]
      kaaerrl[wclr,1] = powlow[0]

  ENDFOR 
  IF nclrs NE 1 THEN multiplot,/reset
;stop

  ;; output fit parameters

  IF (nclrs EQ 4) THEN BEGIN 

      ss=create_struct('norm',fltarr(5),$
                       'pow',fltarr(5),$
                       'normhigh',fltarr(5),$
                       'normlow',fltarr(5),$
                       'powhigh',fltarr(5),$
                       'powlow',fltarr(5), $
                       'wfit', wfit,$
                       'meanr', lum.meanr,$
                       'glumdense', glumdiff,$
                       'glumdenserr',gdifferror,$
                       'rlumdense', rlumdiff,$
                       'rlumdenserr',rdifferror,$
                       'ilumdense', ilumdiff,$
                       'ilumdenserr',idifferror,$
                       'zlumdense', zlumdiff,$
                       'zlumdenserr',zdifferror,$
                       'rmin', lum.rmin,$
                       'rmax', lum.rmax_act,$
                       'gtlumdense', rtlumdiff,$
                       'gtlumdenserr', rtdifferror,$
                       'rtlumdense', rtlumdiff,$
                       'rtlumdenserr', rtdifferror,$
                       'itlumdense', itlumdiff,$
                       'itlumdenserr', itdifferror,$
                       'ztlumdense', ztlumdiff,$
                       'ztlumdenserr', ztdifferror)

      ss.norm[*] = kaa[*,0]
      ss.pow[*] = kaa[*,1]
      ss.normhigh[*] = kaaerrh[*,0]
      ss.normlow[*] = kaaerrl[*,0]
      ss.powhigh[*] = kaaerrh[*,1]
      ss.powlow[*] = kaaerrl[*,1]

      IF NOT keyword_set(nodofits) THEN BEGIN 
          print
          print,'Outputting fits file: ',outfits
          print
          mwrfits, ss, outfits, /create
      ENDIF 

  ENDIF 


;; skip mass stuff except for total sample

  GOTO,jump
  IF type NE '' OR nclrs NE 5 THEN GOTO,jump

  if dtype eq 'X' then key=get_kbrd(1)

  xtitle='Aperture Radius [h!U'+!tsym.minus+'1!N kpc]'
  meanlum = [0.784103, 0.889751, 1.51780, 2.07778, 2.58957]

  mdir='/sdss5/data0/lensout/stripe10/'
  mstr=mrdfits(mdir+'matchrN5_main_zgal_gal_stripe10_r_corr_N1.fit',1)
  w=where(mstr.meanr NE 0.0)

;  minxx = min(xx)
  minxx=0.0
  maxx = 1000./1000.
  xx = arrscl( findgen(1000), minxx, maxx )

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plot integrated luminosity of neighbors
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  lines = [0,1,2,3,4]
  FOR i=4L, 1, -1 DO BEGIN 

      wclr=clrs[i]

      alpha=-kaa[i,1]
      norm=kaa[i,0]
      lumdens = norm*(xx)^(-alpha)
      lumdens_meanrmax = 2./(2.-alpha)*lumdens ;lum/Mpc^2
      ;lumtot = (lumdens_meanrmax*!pi*xx^2 + meanlum[wclr]);*1.e10
      lumtot = lumdens_meanrmax*!pi*xx^2

      ytitle='Luminosity (h!U'+!tsym.minus+'2!N L'+sunsymbol()+')'

      IF i EQ 4 THEN BEGIN 
          aplot,!gratio,xx*1000.,lumtot,line=lines[i],$
            ytitle=ytitle,xtitle=xtitle,title='Neighbors'
      ENDIF ELSE BEGIN 
          oplot,xx*1000.,lumtot,line=lines[i]
      ENDELSE 

;      oplot,xx*1000.,denscontratio,line=2
;      legend,[!tsym.delta_cap+!tsym.sigma_cap+'/'+!tsym.delta_cap+'L', $
;              'M(<R)/L(<R)'], line=[2,0],thick=[!p.thick,!p.thick]

  ENDFOR 
  legend,!colorsp,line=lines,thick=replicate(!p.thick,5)
  if dtype eq 'X' then key=get_kbrd(1)


  ;; now of neighbors + central galaxy
  lines = [0,1,2,3,4]
  FOR i=1L, 4 DO BEGIN 

      wclr=clrs[i]

      alpha=-kaa[i,1]
      norm=kaa[i,0]
      lumdens = norm*(xx)^(-alpha)
      lumdens_meanrmax = 2./(2.-alpha)*lumdens ;lum/Mpc^2
      ;lumtot = (lumdens_meanrmax*!pi*xx^2 + meanlum[wclr]);*1.e10
      lumtot = (lumdens_meanrmax*!pi*xx^2 + meanlum[wclr])/meanlum[wclr]
      ;lumtot = lumdens_meanrmax*!pi*xx^2

      ytitle='Luminosity (h!U'+!tsym.minus+'2!N L'+sunsymbol()+')'

      IF i EQ 0 THEN BEGIN 
          aplot,!gratio,xx*1000.,lumtot,line=lines[i],$
            ytitle=ytitle,xtitle=xtitle,title='L/L!DCent!N';,title='Central Galaxy + neighbors'
      ENDIF ELSE BEGIN 
          oplot,xx*1000.,lumtot,line=lines[i]
      ENDELSE 

;      oplot,xx*1000.,denscontratio,line=2
;      legend,[!tsym.delta_cap+!tsym.sigma_cap+'/'+!tsym.delta_cap+'L', $
;              'M(<R)/L(<R)'], line=[2,0],thick=[!p.thick,!p.thick]

  ENDFOR 
  legend,!colorsp,line=lines,thick=replicate(!p.thick,5)
  if dtype eq 'X' then key=get_kbrd(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now mass-to-light ratio using density contrast
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  lines=[0,1,2,3,4]
  FOR i=0L, nclrs-1 DO BEGIN 

      wclr=clrs[i]

      alpha=-kaa[i,1]
      norm=kaa[i,0]
      userad = mstr.rmax_act[wfit]/1000.
      lumdens = norm*(userad)^(-alpha)

      ;; density contrast
      lumdenscont = (lumdens*alpha/(2.-alpha) + meanlum[wclr]/!pi/userad^2)*1.e10 ;lum/Mpc^2
      ;; mean density contrast within radius
      lumdenscont_meanrmax = 2./(2.-alpha)*lumdenscont

      rratio=mstr.rmin/mstr.rmax_act[wfit]
      lumdenscont_rdiff = $
        lumdenscont_meanrmax*(1.-rratio^(2.-alpha))/(1. - rratio^2 )
      
      ;; convert to mass/Mpc^2
      tsigma = mstr.tsigma[wfit]*1.e12
      denscontratio = tsigma/lumdenscont_rdiff
;      IF i EQ 3 THEN print,denscontratio

      ytitle='M/L (h M'+sunsymbol()+' /L'+sunsymbol()+')'
      message=!colorsp[wclr]

      IF i EQ 0 THEN BEGIN 
          aplot,!gratio,mstr.rmax_act[wfit],denscontratio,line=lines[i],$
            ytitle=ytitle,xtitle=xtitle,title=title
      ENDIF ELSE BEGIN 
          oplot,mstr.rmax_act[wfit],denscontratio,line=lines[i]
      ENDELSE 

  ENDFOR 
  legend,!colorsp,line=lines,thick=replicate(!p.thick,5)

  if dtype eq 'X' then key=get_kbrd(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now mass-to-light ratio using model
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  yrange=!y.crange
  M0 = 2.8
  Malpha = 0.8
  denscont = M0*(xx)^(-Malpha)     ;Msolar/pc^2
  denscont = denscont*1.e12     ;Msolar/Mpc^2
  ;; convert to density
  density = denscont*(2.-Malpha)/Malpha
  ;; convert to density within radius
  mean_density = 2./(2.-Malpha)*density
  masstot = mean_density*!pi*xx^2

  lines = [0,1,2,3,4]
  FOR i=0L, nclrs-1 DO BEGIN 

      wclr=clrs[i]

      alpha=-kaa[i,1]
      norm=kaa[i,0]
      lumdens = norm*(xx)^(-alpha)
      lumdens_meanrmax = 2./(2.-alpha)*lumdens ;lum/Mpc^2
      lumtot = (lumdens_meanrmax*!pi*xx^2 + meanlum[wclr])*1.e10
      
      ratio = masstot/lumtot
      ytitle='M/L (h M'+sunsymbol()+' /L'+sunsymbol()+')'
      message=!colorsp[wclr]

      IF i EQ 0 THEN BEGIN 
          aplot,!gratio,xx*1000.,ratio,line=lines[i],$
            ytitle=ytitle,xtitle=xtitle,title=title,yrange=yrange
      ENDIF ELSE BEGIN 
          oplot,xx*1000.,ratio,line=lines[i]
      ENDELSE 

;      oplot,xx*1000.,denscontratio,line=2
;      legend,[!tsym.delta_cap+!tsym.sigma_cap+'/'+!tsym.delta_cap+'L', $
;              'M(<R)/L(<R)'], line=[2,0],thick=[!p.thick,!p.thick]

  ENDFOR 
  legend,!colorsp,line=lines,thick=replicate(!p.thick,5)

jump:

  IF keyword_set(doplot) THEN endplot

END 
