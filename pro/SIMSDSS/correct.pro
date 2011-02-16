PRO correct, cat, e1, e2, err, $
             dilution=dilution,$
             inputpsf=inputpsf,$
             psfe1=psfe1,psfe2=psfe2, $
             e1mom=e1mom, e2mom=e2mom, $
             olde1mom=olde1mom, olde2mom=olde2mom, $
             wg=wg, ws=ws,$
             norho=norho,$
             plt=plt, drplot=drplot, psfile=psfile, silent=silent


  IF n_params() EQ 0 THEN BEGIN
    print,'-Syntax: correct, cat, e1, e2, err,'
    print,'         dilution=dilition,'
    print,'         inputpsf=inputpsf,'
    print,'         psfe1=psfe1, psfe2=psfe2,'
    print,'         e1mom=e1mom, e2mom=e2mom,'
    print,'         wg=wg, ws=ws,'
    print,'         plt=plt, silent=silent'
    print
    print," use /dilution to apply 1/(1-R)"
    return
  ENDIF 

  IF NOT keyword_set(silent) THEN silent = 0
  pold = !p.multi
  !p.multi = [0,1,2]
  binsize=0.05

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Test the galaxies
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  n=n_elements(cat)
  wg = where(cat.goodflag EQ 1 AND $
             cat.stype EQ 1 AND $
             cat.mag_best LT 22.5 AND $
             cat.sr LT 0.8 AND $        ;;; Play with these cuts
             cat.uncert_ad LT 1.0 AND $
             cat.uncert_ad GT .01 $
             ,ngood)     
  ws = where(cat.goodflag EQ 1 AND $
             cat.stype EQ 2 AND $
             cat.mag_best LT 20.0 AND $
             cat.uncert_ad LT 1.0 $  
             ,ngood2)
  wtot = [wg,ws]

  IF ngood NE 0 AND ngood2 NE 0 THEN BEGIN

    IF NOT silent THEN BEGIN
      print,'Total objects: ',strtrim(string(n),2)
      print,'Using ',strtrim(string(ngood),2),' Good Galaxies'
      print,'Using ',strtrim(string(ngood2),2),' Good stars'
      print
    ENDIF 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; Use R value calculated from knowing which things were stars.
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    e1 = cat.e1_ad
    e2 = cat.e2_ad
    err = cat.uncert_ad

    ;; errors in original ellipticities
    serr = err[ws]  ;; star errors
    gerr = err[wg]

    IF NOT keyword_set(inputpsf) THEN BEGIN
      wmom,e1[ws], serr, psfe1
      wmom,e2[ws], serr, psfe2
      extra = '   Measured PSF'
    ENDIF ELSE BEGIN 
      psfe1 = median(cat[ws].se1)
      psfe2 = median(cat[ws].se2)
      extra='   Input PSF'
    ENDELSE 

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; Corrected e1 and e2 for anisotropy and dilution
    ;; Correct uncertainty for dilution
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    IF keyword_set(norho) THEN BEGIN 
      m=median(cat[ws].x2_ad + cat[ws].y2_ad)
      R=m/(cat[wg].x2_ad + cat[wg].y2_ad)
    ENDIF ELSE BEGIN
      R=cat[wg].sr
    ENDELSE 

    IF NOT keyword_set(dilution) THEN BEGIN 
        e1[wg] = ( e1[wg] - R*psfe1 )
        e2[wg] = ( e2[wg] - R*psfe2 )
    ENDIF ELSE BEGIN 
        e1[wg] = ( e1[wg] - R*psfe1 )/(1-R)
        e2[wg] = ( e2[wg] - R*psfe2 )/(1-R)
        err[wg] = err[wg]/(1-R)
    ENDELSE 

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Moments of corrected measurements.
    ; We have assumed that uncertainty is only changed by dilution
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    mede1corr = median(e1[wg])
    mede2corr = median(e2[wg])

    se1=sqrt( (moment(e1[wg]) )[1] )
    se2=sqrt( (moment(e2[wg]) )[1] )

    wmom, e1[wg], sqrt(se1^2+err[wg]^2), me1corr, sige1corr, we1err
    wmom, e2[wg], sqrt(se2^2+err[wg]^2), me2corr, sige2corr, we2err
    e1mom = [me1corr, we1err, mede1corr]
    e2mom = [me2corr, we2err, mede2corr]

    e1diffcorr = e1[wg] - cat[wg].se1
    e2diffcorr = e2[wg] - cat[wg].se2

    mede1diffcorr = median(e1diffcorr)
    mede2diffcorr = median(e2diffcorr)

    se1=0. && se2=0.
    se1=sqrt( (moment(e1diffcorr[wg]))[1] )
    se2=sqrt( (moment(e2diffcorr[wg]))[1] )

    wmom,e1diffcorr,err[wg],me1diffcorr,sig1,sige1diffcorr
    wmom,e2diffcorr,err[wg],me2diffcorr,sig2,sige2diffcorr

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Moments of original measurements.
    ; Add widths and uncertainty in quadrature
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    mede1 = median(cat[wg].e1_ad)
    mede2 = median(cat[wg].e2_ad)

    se1=0. && se2=0.
    se1=sqrt( (moment(cat[wg].e1_ad))[1] )
    se2=sqrt( (moment(cat[wg].e2_ad))[1] )
  
    wmom, cat[wg].e1_ad, sqrt(se1^2+gerr^2), me1, sige1, we1err
    wmom, cat[wg].e2_ad, sqrt(se2^2+gerr^2), me2, sige2, we2err
    olde1mom = [me1, we1err, mede1]
    olde2mom = [me2, we2err, mede2]
    
    e1diff = cat[wg].e1_ad - cat[wg].se1
    e2diff = cat[wg].e2_ad - cat[wg].se2
    
    ;;; no uncertainty in se1
    mede1diff = median(e1diff)
    mede2diff = median(e2diff)

    se1=0. && se2=0.
    se1=sqrt( (moment(e1diff))[1] )
    se2=sqrt( (moment(e2diff))[1] )

    wmom, e1diff, sqrt(se1^2+gerr^2), me1diff, sige1diff
    wmom, e2diff, sqrt(se2^2+gerr^2), me2diff, sige2diff

    IF NOT silent THEN BEGIN 
      print,'e1 PSF: ',psfe1
      print,'e2 PSF: ',psfe2
      print
      print,'Original mean e1:  ', me1
      print,'Median Orig e1: ', mede1
      print,'Original sigma e1: ',sige1
      print,'Original mean e2:  ', me2
      print,'Median Orig e2: ',mede2
      print,'Original sigma e2: ',sige2
      print
      print,'Corrected mean e1: ', me1corr
      print,'Median Corr e1: ', mede1corr
      print,'Corrected sigma e1: ', sige1corr
      print,'Corrected mean e2: ', me2corr
      print,'Median Corr e2: ',mede2corr
      print,'Corrected sigma e1: ', sige2corr
    ENDIF 

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; done with calculations
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Now do some plots if the user wants them
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    IF keyword_set(plt) THEN BEGIN 
      ;; to get peaks
      plothist, e1[wg], xhist, yhist, bin=binsize,/noplot
      peak = max(yhist)

      plothist, e1diffcorr,xhist,yhist,bin=binsize,/noplot
      peakdiff = max(yhist)

      t1 = 'Mean = '
      t2 = 'Sigma = '
      t3 = 'Median = '
      ytitle='number'
      tit1='Uncorrected Galaxies'
      tit2='Corrected Galaxies'
      continue = 1
      WHILE (continue) DO BEGIN 
        continue=0

        ;; histograms of e1 and e2
        ;; Uncorrected
        plothist, cat[wg].e1_ad,bin=binsize,xtitle='e1',ytitle=ytitle,$
          title = tit1,xrange=[-1,1],peak=peak
        xyouts,-.9,.9*peak,t1+strtrim(string(me1),2)
        xyouts, -.9,.8*peak, t2+strtrim(string(sige1),2)
        xyouts, -.9, .7*peak, t3+strtrim(string(mede1),2)

        plothist, cat[wg].e2_ad,bin=binsize,xtitle='e2',ytitle=ytitle,$
          title = tit1,xrange=[-1,1],peak=peak
        xyouts,-.9,.9*peak,t1+strtrim(string(me2),2)
        xyouts, -.9,.8*peak, t2+strtrim(string(sige2),2)
        xyouts, -.9, .7*peak, t3+strtrim(string(mede2),2)

        key=get_kbrd(20)
        IF key EQ 'q' THEN BEGIN
          !p.multi=pold & return
        ENDIF

        ;; Corrected
        plothist, e1[wg],bin=binsize,xtitle='e1',ytitle=ytitle,$
          title = tit2+extra,xrange=[-1,1],peak=peak
        xyouts,-.9,.9*peak,t1+strtrim(string(me1corr),2)
        xyouts, -.9,.8*peak, t2+strtrim(string(sige1corr),2)
        xyouts, -.9, .7*peak, t3+strtrim(string(mede1corr),2)

        plothist, e2[wg],bin=binsize,xtitle='e2',ytitle=ytitle,$
          title = tit2+extra,xrange=[-1,1],peak=peak
        xyouts,-.9,.9*peak,t1+strtrim(string(me2corr),2)
        xyouts, -.9,.8*peak, t2+strtrim(string(sige2corr),2)
        xyouts, -.9, .7*peak, t3+strtrim(string(mede2corr),2)

        key=get_kbrd(20)
        IF key EQ 'q' THEN BEGIN
          !p.multi=pold & return
        ENDIF 
        IF key EQ 'p' THEN continue = 1
      ENDWHILE

      IF keyword_set(psfile) THEN n=2 ELSE n=1
      FOR i=0, n-1 DO BEGIN
        IF i EQ 1 THEN makeps,psfile

        continue = 1
        WHILE (continue) DO BEGIN 
          continue=0
          ;; Histograms if emeas - e
          ;; Uncorrected
          plothist, e1diff,bin=binsize,xtitle='e1meas - e1',ytitle=ytitle,$
            title = tit1,xrange=[-1,1],peak=peak
          xyouts, -.9,.9*peak, t1+strtrim(string(me1diff),2)
          xyouts, -.9,.8*peak, t2+strtrim(string(sige1diff),2)
          xyouts, -.9, .7*peak, t3+strtrim(string(mede1diff),2)

          plothist, e2diff,bin=binsize,xtitle='e2meas - e2',ytitle=ytitle,$
            title = tit1,xrange=[-1,1],peak=peak
          xyouts, -.9,.9*peak, t1+strtrim(string(me2diff),2)
          xyouts, -.9,.8*peak, t2+strtrim(string(sige2diff),2)
          xyouts, -.9, .7*peak, t3+strtrim(string(mede2diff),2)
        
          key=get_kbrd(20)
          IF key EQ 'q' THEN BEGIN
            !p.multi=pold & return
          ENDIF

          ;; Corrected
          plothist, e1diffcorr,bin=binsize,xtitle='e1meas - e1',ytitle=ytitle,$
            title = tit2+extra,xrange=[-1,1],peak=peak
          xyouts, -.9,.9*peak, t1+strtrim(string(me1diffcorr),2)
          xyouts, -.9,.8*peak, t2+strtrim(string(sige1diffcorr),2)
          xyouts, -.9, .7*peak, t3+strtrim(string(mede1diffcorr),2)

          plothist, e2diffcorr,bin=binsize,xtitle='e2meas - e2',ytitle=ytitle,$
            title = tit2+extra,xrange=[-1,1],peak=peak
          xyouts, -.9,.9*peak, t1+strtrim(string(me2diffcorr),2)
          xyouts, -.9,.8*peak, t2+strtrim(string(sige2diffcorr),2)
          xyouts, -.9, .7*peak, t3+strtrim(string(mede2diffcorr),2)


          key=get_kbrd(20)
          IF key EQ 'q' THEN BEGIN
            !p.multi=pold & return
          ENDIF 
          IF key EQ 'p' THEN continue = 1
        ENDWHILE 
        IF i EQ 1 THEN ep
      ENDFOR 
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; End of plotting stuff
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ENDIF 

    fitlin, cat[wg].sr, e1[wg]-cat[wg].se1, err[wg],a,siga,b,sigb,/silent
    IF keyword_set(drplot) THEN BEGIN
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Make plot of e1corr - e1real vs r
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      x1=0
      x2=1
      yfac1=.9
      yfac2=.7
      yrange=[-1,1]
      xrange=[0,1]
      xtitle='R (polarizability)'

      ytitle='e1 corrected  -  e1 input'
      fitlin, cat[wg].sr, e1[wg]-cat[wg].se1, err[wg],a,siga,b,sigb,/silent
      plot,cat[wg].sr,e1[wg]-cat[wg].se1,psym=7,$
        yrange=yrange, xrange=xrange,$
        xtitle=xtitle,ytitle=ytitle,xstyle=1,ystyle=1
      oploterr, cat[wg].sr,e1[wg]-cat[wg].se1,err[wg]
      oplot,[x1,x2],[0,0]
      oplot,[x1,x2],[a+b*x1,a+b*x2]
      out1='a = '+ntostr(a,6)+' b = '+ntostr(b,6)+$
        ' sig(b) = '+ntostr(sigb,6)
      out2='S/N (a) = '+ntostr(abs(a/siga),6)+$
        ' S/N (b) = '+ntostr(abs(b/sigb),6)
      xyouts,xrange[0]+.1,  yfac1*yrange[1],out1
      xyouts,xrange[0]+.1, yfac2*yrange[1],out2

      ytitle='e2 corrected  -  e2 input'
      fitlin, cat[wg].sr, e2[wg]-cat[wg].se2, err[wg],a,siga,b,sigb,/silent
      plot,cat[wg].sr,e1[wg]-cat[wg].se1,psym=7,$
        yrange=yrange, xrange=xrange,$
        xtitle=xtitle,ytitle=ytitle,xstyle=1,ystyle=1
      oploterr, cat[wg].sr,e1[wg]-cat[wg].se1,err[wg]
      oplot,[x1,x2],[0,0]
      oplot,[x1,x2],[a+b*x1,a+b*x2]
      out1='a = '+ntostr(a,6)+' b = '+ntostr(b,6)+$
        ' sig(b) = '+ntostr(sigb,6)
      out2='S/N (a) = '+ntostr(abs(a/siga),6)+$
        ' S/N (b) = '+ntostr(abs(b/sigb),6)
      xyouts,xrange[0]+.1, yfac1*yrange[1],out1
      xyouts,xrange[0]+.1, yfac2*yrange[1],out2

    ENDIF

  ENDIF ELSE print,'No good objects found.'
  !p.multi = pold
  return


END




