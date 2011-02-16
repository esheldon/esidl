PRO objshear, lcat, scat, rmin, rmax, binsize, $
              run1=run1, run2=run2, clr=clr,$
              step=step, addstr=addstr, $
              outdir=outdir, $
              title=title, wgood=wgood, $
              check=check,  $
              datfile=datfile, psfile=psfile, zcheck=zcheck, $
              maxe=maxe

  IF n_params() LT 5 THEN BEGIN
      print,'-Syntax: objshear_fix, run1, run2, clr, rmin, rmax, binsize,'
      print,'   scat=scat, lcat=lcat, step=step, addstr=addstr, outdir=outdir,'
      print,'   random=random, title=title, wgood=wgood, check=check, '
      print,'   datfile=datfile, psfile=psfile'
      return
  ENDIF 

  colors = ['u','g','r','i','z']

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  COMMON seed,seed

  time = systime(1)
  IF n_elements(maxe) EQ 0 THEN maxe = .2
  vint = .32^2
 

  ;; Size by which to group things  300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = 'gal_gal'
  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss4/data1/esheldon/GAL_GAL/'
  IF NOT keyword_set(check) THEN check=0
  IF NOT keyword_set(zcheck) THEN zcheck=0

  IF n_elements(run1) EQ 0 THEN r1str='' ELSE r1str = '_'+ntostr(run1)
  IF n_elements(run2) EQ 0 THEN r2str='' ELSE r2str = '_'+ntostr(run2)
  IF n_elements(clr) EQ 0 THEN clr_str = '' ELSE clr_str = '_'+colors[clr]

  d2r = !dpi/180.0d0            ;Change degrees to radians.
  r2d = 180./!dpi
  radmin = rmin/3600.             ;degrees
  radmax = rmax/3600.             ;degrees
  bsize=binsize/3600.         ;degrees

  nbin = long( (radmax - radmin)/bsize ) + 1

  print
  print,'Using ',ntostr(nbin),' bins between ',ntostr(radmin*3600.), $
        ' and ',ntostr(radmax*3600.),' arcseconds'

  dfac = 206264.8062471

  stags = tag_names(scat)
  
  momerr = where(stags EQ 'MOMERR',nerr)
  IF nerr EQ 0 THEN BEGIN
      momerr = where(stags EQ 'UNCERT',nerr)
      IF nerr EQ 0 THEN BEGIN 
          print
          print,'No valid moment uncertainty tag'
          print
      ENDIF 
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; declare some arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  npair = lonarr(nbin)
  shear = fltarr(nbin)
  ortho = shear
  shearerr = shear
  orthoerr = shear
  meanr = shear

  etansum    = fltarr(nbin)
  eradsum    = etansum
  etanerrsum = etansum
  eraderrsum = etansum
  wsum       = etansum
  rsum       = etansum
  npsum      = etansum
  Sshsum     = 0.               ;Scalar

  tmpetan = etansum
  tmperad = etansum
  tmpetanerr = etansum
  tmperad = etansum
  tmperaderr = etansum
  tmpwsum = etansum
  tmprsum = etansum
  tmpnpsum = etansum
  tmpSsh = 0L

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output postscript file and datafile
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  prename = outdir+addstr+r1str+r2str+clr_str
  
  psfile = prename + '_N1.ps'
  datfile = prename + '_N1.dat'
  sumfile = prename + '_sum_N1.dat'

  WHILE exist(psfile) OR exist(datfile) DO BEGIN
      psfile = newname(psfile)
      datfile = newname(datfile)
      sumfile = newname(sumfile)
  ENDWHILE 
  print
  print,'PS file: ',psfile
  print,'Dat file: ',datfile
  print,'Sum file: ',sumfile
  print

  logname = outdir+'log.txt'
  openw, lun1, logname, /get_lun

  log_mess = '      npair       tot npair       '
  printf,lun1,log_mess

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; source cat
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nsource = n_elements(scat)
  wsource = lindgen(nsource)

  FIRSTRA = scat[0].ra
  LASTRA  = scat[nsource-1].ra

  ;; check if runs cross ra=0.
  IF FIRSTRA GT LASTRA THEN flag=1 ELSE flag=0
      
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Lens cat
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ninit = n_elements(lcat)

  maxdec = max(scat.dec )
  mindec = min(scat.dec )
 
  wlens = sort(lcat.ra)
  nlens = ninit

  IF zcheck THEN BEGIN 
      lu=lcat[wlens].petrocounts[0]
      lg=lcat[wlens].petrocounts[1]
      lr=lcat[wlens].petrocounts[2]
      li=lcat[wlens].petrocounts[3]
      lz=lcat[wlens].petrocounts[4]

      lzl = .01
      lzu = .4
      print,'Cutting by photoz'
      tred=-3.39+8.65*lr-4.32*lg+0.849*lu-4.85*li-0.121*lr*lr- $
        0.022*lr*lg+0.157*lg*lg-0.582*lr*lu+0.136*lg*lu- $
        0.000593*lu*lu+0.437*lr*li-0.203*lg*li+0.399*lu*li-0.208*li*li

      wg = where(tred GE lzl AND tred LE lzu, nlens)
      IF nlens EQ 0 THEN BEGIN
          print,'No lenses left!'
          return
      ENDIF 
      wlens = wlens[wg]
  ENDIF 

  ;; NEW2   note still treating as x-y coord. system
;  wlens = where( (maxdec - lcat.dec GE radmax) AND $
;                     (lcat.dec - mindec GE radmax), nlens)
;  IF flag THEN BEGIN
;      rmdiff1 = 360. - FIRSTRA
;      IF rmdiff1 LT radmax THEN BEGIN
;          ;; will still only use rectangle with ra > 0  Reset flag
;          flag=0
;          rmd = radmax-rmdiff1
;          wlenstmp = where( (lcat[wlens].ra GE rmd) AND $
;                            (LASTRA - lcat[wlens].ra GE radmax), nlens)
;      ENDIF ELSE IF LASTRA LT radmax THEN BEGIN
;          ;; will still only use rectangle with ra > 0 Reset flag
;          flag=0
;          rmdiff2 = radmax-LASTRA
;          wlenstmp = where( (lcat[wlens].ra - FIRSTRA GE radmax) AND $
;                            (lcat[wlens].ra LE (360.-rmdiff2)), nlens)
;     ENDIF ELSE BEGIN 
;          wlenstmp = where( (lcat[wlens].ra - FIRSTRA GE radmax) OR $
;                            (LASTRA - lcat[wlens].ra GE radmax), nlens)
;      ENDELSE 
;  ENDIF ELSE BEGIN 
;      wlenstmp = where( (lcat[wlens].ra - FIRSTRA GE radmax) AND $
;                        (LASTRA - lcat[wlens].ra GE radmax), nlens)
;  ENDELSE 
;  wlens=wlens[wlenstmp]


  print
  print,'LASTRA = ',LASTRA
  print,'FIRSTRA = ',FIRSTRA
  print,'Finally FLAG = ',FLAG
  print,'Using '+ntostr(nlens)+'/'+ntostr(ninit)+' lenses'
  print
  IF nlens LT 100 THEN sym = 1 ELSE sym = 3
  wlensold = wlens

  IF check THEN BEGIN           ;May just want to find the good lenses.
      wgood = wlens
      return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Estimate shear around all of these lenses
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  lra = lcat.ra
  ldec = lcat.dec
  
  ttt='Lenses'

  nstepOld = nlens/step
  nstep=nstepOld
  left = nlens MOD step

  ;; To account for leftover stuff
  IF nstepOld EQ 0 THEN BEGIN
      nstepOld = -10
      step = left
      nstep = 1
  ENDIF ELSE BEGIN
      IF left NE 0 THEN nstep = nstepOld + 1
  ENDELSE 

  indices = lindgen(nlens)
  FOR group = 0L, nstep-1 DO BEGIN
      IF group EQ nstepOld THEN BEGIN
          ind = indices[ group*step: group*step+left-1  ]
          ii = wlens[ind]
          step = left
      ENDIF ELSE BEGIN
          ind = indices[ group*step : (group+1)*step -1 ]
          ii = wlens[ind]
      ENDELSE 

      ;; Choose sources around this lens group
      ;; they are sorted by ra, but it could go over ra=0
      maxii = wlens[ max(ind) ] & minii = wlens[ min(ind) ]
;      maxii = max(ii) & minii = min(ii)

      tmpmaxra = lra[maxii]+radmax
      tmpminra = lra[minii]-radmax
      
      diff1 = 360. - tmpmaxra
      
      ;; quick fix
      IF flag NE 0 THEN BEGIN 
          CASE 1 OF 
              diff1 LT 0.: wsrc = where(scat.ra GE tmpminra OR $
                                        scat.ra LE abs(diff1), nwsrc )
              tmpminra LT 0: wsrc = where(scat.ra GE (360.-abs(tmpminra)) OR $
                                          scat.ra LE tmpmaxra, nwsrc)
              tmpmaxra LT tmpminra: wsrc = where(scat.ra GE tmpminra OR $
                                                 scat.ra LE tmpmaxra, nwsrc)
              ELSE: wsrc = where(scat.ra LE tmpmaxra AND $
                                 scat.ra GE tmpminra, nwsrc)
          ENDCASE 
      ENDIF ELSE BEGIN 
          
          ra1 = tmpminra
          ra2 = tmpmaxra
          binary_search, scat.ra, ra1, i1
          binary_search, scat.ra, ra2, i2
          IF i1 EQ -1 THEN i1 = 0
          IF i2 EQ -1 THEN i2 = nsource-1
          wsrc = wsource[i1:i2]
          nwsrc = n_elements(wsrc)
      ENDELSE 
      
      IF nwsrc NE 0 THEN BEGIN 
          
          xrel = dblarr(nwsrc)
          yrel = xrel
          print,'Lens Group = ',ntostr(group+1)+'/'+ntostr(nstep)
          FOR gi=0L, step-1 DO BEGIN ;Loop over gal in group
              IF (gi MOD 10) EQ 0 THEN  print,'.',format='(a,$)'
              index = ii[gi]
              cendec = ldec[index]
              cenra  = lra[index]

              ;; choose sources around this lens
              tmpmaxra = cenra+radmax
              tmpminra = cenra-radmax
              diff1 = 360. - tmpmaxra
              IF flag NE 0 THEN BEGIN 
                  ;; Quick Fix
                  wsrc2=wsrc
                  nwsrc2=nwsrc
              ENDIF ELSE BEGIN 
                  CASE 1 OF 
                  diff1 LT 0.: BEGIN & ra1=tmpminra & ra2=abs(diff1) & END 
                  tmpminra LT 0: BEGIN & ra1=360.-abs(tmpminra) & ra2=tmpmaxra & END 
                  ELSE: BEGIN & ra1=tmpminra & ra2=tmpmaxra & END 
                  ENDCASE 
                  binary_search, scat[wsrc].ra, ra1, i1
                  binary_search, scat[wsrc].ra, ra2, i2
                  IF i1 EQ -1 THEN i1 = 0
                  IF i2 EQ -1 THEN i2 = nwsrc-1
                  wsrc2 = wsrc[i1:i2]
                  nwsrc2 = n_elements(wsrc2)
              ENDELSE 
              
              ;; If there are any sources left, measure the shear
              IF nwsrc2 NE 0 THEN BEGIN 
                  radiff = (cenra-scat[wsrc2].ra)*d2r
                  cosradiff = cos(radiff)
                  
                  tcendec = cendec*d2r
                  tcenra = cenra*d2r
                  sincendec=sin(tcendec)
                  coscendec=cos(tcendec)
                  
                  sinscatdec = sin(scat[wsrc2].dec*d2r)
                  cosscatdec = cos(scat[wsrc2].dec*d2r)
                  
                  ;; distance in degrees
                  R=acos( sincendec*sinscatdec + $
                          coscendec*cosscatdec*cosradiff )*r2d
                  
                  hist=histogram(R, binsize=bsize, min=radmin, $
                                 max=radmax,rever=rev_ind)
                  numbin=n_elements(hist)
                  IF numbin NE nbin THEN BEGIN
                      print,'What!!!'
                      print,numbin,nbin
                      return
                  ENDIF 
                  whist = where(hist NE 0, nhist)
                  
                  ;; Check if there are any in this annulus rmin-rmax
                  IF nhist NE 0 THEN BEGIN 
                      
                      ;;nbin is predefined. Last bin will hold few
                      npair[*] = 0L
                      FOR i=0L, nhist-1 DO BEGIN 
                          binnum = whist[i]
                          w=rev_ind( rev_ind(binnum):rev_ind(binnum+1)-1 )
                          
                          npair[binnum] = n_elements(w)
                          
                          theta=atan(sin(radiff[w]), $
                                     (sincendec*cosradiff[w] - $
                                      coscendec*sinscatdec[w]/cosscatdec[w]) )-!dpi
                          xrel[w] = R[w]*cos(theta)
                          yrel[w] = R[w]*sin(theta)
                          diffsq=xrel[w]^2 - yrel[w]^2
                          xy=xrel[w]*yrel[w]
                          
                          e1prime=-(scat[wsrc2[w]].e1*diffsq + $
                                    scat[wsrc2[w]].e2*2.*xy  )/R[w]^2
                          e2prime= (scat[wsrc2[w]].e1*2.*xy - $
                                    scat[wsrc2[w]].e2*diffsq )/R[w]^2
                          wts = 1./( vint + scat[wsrc2[w]].(momerr[0])^2)
                          
                          tmprsum[binnum] = total(R[w])
                          tmpnpsum[binnum] = npair[binnum]
                          
                          tmpetan[binnum] = total(e1prime*wts)
                          tmperad[binnum] = total(e2prime*wts)
                          tmpetanerr[binnum]=total(wts^2*e1prime^2)
                          tmperaderr[binnum]=total(wts^2*e2prime^2)
                          tmpSsh = tmpSsh+total( wts*(1.-vint*wts*e1prime^2) )
                          
                          tmpwsum[binnum] = total(wts)
                          
                          ;; free memory
                          theta = 0 & e1prime=0 & e2prime=0 & wts=0 & w=0
                          diffsq = 0 & xy=0
                      ENDFOR 
                      
                      ;; make sure it has a reasonably symmetric distribution
                      ;; of sources behind it.  (this is what phil does)
                      wg=where(xrel NE 0., ng)
                      ixx = total(xrel[wg]^2)
                      iyy = total(yrel[wg]^2)
                      ixy = total(xrel[wg]*yrel[wg])
                      dd = ixx+iyy
                      ie1 = (ixx-iyy)/dd
                      ie2 = 2.*ixy/dd
                      ie = sqrt( ie1^2 + ie2^2 )
                      mm = 3./sqrt(ng)
                      
                      ;; maxe defined above
                      IF ie LE max([mm, maxe]) THEN BEGIN 
                          npsum[whist] = npsum[whist] + tmpnpsum[whist]
                          rsum[whist] = rsum[whist] + tmprsum[whist]
                          etansum[whist] = etansum[whist] + tmpetan[whist]
                          eradsum[whist] = eradsum[whist] + tmperad[whist]
                          
                          etanerrsum[whist]=etanerrsum[whist]+tmpetanerr[whist]
                          eraderrsum[whist]=eraderrsum[whist]+tmperaderr[whist]
                          
                          Sshsum = Sshsum+tmpSsh
                          
                          wsum[whist] = wsum[whist] + tmpwsum[whist]
                      ENDIF ELSE BEGIN 
                          ng=0
                      ENDELSE 
                      ;; reinitialize the arrays
                      xrel[wg] = 0. & yrel[wg] = 0.
                      tmpetan[whist] = 0.
                      tmperad[whist] = 0.
                      tmpetanerr[whist] = 0.
                      tmperaderr[whist] = 0.
                      tmpwsum[whist] = 0.
                      tmprsum[whist] = 0.
                      tmpnpsum[whist] = 0.
                      tmpSsh = 0L
                  ENDIF 
                  ;; One final check.  If not passed, remember that this lens wasn't used
                  IF ng EQ 0 THEN BEGIN 
                      print,'|',format='(a,$)'
                      indices[ ind[gi] ] = -1
                  ENDIF 

                  printf, lun1, ng, total(npsum),format='(2(F0,:,1X))'
                  flush, lun1
                  R = 0
                  radiff=0      ;Free memory
              ENDIF ELSE BEGIN
                  print,'|',format='(a,$)'
                  indices[ ind[gi] ] = -1
              ENDELSE 
          ENDFOR 
          print
          xrel = 0 & yrel = 0
      ENDIF ELSE BEGIN 
          print,'Two-'
          indices[ ind ] = -1
      ENDELSE 
      mradiff = 0               ;Free memory
      mdist = 0
  ENDFOR 
  ;; Remove unused lenses
  wbad = where(indices EQ -1, nwbad)
  IF nwbad NE 0 THEN remove, wbad, wlens
  lensused = n_elements(wlens)
  
  wgood = wlens
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find averages from sums
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Finally used ',ntostr(lensused),'/',ntostr(nlens),' Lenses'
  
  FOR i=0L, nbin-1 DO BEGIN
      shear[i]    = etansum[i]/wsum[i]/2.
      ortho[i]    = eradsum[i]/wsum[i]/2.
      shearerr[i] = sqrt(etanerrsum[i]/wsum[i]^2)/2.
      orthoerr[i] = sqrt(eraderrsum[i]/wsum[i]^2)/2.
      meanr[i]   = rsum[i]/npsum[i]
  ENDFOR 
  Ssh = Sshsum/total(wsum)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Print the data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  makeps, psfile, /noland
  openw, lun2, datfile, /get_lun
  openw, lun3, sumfile, /get_lun
  
  !textunit = lun2
  message1 = '#        meanr          shear         shearerr'
  message1 = message1 +  '         ortho          orthoerr     npair'
  
  printf, lun2, 'Nlenses:   ',lensused
  printf, lun2, 'TotPairs:   ',total(npsum)
  printf, lun2, 'binwidth(arcsec): ',binsize
  printf, lun2, 'Ssh', Ssh
  printf, lun2
  printf, lun2, message1
  printf, lun2
  fmt='(5E, F15.1)'
  
  colprint, meanr*3600., shear, shearerr, ortho, orthoerr, npsum, $
    lun=lun2, format=fmt
  
  free_lun, lun2
  
  !textunit = lun3
  message2 = '#      rsum(deg)     etansum     etanerrsum         '
  message2 = message2 + 'eradsum       eraderrsum           '
  message2 = message2 + 'wsum         npair'
  
  printf, lun3, 'Nlenses:   ',lensused
  printf, lun3, 'TotPairs:   ',total(npsum)
  printf, lun3, 'binwidth(arcsec): ',binsize
  printf, lun3, 'Sshsum: ',Sshsum ;Later, divide by total(wsum)
  printf, lun3
  printf, lun3, message2
  printf, lun3
  fmt = '(6E, F15.1)'
  colprint, rsum, etansum, etanerrsum, eradsum, eraderrsum, wsum, npsum, $
    lun=lun3, format=fmt

  etansum[*] = 0. & eradsum[*] = 0. & etanerrsum[*] = 0. & eraderrsum[*]=0.
  wsum[*] = 0. & rsum[*] = 0. & npair[*] = 0.
  npsum[*] = 0.
  
  free_lun, lun3
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make some plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; New
  IF flag THEN xrange=[0, 360] ELSE xrange=[FIRSTRA, LASTRA]
  plot,lra[wlens], ldec[wlens], $
    xrange=xrange, yrange=[mindec,maxdec],$
    xtitle='RA',ytitle='DEC', title=ttt, $
    psym=sym
  
  
  IF nbin GT 1 THEN BEGIN 
      pold=!p.multi
      !p.multi = [0,1,2]

      tit = 'Lenses'

      wmax1 = where(shear EQ max(shear))
      wmin1 = where(shear EQ min(shear))
      max1 = shear[wmax1] + shearerr[wmax1]
      min1 = shear[wmin1] - shearerr[wmin1]

      wmax2 = where(ortho EQ max(ortho))
      wmin2 = where(ortho EQ min(ortho))
      max2 = ortho[wmax2] + orthoerr[wmax2]
      min2 = ortho[wmin2] - orthoerr[wmin2]

      yrange = prange(shear/Ssh,ortho/Ssh,shearerr,orthoerr)
  
      xt='Radius'
      yt='"Tangential Shear"'
      ploterror, meanr*3600., shear/Ssh, shearerr, $
        psym=1,yrange=yrange,xtitle=xt,ytitle=yt, title=tit
      oplot,[0,5000],[0,0]
      
      yt='"Radial Shear"'
      ploterror, meanr*3600., ortho/Ssh, orthoerr, $
        psym=1,yrange=yrange,xtitle=xt,ytitle=yt,$
        title=ntostr(lensused)+' lenses'
      oplot,[0,5000],[0,0]
      
      !p.multi=pold
      
  ENDIF 

  step=oldstep
  ep

  close, lun1
  free_lun, lun1

  ptime, systime(1)-time
  return
END 
