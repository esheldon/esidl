PRO regressgal, run1, run2, clr, rmin, rmax, binsize, $
                scat=scat, lcat=lcat, step=step, addstr=addstr, $
                indir=indir, outdir=outdir, $
                random=random,nrand=nrand,title=title, wgood=wgood, $
                check=check, FIRSTRA=FIRSTRA, LASTRA=LASTRA, $
                datfile=datfile, psfile=psfile

  IF n_params() EQ 0 THEN BEGIN
      print,'-Syntax: regressgal, run1, run2, clr, rmin, rmax, binsize,'
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

  vint = .32^2

  ;; Size by which to group things  300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 16000L
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = ''
  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss4/data1/esheldon/REGRESS/'
  IF NOT keyword_set(random) AND (n_elements(nrand) EQ 0) THEN random = 0 $
  ELSE random = 1
  IF random THEN rep = 2 ELSE rep = 1
  IF NOT keyword_set(check) THEN check=0

  r1str = ntostr(run1)
  r2str = ntostr(run2)

  d2r = !dpi/180.0d0            ;Change degrees to radians.
  radmin = rmin/3600.             ;degrees
  radmax = rmax/3600.             ;degrees
  bsize=binsize/3600.         ;degrees

  nbin = long( (radmax - radmin)/bsize ) + 1 ;Histogram always rounds up

  print
  print,'Using ',ntostr(nbin),' bins between ',ntostr(radmin*3600.), $
        ' and ',ntostr(radmax*3600.),' arcseconds'

  dfac = 206264.8062471

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; declare some arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  meanr = fltarr(nbin)

  rsum = fltarr(nbin)
  wsum = fltarr(nbin)
  npsum = fltarr(nbin)
  Sshsum = 0.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output postscript file and datafile
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  endstr = r1str+'_'+r2str+'_'+colors[clr]+'_N1.txt'
  matrix_file= outdir+addstr+'matrix_'+endstr
  evec_file =  outdir+addstr+'evec_'+endstr
  rsum_file =  outdir+addstr+'rsum_'+endstr

  WHILE exist(matrix_file) OR exist(evec_file) DO BEGIN 
      matrix_file=newname(matrix_file)
      evec_file = newname(evec_file)
      rsum_file = newname(rsum_file)
  ENDWHILE

  print
  print,'Matrix file: ',matrix_file
  print,'e vector file: ',evec_file
  print,'rSum file: ',rsum_file
  print

  openw, matlun, matrix_file, /get_lun
  openw, eveclun, evec_file, /get_lun

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Get the files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nnm1 = 'srcgal'
  nnm2 = 'lensgal'


  IF n_elements(indir) EQ 0 THEN dir = '/sdss4/data1/esheldon/CORRECTED/' $
  ELSE dir = indir
  sname=dir+'run'+r1str+'_'+r2str+'_'+nnm1+'_'+colors[clr]+'_overlap.fit'
  lname=dir+'run'+r1str+'_'+r2str+'_'+nnm2+'_'+colors[clr]+'_overlap.fit'
  IF n_elements(scat) EQ 0 THEN BEGIN 
      IF NOT exist(sname) THEN BEGIN 
          sname=dir+'run'+r2str+'_'+r1str+'_'+nnm1+'_'+$
            colors[clr]+'_overlap.fit'
          IF NOT exist(sname) THEN BEGIN
              print,'No source overlap file found for these two runs'
              return
          ENDIF 
      ENDIF 
      scat = mrdfits(sname, 1, hdr1)
      LASTRA=sxpar(hdr1, 'LASTRA')
      FIRSTRA=sxpar(hdr1, 'FIRSTRA')
  ENDIF 
  IF n_elements(LASTRA) EQ 0 THEN BEGIN
      print,'Must give LASTRA if file not read in'
      return
  ENDIF 

  IF n_elements(lcat) EQ 0 THEN BEGIN 
      IF NOT exist(lname) THEN BEGIN 
          lname=dir+'run'+r2str+'_'+r1str+'_'+nnm2+'_'+$
            colors[clr]+'_overlap.fit'
          IF NOT exist(lname) THEN BEGIN
              print,'No lens overlap file found for these two runs'
              return
          ENDIF 
      ENDIF 
      lcat = mrdfits(lname, 1, hdr2)
  ENDIF 

  ;; check if runs cross ra=0.
  IF FIRSTRA GT LASTRA THEN flag=1 ELSE flag=0

  ninit = n_elements(scat)

  maxdec = max(scat.dec )
  mindec = min(scat.dec )

  ;; NEW2   note still treating as x-y coord. system
  wsource = where( (maxdec - scat.dec GE radmax) AND $
                   (scat.dec - mindec GE radmax), nsource)
  IF flag THEN BEGIN
      rmdiff1 = 360. - FIRSTRA
      IF rmdiff1 LT radmax THEN BEGIN
          ;; will still only use rectangle with ra > 0  Reset flag
          flag=0
          rmd = radmax-rmdiff1
          wsourcetmp = where( (scat[wsource].ra GE rmd) AND $
                            (LASTRA - scat[wsource].ra GE radmax), nsource)
      ENDIF ELSE IF LASTRA LT radmax THEN BEGIN
          ;; will still only use rectangle with ra > 0 Reset flag
          flag=0
          rmdiff2 = radmax-LASTRA
          wsourcetmp = where( (scat[wsource].ra - FIRSTRA GE radmax) AND $
                            (scat[wsource].ra LE (360.-rmdiff2)), nsource)
      ENDIF ELSE BEGIN 
          wsourcetmp = where( (scat[wsource].ra - FIRSTRA GE radmax) OR $
                            (LASTRA - scat[wsource].ra GE radmax), nsource)
      ENDELSE 
  ENDIF ELSE BEGIN 
      wsourcetmp = where( (scat[wsource].ra - FIRSTRA GE radmax) AND $
                          (LASTRA - scat[wsource].ra GE radmax), nsource)
  ENDELSE 
  
  wsource=wsource[wsourcetmp]

  print
  print,'LASTRA = ',LASTRA
  print,'FIRSTRA = ',FIRSTRA
  print,'Finally FLAG = ',FLAG
  print,'Using '+ntostr(nsource)+'/'+ntostr(ninit)+' sources'
  print

  nlens = n_elements(lcat)
  IF nlens LT 100 THEN sym = 1 ELSE sym = 3
;  wlensold = wlens

  IF check THEN BEGIN           ;May just want to find the good lenses.
      wgood = wlens
      return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop if random points are to be output as well
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Use lenscat instead of lcat from here on so we can change it.

  ;; New
  lra = lcat.ra
  ldec = lcat.dec
  sra = scat.ra
  sdec = scat.dec

  nsource = n_elements(scat)

  nstepOld = nsource/step
  nstep=nstepOld
  left = nsource MOD step

  ;; To account for leftover stuff
  IF nstepOld EQ 0 THEN BEGIN
      nstepOld = -10
      step = left
      nstep = 1
  ENDIF ELSE BEGIN
      IF left NE 0 THEN nstep = nstepOld + 1
  ENDELSE 

  sourcenum = 0L
  wsource = lindgen(nsource)    ;don't throw any away for now
  indices = lindgen(nsource)
  FOR group = 0L, nstep-1 DO BEGIN
      IF group EQ nstepOld THEN BEGIN
          ind = indices[ group*step: group*step+left-1  ]
          ii = wsource[ind]
          step = left
      ENDIF ELSE BEGIN
          ind = indices[ group*step : (group+1)*step -1 ]
          ii = wsource[ind]
      ENDELSE 

      ;; they are sorted by ra, but it could go over ra=0
      maxii = max(ii) & minii = min(ii)
      tmpmaxra = sra[maxii]+radmax
      tmpminra = sra[minii]-radmax

      diff1 = 360. - tmpmaxra
      
      CASE 1 OF 
          diff1 LT 0.: wlns = where(lra GE tmpminra OR $
                                    lra LE abs(diff1), nwlns )
          tmpminra LT 0: wlns = where(lra GE (360.-abs(tmpminra)) OR $
                                      lra LE tmpmaxra, nwlns)
          tmpmaxra LT tmpminra: wlns = where(lra GE tmpminra OR $
                                             lra LE tmpmaxra, nwlns)
          ELSE: wlns = where(lra LE tmpmaxra AND $
                             lra GE tmpminra, nwlns)
      ENDCASE 

      print,nwlns
      IF nwlns NE 0 THEN BEGIN 

          print,'Source Group = ',ntostr(group+1)+'/'+ntostr(nstep)
          FOR gi=0L, step-1 DO BEGIN ;Loop over gal in group
              IF (gi MOD 1000) EQ 0 THEN  print,'.',format='(a,$)'
              index = ii[gi]
              cendec = sdec[index]
              cenra  = sra[index]
 
              radiff = (cenra-lra[wlns])*d2r
              cosradiff = cos(radiff)

              tcendec = cendec*d2r
              tcenra = cenra*d2r
              sincendec=sin(tcendec)
              coscendec=cos(tcendec)

              sinlcatdec = sin(ldec[wlns]*d2r)
              coslcatdec = cos(ldec[wlns]*d2r)

              ;; distance in degrees
              R=acos( sincendec*sinlcatdec + $
                      coscendec*coslcatdec*cosradiff)/d2r
                  

              hist=histogram(R, binsize=bsize, min=radmin, max=radmax,$
                             rever=rev_ind)
              whist = where(hist NE 0, nhist)

              npair = lonarr(nbin)

              IF nhist NE 0 THEN BEGIN 
                  FOR i=0L, nhist-1 DO BEGIN 
                      
                      binnum = whist[i]
                      w=rev_ind( rev_ind(binnum):rev_ind(binnum+1)-1 )
                      npair[binnum] = n_elements(w)

                      theta2=2.*(atan(sin(radiff[w]), $
                                      (sincendec*cosradiff[w] - $
                                       coscendec*sinlcatdec[w]/coslcatdec[w]))$
                                 -!dpi)

                      npsum[binnum] = npsum[binnum] + npair[binnum]
                      rsum[binnum] = rsum[binnum] + total(R[w])
                      aij = total(-cos(theta2))
                      bij = total(-sin(theta2))

                      printf, matlun, sourcenum, binnum, aij, bij, $
                              format='( 2(I0,:,1X), 2(F0,:,1X) )'


                  ENDFOR 
                  printf, eveclun, sourcenum, scat[index].e1, scat[index].e2,$
                          scat[index].uncert, format='( (I0,:,1X),3(F0,:,1X) )'
                  sourcenum = sourcenum + 1L
              ENDIF 

              totalnpair = total(npair)
              IF totalnpair EQ 0 THEN BEGIN ;If Five
                  ;print,'Two-',ntostr(gi)
                  indices[ ind[gi] ] = -1
              ENDIF 
                  
              radiff=0          ;Free memory
              R = 0
          ENDFOR 
          print
      ENDIF ELSE BEGIN 
          print,'One-'
          indices[ ind ] = -1
      ENDELSE 
      mradiff = 0               ;Free memory
      mdist = 0
  ENDFOR 

  free_lun, matlun
  free_lun, eveclun

  ;; Remove unused sources
  wbad = where(indices EQ -1, nwbad)
  IF nwbad NE 0 THEN remove, wbad, wsource
  sourceused = n_elements(wsource)

  wgood = wsource

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find averages from sums
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Finally used ',ntostr(sourceused),'/',ntostr(nsource),' Sources'

  FOR i=0L, nbin-1 DO BEGIN
      meanr[i]   = rsum[i]/npsum[i]
  ENDFOR 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Print the data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  openw, rsumlun, rsum_file, /get_lun

  !textunit = rsumlun
  message1 = '#        meanr          npair        '
      
  printf, rsumlun, 'Nsouces:   ',sourceused
  printf, rsumlun, 'TotPairs:   ',total(npsum)
  printf, rsumlun, 'binwidth(arcsec): ',binsize
  printf, rsumlun
  printf, rsumlun, message1
  printf, rsumlun
  fmt='(3(F0,:,1X))'                    ;15=#of spaces to use(including decimal. ]
                                ;8=#decimal places
  forprint, meanr*3600., rsum, npsum, $
    TEXT=5, $
    F=fmt, $
    /silent

  free_lun, rsumlun

  ptime, systime(1)-time
  return
END 
