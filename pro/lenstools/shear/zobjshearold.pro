PRO zobjshearold, run1, run2, clr, rmin, rmax, binsize, $
                  scat=scat, lcat=lcat, step=step, addstr=addstr, $
                  outdir=outdir, indir=indir,$
                  random=random, title=title, wgood=wgood, check=check, $
                  datfile=datfile, psfile=psfile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    
;       
; PURPOSE:
;    
;
; CALLING SEQUENCE:
;    
;
; INPUTS: 
;    
;
; OPTIONAL INPUTS:
;    
;
; KEYWORD PARAMETERS:
;    
;       
; OUTPUTS: 
;    
;
; OPTIONAL OUTPUTS:
;    
;
; CALLED ROUTINES:
;    
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF n_params() EQ 0 THEN BEGIN
      print,'-Syntax: zobjshear, run1, run2, clr, rmin, rmax, binsize,'
      print,'   scat=scat, lcat=lcat, step=step, addstr=addstr, outdir=outdir,'
      print,'   random=random, title=title, wgood=wgood, check=check,'
      print,'   datfile=datfile, psfile=psfile'
      return
  ENDIF 

  colors = ['u','g','r','i','z']

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  COMMON seed,seed
  time = systime(1)
  ;; Size by which to group things  300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = 'gal_gal_'
  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss4/data1/esheldon/GAL_GAL/'
  IF n_elements(indir) EQ 0 THEN indir = '/sdss4/data1/esheldon/CORRECTED/'
  IF NOT keyword_set(random) THEN random = 0
  IF random THEN rep = 2 ELSE rep = 1
  IF NOT keyword_set(check) THEN check=0

  r1str = ntostr(run1)
  r2str = ntostr(run2)

  nbin = long( (rmax - rmin)/binsize )

  print
  print,'Using ',ntostr(nbin),' bins between ',ntostr(rmin), $
        ' and ',ntostr(rmax),' kpc'

  zsource = .4


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; declare some arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;  npair = lonarr(nbin)
  shear = fltarr(nbin)
  ortho = shear
  shearerr = shear
  orthoerr = shear

  sigma = shear
  orthosig = shear
  sigmaerr = shear
  orthosigerr = shear

  meanr = shear

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output postscript file and datafile
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  prename = outdir+addstr+r1str+'_'+r2str+'_'+colors[clr]
  
  psfile = prename + '_N1.ps'
  datfile = prename + '_N1.dat'
  WHILE exist(psfile) OR exist(datfile) DO BEGIN
      psfile = newname(psfile)
      datfile = newname(datfile)
  ENDWHILE 
  print
  print,'PS file: ',psfile
  print,'Dat file: ',datfile
  print

  logname = outdir+'log.txt'
  openw, lun1, logname, /get_lun

  log_mess = '      npair       tot npair       '
  printf,lun1,log_mess

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Get the files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; For now, use all source galaxies, but only use lenses with good photoz
  nnm1 = 'srcgal'
  nnm2 = 'zlensgal'

  sname=indir+'run'+r1str+'_'+r2str+'_'+nnm1+'_'+colors[clr]+'_overlap.fit'
  lname=indir+'run'+r1str+'_'+r2str+'_'+nnm2+'_'+colors[clr]+'_overlap.fit'
  IF NOT exist(sname) THEN BEGIN 
      sname=indir+'run'+r2str+'_'+r1str+'_'+nnm1+'_'+colors[clr]+'_overlap.fit'
      lname=indir+'run'+r2str+'_'+r1str+'_'+nnm2+'_'+colors[clr]+'_overlap.fit'
  ENDIF 
  IF NOT exist(sname) THEN BEGIN
      print,'No overlap file exists for these two runs'
      return
  ENDIF 
  IF n_elements(scat) EQ 0 THEN scat = mrdfits(sname, 1, hdr1)
  IF n_elements(lcat) EQ 0 THEN lcat = mrdfits(lname, 1, hdr2)
  lenscat=lcat                  ;Need copy for random

  ninit = n_elements(lcat)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up sigma crit and Dlens
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  tagnames = tag_names(lcat)
  wz = where(tagnames EQ 'Z', nz)
  IF nz EQ 0 THEN BEGIN
      wz = where(tagnames EQ 'PHOTOZ', nz)
      
      IF nz EQ 0 THEN BEGIN
          print,'Lens structure must have "Z" or "PHOTOZ" flag'
          return
      ENDIF 
  ENDIF 

  h=1.
  sigmacrit = sdss_sigma_crit( lcat.(wz[0]), zsource, wgood=wlens, h=h)

  DL = angdist( lcat.(wz[0]), h=h )*1000. ;Convert from Mpc to kpc

  nlens1 = n_elements(wlens)
  print,'Threw out ',ntostr(ninit - nlens1),' lenses because too deep or zero'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Get lenses that aren't too close to edge
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Can get maximum possible angle from each lens
  ;; angmax = rmax/DL

  angmax = rmax/DL*180./!pi

  maxra  = max( scat.ra  )
  minra  = min( scat.ra  )
  maxdec = max( scat.dec )
  mindec = min( scat.dec )

  bad = -1
  FOR i=0, nlens1-1 DO BEGIN 
      ii = wlens[i]
     IF ( ( maxra - lcat[ii].ra LE angmax[ii]) OR $
          ( lcat[ii].ra - minra LE angmax[ii]) OR $
          ( maxdec - lcat[ii].dec LE angmax[ii]) OR $
          ( lcat[ii].dec - mindec LE angmax[ii]) ) THEN BEGIN 
         IF bad[0] EQ -1 THEN bad = i ELSE bad=[bad,i]
     ENDIF 
  ENDFOR 
 
  IF bad[0] NE -1 THEN remove, bad, wlens
  wlensold = wlens
  nlens = n_elements(wlens)

  print,'Threw out ',ntostr(nlens1-nlens),' edge lenses'
  print,'Using '+ntostr(nlens)+' lenses'
  print
  IF nlens LT 100 THEN sym = 1 ELSE sym = 3

  IF check THEN BEGIN           ;May just want to find the good lenses.
      wgood = wlens
      return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop if random points are to be used as well
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Use lenscat instead of lcat from here on so we can change it.

  rand=0
  FOR plt_i=0, rep-1 DO BEGIN 


      IF plt_i EQ 1 THEN rand = 1
      IF rand THEN BEGIN

          wlens = wlensold
          step=oldstep
          print,'Using random RA and DEC points'

          ;; WARNING! THIS WON'T WORK EXCEPT IF SCAN IS ALIGNED WITH RA!
          tmaxra = max(lenscat[wlens].ra)
          tminra = min(lenscat[wlens].ra)
          tmaxdec= max(lenscat[wlens].dec)
          tmindec= min(lenscat[wlens].dec)

          ratmp  = randomu(seed, nlens)
          dectmp = randomu(seed, nlens)
          ratmp = arrscl(ratmp,  tminra,  tmaxra,  arrmin=0.,arrmax=1.)
          dectmp= arrscl(dectmp, tmindec, tmaxdec, arrmin=0.,arrmax=1.)

          lenscat[wlens].ra  = ratmp
          lenscat[wlens].dec = dectmp
          ttt = 'Random Points'

      ENDIF ELSE ttt='Lenses'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Estimate shear around all of these lenses
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


      ;; do foreground galaxies one "step" at a time to save search time
      ;; currently works best if sorted by RA (assuming equatorial scan)
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

      ntot=0.
      ;; WARNING! NEEDS TO BE FIXED FOR OFF EQUATORIAL SCANS
      indices = lindgen(nlens)
      FOR group = 0L, nstep-1 DO BEGIN
          IF group EQ nstepOld THEN BEGIN
              ind = indices[ group*step: group*step+left-1L  ]
              ii = wlens[ind]
              step = left
          ENDIF ELSE BEGIN
              ind = indices[ group*step : (group+1L)*step -1L ]
              ii = wlens[ind]
          ENDELSE 

          angmax2 = max(angmax[ii])

          ;; Choose stuff withing the maximum angular distance first.
          tmpmaxra = max( lenscat[ii].ra )
          tmpminra = min( lenscat[ii].ra )
          tmpmaxdec = max( lenscat[ii].dec )
          tmpmindec = min( lenscat[ii].dec )
          wsrc = where( scat.ra LE (tmpmaxra + angmax2) AND $
                        scat.ra GE (tmpminra - angmax2), nwsrc)
          IF nwsrc NE 0 THEN BEGIN 
              wsrc2 = where( scat[wsrc].dec LE (tmpmaxdec + angmax2) AND $
                             scat[wsrc].dec GE (tmpmindec - angmax2), nsrc)
              IF nsrc NE 0 THEN BEGIN 
                  wsrc = wsrc[ wsrc2 ]
                  print,'Lens Group = ',ntostr(group+1)+'/'+ntostr(nstep)
                  FOR gi=0L, step-1 DO BEGIN ;Loop over gal in group

                      index = ii[gi]
                      cenx = lenscat[index].dec
                      ceny = lenscat[index].ra

                      angmax_i = angmax[index]
    
                      ;;  Note: in the future, x,y not equal dec,ra!!  
                      ;;Will have to do a coordinate transformation
                      xrel = scat[wsrc].dec - cenx
                      yrel = scat[wsrc].ra  - ceny
                      ws1=where( xrel GE -angmax_i AND $
                                 xrel LE angmax_i, nw1)
                      IF nw1 NE 0 THEN BEGIN 
                          ws2=where( yrel[ws1] GE -angmax_i AND $
                                     yrel[ws1] LE angmax_i, nw2)
                          IF nw2 NE 0 THEN BEGIN 
                              ws = wsrc[ ws1[ws2] ]

                              ntot = ntot+nw2

                              ;; Convert xrel, yrel to physical distance 
                              ;;(kpc) from center of lens.
                              xrel = xrel[ ws1[ws2] ]*!pi/180.*DL[index]
                              yrel = yrel[ ws1[ws2] ]*!pi/180.*DL[index]

                              ;; IMPORTANT to get x and y right!
                              sum_sigma, sigmacrit[index], $
                                scat[ws].e1, scat[ws].e2, scat[ws].uncert, $
                                xrel, yrel, $
                                rmin, rmax, binsize, $
                                etansum, eradsum, tansigsum, radsigsum, $
                                etanerrsum, eraderrsum, $
                                tansigerrsum, radsigerrsum, $
                                wsum, npsum, rsum, npair
                      
                              totalnpair = total(npair)
                              IF totalnpair NE 0 THEN BEGIN
                                  ;print,gi+1
                              ENDIF ELSE indices[ ind[gi] ] = -1
                              printf, lun1, totalnpair, total(npsum)
                              flush,lun1
                          ENDIF ELSE indices[ ind[gi] ] = -1
                      ENDIF ELSE indices[ ind[gi] ] = -1
                  ENDFOR 
              ENDIF ELSE indices[ ind ] = -1
          ENDIF ELSE indices[ ind ] = -1
      ENDFOR 

      ;; Remove unused lenses
      wbad = where(indices EQ -1, nwbad)
      IF nwbad NE 0 THEN remove, wbad, wlens
      lensused = n_elements(wlens)

      IF NOT rand THEN wgood = wlens

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find averages from sums
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      print,'Initial count: ',ntot
      print,'Used count: ',total(npsum)
      print
      print,'Finally used ',ntostr(lensused),'/',ntostr(nlens),' Lenses'

      FOR i=0L, nbin-1 DO BEGIN
          shear[i]    = etansum[i]/wsum[i]/2. ; 2's convert to shear.
          ortho[i]    = eradsum[i]/wsum[i]/2.
          shearerr[i] = sqrt( etanerrsum[i]/wsum[i]^2 )/2.
          orthoerr[i] = sqrt( eraderrsum[i]/wsum[i]^2 )/2.

          sigma[i]    = tansigsum[i]/wsum[i]/2.
          orthosig[i] = radsigsum[i]/wsum[i]/2.
          sigmaerr[i] = sqrt( tansigerrsum[i]/wsum[i]^2 )/2.
          orthosigerr[i] = sqrt( radsigerrsum[i]/wsum[i]^2 )/2.

          meanr[i]   = rsum[i]/npsum[i]
      ENDFOR 
      etansum[*] = 0. & eradsum[*] = 0. & etanerrsum[*] = 0. & eraderrsum[*]=0.
      tansigsum[*] = 0. & radsigsum[*] = 0.
      tansigerrsum[*] = 0. & radsigerrsum[*] = 0.
      wsum[*] = 0. & rsum[*] = 0. & npair[*] = 0.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Print the data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF NOT rand THEN BEGIN 
          makeps, psfile, /noland
          openw, lun2, datfile, /get_lun
          !textunit = lun2
      ENDIF ELSE BEGIN
          tmp = str_sep(datfile,addstr)
          randname = tmp[0] + addstr+'rand_' + tmp[1]
          openw, lun2, randname, /get_lun
          !textunit = lun2
      ENDELSE 

      message = '        meanr          shear        shearerr'
      message = message +  '        ortho         orthoerr'
      message = message +  '        sigma         sigmaerr'
      message = message +  '       orthosig      orthosigerr        npair'

      printf, lun2, 'Nlenses:   ',lensused
      printf, lun2, 'TotPairs:   ',total(npsum)
      printf, lun2, 'binwidth(arcsec): ',binsize
      printf, lun2, message
      printf, lun2, ' '
      fmt='(9E, F15.1)'

      colprint, meanr, shear, shearerr, ortho, orthoerr, $
                sigma, sigmaerr, orthosig, orthosigerr, npsum, $
                lun=lun2, format=fmt

      close, lun2
      free_lun, lun2

      npsum[*] = 0.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make some plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      plot,lenscat[wlens].ra, lenscat[wlens].dec, $
           xrange=[minra,maxra], yrange=[mindec,maxdec],$
           xtitle='RA',ytitle='DEC', title=ttt, $
           psym=sym

      IF nbin GT 1 THEN BEGIN 
          pold=!p.multi
          !p.multi = [0,1,2]

          IF rand THEN tit='Random Points' ELSE BEGIN
              IF n_elements(title) EQ 0 THEN BEGIN
                  tit = 'Centered on Galaxies'
              ENDIF ELSE tit=title
          ENDELSE 
          xt='Projected Radius (h^-1 kpc)'

          IF NOT rand THEN BEGIN 
              wmax1 = where(shear EQ max(shear))
              wmin1 = where(shear EQ min(shear))
              max1 = shear[wmax1] + shearerr[wmax1]
              min1 = shear[wmin1] - shearerr[wmin1]

              wmax2 = where(ortho EQ max(ortho))
              wmin2 = where(ortho EQ min(ortho))
              max2 = ortho[wmax2] + orthoerr[wmax2]
              min2 = ortho[wmin2] - orthoerr[wmin2]

              yrange1 = [min([min1,min2]),max([max1,max2])]

              wmax1 = where(sigma EQ max(sigma))
              wmin1 = where(sigma EQ min(sigma))
              max1 = sigma[wmax1] + sigmaerr[wmax1]
              min1 = sigma[wmin1] - sigmaerr[wmin1]

              wmax2 = where(orthosig EQ max(orthosig))
              wmin2 = where(orthosig EQ min(orthosig))
              max2 = orthosig[wmax2] + orthosigerr[wmax2]
              min2 = orthosig[wmin2] - orthosigerr[wmin2]

              yrange2 = [min([min1,min2]),max([max1,max2])]
          ENDIF 
          yt='Tangential Shear'
          ploterr, meanr, shear, shearerr, $
            psym=1,yrange=yrange1,xtitle=xt,ytitle=yt, title=tit
          oplot,[0,5000],[0,0]

          yt='Ortho-tangential Shear'
          ploterr, meanr, ortho, orthoerr, $
            psym=1,yrange=yrange1,xtitle=xt,ytitle=yt,$
            title=ntostr(lensused)+' lenses'
          oplot,[0,5000],[0,0]


          yt='<Sigma(<r)> - Sigma(r)   (h Msun/pc^2)'
          ploterr, meanr, sigma, sigmaerr, $
            psym=1, yrange=yrange2, xtitle=xt, ytitle=yt, title=tit
          oplot, [0,5000],[0,0]

          yt='Ortho'
          ploterr, meanr, orthosig, orthosigerr, $
            psym=1, yrange=yrange2, xtitle=xt, ytitle=yt,$
            title=ntostr(lensused)+' lenses'
          oplot, [0,5000],[0,0]

          !p.multi=pold

      ENDIF 
  ENDFOR ;; Looping over lens/random

  step=oldstep
  ep

  close, lun1
  free_lun, lun1

  ptime, systime(1)-time
  return
END 
