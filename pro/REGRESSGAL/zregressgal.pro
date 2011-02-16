
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    ZOBJSHEAR_LAMBDA
;       
; PURPOSE:
;    Measure the shear and density contrast for a set of lenses. Structure
;    must have a redshift tags.
;
; CALLING SEQUENCE:
;    zobjshear_lambda, stripe, lenscat, scat, clr, rmin, rmax, binsize, $
;                      use_lambda=use_lambda, $
;                      step=step, addstr=addstr, $
;                      outdir=outdir, $
;                      maxe=maxe, $
;                      wgood=wgood, $
;                      usecat=usecat, $
;                      check=check, $
;                      datfile=datfile, sumfile=sumfile, zfile=zfile, $
;                      lensumfile=lensumfile
;
; INPUTS: 
;    stripe: stripe number in integer form.
;    lenscat: lens catalog. Must have either (ra,dec) or (lambda,eta)
;    scat: source catalog. Must contain (lambda,eta), and must be sorted
;          by lambda.
;    clr: the sdss bandpass
;    rmin: minimum radius in kpc
;    rmax: maximum radius in kpc
;    binsize: size of bins in kpc
;
; OPTIONAL INPUTS:
;    /use_lambda: lambda=0.7 else lambda=0 (omega_tot=1)
;    step: how many lenses to use in a group. default is 300
;    addstr: string to add to front of output file names.
;    outdir: output directory
;    maxe: maximum allowed ellipticity of source galaxy distribution.

; KEYWORD PARAMETERS:
;    /check: just see which lenses are in allowd redshift range, set wgood and return.
;       
; OUTPUTS: 
;    datafile, sumfile, zfile, lensumfile
;
; OPTIONAL OUTPUTS:
;    wgood: the indices of used lenses.
;    usecat: the final lenses.
;    datafile, sumfile, zfile, lensumfile: the names of output files
;
; CALLED ROUTINES:
;    NTOSTR
;    TAG_EXIST
;    ZSHSTRUCT
;    ZSUMSTRUCT
;    ZLENSUMSTRUCT
;    EXIST
;    NEWNAME
;    (EQ2SURVEY)
;    (ROTATE_LAMBDA)
;    ANGDIST_LAMBDA
;    COPY_STRUCT
;    SDSS_SIGMA_CRIT
;    ZOBJSHEAR_LAMBDA_GET_SOURCES (in this file)
;    BINARY_SEARCH
;    ZSHHDR
;    ZSUMHDR
;    MWRFITS
;    
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    03-DEC-2000: First documentation. Erin Scott Sheldon.
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO bin_lenses_by_redshift, z, rmax, rev_ind, whist, nhist, zmaxangle

  ;; note that the histogram preserves whatever sort there is

  zbinsize=.01
  rmax=1000.
  hist = histogram(z, binsize=zbinsize, rever=rev_ind, min=0.0)
  whist=where(hist NE 0, nhist)

  zmaxangle = fltarr(nhist)

  IF nhist NE 0 THEN BEGIN 

      FOR i=0L, nhist-1 DO BEGIN 

          binnum=whist[i]
          wz=rev_ind[ rev_ind[binnum]:rev_ind[binnum+1]-1 ]
          nbin=n_elements(wz)

          ;; choosing minimum z in bin will gaurantee all
          ;; are found withing zmaxangle
          minz = min( z[wz] )

          tmpDL = angdist_lambda( minz )*1000.

          zmaxangle[binnum] = rmax/tmpDL*180./!pi ;degrees

          print,nbin,meanz,maxangle[binnum]*3600.
      ENDFOR 

  END 

  message,'testing'

END 

PRO zregressgal_get_neighbors, lambda, minlam, maxlam, $
                               w, num, issouth=issouth

  IF (maxlam GT 180.) OR (minlam LT -180.) OR $
    (minlam GT maxlam) THEN BEGIN                  ; We DO cross -180,180 mark
      print,'Crossed lambda=[180,-180] minlam = '+ntostr(minlam)+' maxlam = '+ntostr(maxlam)
      
      IF maxlam GT 180. THEN BEGIN 
          diff = maxlam - 180.
          w1 = where(lambda LE 180. AND $
                     lambda GE minlam, n1)
          w2 = where(lambda GE -180. AND $
                     lambda LE (-180.+diff), n2)
      ENDIF ELSE IF minlam LT -180 THEN BEGIN 
          diff = -180. - minlam
          w1 = where(lambda LE maxlam AND $
                     lambda GE -180., n1)
          w2 = where(lambda LE 180. AND $
                     lambda GE (180. - diff), n2)
      ENDIF ELSE BEGIN 
          w=where( (lambda GE minlam) OR $
                   (lambda LE maxlam), num)
          return
      ENDELSE 
      CASE 1 OF
          (n1 NE 0) AND (n2 NE 0): BEGIN 
              w=[w1,w2]
              num=n1+n2
          END 
          (n1 NE 0): BEGIN
              w=w1
              num=n1
          END 
          (n2 NE 0): BEGIN
              w=w2
              num=n2
          END 
          ELSE: BEGIN
              w=-1
              num=0
          END 
      ENDCASE 
  ENDIF ELSE BEGIN ; We don't cross -180,180 mark
      ntot=n_elements(lambda)
      w=lindgen(ntot)

      IF NOT keyword_set(issouth) THEN BEGIN 
          binary_search, lambda, minlam, i1, /round, /edgedefault
          binary_search, lambda, maxlam, i2, /round, /edgedefault
      ENDIF ELSE BEGIN 
          IF maxlam GT 0 THEN w=where(lambda GE 0.,ntot) $
          ELSE w=where(lambda LE 0.,ntot)
          IF ntot EQ 0 THEN BEGIN 
              w=-1
              num=0
              return
          ENDIF 
          binary_search, lambda[w], minlam, i1, /round, /edgedefault
          binary_search, lambda[w], maxlam, i2, /round, /edgedefault
      ENDELSE 
      CASE 1 OF
          (i1 NE -1) AND (i2 NE -1): BEGIN
              w=w[i1:i2]
              num=n_elements(w)
          END 
          (i1 EQ -1) AND (i2 NE -1): BEGIN
              w=w[0:i2]
              num=n_elements(w)
          END 
          (i2 EQ -1) AND (i1 NE -1): BEGIN
              w=w[i1:ntot-1]
              num=n_elements(w)
          END 
          ELSE: BEGIN
              w=-1
              num=0
          END 
      ENDCASE 
  ENDELSE 
  return
END 

PRO zregressgal, stripe, lenscat, scat, clr, rmin, rmax, binsize, $
                      use_lambda=use_lambda, $
                      step=step, addstr=addstr, $
                      outdir=outdir, $
                      wgood=wgood, $
                      usecat=usecat, $
                      check=check, $
                      datfile=datfile, sumfile=sumfile, zfile=zfile, $
                      lensumfile=lensumfile, $
                      maxe=maxe, simulation=simulation

  IF n_params() LT 7 THEN BEGIN
      print,'-Syntax: zobjshear, stripe, lenscat, scat, clr, rmin, rmax, binsize,'
      print,'  use_lambda=use_lambda,'
      print,'  step=step, addstr=addstr, '
      print,'  outdir=outdir, '
      print,'  wgood=wgood, '
      print,'  usecat=usecat, '
      print,'  check=check, '
      print,'  datfile=datfile, sumfile=sumfile, zfile=zfile, '
      print,'  lensumfile=lensumfile, '
      print,'  maxe=maxe '
      print,'-> lenscat may have ra,dec or lambda, eta but scat must contain lambda,eta'
      print,'-> scat must be sorted by lambda'
      return
  ENDIF 

  colors = ['u','g','r','i','z']

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;COMMON zregressgal_block, zrev_ind, whist, nhist, zmaxangle


  time = systime(1)

  vint = .32^2
  IF n_elements(maxe) EQ 0 THEN maxe = .2
  IF NOT keyword_set(use_lambda) THEN omegamat=1.0 ELSE omegamat=0.3

  ;; Size by which to group things  300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = 'zregressgal'
  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss4/data1/esheldon/TMP/'
  IF NOT keyword_set(check) THEN check=0

  stripestr = '_stripe'+ntostr(stripe)
  clr_str = '_'+colors[clr]

  d2r = !dpi/180.0d0            ;Change degrees to radians.
  r2d = 180./!dpi               ;radians to degrees

  nbin = long( (rmax - rmin)/binsize ) + 1 

  print
  print,'-------------------------------------------------------------'
  IF NOT keyword_set(use_lambda) THEN print,'Using lambda = 0' $
  ELSE print,'Using lambda = 0.3'
  print,'Rmin: ',rmin
  print,'Rmax: ',rmax
  print,'Binsize: ',binsize
  print,'Step: ',step
  print,'Using ',ntostr(nbin),' bins between ',ntostr(rmin), $
        ' and ',ntostr(rmax),' kpc'

  IF NOT tag_exist(scat,'momerr',index=momerr) THEN BEGIN
      IF NOT tag_exist(scat,'uncert',index=momerr) THEN BEGIN 
          print
          print,'No valid moment uncertainty tag'
          print
          return
      ENDIF 
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; declare some arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  arrval = fltarr(nbin)

  rsum       = arrval
  npsum      = arrval
  npair      = arrval

  lens_npair = arrval
  lens_struct = lenscat[0]
  zero_struct,lens_struct

  lensumstruct = create_struct(lens_struct, $
                               'zindex', 0L,$
                               'totpairs', 0.0, $
                               'npair', arrval, $
                               'rsum', arrval, $
                               'scritinv', 0.0)
  lensum = replicate(lensumstruct, n_elements(lenscat))
  
  av = arrval
  bv = arrval
  ata = fltarr(nbin,nbin)
  btb = ata

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  prename = outdir+addstr+stripestr+clr_str
  
  datfile = prename + '_N1.sav'
  lensumfile = prename + '_lensum_N1.fit'
  source_positionfile = prename+'_position_N1.fit'
  WHILE exist(datfile) DO BEGIN
      datfile = newname(datfile)
      lensumfile = newname(lensumfile)
      source_positionfile = newname(source_positionfile)
  ENDWHILE 
  print
  print,'Dat save file: ',datfile
  print,'Lens sum file: ',lensumfile
  print,'Position file: ',source_positionfile
  print

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; source cat
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; they are already sorted
  nsource = n_elements(scat)
  ninit=nsource
  wsource = lindgen(nsource)

  FIRSTLAMBDA = scat[0].lambda
  LASTLAMBDA  = scat[nsource-1].lambda 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Lens cat: Set up sigma crit and Dlens. 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  lcat = lenscat                ;make copy which will be sorted by lambda
  lninit = n_elements(lcat)

  IF NOT tag_exist(lcat,'z1d',index=wz) THEN BEGIN
      IF NOT tag_exist(lcat,'z',index=wz) THEN BEGIN 
          IF NOT tag_exist(lcat,'photoz',index=wz) THEN BEGIN
              print,'Lens structure must have "Z" or "PHOTOZ" or "Z1D" flag'
              return
          ENDIF 
      ENDIF 
  ENDIF 

  IF tag_exist(lcat, 'lambda') AND tag_exist(lcat, 'eta') THEN BEGIN 
      llambda = lcat.lambda
      leta = lcat.eta
  ENDIF ELSE BEGIN 
      IF tag_exist(lcat, 'ra') AND tag_exist(lcat, 'dec') THEN BEGIN 
          print
          print,'-------------------------------------'
          print,'Converting lens RA,DEC to LAMBDA,ETA'
          print,'-------------------------------------'
          print
          eq2survey, lcat.ra, lcat.dec, llambda, leta
      ENDIF ELSE BEGIN 
          print,'Neither (LAMBDA, ETA) or (RA,DEC) found in lens tags'
          return
      ENDELSE 
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; If it is a southern stripe, deal with fact that it
  ;; crosses [-180,180] mark (possibly)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF stripe GT 45 THEN BEGIN
      print
      print,'This is a Southern Stripe.'
      print
      issouth = 1
      ;; rotate, sort, then rotate back
      rotate_lambda, llambda
      s = sort(llambda)
      lcat = temporary(lcat[s])
      llambda = temporary( llambda[s] )
      leta = temporary( leta[s] )
      ;; rotate back
      rotate_lambda, llambda
  ENDIF ELSE BEGIN 
      ;; sort them
      issouth=0
      s = sort(llambda)
      lcat = temporary(lcat[s])
      llambda = temporary( llambda[s] )
      leta = temporary( leta[s] )
  ENDELSE 

  ;; subscripts for lenses. This will change as we throw
  ;; out lenses
  wlens = lindgen(lninit)

  h=1.
  print,'Using h = ',h
  print,'Calculating 1/sigma_crit'
  print

  sigcritinv = sdss_sigma_crit(stripe, clr, lcat.(wz[0]), $
                               wgood=wgood, use_lambda=use_lambda)
  sigmacrit = 1./sigcritinv

  ;;Convert from Mpc to kpc
  DL = angdist_lambda( lcat.(wz[0]), h=h, omegamat=omegamat)*1000. 
  angmax = rmax/DL*180./!pi     ;angles in degrees

  ;; don't want angle to be larger than this angle
  MAX_ALLOWED_ANGLE = 0.35 ;; degrees (z ~ 0.06)
  ;;MAX_ALLOWED_ANGLE = 3.
  wgood2 = where(angmax[wgood] LE MAX_ALLOWED_ANGLE)
  wgood = wgood[wgood2]
  wlens = wlens[wgood]

  nlens1 = n_elements(wlens)
  print,'Threw out ',ntostr(lninit - nlens1),' lenses because too deep or close'

  ;; copy in the stuff we need into lensum struct
  copy_struct, lcat, lensum
  lensum.scritinv = sigcritinv

  ;; now bin lenses by redshift; return the rev_ind of each bin, and
  ;; the max angle for each bin.

  ;bin_lenses_by_redshift, lcat[wlens].(wz[0]), rmax, $
  ;  zrev_ind, zwhist, znhist, $
  ;  zmaxangle

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; apply mask to sources
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT keyword_set(simulation) THEN BEGIN 
      print,'Applying mask'
      read_stripe_mask, stripe, mask, /regress
      apply_mask, mask, scat.lambda, bad, wsource, edgecut=MAX_ALLOWED_ANGLE
      nsource = n_elements(wsource)

      read_etarange, stripe, etarange

      print,'Interpolating eta'
      maxeta = interpol(etarange.maxeta, etarange.lambda, scat[wsource].lambda)
      mineta = interpol(etarange.mineta, etarange.lambda, scat[wsource].lambda)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Now eta edge cut
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      print
      print,'Making eta edge cut'
      twsource = wsource
      FOR i=0L, nsource-1 DO BEGIN 
          index = wsource[i]
          
          IF ( (abs(maxeta[i]-scat[index].eta) LT MAX_ALLOWED_ANGLE ) OR $
               (abs(mineta[i]-scat[index].eta) LT MAX_ALLOWED_ANGLE ) ) THEN BEGIN 
              twsource[i] = -1
          ENDIF   
          
      ENDFOR 

      wtmp = where(twsource NE -1,nsource)
      IF nsource EQ 0 THEN message,'No objects passed distance cut!'
      print,'Threw out ',ntostr(n_elements(wsource)-nsource),' sources on edge cut'
      print
      IF (nsource NE ninit) THEN wsource = wsource[wtmp] 
  

  ENDIF ELSE BEGIN 
      nsource=n_elements(scat)
      wsource=lindgen(nsource)
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; print some stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'FIRSTLAMBDA = ',FIRSTLAMBDA
  print,'LASTLAMBDA = ',LASTLAMBDA
  print,'Using '+ntostr(nsource)+'/'+ntostr(ninit)+' sources'
  print
  
  ttt='Lenses'

;  simpctable
;  plot,llambda[wlens],leta[wlens],psym=3
;  oplot,scat[wsource].lambda,scat[wsource].eta,psym=3,color=!red
;return
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Estimate shear around all of these lenses
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

      ;; Choose lenses around this source group
      ;; they are sorted by lambda, but it could go over lambda=180      
      angmax2 = max(angmax[ii])
      
      maxii = wsource[ max(ind) ] & minii = wsource[ min(ind) ]

      maxlam1 = scat[maxii].lambda + MAX_ALLOWED_ANGLE
      minlam1 = scat[minii].lambda - MAX_ALLOWED_ANGLE

      zregressgal_get_neighbors, llambda[wlens], minlam1, maxlam1, $
        wlns, nwlns, $
        issouth=issouth

;      wlns = where( (abs(llambda[wlens]-scat[maxii].lambda) LT MAX_ALLOWED_ANGLE) $
;                    OR $
;                    (abs(llambda[wlens]-scat[minii].lambda) LT MAX_ALLOWED_ANGLE),$
;                    nwlns)

;help,wlens,nwlns
;message,'stopping'

      IF nwlns NE 0 THEN BEGIN 
          wlns=wlens[wlns]

          ;xrel = dblarr(nwlns)
          ;yrel = xrel

          print,'Source Group = ',ntostr(group+1)+'/'+ntostr(nstep)
          FOR gi=0L, step-1 DO BEGIN ;Loop over sources in group
              IF (gi MOD 10) EQ 0 THEN  print,'.',format='(a,$)'
              index = ii[gi]
              ceneta = scat[index].eta ;source center
              cenlam  = scat[index].lambda

              ;; choose lenses around this source
              maxlam2 = cenlam+MAX_ALLOWED_ANGLE
              minlam2 = cenlam-MAX_ALLOWED_ANGLE

              zregressgal_get_neighbors, llambda[wlns], minlam2, maxlam2, $
                twlns, nwlns2, $
                issouth=issouth

;              twlns = where( (abs(llambda[wlns]-scat[maxii].lambda) LT MAX_ALLOWED_ANGLE) $
;                            OR $
;                            (abs(llambda[wlns]-scat[minii].lambda) LT MAX_ALLOWED_ANGLE),$
;                            nwlns2)

              ;; If there are any sources left, measure the shear
              IF nwlns2 GT 1 THEN BEGIN 

                  ;; calculate angular distance (dec->lambda in distance calc)

                  wlns2 = wlns[twlns]
                  tlenseta = leta[wlns2]
                  etadiff = (ceneta-tlenseta)*d2r
                  cosetadiff = cos(etadiff)
                  
                  tcenlam = cenlam*d2r ;lens center in radians
                  sincenlam = sin(tcenlam)
                  coscenlam = cos(tcenlam)
                  
                  tlenslam = llambda[wlns2]*d2r
                  sinlenslam = sin(tlenslam)
                  coslenslam = cos(tlenslam)
;                  tscatlam = tscatlam*r2d
                  ;; Find distance in kpc
                  args = sincenlam*sinlenslam + coscenlam*coslenslam*cosetadiff
                  warg = where(args LT -1., narg)
                  IF narg NE 0 THEN args[warg] = -1.
                  warg = where(args GT 1., narg)
                  IF narg NE 0 THEN args[warg] = 1.

                  ;sig_crit = sigmacrit[wlns2]
                  sig_crit = replicate(1.0, nwlns2)
                  R=acos( args )*DL[wlns2]
                  theta = atan( sin(etadiff), $
                                (sincenlam*cosetadiff - $
                                 coscenlam*sinlenslam/coslenslam) )-!dpi

                  tnan=where(finite(theta,/NAN),ntnan)
                  IF ntnan NE 0 THEN BEGIN
                      ;;print,'found a nan'
                      R[tnan] = 0.0
                      theta[tnan] = 0.0
                  ENDIF 
                  npair[*] = 0L
                  tlens_npair = lonarr(nwlns2,nbin)
                  tlens_rsum = fltarr(nwlns2,nbin)
                  errsend = sqrt(scat[index].momerr^2 + 0.32^2)
                  build_system, R, theta, sig_crit, $
                                scat[index].e1,scat[index].e2,errsend,$
                                rmin,rmax,binsize,nbin,$
                                ata,btb,av,bv,npair,rsum,tlens_npair,tlens_rsum
                  ;;IF (gi MOD 100) EQ 0 THEN BEGIN
           
                      ;;print
                      ;;print,'ata',max(ata),min(ata)
                      ;;print,'btb',max(btb),min(btb)
                      ;;print,'av',max(av),min(av)
                      ;;print,'bv',max(bv),min(bv)
                      ;message,'stopping:  gi = '+ntostr(gi)

                  ;;ENDIF 
                  nadd = total(npair)
                  IF nadd EQ 0. THEN BEGIN
                      print,'|',format='(a,$)'
                      indices[ ind[gi] ] = -1
                  ENDIF ELSE BEGIN 
                      lensum[wlns2].npair[*] = lensum[wlns2].npair[*]+tlens_npair
                      lensum[wlns2].rsum[*] = lensum[wlns2].rsum[*]+tlens_rsum
                      npsum[*] = npsum[*] + npair[*]
                  ENDELSE 

                  ;; Free memory
                  setzero, R, theta, etadiff, cosetadiff, sinlenslam,coslenslam,$
                           args, tlens_npair, tlens_rsum
              ENDIF ELSE BEGIN
                  print,'|',format='(a,$)'
                  indices[ ind[gi] ] = -1
              ENDELSE 

          ENDFOR 
          print
          ;xrel = 0 & yrel = 0
      ENDIF ELSE BEGIN 
          print,'Two-'
          indices[ ind ] = -1
      ENDELSE 
      metadiff = 0               ;Free memory
      mdist = 0
  ENDFOR 
  ;; Remove unused lenses
  wbad = where(indices EQ -1, nwbad)
  IF nwbad NE 0 THEN remove, wbad, wsource
  sourceused = n_elements(wsource)
  
  print
  print,'Finally used ',ntostr(sourceused),'/',ntostr(n_elements(scat)),' Sources'
  
  print,'Outputting to : ',datfile
  print,'Outputting lensum file: ',lensumfile
  print,'Outputting position file: ',source_positionfile

  nlens=n_elements(lensum)
  FOR i=0L, nlens-1 DO BEGIN 
      lensum[i].totpairs = total(lensum[i].npair)
  ENDFOR 

  wkeep = where(lensum.totpairs NE 0.0,nkeep)
  lensum = lensum[wkeep]
  lensum.zindex = lindgen(nkeep)
  mwrfits, lensum, lensumfile, /create

  meanr = rsum/npsum
  npair = npsum

  save, ata, btb, av, bv, npsum, rsum, npair, meanr, filename=datfile

  ;; output position file
  pos=create_struct('lambda', 0d,$
                    'eta', 0d)
  pos=replicate(pos, sourceused)
  pos.lambda = scat[wsource].lambda
  pos.eta = scat[wsource].eta

  mwrfits, pos, source_positionfile, /create

  ptime, systime(1)-time
  return
END 
