
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    RANDSHEAR_LAMBDA
;       
; PURPOSE:
;    Measure the shear and density contrast for a set of lenses. Structure
;    must have a redshift tags.
;
; CALLING SEQUENCE:
;    randshear_lambda, stripe, lenscat, scat, clr, rmin, rmax, binsize, nrand, $
;                      step=step, addstr=addstr, $
;                      outdir=outdir, $
;                      maxe=maxe, $
;                      wgood=wgood, $
;                      usecat=usecat, $
;                      datfile=datfile, sumfile=sumfile, $
;                      lensumfile=lensumfile
;
; INPUTS: 
;    stripe: stripe number in integer form.
;    lenscat: lens catalog
;    scat: source catalog. Must contain lambda and eta, and must be sorted
;          by lambda.
;    clr: the sdss bandpass
;    rmin: minimum radius in arcsec
;    rmax: maximum radius in arcsec
;    binsize: size of bins in arcsec
;    nrand: number of random points
;
; OPTIONAL INPUTS:
;    step: how many lenses to use in a group. default is 300
;    addstr: string to add to front of output file names.
;    outdir: output directory
;    maxe: maximum allowed ellipticity of source galaxy distribution.
;
; KEYWORD PARAMETERS:
;       
; OUTPUTS: 
;    datafile, sumfile, lensumfile
;
; OPTIONAL OUTPUTS:
;    wgood: the indices of used lenses.
;    usecat: the final lenses.
;    datafile, sumfile, lensumfile: the names of output files
;
; CALLED ROUTINES:
;    NTOSTR
;    TAG_EXIST
;    SHSTRUCT
;    SUMSTRUCT
;    LENSUMSTRUCT
;    EXIST
;    NEWNAME
;    (EQ2SURVEY)
;    (ROTATE_LAMBDA)
;    ARRSCL
;    COPY_STRUCT
;    RANDSHEAR_LAMBDA_GET_SOURCES (in this file)
;    BINARY_SEARCH
;    SHHDR
;    SUMHDR
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

PRO randshear_lambda_get_sources, lambda, minlam, maxlam, w, num, issouth=issouth

  IF (maxlam GT 180.) OR (minlam LT -180.) OR $
    (minlam GT maxlam) THEN BEGIN                  ; We DO cross -180,180 mark

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
          binary_search, lambda, minlam, i1, /round
          binary_search, lambda, maxlam, i2, /round
      ENDIF ELSE BEGIN 
          IF maxlam GT 0 THEN w=where(lambda GE 0.,ntot) $
          ELSE w=where(lambda LE 0.,ntot)
          IF ntot EQ 0 THEN BEGIN 
              w=-1
              num=0
              return
          ENDIF 
          binary_search, lambda[w], minlam, i1, /round
          binary_search, lambda[w], maxlam, i2, /round
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

PRO randshear_lambda, stripe, scat, clr, rmin, rmax, binsize,nrand, $
                       step=step, addstr=addstr, $
                       outdir=outdir, $
                       wgood=wgood, $
                       datfile=datfile, sumfile=sumfile, $
                       lensumfile=lensumfile, $
                       maxe=maxe, lensra=lensra, lensdec=lensdec

  IF n_params() LT 7 THEN BEGIN
      print,'-Syntax: randshear, stripe, scat, clr, rmin, rmax, binsize, nrand '
      print,'  step=step, addstr=addstr, '
      print,'  outdir=outdir, '
      print,'  title=title, wgood=wgood, '
      print,'  datfile=datfile, sumfile=sumfile, '
      print,'  lensumfile=lensumfile, '
      print,'  maxe=maxe '
      return
  ENDIF 

  colors = ['u','g','r','i','z']

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  COMMON seed,seed

  time = systime(1)

  vint = .32^2
  IF n_elements(maxe) EQ 0 THEN maxe = .2

  ;; Size by which to group things 300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L ELSE step = long(step)
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = ''
  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss4/data1/esheldon/TMP/'
  
  stripestr = '_stripe'+ntostr(stripe)
  clr_str = '_'+colors[clr]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; maxit is number of retries until random point
  ;; passes symmetry check
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  maxit = 20
                  

  d2r = !dpi/180.0d0            ;Change degrees to radians.
  r2d = 180./!dpi               ;radians to degrees

  radmax = rmax/3600.           ;rmax in degrees

  nbin = long( (rmax - rmin)/binsize ) + 1 

  print
  print,'-------------------------------------------------------------'
  print,'Rmin: ',rmin
  print,'Rmax: ',rmax
  print,'Binsize: ',binsize
  print,'Step: ',step
  print,'Using ',ntostr(nbin),' bins between ',ntostr(rmin), $
        ' and ',ntostr(rmax),' arcsec'

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

  shstruct = shstruct(arrval)
  sumstruct = sumstruct(arrval)
  lensumstruct = lensumstruct(arrval)

  lensum = replicate(lensumstruct, nrand )

  etansum    = arrval
  eradsum    = arrval
  etanerrsum = arrval
  eraderrsum = arrval

  wsum       = arrval
  rsum       = arrval
  npsum      = arrval
  npair      = arrval
  Sshsum     = 0.               ;Scalar
  rmax_act   = arrval
  rmax_act_tmp = arrval

  tmpetan = arrval
  tmperad = arrval
  tmpetanerr = arrval
  tmperad = arrval
  tmperaderr = arrval

  tmpwsum = arrval
  tmprsum = arrval
  tmpnpsum = arrval
  tmpSsh = 0.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF addstr EQ '' THEN prename = outdir+addstr+'rand'+stripestr+clr_str $
  ELSE prename = outdir+addstr+'_rand'+stripestr+clr_str

  datfile = prename + '_N1.fit'
  sumfile = prename + '_sum_N1.fit'
  lensumfile = prename + '_lensum_N1.fit'

  WHILE exist(datfile) DO BEGIN
      datfile = newname(datfile)
      sumfile = newname(sumfile)
      lensumfile = newname(lensumfile)
  ENDWHILE 
  print
  print,'Dat file: ',datfile
  print,'Sum file: ',sumfile
  print,'Lens sum file: ',lensumfile
  print


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; source cat
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; they are already sorted
  nsource = n_elements(scat)
  wsource = lindgen(nsource)

  FIRSTLAMBDA = scat[0].lambda
  LASTLAMBDA  = scat[nsource-1].lambda

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Lens cat: get random points
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nlens=nrand
  wlens = lindgen(nrand)

  llambda  = randomu(seed, nrand)
  leta = fltarr(nrand)          ;will generate on the fly later
  s = sort(llambda)
  llambda  = temporary(llambda[s])

  IF stripe GT 45 THEN BEGIN
      print
      print,'This is a Southern Stripe.'
      print
      issouth = 1
      ;; rotate, then rotate back
      rotate_lambda, FIRSTLAMBDA
      rotate_lambda, LASTLAMBDA
      llambda = arrscl(llambda, FIRSTLAMBDA+radmax, LASTLAMBDA-radmax, arrmin=0.,arrmax=1.)

      rotate_lambda, FIRSTLAMBDA
      rotate_lambda, LASTLAMBDA
      rotate_lambda, llambda

  ENDIF ELSE BEGIN 
      llambda = arrscl(llambda, FIRSTLAMBDA+radmax, LASTLAMBDA-radmax, arrmin=0.,arrmax=1.)
  ENDELSE 

  nlens = n_elements(wlens)
  print
  print,'FIRSTLAMBDA = ',FIRSTLAMBDA
  print,'LASTLAMBDA = ',LASTLAMBDA
  print,'Using '+ntostr(nlens)+' random points'
  print

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Estimate shear around all of these lenses
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
            
      maxii = wlens[ max(ind) ] & minii = wlens[ min(ind) ]

      maxlam1 = llambda[maxii]+radmax
      minlam1 = llambda[minii]-radmax

      randshear_lambda_get_sources, scat.lambda, minlam1, maxlam1, wsrc, nwsrc, $
        issouth=issouth

      IF nwsrc NE 0 THEN BEGIN 
          
          ;; generate eta's
          ;; Stripe is not rectangular in lambda-eta, this helps to 
          ;; keep random points from being thrown out

          tmpmaxeta = max(scat[wsrc].eta)
          tmpmineta = min(scat[wsrc].eta)
          leta[ii] = arrscl( randomu(seed, step), tmpmineta, tmpmaxeta, arrmin=0., arrmax=1.)

          xrel = dblarr(nwsrc)
          yrel = xrel
          IF group NE 0 THEN BEGIN 
              wsum_total = total(wsum)
              mean_shear = ntostr( total(etansum)/wsum_total/2. )
              mean_shear_err = ntostr( sqrt( total(etanerrsum)/wsum_total^2 )/2. )
          ENDIF ELSE BEGIN 
              mean_shear = '0'
              mean_shear_err = '0'
          ENDELSE 
          print,'Lens Group = ',ntostr(group+1)+'/'+ntostr(nstep),$
            '  Mean shear = ',mean_shear,' +/- ',mean_shear_err
          FOR gi=0L, step-1 DO BEGIN ;Loop over lenses in group
              
              IF (gi MOD 10) EQ 0 THEN  print,'.',format='(a,$)'
              index = ii[gi]
              ceneta = leta[index]
              cenlam  = llambda[index]
                            
              ;; choose sources around this lens
              maxlam2 = cenlam+radmax
              minlam2 = cenlam-radmax

              randshear_lambda_get_sources, scat[wsrc].lambda, minlam2, maxlam2, $
                twsrc, nwsrc2, $
                issouth=issouth
              
              ;; If there are any sources left, measure the shear
              IF nwsrc2 GT 1 THEN BEGIN 

                  wsrc2 = wsrc[twsrc]

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; Try to redo a few times until passes checks
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  redo=0
                  WHILE (redo NE -1) AND (redo LT maxit) DO BEGIN 
                      redo=redo+1
                      ;; calculate angular distance (dec->lambda in distance calc)
                      tscateta = scat[wsrc2].eta
                      etadiff = (ceneta-tscateta)*d2r
                      cosetadiff = cos(etadiff)
                      
                      tceneta = ceneta*d2r
                      tcenlam = cenlam*d2r
                      sincenlam = sin(tcenlam)
                      coscenlam = cos(tcenlam)
                  
                      tscatlam = scat[wsrc2].lambda*d2r
                      sinscatlam = sin(tscatlam)
                      cosscatlam = cos(tscatlam)
                      tscatlam = tscatlam*r2d

                      ;; Find distance in arcseconds
                      args = sincenlam*sinscatlam + coscenlam*cosscatlam*cosetadiff
                      warg = where(args LT -1., narg)
                      IF narg NE 0 THEN args[warg] = -1.
                      warg = where(args GT 1., narg)
                      IF narg NE 0 THEN args[warg] = 1.

                      R=acos( args )*r2d*3600.
                  
                      hist=histogram(R, binsize=binsize, min=rmin, $
                                     max=rmax,rever=rev_ind)
                      numbin=n_elements(hist)
                      IF numbin NE nbin THEN BEGIN
                          print,'What!!!'
                          print,numbin,nbin
                          return
                      ENDIF 
                      whist = where(hist NE 0, nhist)
                      ng=0
                      ;; Check if there are any in this annulus rmin-rmax
                      IF nhist NE 0 THEN BEGIN 
                          
                          
                          ;;nbin is predefined. Last bin will hold few
                          npair[*] = 0L
                          rmax_act_tmp[*] = 0.
                          FOR i=0L, nhist-1 DO BEGIN 
                              binnum = whist[i]
                              w=rev_ind( rev_ind(binnum):rev_ind(binnum+1)-1 )
                              
                              rmax_act_tmp[binnum] = max(R[w])
                              rmax_act[binnum] = max( [rmax_act[binnum], rmax_act_tmp[binnum]] )
                              npair[binnum] = n_elements(w)
                              
                              ;; Original way
                              theta=atan(sin(abs(etadiff[w])), $
                                         (sincenlam*cosetadiff[w] - $
                                          coscenlam*sinscatlam[w]/cosscatlam[w]) )-!dpi

                              xrel[w] = R[w]*sin(theta)
                              yrel[w] = R[w]*cos(theta)
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
                              diffsq=0 & xy=0 & theta=0 & e1prime=0 & e2prime=0 & wts=0 & w=0
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
                          
                              ;; copy in individual lens stuff
                              lensum[index].totpairs = total(tmpnpsum[whist])
                              lensum[index].npair[whist] = tmpnpsum[whist]
                              lensum[index].ie = ie
                              lensum[index].rsum[whist] = tmprsum[whist]
                              lensum[index].rmax_act[whist] = rmax_act_tmp[whist]
                              lensum[index].etansum[whist] = tmpetan[whist]
                              lensum[index].eradsum[whist] = tmperad[whist]

                              lensum[index].etanerrsum[whist] = tmpetanerr[whist]
                              lensum[index].eraderrsum[whist] = tmperaderr[whist]
                              
                              lensum[index].sshsum = tmpSsh
                              
                              lensum[index].wsum[whist] = tmpwsum[whist]

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
                          tmpnpsum[whist] = 0
                          tmpSsh = 0L
                      ENDIF 
                      ;; One final check.  If not passed, remember that this lens wasn't used
                      IF ng EQ 0 THEN BEGIN 
                          qeta = randomu(seed,10)
                          qeta = (arrscl( qeta, tmpmineta, tmpmaxeta, arrmin=0.,arrmax=1.))[0]
                          leta[index] = qeta
                          ceneta = qeta
                      ENDIF ELSE redo=-1
                      R = 0
                      radiff=0  ;Free memory
                  ENDWHILE 
                  IF redo NE -1 THEN BEGIN
                      print,'/',format='(a,$)'
                      indices[ind[gi]] = -1
                  ENDIF 
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
  
  lensw = lensum[wlens].totpairs
  lenswsum = total( lensw )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find averages from sums
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Finally used ',ntostr(lensused),'/',ntostr(nlens),' Lenses'

  totpairs = 0L

  IF npsum[nbin-1] LT npsum[nbin-2] THEN BEGIN 
      ;; This means last bin was incomplete
      nbin = nbin-1
  ENDIF 

  FOR i=0L, nbin-1 DO BEGIN 
      ;; shear struct for this set of lenses
      shstruct.meanr[i] = rsum[i]/npsum[i]

      shstruct.shear[i] = etansum[i]/wsum[i]/2.
      shstruct.shearerr[i] = sqrt(etanerrsum[i]/wsum[i]^2)/2.
      shstruct.ortho[i] = eradsum[i]/wsum[i]/2.
      shstruct.orthoerr[i] = sqrt(eraderrsum[i]/wsum[i]^2)/2.

      shstruct.npair[i] = npsum[i]

      ;; Now totals within each rmax_act[i]
      shstruct.tshear[i] = total(etansum[0:i])/total(wsum[0:i])/2.
      shstruct.tshearerr[i] = sqrt( total(etanerrsum[0:i])/total(wsum[0:i])^2 )/2.
      shstruct.tortho[i] = total(eradsum[0:i])/total(wsum[0:i])/2.
      shstruct.torthoerr[i] = sqrt( total(eraderrsum[0:i])/total(wsum[0:i])^2 )/2.

      shstruct.tnpair[i] = total(npsum[0:i])

      ;; calculate area, density of background galaxies
      R1 = rmin + i*binsize
      R2 = rmin + (i+1)*binsize
      shstruct.area[i] = !pi*(R2^2 - R1^2)
      shstruct.density[i] = shstruct.npair[i]/shstruct.area[i]/lensused
      
      ; sum struct for adding to other sets of lenses
      sumstruct.rsum[i] = rsum[i]

      sumstruct.etansum[i] = etansum[i]
      sumstruct.etanerrsum[i] = etanerrsum[i]
      sumstruct.eradsum[i] = eradsum[i]
      sumstruct.eraderrsum[i] = eraderrsum[i]

      sumstruct.wsum[i] = wsum[i]
      sumstruct.npair[i] = npsum[i]

      totpairs = totpairs + npsum[i]
  ENDFOR 
  Ssh = Sshsum/total(wsum)

  shstruct.nlenses = lensused
  shstruct.totpairs = totpairs
  shstruct.binsize = binsize
  shstruct.rmin = rmin
  shstruct.rmax = rmax
  shstruct.rmax_act = rmax_act
  shstruct.ssh = ssh

  sumstruct.nlenses = lensused
  sumstruct.totpairs = totpairs
  sumstruct.binsize = binsize
  sumstruct.rmin = rmin
  sumstruct.rmax = rmax
  sumstruct.rmax_act = rmax_act
  sumstruct.sshsum = Sshsum
  sumstruct.lenswsum = lenswsum

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Ouput the data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Shear for these lenses
  shhdr = shhdr(shstruct)   
  mwrfits, shstruct, datfile, shhdr, /create

  ;; Sum file for these lenses
  sumhdr = sumhdr(sumstruct)
  lhdr = sumhdr
  mwrfits, sumstruct, sumfile, sumhdr, /create

  ;; File with structure for each lens used
  lensum = temporary(lensum[wlens])
  mwrfits, lensum, lensumfile, lhdr, /create

  step=oldstep

  ptime, systime(1)-time
  return
END 
