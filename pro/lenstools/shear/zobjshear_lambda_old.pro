
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

PRO zobjshear_lambda, stripe, lenscat, scat, clr, rmin, rmax, binsize, $
                      commonsrc=commonsrc,$
                      use_lambda=use_lambda, $
                      step=step, addstr=addstr, $
                      outdir=outdir, $
                      wgood=wgood, $
                      usecat=usecat, $
                      check=check, $
                      datfile=datfile, sumfile=sumfile, zfile=zfile, $
                      lensumfile=lensumfile, $
                      maxe=maxe

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

  time = systime(1)

  vint = .32^2
  IF n_elements(maxe) EQ 0 THEN maxe = .2
  IF NOT keyword_set(use_lambda) THEN omegamat=1.0 ELSE omegamat=0.3

  ;; Size by which to group things  300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = 'zgal_gal'
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

  shstruct = zshstruct(arrval)
  sumstruct = zsumstruct(arrval)
  lensumstruct = zlensumstruct(arrval)
  lensumstruct = create_struct(lenscat[0], lensumstruct)

  lensum = replicate(lensumstruct, n_elements(lenscat) )

  etansum    = arrval
  eradsum    = arrval
  etanerrsum = arrval
  eraderrsum = arrval

  tansigsum  = arrval
  radsigsum  = arrval
  tansigerrsum = arrval
  radsigerrsum = arrval

  wsum       = arrval
  rsum       = arrval
  npsum      = arrval
  npair      = arrval
  Sshsum     = 0.               ;Scalar
  wsum_ssh   = 0.
  rmax_act   = arrval
  rmax_act_tmp = arrval

  ;; temporary
  tmpetan = arrval
  tmperad = arrval
  tmpetanerr = arrval
  tmperad = arrval
  tmperaderr = arrval

  tmptansig = arrval
  tmpradsig = arrval
  tmptansigerr = arrval
  tmpradsigerr = arrval

  tmpwsum = arrval
  tmprsum = arrval
  tmpnpsum = arrval
  tmpwsum_ssh = 0.
  tmpSsh = 0.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  prename = outdir+addstr+stripestr+clr_str
  
  datfile = prename + '_N1.fit'
  sumfile = prename + '_sum_N1.fit'
  zfile   = prename + '_z_N1.fit'
  lensumfile = prename + '_lensum_N1.fit'
  groupfile = prename +'_groupstats_N1.fit'

  WHILE exist(datfile) DO BEGIN
      datfile = newname(datfile)
      sumfile = newname(sumfile)
      zfile = newname(zfile)
      lensumfile = newname(lensumfile)
      groupfile = newname(groupfile)
  ENDWHILE 
  print
  print,'Dat file: ',datfile
  print,'Sum file: ',sumfile
  print,'Lens sum file: ',lensumfile
  print,'Redshift file: ',zfile
  print,'Group stats file: ',groupfile
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
; Lens cat: Set up sigma crit and Dlens. 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  lcat = lenscat                ;make copy which will be sorted by lambda
  ninit = n_elements(lcat)

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
  wlens = lindgen(ninit)

  h=1.
  print,'Using h = ',h
  print,'Calculating 1/sigma_crit'
  print

  sigcritinv = sdss_sigma_crit(stripe, clr, lcat.(wz[0]), $
                               wgood=wgood, use_lambda=use_lambda, $
                               commonsrc=commonsrc)
  sigmacrit = 1./sigcritinv

  ;;Convert from Mpc to kpc
  DL = angdist_lambda( lcat.(wz[0]), h=h, omegamat=omegamat)*1000. 
  angmax = rmax/DL*180./!pi     ;angles in degrees

  ;; don't want angle to be larger than 1/rfac of a stripe (until we use 
  ;; multiple stripes)
  rfac = 1.0
  max_allowed_angle = 2.5/rfac
  wgood2 = where(angmax[wgood] LE max_allowed_angle)
  wgood = wgood[wgood2]
  wlens = wlens[wgood]

  nlens1 = n_elements(wlens)
  print,'Threw out ',ntostr(ninit - nlens1),' lenses because too deep or close'

  ;; apply mask (rotated frame)

  read_stripe_mask, stripe, mask, tsgals=tsgals, indir=maskdir
  IF issouth THEN BEGIN
      rotate_lambda, llambda
      rotate_lambda, FIRSTLAMBDA
      rotate_lambda, LASTLAMBDA
  ENDIF 

  apply_mask, mask, llambda[wlens], tbad, tgood, edgecut=angmax[wlens]
  IF tgood[0] EQ -1 THEN message,'No objects passed mask cuts' 
  nlens2 = n_elements(tgood)
  print,'Threw out ',ntostr(nlens1-nlens2),' lenses in mask'
  wlens = wlens[tgood]

  ;; apply lambda cut
  wtmp = where( ((llambda[wlens] - angmax[wlens]) GE FIRSTLAMBDA ) AND $
                ((llambda[wlens] + angmax[wlens]) LE LASTLAMBDA ), nlens)
  IF nlens EQ 0 THEN message,'No objects passed lambda cut'
  wlens = wlens[wtmp]
  print,'Threw out ',ntostr(nlens2-nlens),' lenses in lambda cut'

  IF issouth THEN BEGIN
      rotate_lambda, llambda
      rotate_lambda, FIRSTLAMBDA
      rotate_lambda, LASTLAMBDA
  ENDIF 
 
  ;; copy in the stuff we need into lensum struct
  copy_struct, lcat, lensum
  lensum.scritinv = sigcritinv

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; print some stuff
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  maxeta = max(scat.eta )
  mineta = min(scat.eta )

  IF check THEN BEGIN           ;May just want to find the good lenses.
      wgood = wlens
      return
  ENDIF 

  nlens = n_elements(wlens)
  print
  print,'FIRSTLAMBDA = ',FIRSTLAMBDA
  print,'LASTLAMBDA = ',LASTLAMBDA
  print,'Using '+ntostr(nlens)+'/'+ntostr(ninit)+' lenses'
  print
  IF nlens LT 100 THEN sym = 1 ELSE sym = 3
  
  ttt='Lenses'

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
  groupstruct = create_struct('groupid',lindgen(nstep),$
                              'groupsize',step,$
                              'nlens',lonarr(nstep),$
                              'shear',replicate(-1.e10,nstep),$
                              'shearerr',replicate(-1.e10,nstep),$
                              'mean_lambda',replicate(-1.e10,nstep),$
                              'mean_eta',replicate(-1.e10,nstep))

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
      ;; they are sorted by lambda, but it could go over lambda=180      
      angmax2 = max(angmax[ii])
      
      maxii = wlens[ max(ind) ] & minii = wlens[ min(ind) ]

      maxlam1 = llambda[maxii]+angmax2
      minlam1 = llambda[minii]-angmax2

      zobjshear_lambda_get_sources, scat.lambda, minlam1, maxlam1, wsrc, nwsrc, $
        issouth=issouth

      ;; for group statistics
      groupused=0L & lamsum=0d & etasum=0d
      IF nwsrc NE 0 THEN BEGIN 
          
          xrel = dblarr(nwsrc)
          yrel = xrel
          wsum_total = total(wsum)
          IF wsum_total GT 0.0 THEN BEGIN 
              mean_shear_p = ntostr( total(etansum)/wsum_total/2. )
              mean_shear_err_p = ntostr( sqrt( total(etanerrsum)/wsum_total^2 )/2. )
          ENDIF ELSE BEGIN 
              mean_shear_p = '0' 
              mean_shear_err_p = '0'
          ENDELSE 
          print,'Lens Group = ',ntostr(group+1)+'/'+ntostr(nstep),$
            '  Mean shear = ',mean_shear_p,' +/- ',mean_shear_err_p
          FOR gi=0L, step-1 DO BEGIN ;Loop over lenses in group
              IF (gi MOD 10) EQ 0 THEN  print,'.',format='(a,$)'
              index = ii[gi]
              ceneta = leta[index] ;lens center
              cenlam  = llambda[index]

              sig_crit = sigmacrit[index]
              sig_inv = sigcritinv[index]
              IF sig_inv EQ -1000. THEN BEGIN 
                  print,'What!'
                  return
              ENDIF 
              angmax_i = angmax[index]

              ;; choose sources around this lens
              maxlam2 = cenlam+angmax_i
              minlam2 = cenlam-angmax_i

              zobjshear_lambda_get_sources, scat[wsrc].lambda, minlam2, maxlam2, $
                twsrc, nwsrc2, $
                issouth=issouth

              ;; If there are any sources left, measure the shear
              IF nwsrc2 GT 1 THEN BEGIN 

                  ;; calculate angular distance (dec->lambda in distance calc)

                  wsrc2 = wsrc[twsrc]
                  tscateta = scat[wsrc2].eta
                  etadiff = (ceneta-tscateta)*d2r
                  cosetadiff = cos(etadiff)
                  
                  tcenlam = cenlam*d2r ;lens center in radians
                  sincenlam = sin(tcenlam)
                  coscenlam = cos(tcenlam)
                  
                  tscatlam = scat[wsrc2].lambda*d2r
                  sinscatlam = sin(tscatlam)
                  cosscatlam = cos(tscatlam)
;                  tscatlam = tscatlam*r2d
                  ;; Find distance in kpc
                  args = sincenlam*sinscatlam + coscenlam*cosscatlam*cosetadiff
                  warg = where(args LT -1., narg)
                  IF narg NE 0 THEN args[warg] = -1.
                  warg = where(args GT 1., narg)
                  IF narg NE 0 THEN args[warg] = 1.

                  R=acos( args )*DL[index]
                  
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
                          theta=atan(sin(etadiff[w]), $
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

                          ;; also weighting by 1/sigmacrit
                          wts_ssh = 1./( vint + scat[wsrc2[w]].(momerr[0])^2)
                          wts = wts_ssh*sig_inv^2
                              
                          tmprsum[binnum] = total(R[w])
                          tmpnpsum[binnum] = npair[binnum]
                              
                          tmpetan[binnum] = total(e1prime*wts)
                          tmperad[binnum] = total(e2prime*wts)
                          tmptansig[binnum] = tmpetan[binnum]*sig_crit
                          tmpradsig[binnum] = tmperad[binnum]*sig_crit

                          tmpetanerr[binnum]=total(wts^2*e1prime^2)
                          tmperaderr[binnum]=total(wts^2*e2prime^2)
                          tmptansigerr[binnum] = tmpetanerr[binnum]*sig_crit^2
                          tmpradsigerr[binnum] = tmperaderr[binnum]*sig_crit^2

                          ;; Ssh weighted differently: no weight by sig_crit
                          tmpSsh = tmpSsh+total( wts_ssh*(1.-vint*wts_ssh*e1prime^2) )
                              
                          tmpwsum[binnum] = total(wts)
                          tmpwsum_ssh = tmpwsum_ssh + total(wts_ssh)

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
                          tansigsum[whist] = tansigsum[whist] + tmptansig[whist]
                          radsigsum[whist] = radsigsum[whist] + tmpradsig[whist]

                          etanerrsum[whist]=etanerrsum[whist]+tmpetanerr[whist]
                          eraderrsum[whist]=eraderrsum[whist]+tmperaderr[whist]
                          tansigerrsum[whist]=tansigerrsum[whist]+tmptansigerr[whist]
                          radsigerrsum[whist]=radsigerrsum[whist]+tmpradsigerr[whist]
                              
                          Sshsum = Sshsum+tmpSsh
                             
                          wsum[whist] = wsum[whist] + tmpwsum[whist]
                          wsum_ssh = wsum_ssh + tmpwsum_ssh

                          ;; copy in individual lens stuff
                          lensum[index].totpairs = total(tmpnpsum[whist])
                          lensum[index].npair[whist] = tmpnpsum[whist]
                          lensum[index].ie = ie
                          lensum[index].rmax_act[whist] = rmax_act_tmp[whist]
                          lensum[index].rsum[whist] = tmprsum[whist]
                          lensum[index].etansum[whist] = tmpetan[whist]
                          lensum[index].eradsum[whist] = tmperad[whist]
                          lensum[index].tansigsum[whist] = tmptansig[whist]
                          lensum[index].radsigsum[whist] =  tmpradsig[whist]

                          lensum[index].etanerrsum[whist] = tmpetanerr[whist]
                          lensum[index].eraderrsum[whist] = tmperaderr[whist]
                          lensum[index].tansigerrsum[whist] = tmptansigerr[whist]
                          lensum[index].radsigerrsum[whist] = tmpradsigerr[whist]

                          lensum[index].sshsum = tmpSsh

                          lensum[index].wsum[whist] = tmpwsum[whist]
                          lensum[index].wsum_ssh = tmpwsum_ssh

                      ENDIF ELSE BEGIN 
                          ng=0
                      ENDELSE 
                      ;; reinitialize the arrays
                      xrel[wg] = 0. & yrel[wg] = 0.
                      tmpetan[whist] = 0.
                      tmperad[whist] = 0.
                      tmptansig[whist] = 0.
                      tmpradsig[whist] = 0.

                      tmpetanerr[whist] = 0.
                      tmperaderr[whist] = 0.
                      tmptansigerr[whist] = 0.
                      tmpradsigerr[whist] = 0.

                      tmpwsum[whist] = 0.
                      tmprsum[whist] = 0.
                      tmpnpsum[whist] = 0
                      tmpwsum_ssh = 0.
                      tmpSsh = 0L
                  ENDIF ;ELSE print,'Dohh!',format='(a,$)'
                  ;; One final check.  If not passed, remember that this 
                  ;; lens wasn't used
                  IF ng EQ 0 THEN BEGIN 
                      print,'/',format='(a,$)'
                      indices[ ind[gi] ] = -1
                  ENDIF ELSE BEGIN
                      groupused=groupused+1 ;for group statistics
                      lamsum=lamsum+cenlam
                      etasum=etasum+ceneta
                  ENDELSE 
                  R = 0
                  etadiff=0      ;Free memory
              ENDIF ELSE BEGIN
                  print,'|',format='(a,$)'
                  indices[ ind[gi] ] = -1
              ENDELSE 

          ENDFOR 
          print
          xrel = 0 & yrel = 0
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; fill in group statistics
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          wsum_total = total(wsum)
          IF wsum_total GT 0.0 THEN BEGIN 
              mean_shear = total(etansum)/wsum_total/2.
              mean_shear_err = sqrt( total(etanerrsum)/wsum_total^2 )/2.
          endif else begin 
              mean_shear=0.0
              mean_shear_err = 0.0
          endelse 
          
          groupstruct.nlens[group] = groupused
          groupstruct.shear[group] = mean_shear
          groupstruct.shearerr[group] = mean_shear_err
          groupstruct.mean_lambda[group] = lamsum/groupused
          groupstruct.mean_eta[group] = etasum/groupused
      ENDIF ELSE BEGIN 
          print,'Two-'
          indices[ ind ] = -1
      ENDELSE 
      metadiff = 0               ;Free memory
      mdist = 0
  ENDFOR 
  ;; Remove unused lenses
  wbad = where(indices EQ -1, nwbad)
  IF nwbad NE 0 THEN remove, wbad, wlens
  lensused = n_elements(wlens)
  
  wgood = wlens
  
  lensw = sigcritinv[wlens]^2*lensum[wlens].totpairs
  lenswsum = total( lensw )
  zsum = total(lcat[wlens].(wz[0])*lensw)
  zmean = zsum/lenswsum
  print
  print,'zmean = ',zmean
  print

  scritinvsum = total( sigcritinv[wlens]*lensw )
  meanscritinv = scritinvsum/lenswsum
  print
  print,'meanscritinv: ',meanscritinv
  print

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

      shstruct.sigma[i] = tansigsum[i]/wsum[i]/2.
      shstruct.sigmaerr[i] = sqrt( tansigerrsum[i]/wsum[i]^2 )/2.
      shstruct.orthosig[i] = radsigsum[i]/wsum[i]/2.
      shstruct.orthosigerr[i] = sqrt( radsigerrsum[i]/wsum[i]^2 )/2.
      shstruct.npair[i] = npsum[i]

      ;; Now totals within each rmax_act[i]
      shstruct.tshear[i] = total(etansum[0:i])/total(wsum[0:i])/2.
      shstruct.tshearerr[i] = sqrt( total(etanerrsum[0:i])/total(wsum[0:i])^2 )/2.
      shstruct.tortho[i] = total(eradsum[0:i])/total(wsum[0:i])/2.
      shstruct.torthoerr[i] = sqrt( total(eraderrsum[0:i])/total(wsum[0:i])^2 )/2.

      shstruct.tsigma[i] = total(tansigsum[0:i])/total(wsum[0:i])/2.
      shstruct.tsigmaerr[i] = sqrt( total(tansigerrsum[0:i])/total(wsum[0:i])^2 )/2.
      shstruct.torthosig[i] = total(radsigsum[0:i])/total(wsum[0:i])/2.
      shstruct.torthosigerr[i] = sqrt( total(radsigerrsum[0:i])/total(wsum[0:i])^2 )/2.
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

      sumstruct.tansigsum[i] = tansigsum[i]
      sumstruct.tansigerrsum[i] = tansigerrsum[i]
      sumstruct.radsigsum[i] = radsigsum[i]
      sumstruct.radsigerrsum[i] = radsigerrsum[i]

      sumstruct.wsum[i] = wsum[i]
      sumstruct.wsum_ssh = wsum_ssh
      sumstruct.npair[i] = npsum[i]

      totpairs = totpairs + npsum[i]
  ENDFOR 
  Ssh = Sshsum/wsum_ssh

  shstruct.nlenses = lensused
  shstruct.totpairs = totpairs
  shstruct.binsize = binsize
  shstruct.rmin = rmin
  shstruct.rmax = rmax
  shstruct.rmax_act = rmax_act
  shstruct.ssh = ssh
  shstruct.zmean = zmean
  shstruct.meanscritinv = meanscritinv
  shstruct.h = h

  sumstruct.nlenses = lensused
  sumstruct.totpairs = totpairs
  sumstruct.binsize = binsize
  sumstruct.rmin = rmin
  sumstruct.rmax = rmax
  sumstruct.rmax_act = rmax_act
  sumstruct.sshsum = Sshsum
  sumstruct.lenswsum = lenswsum
  sumstruct.zsum = zsum
  sumstruct.scritinvsum = scritinvsum
  sumstruct.h = h

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Ouput the data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Shear file for these lenses
  shhdr = zshhdr(shstruct)      
  mwrfits2, shstruct, datfile, shhdr, /create, /destroy

  ;; Sum file for these lenses
  sumhdr = zsumhdr(sumstruct)
  lhdr = sumhdr
  mwrfits2, sumstruct, sumfile, sumhdr, /create, /destroy

  ;; File with structure for each lens used
  lensum = temporary(lensum[wlens])
  lensum.zindex = lindgen(lensused)
  mwrfits2, lensum, lensumfile, lhdr, /create, /destroy

  ;; File containing ra,dec,redshift for each used lens
  zs = create_struct('z', 0., $
                     'sigcritinv', 0., $
                     'lambda', double(0.), $
                     'eta', double(0.) )
  zstruct = replicate(zs, lensused)
  zstruct.z = lcat[wlens].(wz[0])
  zstruct.sigcritinv = sigcritinv[wlens]
  zstruct.lambda = llambda[wlens]
  zstruct.eta = leta[wlens]

  mwrfits2, zstruct, zfile, /create, /destroy

  ;; file containing group stats
  mwrfits2, groupstruct, groupfile, /create, /destroy

  step=oldstep
  usecat = lcat[wlens]

  ptime, systime(1)-time
  return
END 
