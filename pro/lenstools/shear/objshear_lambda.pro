
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    OBJSHEAR_LAMBDA
;       
; PURPOSE:
;    Measure the shear and density contrast for a set of lenses. Structure
;    must have a redshift tags.
;
; CALLING SEQUENCE:
;    objshear_lambda, stripe, lenscat, scat, clr, rmin, rmax, binsize, $
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
;    lenscat: lens catalog. Must have either (ra,dec) or (lambda,eta)
;    scat: source catalog. Must contain (lambda,eta), and must be sorted
;          by lambda.
;    clr: the sdss bandpass
;    rmin: minimum radius in arcsec
;    rmax: maximum radius in arcsec
;    binsize: size of bins in arcsec
;
; OPTIONAL INPUTS:
;    step: how many lenses to use in a group. default is 300
;    addstr: string to add to front of output file names.
;    outdir: output directory
;    maxe: maximum allowed ellipticity of source galaxy distribution.

; KEYWORD PARAMETERS:
;       
; OUTPUTS: 
;    datafile, sumfile, lensumfile
;
; OPTIONAL OUTPUTS:
;    wgood: the indices of used lenses.
;    usecat: the final lenses.
;    datafile, sumfile, zfile, lensumfile: the names of output files
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
;    COPY_STRUCT
;    OBJSHEAR_LAMBDA_GET_SOURCES (in this file)
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

PRO objshear_lambda_get_sources, lambda, minlam, maxlam, w, num, issouth=issouth

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


PRO objshear_lambda, stripe, lenscat, scat, clr, rmin, rmax, binsize, $
                     step=step, addstr=addstr, $
                     outdir=outdir, $
                     wgood=wgood, $
                     usecat=usecat, $
                     datfile=datfile, sumfile=sumfile, $
                     lensumfile=lensumfile, lensew=lensew, $
                     maxe=maxe

  IF n_params() LT 7 THEN BEGIN
      print,'-Syntax: objshear, stripe, lenscat, scat, clr, rmin, rmax, binsize,'
      print,'  step=step, addstr=addstr, '
      print,'  outdir=outdir, '
      print,'  wgood=wgood, '
      print,'  usecat=usecat, '
      print,'  datfile=datfile, sumfile=sumfile,'
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

  ;; Size by which to group things  300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = 'gal_gal'
  IF keyword_set(lensew) THEN addstr = 'lensew_'+addstr

  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss4/data1/esheldon/TMP/'

  stripestr = '_stripe'+ntostr(stripe)
  clr_str = '_'+colors[clr]

  d2r = !dpi/180.0d0            ;Change degrees to radians.
  r2d = 180./!dpi               ;radians to degrees

  radmax = rmax/3600.           ;rmax in degrees
  
  nbin = long( (rmax - rmin)/binsize ) + 1 

  tol = 1.e-6
  angtol_u = 1.- tol
  angtol_l = -1.+tol

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
          message,'No valid moment uncertainty tag'
      ENDIF 
  ENDIF 
  
  IF keyword_set(lensew) THEN BEGIN 
      IF NOT tag_exist(lenscat,'momerr',index=lmomerr) THEN BEGIN
          IF NOT tag_exist(lenscat,'uncert',index=lmomerr) THEN BEGIN 
              print
              message,'/lensew Set But No Valid Moment Uncertainty Tag in lenscat'
          ENDIF 
      ENDIF 
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; declare some arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  arrval = fltarr(nbin)

  shstruct = shstruct(arrval)
  sumstruct = sumstruct(arrval)
  lensumstruct = lensumstruct(arrval)
  lensumstruct = create_struct(lenscat[0], lensumstruct)

  lensum = replicate(lensumstruct, n_elements(lenscat) )

  ;;;;;;;;;;;;;;;;;;;;
  ;; temporary
  ;;;;;;;;;;;;;;;;;;;;

  etansum    = arrval
  etanerrsum = arrval
  wsum       = arrval

  npair      = arrval
  Sshsum     = 0.               ;Scalar
  rmax_act   = arrval
  rmax_act_tmp = arrval

  tmpetan = arrval
  tmpetanerr = arrval
  tmperad = arrval
  tmperaderr = arrval

  tmpwsum = arrval
  tmprsum = arrval
  tmpnpsum = arrval
  tmpSsh = 0.

  ;; halo shape stuff
  tqrsum = arrval
  tqisum = arrval
  tqrerrsum = arrval
  tqierrsum = arrval

  tqradrsum = arrval
  tqradisum = arrval
  tqradrerrsum = arrval
  tqradierrsum = arrval

  tetan1sum = arrval
  tetan1errsum = arrval
  tetan2sum = arrval
  tetan2errsum = arrval

  terad1sum = arrval
  terad1errsum = arrval
  terad2sum = arrval
  terad2errsum = arrval

  tmpw1sum = arrval
  tmpw2sum = arrval



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  prename = outdir+addstr+stripestr+clr_str
  
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

  lcat = lenscat                ;make copy which will be sorted by lambda
  ninit = n_elements(lcat)
  wlens = lindgen(ninit)

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

  ;; copy in the stuff we need into lensum struct
  copy_struct, lcat, lensum

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; print some stuff
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nlens = n_elements(wlens)
  print
  print,'FIRSTLAMBDA = ',FIRSTLAMBDA
  print,'LASTLAMBDA = ',LASTLAMBDA
  print,'Using '+ntostr(nlens)+' lenses'
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

      objshear_lambda_get_sources, scat.lambda, minlam1, maxlam1, wsrc, nwsrc, $
        issouth=issouth

      IF nwsrc NE 0 THEN BEGIN 
          
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
              ceneta = leta[index] ;lens center
              cenlam  = llambda[index]

              ;; choose sources around this lens
              maxlam2 = cenlam+radmax
              minlam2 = cenlam-radmax

              objshear_lambda_get_sources, scat[wsrc].lambda, minlam2, maxlam2, $
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
                  ;; Find distance in arcseconds
                  args = sincenlam*sinscatlam + coscenlam*cosscatlam*cosetadiff
                  warg = where(args LT -1., narg)
                  IF narg NE 0 THEN args[warg] = -1.
                  warg = where(args GT 1., narg)
                  IF narg NE 0 THEN args[warg] = 1.

                  R=acos( args )*r2d*3600.

                  blah=where(R NE R,nblah)
                  IF nblah NE 0 THEN BEGIN
                      ;; don't use these errors
                      print
                      message,'Fixing!!',/inform
                      print
                      R[blah] = 10.^10
                  ENDIF 

                  hist=histogram(R, binsize=binsize, min=rmin, $
                                 max=rmax,rever=rev_ind)
                  
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

                          xrel[w]=R[w]*sin(theta)
                          yrel[w]=R[w]*cos(theta)

                          diffsq=xrel[w]^2 - yrel[w]^2
                          xy=xrel[w]*yrel[w]

                          e1prime=-(scat[wsrc2[w]].e1*diffsq + $
                                    scat[wsrc2[w]].e2*2.*xy  )/R[w]^2
                          e2prime= (scat[wsrc2[w]].e1*2.*xy - $
                                    scat[wsrc2[w]].e2*diffsq )/R[w]^2

                          IF keyword_set(lensew) THEN BEGIN 
                              wts = 1./( vint + $
                                         lcat[index].(lmomerr[0])^2 + $
                                         scat[wsrc2[w]].(momerr[0])^2 )
                          ENDIF ELSE BEGIN 
                              wts = 1./( vint + scat[wsrc2[w]].(momerr[0])^2)
                          ENDELSE 

                          tmprsum[binnum] = total(R[w])
                          tmpnpsum[binnum] = npair[binnum]
                              
                          tmpetan[binnum] = total(e1prime*wts)
                          tmperad[binnum] = total(e2prime*wts)

                          tmpetanerr[binnum]=total(wts^2*e1prime^2)
                          tmperaderr[binnum]=total(wts^2*e2prime^2)

                          tmpSsh = tmpSsh+total( wts*(1.-vint*wts*e1prime^2) )
                          tmpwsum[binnum] = total(wts)

                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                          ;; halo shape stuff
                          ;; Must be lcat chosen to have e1/e2 in this band!
                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                          ltheta=0.5*atan(lcat[index].e2,lcat[index].e1)
                          newx =  xrel[w]*cos(ltheta) + yrel[w]*sin(ltheta)
                          newy = -xrel[w]*sin(ltheta) + yrel[w]*cos(ltheta)

                          newe1 =  scat[wsrc2[w]].e1*cos(2*ltheta)+$
                            scat[wsrc2[w]].e2*sin(2*ltheta)
                          newe2 = -scat[wsrc2[w]].e1*sin(2*ltheta)+$
                            scat[wsrc2[w]].e2*cos(2*ltheta)
                          e1prime = -( newe1*(newx^2-newy^2)+$
                                       newe2*2.*newx*newy )/R[w]^2
                          e2prime =  ( newe1*2.*newx*newy  -  $
                                       newe2*(newx^2 - newy^2) )/R[w]^2
                          tqr = (newx^2-newy^2)/R[w]^2*e1prime 
                          tqi = 2*newx*newy/R[w]^2*e1prime
                          tqradr = -2*newx*newy/R[w]^2*e2prime
                          tqradi =  (newx^2-newy^2)/R[w]^2*e2prime 

                          tqrsum[binnum] = total(tqr*wts)
                          tqisum[binnum] = total(tqi*wts)
                          tqrerrsum[binnum] = total(tqr^2*wts^2)
                          tqierrsum[binnum] = total(tqi^2*wts^2)

                          
                          tqradrsum[binnum] = total(tqradr*wts)
                          tqradisum[binnum] = total(tqradi*wts)
                          tqradrerrsum[binnum] = total(tqradr^2*wts^2)
                          tqradierrsum[binnum] = total(tqradi^2*wts^2)

                          ;; now in 2 half spaces
                          half1 = where( abs(newy) GT abs(newx), num1)
                          half2 = where( abs(newx) GE abs(newy), num2)

                          IF num1 NE 0 THEN BEGIN 
                              tetan1sum[binnum]=total(e1prime[half1]*wts[half1])
                              tetan1errsum[binnum]=$
                                total(e1prime[half1]^2*wts[half1]^2)

                              terad1sum[binnum]=total(e2prime[half1]*wts[half1])
                              terad1errsum[binnum]=$
                                total(e2prime[half1]^2*wts[half1]^2)
                              tmpw1sum[binnum] = total(wts[half1])
                          ENDIF 
                          IF num2 NE 0 THEN BEGIN 
                              tetan2sum[binnum]=total(e1prime[half2]*wts[half2])
                              tetan2errsum[binnum]=$
                                total(e1prime[half2]^2*wts[half2]^2)

                              terad2sum[binnum]=total(e2prime[half2]*wts[half2])
                              terad2errsum[binnum]=$
                                total(e2prime[half2]^2*wts[half2]^2)
                              tmpw2sum[binnum] = total(wts[half2])
                          ENDIF 

                          ;; free memory
                          setzero,diffsq,xy,theta,ltheta,$
                            e1prime,e2prime,wts,w,tqr,tqi

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

                          etansum[whist] = etansum[whist] + tmpetan[whist]
                          etanerrsum[whist] = etanerrsum[whist]+tmpetanerr[whist]
                          wsum[whist] = wsum[whist] + tmpwsum[whist]

                          ;; copy in individual lens stuff
                          lensum[index].totpairs = total(tmpnpsum[whist])
                          lensum[index].npair[whist] = tmpnpsum[whist]
                          lensum[index].ie = ie
                          lensum[index].rmax_act[whist] = rmax_act_tmp[whist]
                          lensum[index].rsum[whist] = tmprsum[whist]

                          lensum[index].etansum[whist] = tmpetan[whist]
                          lensum[index].eradsum[whist] = tmperad[whist]
                          lensum[index].etanerrsum[whist] = tmpetanerr[whist]
                          lensum[index].eraderrsum[whist] = tmperaderr[whist]

                          lensum[index].qrsum[whist] = tqrsum[whist]
                          lensum[index].qisum[whist] = tqisum[whist]
                          lensum[index].qrerrsum[whist] = tqrerrsum[whist]
                          lensum[index].qierrsum[whist] = tqierrsum[whist]

                          lensum[index].qradrsum[whist] = tqradrsum[whist]
                          lensum[index].qradisum[whist] = tqradisum[whist]
                          lensum[index].qradrerrsum[whist] = tqradrerrsum[whist]
                          lensum[index].qradierrsum[whist] = tqradierrsum[whist]

                          lensum[index].etan1sum[whist] = tetan1sum[whist]
                          lensum[index].etan1errsum[whist] = tetan1errsum[whist]
                          lensum[index].etan2sum[whist] = tetan2sum[whist]
                          lensum[index].etan2errsum[whist] = tetan2errsum[whist]
  
                          lensum[index].erad1sum[whist] = terad1sum[whist]
                          lensum[index].erad1errsum[whist] = terad1errsum[whist]
                          lensum[index].erad2sum[whist] = terad2sum[whist]
                          lensum[index].erad2errsum[whist] = terad2errsum[whist]

                          lensum[index].w1sum[whist] = tmpw1sum[whist]
                          lensum[index].w2sum[whist] = tmpw2sum[whist]

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

                      tqrsum[whist] = 0.
                      tqisum[whist] = 0.
                      tqrerrsum[whist] = 0.
                      tqierrsum[whist] = 0.

                      tqradrsum[whist] = 0.
                      tqradisum[whist] = 0.
                      tqradrerrsum[whist] = 0.
                      tqradierrsum[whist] = 0.

                      tmpwsum[whist] = 0.
                      tmprsum[whist] = 0.
                      tmpnpsum[whist] = 0
                      tmpSsh = 0L
                  ENDIF ;ELSE print,'Dohh!',format='(a,$)'
                  ;; One final check.  If not passed, remember that this 
                  ;; lens wasn't used
                  IF ng EQ 0 THEN BEGIN 
                      print,'/',format='(a,$)'
                      indices[ ind[gi] ] = -1
                  ENDIF 

                  R = 0
                  etadiff=0      ;Free memory
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
      metadiff = 0               ;Free memory
      mdist = 0
  ENDFOR 
  ;; Remove unused lenses
  wbad = where(indices EQ -1, nwbad)
  IF nwbad NE 0 THEN remove, wbad, wlens
  lensused = n_elements(wlens)
  
  lensw = lensum[wlens].totpairs
  lenswsum = total( lensw )

  wgood = wlens
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find averages from sums
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Finally used ',ntostr(lensused),'/',ntostr(nlens),' Lenses'

  lensum = temporary(lensum[wlens])
  combine_lensum, lensum, binsize, rmin, rmax, sumstruct, shstruct

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Ouput the data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Shear file for these lenses
  shhdr = shhdr(shstruct)      
  mwrfits, shstruct, datfile, shhdr, /create

  ;; Sum file for these lenses
  sumhdr = sumhdr(sumstruct)
  lhdr = sumhdr
  mwrfits, sumstruct, sumfile, sumhdr, /create

  ;; File with structure for each lens used
  mwrfits2, lensum, lensumfile, lhdr, /create, /destroy

  step=oldstep
  usecat = lcat[wlens]

  ptime, systime(1)-time
  return
END 
