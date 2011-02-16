PRO rand_rtag, scat, clr, rmin, zrand_in, radin, ndup, $
                cosmo=cosmo, $
                run1=run1, run2=run2, $
                step=step, addstr=addstr, $
                outdir=outdir, $
                wgood=wgood, $
                check=check, $
                datfile=datfile, sumfile=sumfile, zfile=zfile, $
                lensumfile=lensumfile, $
                maxe=maxe

  IF n_params() LT 6 THEN BEGIN
      print,'-Syntax: rand_rtag, scat, clr, rmin, rmax, binsize, zrand,, ndup'
      print,'  run1=run1, run2=run2,'
      print,'  step=step, addstr=addstr, '
      print,'  outdir=outdir, '
      print,'  title=title, wgood=wgood, '
      print,'  check=check, '
      print,'  datfile=datfile, sumfile=sumfile, zfile=zfile,'
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
  IF n_elements(cosmo) EQ 0 THEN cosmo=1
  IF cosmo EQ 1 THEN omegamat=1.0 ELSE omegamat=0.3

  ;; Size by which to group things  300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L ELSE step = long(step)
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = ''
  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss4/data1/esheldon/TMP/'
  IF NOT keyword_set(check) THEN check=0

  IF n_elements(run1) EQ 0 THEN r1str='' ELSE r1str = '_'+ntostr(run1)
  IF n_elements(run2) EQ 0 THEN r2str='' ELSE r2str = '_'+ntostr(run2)
  clr_str = '_'+colors[clr]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; maxit is number of retries until random point
  ;; passes symmetry check
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  maxit = 20

  d2r = !dpi/180.0d0            ;Change degrees to radians.
  r2d = 180./!dpi               ;radians to degrees

;  nbin = long( (rmax - rmin)/binsize ) + 1 

  print
  print,'-------------------------------------------------------------'
  print,'Rmin: ',rmin
;  print,'Rmax: ',rmax
;  print,'Binsize: ',binsize
  print,'Step: ',step
;  print,'Using ',ntostr(nbin),' bins between ',ntostr(rmin), $
;        ' and ',ntostr(rmax),' kpc'

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

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; duplicate the input z distribution
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Duplicating: ',ntostr(ndup)

  nrand_in = n_elements(zrand_in)
  nrand = nrand_in*ndup

  FOR dd=0L, ndup-1 DO BEGIN 
      IF dd EQ 0 THEN BEGIN 
          zrand = zrand_in 
          zindex = lindgen(nrand_in)
          rad = radin
      ENDIF ELSE BEGIN 
          zrand = [zrand, zrand_in]
          zindex = [zindex, lindgen(nrand_in)]
          rad = [rad, radin]
      ENDELSE 
  ENDFOR 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; declare some arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  arrval = fltarr(nbin)
  arrval = 0.0

  shstruct = zshstruct(arrval)
  sumstruct = zsumstruct(arrval)
  lensumstruct = zlensumstruct(arrval)
  lensumstruct = create_struct('z', 0., lensumstruct)

  lensum = replicate(lensumstruct, nrand )
  lensum.zindex = zindex

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
  rmax_act   = arrval
  rmax_act_tmp = arrval
  Sshsum     = 0.               ;Scalar
  wsum_ssh   = 0.

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
  tmpSsh = 0.
  tmpwsum_ssh = 0.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF addstr EQ '' THEN prename = outdir+addstr+'rand_rtag'+r1str+r2str+clr_str $
  ELSE prename = outdir+addstr+'_rand_rtag'+r1str+r2str+clr_str

  datfile = prename + '_N1.fit'
  sumfile = prename + '_sum_N1.fit'
  zfile   = prename + '_z_N1.fit'
  lensumfile = prename + '_lensum_N1.fit'

  WHILE exist(datfile) DO BEGIN
      datfile = newname(datfile)
      sumfile = newname(sumfile)
      zfile = newname(zfile)
      lensumfile = newname(lensumfile)
  ENDWHILE 
  print
  print,'Dat file: ',datfile
  print,'Sum file: ',sumfile
  print,'Lens sum file: ',lensumfile
  print,'Redshift file: ',zfile
  print


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; source cat
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; they are already sorted
  nsource = n_elements(scat)
  wsource = lindgen(nsource)

  FIRSTRA = scat[0].ra
  LASTRA  = scat[nsource-1].ra

  ;; check if runs cross ra=0.
  IF FIRSTRA GT LASTRA THEN flag=1 ELSE flag=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Lens cat: Set up sigma crit and Dlens. 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  h=1.
  print,'Using h = ',h
  print,'Calculating 1/sigma_crit'
  print

  wlens = lindgen(nrand)

  sigcritinv = sdss_sigma_crit(clr, zrand, wgood=wgood, cosmo=cosmo)
  sigmacrit = 1./sigcritinv
  wlens = wlens[wgood]
  lensum.scritinv = sigcritinv

  nlens = n_elements(wlens)

  DL = angdist_lambda( zrand, h=h, omegamat=omegamat)*1000. ;Convert from Mpc to kpc
  angmax = rad*1000./DL*180./!pi

  maxdec = max(scat.dec )
  mindec = min(scat.dec )

  lra  = randomu(seed, nrand)
;  ldec = randomu(seed, nrand)
  ldec = fltarr(nrand)
  s = sort(lra)
  lra  = temporary(lra[s])
;  ldec = temporary(ldec[s])

  
  IF flag THEN BEGIN
      bbra = FIRSTRA
      llra = 360. + LASTRA

      lra = arrscl(lra, bbra, llra, arrmin=0., arrmax=1.)
;      ldec= arrscl(ldec, mindec, maxdec, arrmin=0.,arrmax=1.)

      wl = where(lra GT 360., nwl)
      IF nwl NE 0 THEN BEGIN
          lra[wl] = lra[wl] - 360.
      ENDIF 

  ENDIF ELSE BEGIN 
      lra = arrscl(lra,  FIRSTRA,  LASTRA,  arrmin=0.,arrmax=1.)
;      ldec= arrscl(ldec, mindec, maxdec, arrmin=0.,arrmax=1.)
  ENDELSE 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Throw out by _edge_ if requested (rather than just symmety of background dist.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF check THEN BEGIN           ;May just want to find the good lenses.
      wgood = wlens
      return
  ENDIF 

  nlens = n_elements(wlens)
  print,'Threw out ',ntostr(nlens-nrand),' edge lenses'
  print
  print,'LASTRA = ',LASTRA
  print,'FIRSTRA = ',FIRSTRA
  print,'Finally FLAG = ',FLAG
  print,'Using '+ntostr(nlens)+'/'+ntostr(nrand)+' lenses'
  print
  print,'-------------------------------------------------------------'
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
      ;; they are sorted by ra, but it could go over ra=0
      
      angmax2 = max(angmax[ii])
      
      maxii = wlens[ max(ind) ] & minii = wlens[ min(ind) ]
;      maxii = max(ii) & minii = min(ii)

      w=where(lra[ii] GT lra[maxii] OR lra[ii] LT lra[minii], nw)
      IF nw NE 0 THEN BEGIN 
          print,'What the hell!'
          return
      ENDIF 

      tmpmaxra = lra[maxii]+angmax2
      tmpminra = lra[minii]-angmax2

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

          binary_search, scat.ra, ra1, i1, /round
          binary_search, scat.ra, ra2, i2, /round
          IF i1 EQ -1 THEN i1 = 0
          IF i2 EQ -1 THEN i2 = nsource-1
          wsrc = wsource[i1:i2]
          nwsrc = n_elements(wsrc)

      ENDELSE 

      IF nwsrc NE 0 THEN BEGIN 
          
          ;; generate dec's
          tmpmaxdec = max(scat[wsrc].dec)
          tmpmindec = min(scat[wsrc].dec)
          ldec[ii] = arrscl( randomu(seed, step), tmpmindec, tmpmaxdec, arrmin=0., arrmax=1.)

          xrel = dblarr(nwsrc)
          yrel = xrel
          print,'Lens Group = ',ntostr(group+1)+'/'+ntostr(nstep)
          FOR gi=0L, step-1 DO BEGIN ;Loop over lenses in group
              IF (gi MOD 10) EQ 0 THEN  print,'.',format='(a,$)'
              index = ii[gi]
              cendec = ldec[index]
              cenra  = lra[index]

              sig_crit = sigmacrit[index]
              sig_inv = sigcritinv[index]
              IF sig_inv EQ -1000. THEN BEGIN 
                  print,'What!'
                  return
              ENDIF 
             
              angmax_i = angmax[index]
              rmax = rad[index]*1000.

              ;; choose sources around this lens
              tmpmaxra = cenra+angmax_i
              tmpminra = cenra-angmax_i
              diff1 = 360. - tmpmaxra
              IF flag NE 0 THEN BEGIN 
                  ;; Quick Fix
                  wsrc2=wsrc
                  nwsrc2=nwsrc
              ENDIF ELSE BEGIN 
                  ra1 = tmpminra
                  ra2 = tmpmaxra
                  binary_search, scat[wsrc].ra, ra1, i1, /round
                  binary_search, scat[wsrc].ra, ra2, i2, /round
                  IF i1 EQ -1 THEN i1 = 0
                  IF i2 EQ -1 THEN i2 = nwsrc-1
                  wsrc2 = wsrc[i1:i2]
                  nwsrc2 = n_elements(wsrc2)
              ENDELSE 

              ;; If there are any sources left, measure the shear
              IF nwsrc2 NE 0 THEN BEGIN 

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; Try to redo a few times until passes checks
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  redo=0
                  WHILE (redo NE -1) AND (redo LT maxit) DO BEGIN 
                      redo = redo+1

                      radiff = (cenra-scat[wsrc2].ra)*d2r
                      cosradiff = cos(radiff)
                  
                      tcendec = cendec*d2r
                      tcenra = cenra*d2r
                      sincendec=sin(tcendec)
                      coscendec=cos(tcendec)
                  
                      sinscatdec = sin(scat[wsrc2].dec*d2r)
                      cosscatdec = cos(scat[wsrc2].dec*d2r)
                      
                      ;; Find distance in kpc
                      args = sincendec*sinscatdec +coscendec*cosscatdec*cosradiff
                      warg = where(args LT -1., narg)
                      IF narg NE 0 THEN args[warg] = -1.
                      warg = where(args GT 1., narg)
                      IF narg NE 0 THEN args[warg] = 1.

                      R=acos( args )*DL[index]
                      wrf = where(R LE rmax AND R GE rmin, npair)
                      ng=0

                      ;; Check if there are any in this annulus rmin-rmax
                      IF npair NE 0 THEN BEGIN 
                          
                          rmax_act_tmp = max(R[wrf])
                          rmax_act = max( [rmax_act, rmax_act_tmp] )
                              
                          theta=atan(sin(radiff[wrf]), $
                                     (sincendec*cosradiff[wrf] - $
                                      coscendec*sinscatdec[wrf]/cosscatdec[wrf]) )-!dpi
                          xrel[wrf] = R[wrf]*cos(theta)
                          yrel[wrf] = R[wrf]*sin(theta)
                          diffsq=xrel[wrf]^2 - yrel[wrf]^2
                          xy=xrel[wrf]*yrel[wrf]
                              
                          e1prime=-(scat[wsrc2[wrf]].e1*diffsq + $
                                    scat[wsrc2[wrf]].e2*2.*xy  )/R[wrf]^2
                          e2prime= (scat[wsrc2[wrf]].e1*2.*xy - $
                                    scat[wsrc2[wrf]].e2*diffsq )/R[wrf]^2
                          
                          ;; also weighting by 1/sigmacrit
                          wts_ssh = 1./( vint + scat[wsrc2[wrf]].(momerr[0])^2)
                          wts = wts_ssh*sig_inv^2
                          
                          tmprsum = total(R[wrf])
                          tmpnpsum = npair
                          
                          tmpetan = total(e1prime*wts)
                          tmperad = total(e2prime*wts)
                          tmptansig = tmpetan*sig_crit
                          tmpradsig = tmperad*sig_crit

                          tmpetanerr=total(wts^2*e1prime^2)
                          tmperaderr=total(wts^2*e2prime^2)
                          tmptansigerr = tmpetanerr*sig_crit^2
                          tmpradsigerr = tmperaderr*sig_crit^2
                          
                          ;; Ssh weighted differently: no weight by sig_crit
                          tmpSsh = tmpSsh+total( wts_ssh*(1.-vint*wts_ssh*e1prime^2) )
                              
                          tmpwsum = total(wts)
                          tmpwsum_ssh = tmpwsum_ssh + total(wts_ssh)

                          ;; free memory
                          theta = 0 & e1prime=0 & e2prime=0 & wts=0 & w=0
                          diffsq = 0 & xy=0

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
                              npsum = npsum + tmpnpsum
                              rsum = rsum + tmprsum
                              etansum = etansum + tmpetan
                              eradsum = eradsum + tmperad
                              tansigsum = tansigsum + tmptansig
                              radsigsum = radsigsum + tmpradsig

                              etanerrsum=etanerrsum+tmpetanerr
                              eraderrsum=eraderrsum+tmperaderr
                              tansigerrsum=tansigerrsum+tmptansigerr
                              radsigerrsum=radsigerrsum+tmpradsigerr
                              
                              Sshsum = Sshsum+tmpSsh
                              
                              wsum = wsum + tmpwsum
                              wsum_ssh = wsum_ssh + tmpwsum_ssh

                              ;; copy in individual lens stuff
                              lensum[index].totpairs = total(tmpnpsum)
                              lensum[index].npair = tmpnpsum
                              lensum[index].ie = ie
                              lensum[index].rsum = tmprsum
                              lensum[index].rmax_act = rmax_act_tmp
                              lensum[index].etansum = tmpetan
                              lensum[index].eradsum = tmperad
                              lensum[index].tansigsum = tmptansig
                              lensum[index].radsigsum =  tmpradsig

                              lensum[index].etanerrsum = tmpetanerr
                              lensum[index].eraderrsum = tmperaderr
                              lensum[index].tansigerrsum = tmptansigerr
                              lensum[index].radsigerrsum = tmpradsigerr

                              lensum[index].sshsum = tmpSsh

                              lensum[index].wsum = tmpwsum
                              lensum[index].wsum_ssh = tmpwsum_ssh

                              
                          ENDIF ELSE BEGIN 
                              ng=0
                          ENDELSE 
                          ;; reinitialize the arrays
                          xrel[wg] = 0. & yrel[wg] = 0.
                          tmpetan = 0.
                          tmperad = 0.
                          tmptansig = 0.
                          tmpradsig = 0.

                          tmpetanerr = 0.
                          tmperaderr = 0.
                          tmptansigerr = 0.
                          tmpradsigerr = 0.

                          tmpwsum = 0.
                          tmprsum = 0.
                          tmpnpsum = 0
                          tmpwsum_ssh = 0.
                          tmpSsh = 0L
                      ENDIF 
                      ;; One final check.  If not passed, remember that this lens wasn't used
                      IF ng EQ 0 THEN BEGIN 
                          qdec = randomu(seed, 10)
                          qdec = (arrscl( qdec, tmpmindec, tmpmaxdec, arrmin=0., arrmax=1.))[0]
                          ldec[index] = qdec
                          cendec = qdec
                      ENDIF ELSE redo=-1
                      R = 0
                      radiff=0  ;Free memory
                  ENDWHILE 
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

  lensw = sigcritinv[wlens]^2*lensum[wlens].totpairs
  lenswsum = total( lensw )
  zsum = total(zrand[wlens]*lensw)
  zmean = zsum/lenswsum
  print
  print,'zmean = ',zmean
  print

  scritinvsum = total( sigcritinv[wlens]*lensw )
  meanscritinv = scritinvsum/lenswsum
  print
  print,'meanscritinv: ',meanscritinv
  print

  ;; for this, we set rmax_act to the average rtag value
  rtagsum = total( rad[wlens]*1000.*lensw )
  rmax_act = rtagsum/lenswsum
  rmaxerr = sqrt( total( lensw^2*(rad[wlens]*1000. - rmax_act)^2 )/lenswsum^2 )
  rmax = rmaxerr
  print
  print,'mean rtag: ',rmax_act
  print,'rtag error: ',rmax
  print

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find averages from sums
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Finally used ',ntostr(lensused),'/',ntostr(nlens),' Lenses'

  totpairs = 0.0

  ;; shear struct for this set of lenses
  shstruct.meanr = rsum/npsum

  shstruct.shear = etansum/wsum/2.
  shstruct.shearerr = sqrt(etanerrsum/wsum^2)/2.
  shstruct.ortho = eradsum/wsum/2.
  shstruct.orthoerr = sqrt(eraderrsum/wsum^2)/2.

  shstruct.sigma = tansigsum/wsum/2.
  shstruct.sigmaerr = sqrt( tansigerrsum/wsum^2 )/2.
  shstruct.orthosig = radsigsum/wsum/2.
  shstruct.orthosigerr = sqrt( radsigerrsum/wsum^2 )/2.
  shstruct.npair = npsum

  ;; Now totals within each rmax_act
  shstruct.tshear = etansum/wsum/2.
  shstruct.tshearerr = sqrt( etanerrsum/wsum^2 )/2.
  shstruct.tortho = eradsum/wsum/2.
  shstruct.torthoerr = sqrt( eraderrsum/wsum^2 )/2.

  shstruct.tsigma = tansigsum/wsum/2.
  shstruct.tsigmaerr = sqrt( tansigerrsum/wsum^2 )/2.
  shstruct.torthosig = radsigsum/wsum/2.
  shstruct.torthosigerr = sqrt( radsigerrsum/wsum^2 )/2.
  shstruct.tnpair = npsum

  ;; calculate area, density of background galaxies
  R1 = rmin
  R2 = rmin + rmax_act
  shstruct.area = !pi*(R2^2 - R1^2)
  shstruct.density = shstruct.npair/shstruct.area/lensused
      
  ;;sum struct for adding to other sets of lenses
  sumstruct.rsum = rsum

  sumstruct.etansum = etansum
  sumstruct.etanerrsum = etanerrsum
  sumstruct.eradsum = eradsum
  sumstruct.eraderrsum = eraderrsum

  sumstruct.tansigsum = tansigsum
  sumstruct.tansigerrsum = tansigerrsum
  sumstruct.radsigsum = radsigsum
  sumstruct.radsigerrsum = radsigerrsum

  sumstruct.wsum = wsum
  sumstruct.wsum_ssh = wsum_ssh
  sumstruct.npair = npsum

  totpairs = totpairs + npsum

  Ssh = Sshsum/wsum_ssh

  shstruct.nlenses = lensused
  shstruct.totpairs = totpairs
;  shstruct.binsize = binsize
  shstruct.rmin = rmin
  shstruct.rmax = rmax
  shstruct.rmax_act = rmax_act
  shstruct.ssh = ssh
  shstruct.zmean = zmean
  shstruct.meanscritinv = meanscritinv
  shstruct.h = h

  sumstruct.nlenses = lensused
  sumstruct.totpairs = totpairs
;  sumstruct.binsize = binsize
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

  ;; Shear for these lenses
  shhdr = zshhdr(shstruct)    
  mwrfits, shstruct, datfile, shhdr, /create

  ;; Sum file for these lenses
  sumhdr = zsumhdr(sumstruct)
  lhdr = sumhdr
  mwrfits, sumstruct, sumfile, sumhdr, /create

  ;; File with structure for each lens used
  lensum.z = zrand
  lensum = temporary(lensum[wlens])
  mwrfits, lensum, lensumfile, lhdr, /create

  ;; File containing ra,dec,redshift for each used lens
  zs = create_struct('z', 0., $
                     'sigcritinv', 0., $
                     'ra', double(0.), $
                     'dec', double(0.) )
  zstruct = replicate(zs, lensused)
  zstruct.z = zrand[wlens]
  zstruct.sigcritinv = sigcritinv[wlens]
  zstruct.ra = lra[wlens]
  zstruct.dec = ldec[wlens]

  mwrfits, zstruct, zfile, /create


  step=oldstep
;  ep

  ptime, systime(1)-time
  return
END 
