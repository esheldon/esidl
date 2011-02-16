PRO zrandshear, run1, run2, scat, clr, rmin, rmax, binsize, zrand_in, ndup, $
                use_lambda=use_lambda, $
                step=step, addstr=addstr, $
                outdir=outdir, $
                wgood=wgood, $
                check=check, $
                datfile=datfile, sumfile=sumfile, zfile=zfile, $
                lensumfile=lensumfile, $
                maxe=maxe, $
                commonsrc=commonsrc

  IF n_params() LT 9 THEN BEGIN
      print,'-Syntax: zrandshear, run1, run2, scat, clr, rmin, rmax, binsize, zrand,, ndup'
      print,'  use_lambda=use_lambda, '
      print,'  step=step, addstr=addstr, '
      print,'  outdir=outdir, '
      print,'  title=title, wgood=wgood, '
      print,'  check=check, '
      print,'  datfile=datfile, sumfile=sumfile, zfile=zfile,'
      print,'  lensumfile=lensumfile, '
      print,'  maxe=maxe '
      print,'cosmo=1 for omegamat=1 (default)  cosmo=2 for lambda'
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
  IF NOT keyword_set(use_lambda) THEN omegamat=1.0 ELSE omegamat=0.3

  ;; Size by which to group things  300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L ELSE step = long(step)
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = ''
  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss4/data1/esheldon/TMP/'
  IF NOT keyword_set(check) THEN check=0

  r1str = '_'+ntostr(run1)
  r2str = '_'+ntostr(run2)
  clr_str = '_'+colors[clr]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; maxit is number of retries until random point
  ;; passes symmetry check
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  maxit = 20

  d2r = !dpi/180.0d0            ;Change degrees to radians.
  r2d = 180./!dpi               ;radians to degrees

  nbin = long( (rmax - rmin)/binsize ) + 1 

  print
  print,'-------------------------------------------------------------'
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
      ENDIF ELSE BEGIN 
          zrand = [zrand, zrand_in]
          zindex = [zindex, lindgen(nrand_in)]
      ENDELSE 
  ENDFOR 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; declare some arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  arrval = fltarr(nbin)

  shstruct = zshstruct(arrval)
  sumstruct = zsumstruct(arrval)
  lensumstruct = zlensumstruct(arrval)
  lensumstruct = create_struct('z', 0., 'ra', 0d, 'dec', 0d, lensumstruct)

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

  IF addstr EQ '' THEN prename = outdir+addstr+'zrand'+r1str+r2str+clr_str $
  ELSE prename = outdir+addstr+'_zrand'+r1str+r2str+clr_str

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

  CASE 1 OF
      (run1 EQ 1140) OR (run2 EQ 1140): stripe=9
      (run1 EQ 752) OR (run2 EQ 752): stripe=10
      (run1 EQ 94) OR (run2 EQ 94): stripe=82
      ELSE: message,'No info on these runs'
  ENDCASE 

  sigcritinv = sdss_sigma_crit(stripe=stripe, clr, zrand, $
                               wgood=wgood, use_lambda=use_lambda, $
                               commonsrc=commonsrc)
  IF (n_elements(sigcritinv) EQ 1) AND (sigcritinv[0] EQ -1) THEN return
  sigmacrit = 1./sigcritinv
  wlens = wlens[wgood]
  lensum.scritinv = sigcritinv

  nlens = n_elements(wlens)

  DL = angdist_lambda( zrand, h=h, omegamat=omegamat)*1000. ;Convert from Mpc to kpc
  angmax = rmax/DL*180./!pi

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

;      w=where(lra[ii] GT lra[maxii] OR lra[ii] LT lra[minii], nw)
;      IF nw NE 0 THEN BEGIN 
;          print,'What the hell!'
;          return
;      ENDIF 

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

          binary_search, scat.ra, ra1, i1, /round, /edgedefault
          binary_search, scat.ra, ra2, i2, /round, /edgedefault
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
              cendec = ldec[index]
              cenra  = lra[index]

              sig_crit = sigmacrit[index]
              sig_inv = sigcritinv[index]
              IF sig_inv EQ -1000. THEN BEGIN 
                  print,'What!'
                  return
              ENDIF 
              angmax_i = angmax[index]

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
                  binary_search, scat[wsrc].ra, ra1, i1, /round, /edgedefault
                  binary_search, scat[wsrc].ra, ra2, i2, /round, /edgedefault
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
                              lensum[index].rsum[whist] = tmprsum[whist]
                              lensum[index].rmax_act[whist] = rmax_act_tmp[whist]
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

  lensum.ra = lra[wlens]
  lensum.dec = ldec[wlens]

  mwrfits, zstruct, zfile, /create


  step=oldstep
;  ep

  ptime, systime(1)-time
  return
END 
