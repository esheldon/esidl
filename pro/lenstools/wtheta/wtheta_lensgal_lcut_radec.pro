PRO wtheta_lensgal_lcut_radec, wclr, stripe, lenscat, allcat, $
                               rmin, rmax, nbin_OR_binsize, $
                               use_lambda=use_lambda, $
                               step=step, addstr=addstr, $
                               outdir=outdir, $
                               wgood=wgood, $
                               usecat=usecat, $
                               datfile=datfile, sumfile=sumfile, zfile=zfile, $
                               lensumfile=lensumfile, $
                               maxe=maxe,fraclstar=fraclstar,$
                               numfile=numfile,nbeta=nbeta, $
                               tsgals=tsgals, edgecut=edgecut, racut=racut, $
                               maskdir=maskdir, etarangedir=etarangedir,$
                               logbin=logbin

  IF n_params() LT 7 THEN BEGIN
      print,'-Syntax: wtheta_lensgal_lcut_radec, wclr, stripe, lenscat, allcat, rmin, rmax, binsize,'
      print,'  use_lambda=use_lambda,'
      print,'  step=step, addstr=addstr, '
      print,'  outdir=outdir, '
      print,'  wgood=wgood, '
      print,'  usecat=usecat, '
      print,'  datfile=datfile, sumfile=sumfile, zfile=zfile, '
      print,'  lensumfile=lensumfile, '
      print,'  maxe=maxe'
      print,'-> lenscat may have ra,dec or lambda, eta but allcat must contain lambda,eta'
      print,'-> allcat must be sorted by lambda'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; WARNING!!
  ;; Both catalogs must be sorted by ra and transformed to equator and
  ;; not cross ra=0(360)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Common blocks
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  COMMON wtheta_block, Mlimit

  Mstar = [-18.34,-20.04,-20.83,-21.26,-21.55]
  Msun = [6.39,5.07,4.62,4.52,4.48]
  Lstar = 10.0^((Mstar-Msun)/(-2.5))
  Lstar = Lstar/1.e10           ;in units of 10^10 solar lum

  IF n_elements(fraclstar) EQ 0 THEN fraclstar = 0.1
  lcut = lstar[wclr]*fraclstar

;  Mlimit = Mstar[2]-2.5*alog10(fraclstar)

  print
  print,'Going to '+ntostr(fraclstar)+' of L*'
  print

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  colors = ['u','g','r','i','z']
  time = systime(1)

  IF NOT keyword_set(use_lambda) THEN omegamat=1.0 ELSE omegamat=0.3

  ;; Size by which to group things  300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = 'wthetalumweq'
  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss4/data1/esheldon/TMP/'

  IF wclr LT 0 OR wclr GT 4 THEN message,'wclr must be in [0,4]!'

  IF n_elements(numfile) EQ 0 THEN fend=colors[wclr]+'w_N1' $
  ELSE fend=colors[wclr]+'w_N'+ntostr(long(numfile))

  stripestr = '_stripe'+ntostr(stripe)

  d2r = !dpi/180.0d0            ;Change degrees to radians.
  r2d = 180./!dpi               ;radians to degrees

  ;;;;;;;;;;;;;;;;;;;;;;
  ;; log binning?
  ;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(logbin) THEN BEGIN 
      print,'--------------------------'
      print,'Using Logarithmic Binning'
      print,'--------------------------'
      binsize = -1.
      nbin = nbin_or_binsize
  ENDIF ELSE BEGIN 
      binsize = nbin_or_binsize
      nbin = long( (rmax - rmin)/binsize ) ;+ 1 
  ENDELSE 

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

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Beta values
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;IF n_elements(nbeta) EQ 0 THEN nbeta = 41 ;this does beta in [0,2], steps of 0.05
  IF n_elements(nbeta) EQ 0 THEN nbeta = 1
  betamin = 0.0
  betastep = 0.05
  beta = fltarr(nbeta)
  FOR i=0L, nbeta-1 DO beta[i] = betamin+i*betastep
  print
  print,'n_beta = ',nbeta
  print,'max beta = ',max(beta)
  print,'min beta = ',betamin
  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; default area in bins
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  darea = fltarr(nbin)
  R1 = darea & R2 = R1
  IF keyword_set(logbin) THEN BEGIN 
      logrmin = alog10(rmin)
      logrmax = alog10(rmax)
      binsize_log = ( logrmax - logrmin )/nbin
      r_ratio = 10.^binsize_log

      FOR i=0L, nbin-1 DO BEGIN
          R1[i] = rmin*r_ratio^(i)
          R2[i] = rmin*r_ratio^(i+1)
          darea[i] = !pi*(R2[i]^2 - R1[i]^2)/1.e6 ;Mpc^2
      ENDFOR 
  ENDIF ELSE BEGIN 

      FOR i=0L, nbin-1 DO BEGIN 
          R1[i] = rmin + i*binsize
          R2[i] = rmin + (i+1)*binsize
          darea[i] = !pi*(R2[i]^2 - R1[i]^2)/1.e6 ;Mpc^2
      ENDFOR 
  ENDELSE 

  totarea = total(darea)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; declare some arrays
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  arrval = fltarr(nbin)
  arrval2 = fltarr(nbin, nbeta)

  wthetastr=wthetaLumwStruct(arrval,beta)
  wthetasumstr=wthetaLumwSumstruct(arrval,beta)
  wthetalensum=create_struct(lenscat[0], wthetaLumwLensumStruct(arrval,beta) )
  wthetalensum=replicate(wthetalensum, n_elements(lenscat))

  wsum       = arrval
  rsum       = arrval
  lsum       = arrval
  lerrsum1   = arrval
  lerrsum2   = arrval
  lerrsum3   = arrval
  lwsum      = arrval
  npsum      = arrval2
  npair      = arrval

  rmax_act   = arrval
  rmax_act_tmp = arrval

  ;; temporary
  tmprsum = arrval
  tmplsum = arrval
  tmplerrsum1   = arrval
  tmplerrsum2   = arrval
  tmplerrsum3   = arrval
  tmpnpsum = arrval2

  hits=0L

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up output files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  prename = outdir+addstr+stripestr
  
  datfile = prename + '_'+fend+'.fit'
  sumfile = prename + '_sum_'+fend+'.fit'
  zfile   = prename + '_z_'+fend+'.fit'
  lensumfile = prename + '_lumlensum_'+fend+'.fit'

  WHILE fexist(datfile) DO BEGIN
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
  nall = n_elements(allcat)
  wall = lindgen(nall)

  FIRSTRA = allcat[0].ra
  LASTRA = allcat[nall-1].ra

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

  IF tag_exist(lcat, 'ra') AND tag_exist(lcat, 'dec') THEN BEGIN 
      lra = lcat.ra
      ldec = lcat.dec
  ENDIF ELSE BEGIN 
      IF tag_exist(lcat, 'lambda') AND tag_exist(lcat, 'eta') THEN BEGIN 
          print
          print,'-------------------------------------'
          print,'Converting lens RA,DEC to LAMBDA,ETA'
          print,'-------------------------------------'
          print
          survey2eq, lcat.lambda, lcat.eta, lra, ldec
      ENDIF ELSE BEGIN 
          print,'Neither (LAMBDA, ETA) or (RA,DEC) found in lens tags'
          return
      ENDELSE 
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; read etarange and mask files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; not going to use mask for now, just testing
  read_etarange, stripe, etarange, indir=etarangedir
  read_stripe_mask, stripe, mask, tsgals=tsgals, indir=maskdir

  IF stripe GT 45 THEN BEGIN
      print
      print,'This is a Southern Stripe.'
      print
  ENDIF

  ;; subscripts for lenses. This will change as we throw
  ;; out lenses
  wlens = lindgen(ninit)

  h=1.
  print,'Using h = ',h
  print,'Calculating 1/sigma_crit'
  print

  clr=2
  sigcritinv = sdss_sigma_crit(stripe, clr, lcat.(wz[0]), $
                               wgood=wgood, use_lambda=use_lambda)

  nlens1=n_elements(wgood)
  print,'Threw out ',ntostr(ninit - nlens1),' lenses in sdss_sigma_crit'
  sigmacrit = 1./sigcritinv

  ;;Convert from Mpc to kpc
  DL = angdist_lambda( lcat.(wz[0]), h=h, omegamat=omegamat)*1000. 
  angmax = rmax/DL*180./!pi     ;angles in degrees

  ;; don't want angle to be larger than 1/rfac of a stripe (until we use 
  ;; multiple stripes)
  rfac = 1.0
  max_allowed_angle = 2.5/rfac
  wgood2 = where(angmax[wgood] LE max_allowed_angle, nlens2)
  wgood = wgood[wgood2]
  wlens = wlens[wgood]
  print,'Threw out ',ntostr(nlens1 - nlens2),' lenses because maxangle > '+$
    ntostr(max_allowed_angle)+' degrees'

  ;; weights for each lens
  ;; 1.e19 for order 1 numbers
  lensweights = fltarr(ninit)
  lensweights[wlens] = sigcritinv[wlens]^2/DL[wlens]^2*1.e19 

;  wgood2 = where(lcat[wgood].(wz[0]) GT 0.02)
;  wgood = wgood[wgood2]
;  wlens = wlens[wgood]

  ;; cut on redshift for lcut
  ;; these are for mag limits of ??,21,21,21,19.8
  CASE wclr OF
      0: zcut = 0.09
      1: zcut = 0.15
      2: zcut = 0.218
      3: zcut = 0.27
      4: zcut = 0.18
  ENDCASE 
  print
  print,'Making zcut = ',zcut
  print

  wgood3 = where( lcat[wgood].(wz[0]) LE zcut, nlens3)
  wgood = wgood[wgood3]
  wlens = wlens[wgood3]

  print,'Threw out ',ntostr(nlens2 - nlens3),' lenses because of z completeness cut'

  ;; copy in the stuff we need into lensum struct
  copy_struct, lcat, wthetalensum
  ;lensum.scritinv = sigcritinv

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; apply mask. If we are going to use mask on random points, 
  ;; we have to apply it to lenses as well
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;apply_mask, mask, llambda[wlens], tbad, tgood, edgecut=angmax[wlens]
  ;;IF tgood[0] EQ -1 THEN message,'No objects passed mask cuts'
  ;;wlens = wlens[tgood]

  ;;print,'Threw out ',ntostr(nlens1-n_elements(tgood)),' lenses in mask'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; check edge
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  cldec = cos(ldec*d2r)
  dfromend = (LASTRA - lra)/cldec
  dfrombeg = (lra - FIRSTRA)/cldec
  
  IF keyword_set(racut) OR keyword_set(edgecut) THEN BEGIN
      print,'Making RA edge cut'
      wtmp = where( (dfromend[wlens] GT angmax[wlens]) AND $
                    (dfrombeg[wlens] GT angmax[wlens]), nlens)

      IF nlens EQ 0 THEN message,'No lenses left!'
      print,'Threw out ',ntostr(nlens3-nlens),' lenses on RA cut'

      wlens = wlens[wtmp]
  ENDIF 

;  wtmp = where( ((lra[wlens] - angmax[wlens]/cldec[wlens]) GT FIRSTRA ) AND $
;                ((lra[wlens] + angmax[wlens]/cldec[wlens]) LT LASTRA ), nlens)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Check dec
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nlens = n_elements(wlens)
  
  maxdec = max(allcat.dec)
  mindec = min(allcat.dec)
  
  IF keyword_set(edgecut) THEN BEGIN 
      print
      print,'Making edge cut'
      twlens = wlens

      ;; on equator (all are rotated to equator for convenience)
      twlens=where( ( (ldec[wlens]-angmax[wlens]) GT mindec) AND $
                    ( (ldec[wlens]+angmax[wlens]) LT maxdec), nlens2 )

      IF nlens2 EQ 0 THEN message,'No objects passed edge cut!'
      print,'Threw out ',ntostr(nlens-nlens2),' lenses on edge cut'
      wlens = wlens[twlens]
      nlens=nlens2

  ENDIF 
  
  ;; distances in kpc
  dfrommax = (maxdec - ldec)*d2r*DL
  dfrommin = (ldec - mindec)*d2r*DL
  dfromend = dfromend*d2r*DL
  dfrombeg = dfrombeg*d2r*DL

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make sure positive distance from edge
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wtmp=where( (dfrommax[wlens] GT 0.0) AND (dfrommin[wlens] GT 0.0) AND $
              (dfromend[wlens] GT 0.0) AND (dfrombeg[wlens] GT 0.0), ntmp)
  wlens = wlens[wtmp]
  print,nlens-ntmp,' negatives'
  nlens=ntmp

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; print some stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'FIRSTRA = ',FIRSTRA
  print,'LASTRA = ',LASTRA
  print,'Using '+ntostr(nlens)+'/'+ntostr(ninit)+' lenses'
  print
  IF nlens LT 100 THEN sym = 1 ELSE sym = 3
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; count neighbors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  totlwsum = 0.0
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
      angmax2 = max(angmax[ii])
      
      maxii = wlens[ max(ind) ] & minii = wlens[ min(ind) ]

      maxra1 = lra[maxii] + angmax2/cos(ldec[maxii]*d2r)
      minra1 = lra[minii] - angmax2/cos(ldec[minii]*d2r)

      wtheta_get_neighbors_ra, allcat.ra, minra1, maxra1, wneigh, nneigh

      IF nneigh NE 0 THEN BEGIN 

          IF group NE 0 THEN BEGIN 
              IF totlwsum NE 0. THEN BEGIN 
                  darea2=darea^2
                  mmll = total(lsum*darea)/totlwsum/totarea
                  mmllerr = sqrt( total(lerrsum1*darea2) - 2.*total(lerrsum2*darea2)*mmll + total(lerrsum3*darea2)*mmll^2 )/totlwsum/totarea
                  print,'Mean lum dens: '+ntostr(mmll)+' +/- '+ntostr(mmllerr)
              ENDIF 
          ENDIF 

          print,'Lens Group = ',ntostr(group+1)+'/'+ntostr(nstep)

          FOR gi=0L, step-1 DO BEGIN ;Loop over lenses in group
              IF (gi MOD 10) EQ 0 THEN  print,'.',format='(a,$)'
              index = ii[gi]
              cendec = ldec[index] ;lens center
              cenra  = lra[index]

              lensw = lensweights[index]
              totlwsum=totlwsum+lensw

              angmax_i = angmax[index]

              ;; choose sources around this lens brighter
              ;; than certain absolute magnitude

              ccendec = cos(cendec*d2r)
              maxra2 = cenra + angmax_i/ccendec
              minra2 = cenra - angmax_i/ccendec

              wtheta_get_neighbors_ra, allcat[wneigh].ra, minra2, maxra2, $
                twneigh, nneigh2

              ;; If there are any sources left, measure the shear
              IF nneigh2 GT 0 THEN BEGIN 

                  wneigh2 = wneigh[twneigh]

                  wtheta_absmag_meangr, lcat[index].(wz[0]), wclr, $
                    allcat[wneigh2].petrocounts[wclr], $
                    allcat[wneigh2].gr, absmag, lumsolar

                  wneigh3 = where(lumsolar GE lcut, nneigh2)
                  IF nneigh2 EQ 0 THEN GOTO,jump
                  wneigh2 = wneigh2[wneigh3]
                  lumsolar=lumsolar[wneigh3]

                  ;; luminosity weight will be (lum/lstar)^beta
                  ;; must both be in units of 10^10 solar
                  lumw=lumsolar/lstar[wclr]
                  
                  mygcirc, cenra, cendec, $
                    allcat[wneigh2].ra, allcat[wneigh2].dec, $
                    dis, /radians_out

                  R = dis*DL[index]
                  IF nneigh2 EQ 1 THEN R=[R]

                  IF NOT keyword_set(logbin) THEN BEGIN 
                      hist=histogram(R, binsize=nbin_or_binsize, min=rmin, $
                                     max=rmax,rever=rev_ind)
                  ENDIF ELSE BEGIN 
                      logbin, R, rmin, rmax, nbin_or_binsize, hist, rev_ind
                  ENDELSE 

                  whist = where(hist NE 0, nhist)
                  ng=0
                  ;; Check if there are any in this annulus rmin-rmax
                  thit=0L
                  IF nhist NE 0 THEN BEGIN 
                          
                      ;;nbin is predefined. Last bin will hold few
                      npair[*] = 0L
                      rmax_act_tmp[*] = 0.
                      FOR binnum=0L, nbin-1 DO BEGIN 
                          
                          IF rev_ind[binnum] NE rev_ind[binnum+1] THEN BEGIN 

                              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                              ;; area Mpc^2  Figure out if we are
                              ;; missing any of the bin due to edge
                              ;; effects
                              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                              
                              
;                              tbarea = binarea(dfrommax[index], $
;                                              R1[binnum], $
;                                              R2[binnum],$
;                                              d_edge2=dfrommin[index],hit=tthit)/1.e6
                              circ_rect_intersect, $
                                dfrommax[index], dfrommin[index],$
                                dfromend[index], dfrombeg[index],$
                                R2[binnum], area2, hit=tthit
                              circ_rect_intersect, $
                                dfrommax[index], dfrommin[index],$
                                dfromend[index], dfrombeg[index],$
                                R1[binnum], area1
                              barea = (area2 - area1)/1.e6

                              thit=thit+tthit
                              
                              w=rev_ind( rev_ind(binnum):rev_ind(binnum+1)-1 )
                              
                              rmax_act_tmp[binnum] = max(R[w])
                              rmax_act[binnum] = max( [rmax_act[binnum], rmax_act_tmp[binnum]] )
                              npair[binnum] = n_elements(w)
                              
                              ;; loop over the different beta values
                              tmprsum[binnum]=total(R[w])/barea*lensw
                              FOR ib=0L, nbeta-1 DO BEGIN 
                                  ;; luminosity weighted counts
                                  tmpnpsum[binnum,ib]=total(lumw[w]^beta[ib])/barea*lensw
                              ENDFOR 
                              xi = total(lumsolar[w])/barea
                              tmplsum[binnum] = xi*lensw
                              tmplerrsum1[binnum] = xi^2*lensw^2
                              tmplerrsum2[binnum] = xi*lensw^2
                              tmplerrsum3[binnum] = lensw^2
                          ENDIF 
                      ENDFOR 

                      rsum[whist] = rsum[whist] + tmprsum[whist]
                      npsum[whist,*] = npsum[whist,*] + tmpnpsum[whist,*]
                      lsum[whist] = lsum[whist] + tmplsum[whist]

                      lerrsum1[whist] = lerrsum1[whist] + tmplerrsum1[whist]
                      lerrsum2[whist] = lerrsum2[whist] + tmplerrsum2[whist]
                      lerrsum3[whist] = lerrsum3[whist] + tmplerrsum3[whist]

;                      lwsum[whist] = lwsum[whist] + tmplwsum[whist]
;                      wsum[whist] = wsum[whist] + tmpwsum[whist]

                      ;; copy in individual lens stuff
                      ;; Note totpairs is unweighted
                      wthetalensum[index].totpairs = total(npair[whist])
                      wthetalensum[index].npsum[whist,*] = tmpnpsum[whist,*]
                      wthetalensum[index].rmax_act = R2
                      wthetalensum[index].rmin_act = R1
                      wthetalensum[index].rsum[whist] = tmprsum[whist]
                      wthetalensum[index].lsum[whist] = tmplsum[whist]

                      wthetalensum[index].lerrsum1[whist] = tmplerrsum1[whist]
                      wthetalensum[index].lerrsum2[whist] = tmplerrsum2[whist]
                      wthetalensum[index].lerrsum3[whist] = tmplerrsum3[whist]

;                      wthetalensum[index].lwsum[whist] = tmplwsum[whist]
;                      wthetalensum[index].wsum[whist] = tmpwsum[whist]

                      ;; reinitialize the arrays
                      setarrzero,tmprsum,tmpnpsum,tmplsum,$
                        tmplerrsum1,tmplerrsum2,tmplerrsum3,$
                        npair

                      hits=hits+(thit GT 0)

                  ENDIF ELSE BEGIN 
                      print,'/',format='(a,$)'
                  ENDELSE 

                  R = 0
                  etadiff=0      ;Free memory
              ENDIF ELSE BEGIN
jump:
                  ;; No neighbors of this point
                  print,'|',format='(a,$)'
              ENDELSE 

          ENDFOR 
          print
      ENDIF ELSE BEGIN 
          print,'Two-'
      ENDELSE 
      lastind=max(ind)
      metadiff = 0               ;Free memory
      mdist = 0
  ENDFOR 

  lensused = n_elements(wlens)
  wgood = wlens
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find averages from sums
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Finally used ',ntostr(lensused),'/',ntostr(nlens),' Lenses'
  print,'Number of edge hits: ',hits

  ;; calculate sum of weights
  FOR i=0L, nbin-1 DO BEGIN 

      ;; same weight in each radial bin
      lwsum[i] = total(lensweights[wlens])
      wsum[i] = lwsum[i]

      wthetalensum[wlens].lwsum[i] = lensweights[wlens]
      wthetalensum[wlens].wsum[i] = lensweights[wlens]

  ENDFOR 

  FOR i=0L, nbin-1 DO BEGIN 

      ;; calculate area, density of background galaxies
      R1 = rmin + i*binsize
      R2 = rmin + (i+1)*binsize
      wthetastr.area[i] = !pi*(R2^2 - R1^2)

      ;; shear struct for this set of lenses
      wthetastr.meanr[i] = rsum[i]/npsum[i,0]
      ;; does not include central galaxy
      ml = lsum[i]/lwsum[i]
      wthetastr.meanlum[i] = ml

      wthetastr.meanlumerr[i]=sqrt($
                                    (lerrsum1[i] - $
                                     2.*lerrsum2[i]*ml + $
                                     lerrsum3[i]*ml^2 $
                                    ) > 0. $
                                  )/lwsum[i]

      FOR ib=0L, nbeta-1 DO BEGIN 
          ;; now an average per lens
          wthetastr.npair[i,ib] = npsum[i,ib]/wsum[i]

          ;; Now totals within each rmax_act[i]
          wthetastr.tnpair[i,ib] = total(npsum[0:i,ib])/total(wsum[0:i])

          ;; npair is average per lens, no need to divide by lensused
          wthetastr.density[i,ib]=wthetastr.npair[i,ib]
      ENDFOR 
      ; sum struct for adding to other sets of lenses
      wthetasumstr.rsum[i] = rsum[i]
      wthetasumstr.npsum[i,*] = npsum[i,*]
      wthetasumstr.lsum[i] = lsum[i]

      wthetasumstr.lerrsum1[i] = lerrsum1[i]
      wthetasumstr.lerrsum2[i] = lerrsum2[i]
      wthetasumstr.lerrsum3[i] = lerrsum3[i]

      wthetasumstr.lwsum[i] = lwsum[i]
      wthetasumstr.wsum[i] = wsum[i]

  ENDFOR 
  
  totpairs = total(wthetalensum[wlens].totpairs)

  wthetastr.nlenses = lensused
  wthetastr.totpairs = totpairs
  wthetastr.binsize = binsize
  wthetastr.rmin = rmin
  wthetastr.rmax = rmax
  wthetastr.rmax_act = rmax_act
  wthetastr.h = h

  wthetasumstr.nlenses = lensused
  wthetasumstr.totpairs = totpairs
  wthetasumstr.binsize = binsize
  wthetasumstr.rmin = rmin
  wthetasumstr.rmax = rmax
  wthetasumstr.rmax_act = rmax_act
  wthetasumstr.h = h

  wthetalensum.binsize = binsize
  wthetalensum.rmin = rmin
  wthetalensum.rmax = rmax
  wthetalensum.h = h

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Ouput the data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; headers
  hdr = whdr(wthetastr)      
  sumhdr=hdr
  lhdr=hdr

  ;; Shear file for these lenses
  mwrfits, wthetastr, datfile, hdr, /create

  ;; Sum file for these lenses
  mwrfits, wthetasumstr, sumfile, sumhdr, /create

  ;; File with structure for each lens used
  IF n_elements(wlens) LT n_elements(wthetalensum) THEN $
    wthetalensum = temporary(wthetalensum[wlens])
  wthetalensum.zindex = lindgen(lensused)
  mwrfits2, wthetalensum, lensumfile, lhdr, /create, /destroy

  ;; File containing ra,dec,redshift for each used lens
  zs = create_struct('z', 0., $
                     'lensw',0.0,$
                     'ra', double(0.), $
                     'dec', double(0.) )
  zstruct = replicate(zs, lensused)
  zstruct.z = lcat[wlens].(wz[0])
  zstruct.ra = lra[wlens]
  zstruct.dec = ldec[wlens]

  mwrfits, zstruct, zfile, /create

  step=oldstep
  usecat = lcat[wlens]

  ptime, systime(1)-time
  return
END 
