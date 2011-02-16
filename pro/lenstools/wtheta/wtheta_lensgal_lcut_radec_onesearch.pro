PRO wtheta_lensgal_lcut_radec_onesearch, wclr, stripe, lenscat, allcat, rmin, rmax, binsize, $
                         use_lambda=use_lambda, $
                         step=step, addstr=addstr, $
                         outdir=outdir, $
                         wgood=wgood, $
                         usecat=usecat, $
                         datfile=datfile, sumfile=sumfile, zfile=zfile, $
                         lensumfile=lensumfile, $
                         maxe=maxe,fraclstar=fraclstar,$
                         numfile=numfile,nbeta=nbeta, $
                         tsgals=tsgals, edgecut=edgecut, $
                         maskdir=maskdir, etarangedir=etarangedir

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

  ;; Both catalogs must be sorted by ra
  
  ;; Only one get neighbors is used, checking for errors in searching

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
  IF n_elements(addstr) EQ 0 THEN addstr = 'wthetalumw'
  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss4/data1/esheldon/TMP/'

  IF wclr LT 0 OR wclr GT 4 THEN message,'wclr must be in [0,4]!'

  IF n_elements(numfile) EQ 0 THEN fend=colors[wclr]+'w_N1' $
  ELSE fend=colors[wclr]+'w_N'+ntostr(long(numfile))

  stripestr = '_stripe'+ntostr(stripe)

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
  FOR i=0L, nbin-1 DO BEGIN 
      R1[i] = rmin + i*binsize
      R2[i] = rmin + (i+1)*binsize
      darea[i] = !pi*(R2[i]^2 - R1[i]^2)/1.e6 ;Mpc^2
  ENDFOR 

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
  tmpwsum = arrval
  tmprsum = arrval
  tmplsum = arrval
  tmplerrsum1   = arrval
  tmplerrsum2   = arrval
  tmplerrsum3   = arrval
  tmplwsum = arrval
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
  sigmacrit = 1./sigcritinv

  ;;Convert from Mpc to kpc
  DL = angdist_lambda( lcat.(wz[0]), h=h, omegamat=omegamat)*1000. 
  angmax = rmax/DL*180./!pi     ;angles in degrees

  ;; don't want angle to be larger than 1/rfac of a stripe (until we use 
  ;; multiple stripes)
  rfac = 1.5
  max_allowed_angle = 2.5/rfac
  wgood2 = where(angmax[wgood] LE max_allowed_angle)
  wgood = wgood[wgood2]
  wlens = wlens[wgood]

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

  wgood3 = where( lcat[wgood].(wz[0]) LE zcut )
  wgood = wgood[wgood3]
  wlens = wlens[wgood3]

  nlens1 = n_elements(wlens)
  print,'Threw out ',ntostr(ninit - nlens1),' lenses because too deep or close'

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

  print,'Making lambda cut'
  cldec = cos(ldec*d2r)
  wtmp = where( ((lra[wlens] - angmax[wlens])/cldec[wlens] GE FIRSTRA ) AND $
                ((lra[wlens] + angmax[wlens])/cldec[wlens] LE LASTRA ), nlens)

  IF nlens EQ 0 THEN message,'No lenses left!'
  print,'Threw out ',ntostr(n_elements(wlens)-nlens),' lenses on RA cut'
  wlens = wlens[wtmp]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Check eta
  ;; Make etarange structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nlens = n_elements(wlens)
  
  maxdec = max(allcat.dec)
  mindec = min(allcat.dec)

  dfrommax = (maxdec - ldec)*d2r*DL
  dfrommin = (ldec - mindec)*d2r*DL
  
  IF keyword_set(edgecut) THEN BEGIN 
      print
      print,'Making edge cut'
      twlens = wlens

      ;; on equator (all are rotated to equator for convenience)
      twlens=where( ( (ldec[wlens]-mindec) GT 0.) AND $
                    ( (maxdec - ldec[wlens]) GT 0.), nlens2 )

      IF nlens EQ 0 THEN message,'No objects passed edge cut!'
      print,'Threw out ',ntostr(n_elements(wlens)-nlens),' lenses on edge cut'
      wlens = wlens[twlens]
      nlens=nlens2

  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Calculate distance from edge for each
  ;; lens. Currently only deals with being
  ;; close to ONE edge
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wtmp=where( (dfrommax[wlens] LT 0.0) OR (dfrommin[wlens] LT 0.0), ntmp)
  IF ntmp NE 0 THEN BEGIN 
      print,ntmp,'  negatives'
      ;;forprint,d_edge[wlens[wtmp]]
      ;;stop
      remove, wtmp, wlens
      nlens = nlens - ntmp
  ENDIF 

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


  indices = lindgen(nlens)
  print,'#1/'+ntostr(nlens)
  FOR ii=0L, nlens-1 DO BEGIN

      IF ( ( (ii+1) MOD 300 ) EQ 0) AND (ii NE 0) THEN BEGIN 
          print,'#'+ntostr(ii+1)+'/'+ntostr(nlens)
          IF total(lwsum) NE 0 THEN BEGIN 
              mmll = total(lsum)/total(lwsum)
              print,'  Mean Lum dens: '+ntostr(mmll),' +/- ',$
                ntostr(sqrt( total(lerrsum1) - 2.*total(lerrsum2)*mmll + total(lerrsum3)*mmll^2 )/total(lwsum))
          ENDIF 
      ENDIF 

      index = wlens[ii]
      cendec = ldec[index]      ;lens center
      cenra  = lra[index]

      sig_inv = sigcritinv[index]
      lensw = sig_inv^2/DL[index]^2*1.e19 ;1.e19 for order 1 numbers

      angmax_i = angmax[index]
      
      ;; choose sources around this lens brighter
      ;; than certain absolute magnitude

      ccendec = cos(cendec*d2r)
      maxra2 = (cenra+angmax_i)/ccendec
      minra2 = (cenra-angmax_i)/ccendec
      
      wtheta_get_neighbors_ra, allcat.ra, minra2, maxra2, $
        wneigh2, nneigh2, issouth=issouth

      ;; If there are any sources left, measure the shear
      IF nneigh2 GT 0 THEN BEGIN 
          
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
                  
          gcirc,0,cenra*d2r,cendec*d2r,$
            allcat[wneigh2].ra*d2r,allcat[wneigh2].dec*d2r,dis
          R = dis*DL[index]
          IF nneigh2 EQ 1 THEN R=[R]

          hist=histogram(R, binsize=binsize, min=rmin, $
                         max=rmax,rever=rev_ind)
          numbin=n_elements(hist)
          IF numbin NE nbin THEN BEGIN
              print,'What!!!'
              print,numbin,nbin
              message,' '
          ENDIF 
          whist = where(hist NE 0, nhist)
          ng=0
          ;; Check if there are any in this annulus rmin-rmax
          thit=0L
          IF nhist NE 0 THEN BEGIN 
              
              ;;nbin is predefined. Last bin will hold few
              npair[*] = 0L
              rmax_act_tmp[*] = 0.
              FOR i=0L, nhist-1 DO BEGIN 
                  binnum = whist[i]
                  
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; area Mpc^2  Figure out if we are
                  ;; missing any of the bin due to edge
                  ;; effects
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  
                  barea = binarea(dfrommax[index], $
                                  R1[binnum], $
                                  R2[binnum],$
                                  d_edge2=dfrommin[index],hit=tthit)/1.e6
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
                  tmplwsum[binnum] = lensw
                  tmpwsum[binnum] = lensw
                  
              ENDFOR 
              
              rsum[whist] = rsum[whist] + tmprsum[whist]
              npsum[whist,*] = npsum[whist,*] + tmpnpsum[whist,*]
              lsum[whist] = lsum[whist] + tmplsum[whist]
              
              lerrsum1[whist] = lerrsum1[whist] + tmplerrsum1[whist]
              lerrsum2[whist] = lerrsum2[whist] + tmplerrsum2[whist]
              lerrsum3[whist] = lerrsum3[whist] + tmplerrsum3[whist]
              
              lwsum[whist] = lwsum[whist] + tmplwsum[whist]
              wsum[whist] = wsum[whist] + tmpwsum[whist]
              
              ;; copy in individual lens stuff
              ;; Note totpairs is unweighted
              wthetalensum[index].totpairs = total(npair[whist])
              wthetalensum[index].npsum[whist,*] = tmpnpsum[whist,*]
              wthetalensum[index].rmax_act[whist] = rmax_act_tmp[whist]
              wthetalensum[index].rsum[whist] = tmprsum[whist]
              wthetalensum[index].lsum[whist] = tmplsum[whist]
              
              wthetalensum[index].lerrsum1[whist] = tmplerrsum1[whist]
              wthetalensum[index].lerrsum2[whist] = tmplerrsum2[whist]
              wthetalensum[index].lerrsum3[whist] = tmplerrsum3[whist]
              
              wthetalensum[index].lwsum[whist] = tmplwsum[whist]
              wthetalensum[index].wsum[whist] = tmpwsum[whist]
              
              ;; reinitialize the arrays
              setarrzero,tmprsum,tmpnpsum,tmplsum,$
                tmplerrsum1,tmplerrsum2,tmplerrsum3,$
                tmplwsum,tmpwsum,npair
              
              hits=hits+(thit GT 0)
              
          ENDIF ELSE BEGIN 
              print,'/',format='(a,$)'
          ENDELSE 

          R = 0

      ENDIF ELSE BEGIN 
jump:
          print,'Two-'
      ENDELSE 

  ENDFOR 

  lensused = n_elements(wlens)
  
  wgood = wlens
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find averages from sums
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Finally used ',ntostr(lensused),'/',ntostr(nlens),' Lenses'
  print,'Number of edge hits: ',hits

  IF npsum[nbin-1,0] LT npsum[nbin-2,0] THEN BEGIN 
      ;; This means last bin was incomplete
      nbin = nbin-1
  ENDIF 

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
  
  totpairs = total(wthetalensum.totpairs)

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
