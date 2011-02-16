;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    ZOBJSHEAR_HTM_LOOPLENS
;       
; PURPOSE:
;    Tool for looping over lenses and measureing the shear
;
; CALLING SEQUENCE:
;    zobjshear_htm_looplens
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
;    13-NOV-2002
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO zobjshear_htm_looplens_printf, outlun, $
                                   lensum, index, zlens, $
                                   format

  IF n_elements(format) EQ 0 THEN BEGIN 
      nbin = n_elements(lensum[index].rsum)
      nbinstr = ntostr(nbin)
      format = $
        '('+$
        '1(I0,:,1X)'+$          ;index
        '1(I0,:,1X)'+$          ;zindex (only set for random)
        '1(F0,:,1x)'+$          ;redshift
        '2(F15.10,:,1X)'+$      ;clambda,ceta
        '1(I0,:,1X)'+$          ;pixelmaskflags
        '1(e0,:,1X)'+$          ;scritinv
        '1(F0,:,1X)'+$          ;totpairs
        '1(e0,:,1X)'+$          ;sshsum
        '1(e0,:,1X)'+$          ;wsum_ssh
        '1(e0,:,1X)'+$          ;weight
        '1(F0,:,1X)'+$          ;ie
        nbinstr+'(F0,:,1X)'+$   ;npair
        nbinstr+'(F0,:,1X)'+$   ;rmax_act
        nbinstr+'(F0,:,1X)'+$   ;rmin_act
        nbinstr+'(F0,:,1X)'+$   ;rsum
        nbinstr+'(F0,:,1X)'+$   ;sigma
        nbinstr+'(F0,:,1X)'+$   ;sigmaerr
        nbinstr+'(F0,:,1X)'+$   ;orthosig
        nbinstr+'(F0,:,1X)'+$   ;orthosigerr
        nbinstr+'(e0,:,1X)'+$   ;sigerrsum
        nbinstr+'(e0,:,1X)'+$   ;orthosigerrsum
        nbinstr+'(e0,:,1X)'+$   ;wsum
        nbinstr+'(e0,:,1X)'+$   ;owsum
        ')'
  ENDIF 

  printf, outlun, $
          index, $
          lensum[index].zindex, $
          zlens, $
          lensum[index].clambda, lensum[index].ceta, $
          lensum[index].pixelmaskflags, $
          lensum[index].scritinv, $
          lensum[index].totpairs, $
          lensum[index].sshsum, $
          lensum[index].wsum_ssh, $
          lensum[index].weight, $
          lensum[index].ie, $
          lensum[index].npair, $
          lensum[index].rmax_act, $
          lensum[index].rmin_act, $
          lensum[index].rsum, $ 
          lensum[index].sigma, $
          lensum[index].sigmaerr, $
          lensum[index].orthosig, $
          lensum[index].orthosigerr, $
          lensum[index].sigerrsum, $
          lensum[index].orthosigerrsum, $
          lensum[index].wsum, $
          lensum[index].owsum, format=format

END 

PRO zobjshear_htm_looplens_testquad, maskflags, theta, keep, nkeep

  ;; we assume that the lenses passed at least two of the
  ;; quadrants.  We will keep sources from the first two adjacent
  ;; quadrants that are good.  If there are quadrants left, we will
  ;; check them also

  test12 = !FLAGS_QUAD1_MASKED+!FLAGS_QUAD2_MASKED
  test23 = !FLAGS_QUAD2_MASKED+!FLAGS_QUAD3_MASKED
  test34 = !FLAGS_QUAD3_MASKED+!FLAGS_QUAD4_MASKED
  test41 = !FLAGS_QUAD4_MASKED+!FLAGS_QUAD1_MASKED

  bad12 = maskflags AND test12
  bad23 = maskflags AND test23
  bad34 = maskflags AND test34
  bad41 = maskflags AND test41

  ;; 1+2 or 3+4 are not masked
  IF (bad12 EQ 0) OR (bad34 EQ 0) THEN BEGIN 

      ;; keeping both quadrants
      IF (bad12 EQ 0) AND (bad34 EQ 0) THEN BEGIN 
          nkeep = n_elements(theta)
          keep = lindgen(nkeep)
          return
      ENDIF 

      ;; only keeping one set of quadrants
      IF bad12 EQ 0 THEN BEGIN 
          ;; quads 1+2
          keep = where(theta GE 0.0 AND theta LE !pi, nkeep)
          return
      ENDIF ELSE BEGIN
          ;; quads 3+4
          keep = where(theta GE !pi AND theta LE 2*!pi, nkeep)
          return
      ENDELSE 
  ENDIF 

  ;; 2+3 or 4+1
  IF (bad23 EQ 0) OR (bad41 EQ 0) THEN BEGIN 

      ;; keeping both quadrants
      IF (bad23 EQ 0) AND (bad41 EQ 0) THEN BEGIN 
          nkeep = n_elements(theta)
          keep = lindgen(nkeep)
          return
      ENDIF 

      ;; only keeping one set of quadrants
      IF bad23 EQ 0 THEN BEGIN 
          ;; quads 2+3
          keep = where(theta GE !pi/2. AND theta LE 3.*!pi/2., nkeep)
          return
      ENDIF ELSE BEGIN
          ;; quads 4+1
          ;; watch the boundary at 2.*pi
          keep = where( ( theta GE 3.*!pi/2.0 AND theta LE 2.*!pi) OR $
                        ( theta GE 0.0        AND theta LE !pi/2.), nkeep )
          return
      ENDELSE 

  ENDIF 

  print
  message,'What!?'

END 

PRO zobjshear_htm_looplens, depth, lensum, scat, $
                            step, maxe,$
                            rmin, rmax, nbin_or_binsize, $
                            groupstruct, indices, $
                            logbin=logbin, issouth=issouth, $
                            scinvStruct=scinvStruct, $
                            recorr=recorr, hirata=hirata, outlun=outlun, $
                            deltaFuncPhotoZ=deltaFuncPhotoZ, $
                            noTestQuad=noTestQuad, $
                            checkRidgeLine=checkRidgeLine

  IF n_params() LT 8 THEN BEGIN 
      print,'-Syntax: zobjshear_htm_looplens, depth, lensum, scat, '+$
      print,'                    step, maxe,$'
      print,'                    rmin, rmax, nbin_or_binsize, $'
      print,'                    groupstruct, indices, $'
      print,'                    logbin=logbin, issouth=issouth,$'
      print,'                    scinvStruct=scinvStruct, $'
      print,'                    /recorr, /hirata, /deltaFuncPhotoZ, $'
      print,'                    /noTestQuad, /checkRidgeLine'
      return
  ENDIF 

  nz_photoz_test = 100

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;

  ;; shapenoise
  shapenoise = 0.32
  shapenoise2 = shapenoise^2

  ;; number of points in the mean sigmacrit 
  ;; calculation for photoz's
  npts = 70

  ;; default number in leaflist
  numdefault = 100000L

  IF keyword_set(checkRidgeLine) THEN BEGIN 
      print
      maxbcg_read_colors, bcg_colors, bcg_color_spread
  ENDIF 

  IF n_elements(scinvStruct) EQ 0 THEN mean_scinv = 0 ELSE mean_scinv = 1

  IF keyword_set(deltaFuncPhotoZ) THEN BEGIN 
      print
      print,'-----------------------------------------------------'
      print,'Using photozs as delta function distributions'
      print,'-----------------------------------------------------'
      print
  ENDIF 

  IF keyword_set(deltaFuncPhotoZ) OR mean_scinv THEN BEGIN 
      use_photoz = 1
  ENDIF ELSE BEGIN 
      use_photoz = 0
  ENDELSE 

  wz=getztag(lensum[0])

  IF keyword_set(logbin) THEN BEGIN
      nbin = nbin_or_binsize
  ENDIF ELSE BEGIN
      binsize = nbin_or_binsize
      nbin = long( (rmax - rmin)/binsize ) ;+ 1 
  ENDELSE 

  IF keyword_set(recorr) THEN BEGIN 
      IF NOT tag_exist(scat,'e1_recorr',index=e1tag) THEN $
        message,'No e1_recorr tag'
      IF NOT tag_exist(scat,'e2_recorr',index=e2tag) THEN $
        message,'No e2_recorr tag'
  ENDIF ELSE BEGIN 
      IF NOT tag_exist(scat,'e1',index=e1tag) THEN message,'No e1 tag'
      IF NOT tag_exist(scat,'e2',index=e2tag) THEN message,'No e2 tag'
  ENDELSE 

  IF NOT tag_exist(scat,'clambda') THEN $
    message,'scat must contain CLAMBDA and CETA tags'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; should we print the results to a file?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; print out
  print,'In ZOBJSHEAR_HTM_LOOPLENS'
  IF keyword_set(logbin) THEN print,'Doing logbin'
  print

  nlens = n_elements(lensum)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; here is where zobjshear_lambda_looplens and 
; zobjshear_htm_looplens differ
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; find leafids
  ;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Looking up leafids'
  csurvey2eq, scat.clambda, scat.ceta, sra, sdec
  leafids = htm_index(sra, sdec, depth)

  ;;htmLookupRaDec, sra, sdec, depth, leafids
  delvarx, sra, sdec

  ;; histogram leafids
  print
  print,'Histograming leafids'

  minid = min(leafids, max=maxid)
  leaf_hist = histogram(leafids, min=minid, max=maxid, $
                   reverse_indices=leaf_rev_ind)
  leaf_nhist = n_elements(leaf_hist)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Keep track of the shear in groups
  ;; Expect clambda to be sorted, so this is
  ;; basically groups in clambda
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nstep = nlens/step
  left = nlens MOD step
  IF left NE 0 THEN nstep = nstep+1

  groupstruct = create_struct('groupid',lindgen(nstep),$
                              'groupsize',step,$
                              'nlens',lonarr(nstep),$
                              'sigma',replicate(-1.e10,nstep),$
                              'sigmaerr',replicate(-1.e10,nstep),$
                              'mean_lambda',replicate(-1.e10,nstep),$
                              'mean_eta',replicate(-1.e10,nstep))

  ;; should subscript indices with l_index, just
  indices = lindgen(nlens)
  sigwsum=0.
  sigsum=0.
  lamsum = 0d
  etasum = 0d
  numsum = 0.0
  sshsum = 0.0

  group=0
  FOR index=0L, nlens-1 DO BEGIN 
      
      zlens = lensum[index].(wz)

      ceneta  = lensum[index].ceta
      cenlam  = lensum[index].clambda

      csurvey2eq, lensum[index].clambda, lensum[index].ceta, cenra, cendec
      
      searchrad = lensum[index].angMax*!d2r
          


      leaflist = htm_intersect(cenra, cendec, depth, searchrad)
      nlist=n_elements(leaflist)

      ;;htmIntersectRaDec, cenra, cendec, searchrad, depth, leaflist, $
      ;;  numdefault=numdefault
          
      ;; don't know how many objects match
      ;; need temporary pointers to deal with unknown
      nwsrc=0L
      tmp_ptrlist = ptrarr(nlist)
      tmp_numlist = lonarr(nlist)

      ;; see which leafids match
      FOR leaf=0L, nlist-1 DO BEGIN 
          
          binnum = leaflist[leaf] - minid
          IF (binnum LT leaf_nhist) AND (binnum GE 0) THEN BEGIN 
              
              ;; any matches?
              IF leaf_rev_ind[binnum] NE leaf_rev_ind[binnum+1] THEN BEGIN 
                  
                  keep = $
                    leaf_rev_ind[leaf_rev_ind[binnum]:leaf_rev_ind[binnum+1]-1]
                  nkeep = n_elements(keep)
                  
                  nwsrc = nwsrc+nkeep
                  tmp_numlist[leaf] = nkeep
                  tmp_ptrlist[leaf] = ptr_new(keep,/no_copy)
                  
              ENDIF 
          ENDIF 
          
      ENDFOR 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; combine ptrlist and do further checks
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF nwsrc GT 1 THEN BEGIN 

          beg=0L
          wsrc = lonarr(nwsrc)
          FOR leaf=0L, nlist-1 DO BEGIN 
              IF tmp_numlist[leaf] NE 0 THEN BEGIN 
                  wsrc[beg:beg+tmp_numlist[leaf]-1] = *tmp_ptrlist[leaf]
              ENDIF 
              beg = beg+tmp_numlist[leaf]
          ENDFOR 
          ptrarr_free, tmp_ptrlist
          
          tmp_numlist=0

          mygcirc_survey, cenlam, ceneta, $
            scat[wsrc].clambda, scat[wsrc].ceta, $
            R, theta, /radians_out
                            
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; deal with edges: only choose sources in the 
          ;; appropriate quadrants
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; quadrants calculated in a coordinate system where
          ;; x is lambda, y is eta
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


          IF NOT keyword_set(noTestQuad) THEN BEGIN 

              theta2 = !pi/2. - theta
              maskflags = lensum[index].pixelmaskflags
              zobjshear_htm_looplens_testquad, maskflags, theta2, keep, nwsrc

          
              IF nwsrc NE 0 THEN BEGIN 
                  wsrc = wsrc[keep]
                  R = R[keep]
                  theta = theta[keep]
              ENDIF ELSE BEGIN 
                  wsrc = -1L
              ENDELSE 
          ENDIF

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; For clusters, exclude objects in the ridgeline
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF keyword_set(checkRidgeLine) AND nwsrc NE 0 THEN BEGIN 

              ;; This may not be correct for future catalogs,
              ;; make sure to check

              bcgIndex = round(zlens*1000)

              bcg_gr = bcg_colors[bcgIndex].kgr
              bcg_ri = bcg_colors[bcgIndex].kri

              max_gr = bcg_gr + 2*bcg_color_spread.gr_sig
              min_gr = bcg_gr - 2*bcg_color_spread.gr_sig
              max_ri = bcg_ri + 2*bcg_color_spread.ri_sig
              min_ri = bcg_ri - 2*bcg_color_spread.ri_sig
                  
              bad = where(scat[wsrc].grModel LE max_gr AND $
                          scat[wsrc].grModel GE min_gr AND $
                          scat[wsrc].riModel LE max_ri AND $
                          scat[wsrc].riModel GE min_ri, nbad, $
                          comp=keep, ncomp=nwsrc)

              IF nwsrc NE 0 THEN BEGIN 


;                   !p.multi=[0,0,2]

;                   plot, scat[wsrc].grModel, scat[wsrc].rimodel, $
;                     psym=3,xrange=[0,3],yrange=[-1,2], $
;                     xtitle='g-r',ytitle='r-i'
;                   plot_box,min_gr,max_gr,min_ri,max_ri,color=!red
                  
;                   oplot,scat[wsrc[bad]].grModel,scat[wsrc[bad]].rimodel,$
;                     psym=3,color=!red
;                   oplot,bcg_colors.kgr,bcg_colors.kri,color=!green
;                   oplot,[bcg_colors[bcgIndex].kgr],[bcg_colors[bcgIndex].kri],$
;                     psym=8,color=!blue
                      
                      
                      
;                   photoz_dist, $
;                     scat[wsrc].photoz_z, scat[wsrc].photoz_zerr, $
;                     nz_photoz_test, zvals, pofz, minz=0,maxz=0.5
;                   photoz_dist, $
;                     scat[wsrc[bad]].photoz_z, scat[wsrc[bad]].photoz_zerr, $
;                     nz_photoz_test, bzvals, bpofz
;                   plot, bzvals, bpofz, xrange=[0,0.5]
;                   oplot, zvals, pofz, color=!red
;                   oplot,[zlens,zlens],[0,1000],color=!blue
;                   !p.multi=0
                  
;                   key=prompt_kbrd('hit a key')


                  wsrc = wsrc[keep]
                  R = R[keep]
                  theta = theta[keep]
              ENDIF ELSE BEGIN 
                  wsrc = -1L
              ENDELSE 

          ENDIF 
          
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; passed checks?
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF nwsrc GT 1 THEN BEGIN 

          sig_inv = lensum[index].scritinv
          sig_crit = 1./sig_inv
          IF sig_inv EQ -1000. THEN message,'Bad sigcrit'
          
          R = R*lensum[index].DL

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Can either integrate over the distributions (tabulates)
          ;; or treat them as delta functions
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF keyword_set(deltaFuncPhotoZ) THEN BEGIN 
              sigcrit_inv = $
                sigmacritinv(zlens, scat[wsrc].photoz_z, zlens)
          ENDIF ELSE IF mean_scinv THEN BEGIN 
              zlens_send = replicate(zlens,nwsrc)
              photoz_z = scat[wsrc].photoz_z
              photoz_zerr = scat[wsrc].photoz_zerr
              
              sigcrit_inv = interp_meanscinv(zlens_send, $
                                             photoz_z, $
                                             photoz_zerr, $
                                             scinvStruct)
          ENDIF 
              
          IF NOT keyword_set(logbin) THEN BEGIN 
              hist=histogram(R, binsize=nbin_or_binsize, min=rmin, $
                             max=rmax,rever=rev_ind)
          ENDIF ELSE BEGIN 
              logbin, R, rmin, rmax, nbin_or_binsize, hist, rev_ind
          ENDELSE 

          whist = where(hist[0:nbin-1] NE 0, nhist)
          ng=0L
          xysum = 0.
          xpysum = 0.
          xmysum = 0.
          ;; Check if there are any in this annulus rmin-rmax
          IF nhist NE 0 THEN BEGIN 
              
              FOR binnum=0L, nbin-1 DO BEGIN 
                  
                  IF rev_ind[binnum] NE rev_ind[binnum+1] THEN BEGIN 
                      
                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                      ;; Don't call a function here: too slow in IDL
                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                      
                      ;; Get neighbors in this bin
                      w=rev_ind[ rev_ind[binnum]:rev_ind[binnum+1]-1 ]
                      nw = n_elements(w)
                      
                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                      ;; estimate of density contrast for each
                      ;; source galaxy. If photoz, find mean for
                      ;; each othersize use sig_inv already set above
                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                      IF use_photoz THEN BEGIN
                          ;; get mean density contrast
                          
                          sig_inv = sigcrit_inv[w]
                          
                          ;; sometimes the interpolation fails
                          w2 = where(sig_inv NE 0.0 AND finite(sig_inv), nw)
                          IF nw EQ 0 THEN BEGIN
                              print,'*',format='(a,$)'
                              GOTO,jump 
                          ENDIF ELSE BEGIN
                              sig_inv = sig_inv[w2]
                              w = w[w2]
                          ENDELSE 
                          
                      ENDIF
                      sig_inv2 = sig_inv^2
                      
                      ng = ng + nw
                      
                      ;; how many?
                      lensum[index].npair[binnum]=nw
                      lensum[index].totpairs=lensum[index].totpairs+$
                        lensum[index].npair[binnum]

                      ;; radial statistics
                      lensum[index].rmax_act[binnum] = max(R[w], min=minrr)
                      lensum[index].rmin_act[binnum] = minrr
                      lensum[index].rsum[binnum] = total(R[w])
                      
                      ;; X/Y positions
                      xrel =  R[w]*cos(theta[w])
                      yrel =  R[w]*sin(theta[w])
                                  
                      ;; eta is flipped
                      diffsq=xrel^2 - yrel^2
                      xy=xrel*yrel

                      cos2theta = diffsq/R[w]^2
                      sin2theta = 2.*xy/R[w]^2
                          
                      ;; Tangential/45-degree rotated ellipticities
                      e1prime = -(scat[wsrc[w]].(e1tag)*cos2theta + scat[wsrc[w]].(e2tag)*sin2theta)
                      e2prime =  (scat[wsrc[w]].(e1tag)*sin2theta - scat[wsrc[w]].(e2tag)*cos2theta)
                      
                      e1e2err2 = sign(scat[wsrc[w]].e1e2err)*scat[wsrc[w]].e1e2err^2
                      etan_err2 = (scat[wsrc[w]].e1e1err^2*cos2theta^2 + $
                                   scat[wsrc[w]].e2e2err^2*sin2theta^2 - $
                                   2.*cos2theta*sin2theta*e1e2err2 )
                      shear_err2 = (etan_err2 + shapenoise2)/4.0
                      
                      ortho_err2 = (scat[wsrc[w]].e1e1err^2*sin2theta^2 + $
                                    scat[wsrc[w]].e2e2err^2*cos2theta^2 - $
                                    2.*cos2theta*sin2theta*e1e2err2 )
                      orthoshear_err2 = (ortho_err2 + shapenoise2)/4.0
                      
                      ;; This is an array operation for photoz's
                      denscont = e1prime/2./sig_inv
                      densconterr2 = shear_err2/sig_inv2
                      
                      orthodenscont = e2prime/2./sig_inv
                      orthodensconterr2 = orthoshear_err2/sig_inv2
                      
                      wts = 1./densconterr2
                      wsum = total(wts)
                      lensum[index].sigma[binnum] = total(wts*denscont)/wsum
                      lensum[index].sigmaerr[binnum] = sqrt( 1./wsum )

                      ;; mean rotated by 45 degrees
                      owts = 1./orthodensconterr2
                      owsum = total(owts)
                      lensum[index].orthosig[binnum] = total(owts*orthodenscont)/owsum
                      lensum[index].orthosigerr[binnum] = sqrt( 1./owsum )
                      
                      ;; alternative error estimate
                      lensum[index].sigerrsum[binnum] = total(wts^2*denscont^2)
                      lensum[index].orthosigerrsum[binnum] = total(owts^2*orthodenscont^2)
                      
                      ;; sums of weights
                      lensum[index].wsum[binnum] = wsum
                      lensum[index].owsum[binnum] = owsum
                      
                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                      ;; shear polarizability
                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                      wts_ssh = wts
                          ;; old ways to calculate shear polarizability
;                          lensum[index].sshsum = $
;                            lensum[index].sshsum + $
;                            total( wts_ssh*(1.-shapenoise2*wts_ssh*e1prime^2) )
;                          lensum[index].sshsum = lensum[index].sshsum + total(wts_ssh*(1.-e1prime^2))

                      ;; fractional error contributions
                      tw = 1./(etan_err2 + shapenoise2)
                      f_e = etan_err2*tw
                      f_sn = shapenoise2*tw
                      
                      ;; coefficients (p 596 Bern02) 
                      ;; *** NOTE!!!!:  there is a k1*e^2/2 in Bern02 because
                      ;; its the total ellipticity he is using

                      k0 = f_e*shapenoise2
                      k1 = f_sn^2
                      F = 1. - k0 - k1*e1prime^2
;                          F = 1. - k0 - k1*etot2/2.0
;                          G = etot/2.0*(1. - k0 - k1*etot2)
                      
                      lensum[index].sshsum   = lensum[index].sshsum + total(wts_ssh*F)
;                          lensum[index].sshsum   = lensum[index].sshsum + total(wts_ssh*F + wts_prime*G)
                      lensum[index].wsum_ssh = lensum[index].wsum_ssh + total(wts_ssh)
                        
                      ;; This will be the overall weight for each
                      ;; lens, for use in averaging luminisoty, etc
                      lensum[index].weight = lensum[index].weight + wsum
                      
                      ;; distribution of neighbor positions
                      xpysum = xpysum + total(xrel^2 + yrel^2)
                      xmysum = xmysum + total(diffsq)
                      xysum  = xysum  + total(xy)
                      ;;stop
                      ;; free memory
                      diffsq=0 & xy=0 & e1prime=0 & e2prime=0 & wts=0 & w=0 & xrel=0 & yrel=0
                  ENDIF 
                  jump:
                  
              ENDFOR 
              
              ;; make sure it has a reasonably symmetric distribution
              ;; of sources behind it.  (this is what phil does)
              
              ie1 = xmysum/xpysum
              ie2 = 2.*xysum/xpysum
              ie = sqrt( ie1^2 + ie2^2 )
              lensum[index].ie = ie
              
              mm = 3./sqrt(ng)
              
              IF ie GT max([mm, maxe]) THEN BEGIN
                  ;; Not symmetric enough
                  ng=0 
              ENDIF ELSE BEGIN 
                  ;; OK
                  sigw_w = where(lensum[index].npair GT 0,nsigw_w)
                  IF nsigw_w NE 0 THEN BEGIN 
                      sigw = 1./lensum[index].sigmaerr[sigw_w]^2
                      sigwsum = sigwsum + total(sigw)
                      sigsum = sigsum + total(lensum[index].sigma[sigw_w]*sigw)
                      lamsum = lamsum + lensum[index].clambda
                      etasum = etasum + lensum[index].ceta
                      numsum = numsum + 1.
                      
                      sshsum = sshsum + lensum[index].sshsum
                  ENDIF 
                  
                  ;; print to ascii file
                  IF n_elements(outlun) NE 0 THEN BEGIN 
                      zobjshear_htm_looplens_printf, outlun, $
                        lensum, index, zlens, $
                        format
                  ENDIF 
                  
              ENDELSE 
              
          ENDIF 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; One final check.  If not passed, remember that this 
          ;; lens wasn't used 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF (ng EQ 0) THEN BEGIN
                            
              ;; none passed final cuts
              ;;echo,'/',color='green',/bold,/nonewline
              print,'/',format='(a,$)'
              indices[index] = -1
              
          ENDIF 
          
          
          ;; free some memory
          R = 0
          radiff=0
          
      ENDIF ELSE BEGIN

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; no sources found within leaflist
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          ;; No neighbors found
          ;;echo,'|',color='red',/bold,/nonewline
          print,'|',format='(a,$)'
          indices[ index ] = -1
      ENDELSE 

      ;;;;;;;;;;;;;;;;;;;;;;;
      ;; print info
      ;;;;;;;;;;;;;;;;;;;;;;;

      IF (( index MOD step ) EQ 0) OR (index EQ (nlens-1)) THEN BEGIN
          print
          IF sigwsum GT 0.0 THEN BEGIN 
              mean_denscont = sigsum/sigwsum
              err_denscont = sqrt(1./sigwsum)
              groupstruct.sigma = mean_denscont
              groupstruct.sigmaerr = err_denscont
              groupstruct.mean_lambda = lamsum/numsum
              groupstruct.mean_eta = etasum/numsum

              mean_denscont_str = ntostr(mean_denscont)
              err_denscont_str = ntostr(err_denscont)
              ssh_str = ntostr(sshsum/sigwsum)
          ENDIF ELSE BEGIN 
              mean_denscont_str = '0'
              err_denscont_str = '0'
              ssh_str = '0'
          ENDELSE 
          print,'Lens Group = ',ntostr(group+1)+'/'+ntostr(nstep),$
                ' Mean dens. cont. = ' + $
                mean_denscont_str+' '+!plusminus+' '+err_denscont_str+$
                '   Ssh = '+ssh_str
          group=group+1
      ENDIF 

      IF (index MOD 10) EQ 0 THEN  print,'.',format='(a,$)'

  ENDFOR  

END 
