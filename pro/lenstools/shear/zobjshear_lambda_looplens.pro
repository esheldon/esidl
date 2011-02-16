;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    ZOBJSHEAR_LAMBDA_LOOPLENS
;       
; PURPOSE:
;    Tool for looping over lenses and measureing the shear
;
; CALLING SEQUENCE:
;    zobjshear_looplens, scat, lclambda, lceta, angmax, DL, $
;                       sigcritinv, step, maxe, wlens,$
;                       rmin, rmax, nbin, $
;                       indices, $
;                       random=random, mineta=mineta, maxeta=maxeta, $
;                       logbin=logbin, issouth=issouth, maxit=maxit
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
;    17-JAN-2002
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO zobjshear_lambda_looplens, lensum, scat, lclambda, lceta, angmax, DL, $
                               step, maxe, wlens,$
                               rmin, rmax, nbin_or_binsize, $
                               groupstruct, indices, $
                               random=random, mineta=mineta, maxeta=maxeta, $
                               logbin=logbin, issouth=issouth, maxit=maxit, $
                               prior_zs=prior_zs,  prior_pofzs=prior_pofzs, $
                               scinv_struct=scinv_struct, $
                               recorr=recorr

  IF n_params() LT 12 THEN BEGIN 
      print,'-Syntax: zobjshear_looplens, scat, lclambda, lceta, step, wlens, indices, $'
      print,'         random=random, mineta=mineta, maxeta=maxeta,$'
      print,'         logbin=logbin, issouth=issouth, maxit=maxit, prior_zs=prior_zs,  prior_pofzs=prior_pofzs, scinv_struct=scinv_struct, $'
      print,'         recorr=recorr'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;
  
  ;; shapenoise
  shapenoise = 0.32
  shapenoise2 = shapenoise^2

  ;; number of points in the mean sigmacrit 
  ;; calculation for photoz's
  npts = 200

  IF n_elements(prior_zs) EQ 0 THEN use_photoz=0 ELSE use_photoz=1
  IF n_elements(scinv_struct) EQ 0 THEN use_photoz=0 ELSE use_photoz=1

  wz=getztag(lensum[0])

  IF keyword_set(logbin) THEN BEGIN
      nbin = nbin_or_binsize
  ENDIF ELSE BEGIN
      binsize = nbin_or_binsize
      nbin = long( (rmax - rmin)/binsize ) ;+ 1 
  ENDELSE 

  IF keyword_set(recorr) THEN BEGIN 
      IF NOT tag_exist(scat,'e1_recorr',index=e1tag) THEN message,'No e1_recorr tag'
      IF NOT tag_exist(scat,'e2_recorr',index=e2tag) THEN message,'No e2_recorr tag'
  ENDIF ELSE BEGIN 
      IF NOT tag_exist(scat,'e1',index=e1tag) THEN message,'No e1 tag'
      IF NOT tag_exist(scat,'e2',index=e2tag) THEN message,'No e2 tag'
  ENDELSE 

  IF NOT tag_exist(scat,'clambda') THEN message,'scat must contain CLAMBDA and CETA tags'

  ;; print out
  print,'In ZOBJSHEAR_LOOPLENS'
  IF keyword_set(random) THEN print,'Doing random'
  IF keyword_set(logbin) THEN print,'Doing logbin'
  print

  IF n_elements(maxit) EQ 0 THEN maxit=100
  nlens = n_elements(wlens)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; here is where zobjshear_lambda_looplens and 
; zobjshear_htm_looplens differ
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
                              'sigma',replicate(-1.e10,nstep),$
                              'sigmaerr',replicate(-1.e10,nstep),$
                              'mean_lambda',replicate(-1.e10,nstep),$
                              'mean_eta',replicate(-1.e10,nstep))

  indices = lindgen(nlens)
  sigwsum=0.
  sigsum=0.

  FOR group = 0L, nstep-1 DO BEGIN

      ;;tgrp=systime(1)
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

      maxlam1 = lclambda[maxii]+angmax2
      minlam1 = lclambda[minii]-angmax2

      zobjshear_lambda_get_sources, scat.clambda, $
                                    minlam1, maxlam1, wsrc, nwsrc

      ;; for group statistics
      groupused=0L & lamsum=0d & etasum=0d
      IF nwsrc NE 0 THEN BEGIN 
          
          ;; generate eta's
          ;; Stripe is not rectangular in lambda-eta, this helps to 
          ;; keep random points from being thrown out

          IF keyword_set(random) THEN BEGIN 
              FOR iii=0L, step-1 DO BEGIN 
                  lceta[ii[iii]] = arrscl( randomu(seed, /double), $
                                          mineta[ii[iii]], maxeta[ii[iii]],$
                                          arrmin=0.0, arrmax=1.0)
              ENDFOR 
          ENDIF 
          
          zobjshear_looplens_print_status, group, nstep, sigsum, sigwsum
          FOR gi=0L, step-1 DO BEGIN ;Loop over lenses in group

              ;;tt=systime(1)

              IF (gi MOD 10) EQ 0 THEN  print,'.',format='(a,$)'
              index = ii[gi]

              zlens = lensum[index].(wz)

              ceneta = lceta[index]
              cenlam  = lclambda[index]
              
              sig_inv = lensum[index].scritinv
              sig_crit = 1./sig_inv
              IF sig_inv EQ -1000. THEN BEGIN 
                  print,'What!'
                  return
              ENDIF 
              angmax_i = angmax[index]
              
              ;; choose sources around this lens
              maxlam2 = cenlam+angmax_i
              minlam2 = cenlam-angmax_i

              zobjshear_lambda_get_sources, scat[wsrc].clambda, $
                                            minlam2, maxlam2, $
                                            twsrc, nwsrc2
              
              ;; If there are any sources left, measure the shear
              IF nwsrc2 GT 1 THEN BEGIN 

                  wsrc2 = wsrc[twsrc]

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; Try to redo a few times until passes checks
                  ;; ** if we are doing random **
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  redo=0
                  WHILE (redo NE -1) AND (redo LT maxit) DO BEGIN 
                      redo=redo+1

                      mygcirc_survey, cenlam, ceneta, $
                        scat[wsrc2].clambda, scat[wsrc2].ceta, $
                        R, theta, /radians_out
                      R = R*DL[index]

                      IF use_photoz THEN BEGIN 
                          zlens_send = replicate(zlens,nwsrc2)
                          photoz_z = scat[wsrc2].photoz_z
                          photoz_zerr = scat[wsrc2].photoz_zerr

                          sigcrit_inv = interp_meanscinv(zlens_send, $
                                                         photoz_z, $
                                                         photoz_zerr, $
                                                         scinv_struct)
                      ENDIF 

                      IF NOT keyword_set(logbin) THEN BEGIN 
                          hist=histogram(R, binsize=nbin_or_binsize, min=rmin, $
                                         max=rmax,rever=rev_ind)
                      ENDIF ELSE BEGIN 
                          logbin, R, rmin, rmax, nbin_or_binsize, hist, rev_ind
                      ENDELSE 

                      whist = where(hist[0:nbin-1] NE 0, nhist)
                      ng=0
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
                                      
                                      w2 = where(sig_inv NE 0.0 AND finite(sig_inv), nw)
                                      IF nw EQ 0 THEN BEGIN
                                          GOTO,jump 
                                      ENDIF ELSE BEGIN
                                          sig_inv = sig_inv[w2]
                                          w = w[w2]
                                      ENDELSE 
                                      
                                  ENDIF

                                  ;; how many?
                                  lensum[index].npair[binnum]=n_elements(w)
                                  lensum[index].totpairs=lensum[index].totpairs+$
                                    lensum[index].npair[binnum]
                                  
                                  ;; radial statistics
                                  lensum[index].rmax_act[binnum] = $
                                    max(R[w], min=minrr)
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
                                  e1prime = -(scat[wsrc2[w]].(e1tag)*cos2theta + scat[wsrc2[w]].(e2tag)*sin2theta)
                                  e2prime =  (scat[wsrc2[w]].(e1tag)*sin2theta - scat[wsrc2[w]].(e2tag)*cos2theta)

                                  e1e2err2 = sign(scat[wsrc2[w]].e1e2err)*scat[wsrc2[w]].e1e2err^2
                                  etan_err2 = (scat[wsrc2[w]].e1e1err^2*cos2theta^2 + $
                                               scat[wsrc2[w]].e2e2err^2*sin2theta^2 - $
                                               2.*cos2theta*sin2theta*e1e2err2 ) + shapenoise2

                                  ortho_err2 = (scat[wsrc[w]].e1e1err^2*sin2theta^2 + $
                                                scat[wsrc[w]].e2e2err^2*cos2theta^2 - $
                                                2.*cos2theta*sin2theta*e1e2err2 ) + shapenoise2

                                  ;; This is an array operation for photoz's
                                  denscont = e1prime/2./sig_inv
                                  densconterr2 = etan_err2/(4.*sig_inv^2)
                                  
                                  orthodenscont = e2prime/2./sig_inv
                                  orthodensconterr2 = ortho_err2/(4.*sig_inv^2)

                                  ;; mean density contrast over sources
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
                                  lensum[index].sigerrsum = total(wts^2*denscont^2)
                                  lensum[index].orthosigerrsum = total(owts^2*orthodenscont^2)

;                                  lensum[index].sshsum = $
;                                    lensum[index].sshsum + $
;                                    total( wts_ssh*(1.-shapenoise2*wts_ssh*e1prime^2) )
                                  
                                  wts_ssh = wts
                                  lensum[index].sshsum = lensum[index].sshsum + total(wts_ssh*(1.-e1prime^2))

                                  lensum[index].wsum[binnum] = wsum
                                  lensum[index].wsum_ssh = lensum[index].wsum_ssh + total(wts_ssh)
                                  
                                  ;; This will be the overall weight for each
                                  ;; lens, for use in averaging luminisoty, etc
                                  lensum[index].weight = lensum[index].weight + wsum

                                  ;; distribution of neighbor positions
                                  xpysum = xpysum + total(xrel^2 + yrel^2)
                                  xmysum = xmysum + total(diffsq)
                                  xysum  = xysum  + total(xy)

                                  ;; free memory
                                  diffsq=0 & xy=0 & e1prime=0 & e2prime=0 & wts=0 & w=0 & xrel=0 & yrel=0
                              ENDIF 
                              jump:
                          ENDFOR 

                          ;; make sure it has a reasonably symmetric distribution
                          ;; of sources behind it.  (this is what phil does)

                          ng=total(hist[0:nbin-1])
                          ie1 = xmysum/xpysum
                          ie2 = 2.*xysum/xpysum
                          ie = sqrt( ie1^2 + ie2^2 )
                          lensum[index].ie = ie

                          mm = 3./sqrt(ng)

                          IF ie GT max([mm, maxe]) THEN ng=0 ELSE BEGIN 
                              sigw_w = where(lensum[index].npair GT 0,nsigw_w)
                              IF nsigw_w NE 0 THEN BEGIN 
                                  sigw = 1./lensum[index].sigmaerr[sigw_w]^2
                                  sigwsum = sigwsum + total(sigw)
                                  sigsum = sigsum + total(lensum[index].sigma[sigw_w]*sigw)
                              ENDIF 
                          ENDELSE 

                      ENDIF 
                      ;; One final check.  If not passed, remember that this $
                      ;; lens wasn't used ** if doing random then redo **
                      IF (ng EQ 0) AND keyword_set(random) THEN BEGIN 
                          ;; reset this lens
                          zobjshear_setlensumzero, lensum, index

                          qeta = arrscl( randomu(seed, /double), $
                                         mineta[index], maxeta[index], $
                                          arrmin=0.0,arrmax=1.0 )
                          lceta[index] = qeta
                          ceneta = qeta
                      ENDIF ELSE BEGIN 
                          redo=-1
                      ENDELSE 
                      R = 0
                      radiff=0  ;Free memory
                  ENDWHILE ;; redo while for random
                  IF ng EQ 0 THEN BEGIN
                      ;;echo,'/',color='green',/bold,/nonewline
                      print,'/',format='(a,$)'
                      indices[ind[gi]] = -1
                  ENDIF ELSE BEGIN 
                      groupused=groupused+1 ;for group statistics
                      lamsum=lamsum+cenlam
                      etasum=etasum+ceneta
                  ENDELSE 
              ENDIF ELSE BEGIN
                  ;;echo,'|',color='red',/bold,/nonewline
                  print,'|',format='(a,$)'
                  indices[ ind[gi] ] = -1
              ENDELSE 
              ;;print,systime(1)-tt
              ;;stop

          ENDFOR 
          print
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; fill in group statistics
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF sigwsum GT 0.0 THEN BEGIN 
              mean_denscont = sigsum/sigwsum 
              err_denscont = sqrt(1./sigwsum)
          ENDIF ELSE BEGIN 
              mean_denscont = 0.0
              err_denscont = 0.0
          ENDELSE 
          
          groupstruct.nlens[group] = groupused
          groupstruct.sigma[group] = mean_denscont
          groupstruct.sigmaerr[group] = err_denscont
          groupstruct.mean_lambda[group] = lamsum/groupused
          groupstruct.mean_eta[group] = etasum/groupused
      ENDIF ELSE BEGIN 
          print,'Two-'
          indices[ ind ] = -1
      ENDELSE 

  ENDFOR 

END 
