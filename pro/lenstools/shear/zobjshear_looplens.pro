;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    ZOBJSHEAR_LOOPLENS
;       
; PURPOSE:
;    Tool for looping over lenses and measureing the shear
;
; CALLING SEQUENCE:
;    zobjshear_looplens, scat, llambda, leta, angmax, DL, $
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


PRO zobjshear_looplens, lensum, scat, lra, ldec, angmax, DL, $
                              step, maxe, wlens,$
                              rmin, rmax, nbin_or_binsize, $
                              groupsruct, indices, $
                              random=random, $
                              logbin=logbin, issouth=issouth, maxit=maxit
  
  IF n_params() LT 12 THEN BEGIN 
      print,'-Syntax: zobjshear_looplens, scat, llambda, leta, step, wlens, indices, $'
      print,'         random=random, mineta=mineta, maxeta=maxeta,$'
      print,'         logbin=logbin, issouth=issouth, maxit=maxit'
      return
  ENDIF 


  d2r = !dpi/180.0d0            ;Change degrees to radians.
  r2d = 180./!dpi               ;radians to degrees
  vint = 0.32^2

  IF keyword_set(logbin) THEN BEGIN
      nbin = nbin_or_binsize
  ENDIF ELSE BEGIN
      binsize = nbin_or_binsize
      nbin = long( (rmax - rmin)/binsize ) ;+ 1 
  ENDELSE 

  IF NOT tag_exist(scat,'momerr',index=momerr) THEN BEGIN
      IF NOT tag_exist(scat,'uncert',index=momerr) THEN BEGIN 
          print
          print,'No valid moment uncertainty tag'
          print
          return
      ENDIF 
  ENDIF 

  nscat=n_elements(scat)
  LASTRA = scat[nscat-1].ra

  print,'In ZOBJSHEAR_LOOPLENS'
  IF keyword_set(random) THEN print,'Doing random'
  IF keyword_set(logbin) THEN print,'Doing logbin'
  print
  IF n_elements(maxit) EQ 0 THEN maxit=100

;  maxdec = max(scat.dec, min=mindec)

  maxdec = max(ldec, min=mindec)

  nlens = n_elements(wlens)

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
                              'mean_ra',replicate(-1.e10,nstep),$
                              'mean_dec',replicate(-1.e10,nstep))

  indices = lindgen(nlens)
  etansum=0.
  etanerrsum=0.
  wsum=0.

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

      maxra1 = lra[maxii]+angmax2
      minra1 = lra[minii]-angmax2

      zobjshear_get_sources, scat.ra, minra1, maxra1, LASTRA, wsrc, nwsrc, $
        issouth=issouth

      ;; for group statistics
      groupused=0L & rasum=0d & decsum=0d
      IF nwsrc NE 0 THEN BEGIN 
                    
          zobjshear_looplens_print_status, group, nstep, wsum, etansum, etanerrsum
          FOR gi=0L, step-1 DO BEGIN ;Loop over lenses in group

              ;;tt=systime(1)

              IF (gi MOD 10) EQ 0 THEN  print,'.',format='(a,$)'
              index = ii[gi]
              cendec = ldec[index]
              cenra  = lra[index]
              
              sig_inv = lensum[index].scritinv
              sig_crit = 1./sig_inv
              IF sig_inv EQ -1000. THEN BEGIN 
                  print,'What!'
                  return
              ENDIF 
              angmax_i = angmax[index]
              
              ;; choose sources around this lens
              maxra2 = cenra+angmax_i
              minra2 = cenra-angmax_i

              zobjshear_get_sources, scat[wsrc].ra, minra2, maxra2, LASTRA, $
                twsrc, nwsrc2, $
                issouth=issouth
              
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

                      mygcirc, cenra, cendec, $
                        scat[wsrc2].ra, scat[wsrc2].dec, $
                        R, theta, /radians_out
                      R = R*DL[index]

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
                                  
                                  w=rev_ind[ rev_ind[binnum]:rev_ind[binnum+1]-1 ]
                                  
                                  lensum[index].npair[binnum]=n_elements(w)
                                  lensum[index].totpairs=lensum[index].totpairs+$
                                    lensum[index].npair[binnum]
                                  
                                  lensum[index].rmax_act[binnum] = $
                                    max(R[w], min=minrr)
                                  lensum[index].rmin_act[binnum] = minrr

                                  ;;IF NOT issouth THEN BEGIN 
                                      xrel =  R[w]*cos(theta[w])
                                  ;;ENDIF ELSE BEGIN 
                                      ;; flipped in the south
                                  ;;    xrel = -R[w]*cos(theta[w])
                                  ;;ENDELSE 
                                  yrel =  R[w]*sin(theta[w])
                                  
                                  diffsq=xrel^2 - yrel^2
                                  xy=xrel*yrel

                                  e1prime=-(scat[wsrc2[w]].e1*diffsq + $
                                            scat[wsrc2[w]].e2*2.*xy  )/R[w]^2
                                  e2prime= (scat[wsrc2[w]].e1*2.*xy - $
                                            scat[wsrc2[w]].e2*diffsq )/R[w]^2
                                  
                                  ;; also weighting by 1/sigmacrit
                                  wts_ssh = 1./( vint + scat[wsrc2[w]].(momerr[0])^2)
                                  wts = wts_ssh*sig_inv^2
                                  
                                  lensum[index].rsum[binnum] = total(R[w])

                                  lensum[index].etansum[binnum] = $
                                    total(e1prime*wts)
                                  lensum[index].eradsum[binnum] = $
                                    total(e2prime*wts)

                                  lensum[index].etanerrsum[binnum] = $
                                    total(wts^2*e1prime^2)
                                  lensum[index].eraderrsum[binnum] = $
                                    total(wts^2*e2prime^2)
                                  
                                  ;; Ssh weighted differently: no weight by sig_crit
                                  lensum[index].sshsum = $
                                    lensum[index].sshsum + $
                                    total( wts_ssh*(1.-vint*wts_ssh*e1prime^2) )

                                  lensum[index].wsum[binnum] = total(wts)
                                  lensum[index].wsum_ssh = $
                                    lensum[index].wsum_ssh + $
                                    total(wts_ssh)
                                  
                                  xpysum = xpysum + total(xrel^2 + yrel^2)
                                  xmysum = xmysum + total(diffsq)
                                  xysum  = xysum  + total(xy)

                                  ;; free memory
                                  diffsq=0 & xy=0 & e1prime=0 & e2prime=0 & wts=0 & w=0 & xrel=0 & yrel=0
                              ENDIF 
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
                              etansum = etansum + $
                                total(lensum[index].etansum)
                              etanerrsum = etanerrsum + $
                                total(lensum[index].etanerrsum)
                              wsum = wsum + total(lensum[index].wsum)
                          ENDELSE 

                      ENDIF 
                      ;; One final check.  If not passed, remember that this 
                      ;; lens wasn't used ** if doing random then redo **
                      IF (ng EQ 0) AND keyword_set(random) THEN BEGIN 
                          ;;!!!!!!!!!!!!!!!!!!!! FIX THIS !!!!!!!!!!!!
                          ;; reset these sums
                          lensum[index].sshsum = 0.
                          lensum[index].wsum_ssh = 0.

                          qdec = arrscl( randomu(seed), $
                                         mindec, maxdec, $
                                         arrmin=0.0,arrmax=1.0 )
                          ldec[index] = qdec
                          cendec = qdec
                      ENDIF ELSE BEGIN 
                          redo=-1
                      ENDELSE 
                      R = 0
                      radiff=0  ;Free memory
                  ENDWHILE ;; redo while for random
                  IF ng EQ 0 THEN BEGIN
                      echo,'/',color='green',/bold,/nonewline
                      indices[ind[gi]] = -1
                  ENDIF ELSE BEGIN 
                      groupused=groupused+1 ;for group statistics
                      rasum=rasum+cenra
                      decsum=decsum+cendec
                  ENDELSE 
              ENDIF ELSE BEGIN
                  echo,'|',color='red',/bold,/nonewline
                  indices[ ind[gi] ] = -1
              ENDELSE 
              ;;print,systime(1)-tt
              ;;stop

          ENDFOR 
          print
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; fill in group statistics
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF wsum GT 0. THEN BEGIN 
              mean_shear = etansum/wsum/2. 
              mean_shear_err = sqrt(etanerrsum)/wsum/2.
          endif else begin 
              mean_shear=0.0
              mean_shear_err = 0.0
          endelse 
          
          groupstruct.nlens[group] = groupused
          groupstruct.shear[group] = mean_shear
          groupstruct.shearerr[group] = mean_shear_err
          groupstruct.mean_ra[group] = rasum/groupused
          groupstruct.mean_dec[group] = decsum/groupused
      ENDIF ELSE BEGIN 
          print,'Two-'
          indices[ ind ] = -1
      ENDELSE 
      mradiff = 0               ;Free memory
      mdist = 0
      ;;ptime,systime(1)-tgrp
  ENDFOR 


END 
