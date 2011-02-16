;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
; THIS CODE IS BOGUS: YOU SHOULD NOT INCLUDE THE COVARIANCE BETWEEN
; BANDPASSES. Use combine_ellip_cove1e2.pro
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; combine measurements from different bandpasses, accounting for 
; the covariance between bandpasses.  Returns the combined values plus the
; indices of objects that had at least one good ellipticity
; measurement among the clrs input
;
; objects without a good measurement will have -9999. for all
; values, and will have combine_flag == 0
;
; recommend sending a catalog with good measurements in r-band
; for every object
;
; NOTE: Unlike combine_ellip_cove1e2, this is specifically designed
;       to deal with combining g,r,i exposures


PRO ce_calc_cov_3band, corr, cove1e2, cov3

  FOR i=0,1 DO BEGIN 
      FOR j=0,1 DO BEGIN 
          FOR k=0,2 DO BEGIN 
              FOR l=0,2 DO BEGIN 
                  
                  ind1 = 3*i+k
                  ind2 = 3*j+l
                  ;;IF k NE l THEN BEGIN 
                      ;; k ne l
                      cov3[ind1,ind2] = $
                        corr[ind1,ind2]*sqrt(cove1e2[k,i,i]*cove1e2[l,j,j])
                      ;;cov3[ind1,ind2] = 0.0

                      ;;cov3[ind1,ind2] = $
                      ;;  corr[ind1,ind2]*min([cove1e2[k,i,i],cove1e2[l,j,j]])
                  ;;ENDIF ELSE BEGIN 
                      ;; k=l
                  ;;    cov3[ind1,ind2] = cove1e2[k,i,j]
                  ;;ENDELSE 
                  ;;IF (i NE j) AND (k NE l) THEN cov3[ind1,ind2]=0.0
                  ;;IF (i NE j) THEN cov3[ind1,ind2]=0.0
              ENDFOR 
          ENDFOR 
      ENDFOR 
  ENDFOR 

END 


PRO ce_calc_cov, corr, cove1e2, wclr, cov

  nclr = n_elements(wclr)
  IF nclr LT 2 THEN message,"Don't call ce_calc_cov unless nclr = 2 or 3!"

  FOR i=0,1 DO BEGIN 
      FOR j=0,1 DO BEGIN 
          FOR tk=0,nclr-1 DO BEGIN 
              k = wclr[tk]

              ind1 = nclr*i+tk
              cind1 = 3*i+k

              FOR tl=0,nclr-1 DO BEGIN 
                  l = wclr[tl]

                  ind2 = nclr*j+tl
                  cind2 = 3*j+l

                  ;; for nclr=3, ind1 = cind1
                  ;; for nclr=2, must subscript corr and cove1e2, which 
                  ;; have info for 3 bands, with the absolute color k,l and 
                  ;; subscript the cov with the relative values

                  ;;IF (ind1 NE cind1) OR (ind2 NE cind2) THEN message,'WHAT!!!'
                  IF k NE l THEN BEGIN 
                      ;; k ne l
                      cov[ind1,ind2] = $
                        corr[cind1,cind2]*sqrt(cove1e2[k,i,i]*cove1e2[l,j,j])
                      ;;cov[ind1,ind2] = 0.0

                      ;;cov[ind1,ind2] = $
                      ;;  corr[cind1,cind2]*min([cove1e2[k,i,i],cove1e2[l,j,j]])
                  ENDIF ELSE BEGIN 
                      ;; case of k=l: copy in covariance term
                      cov[ind1,ind2] = cove1e2[k,i,j]
                  ENDELSE 
              ENDFOR 
          ENDFOR 
      ENDFOR 
  ENDFOR 

  ;; make sure all off-diagonal terms are less than smallest
  ;; diagonal term

  ind = lindgen(nclr*nclr)
  ix = ind MOD nclr
  iy = ind/nclr

  offdiag = where(ix NE iy)
  diag = where(ix EQ iy)

  mindiag = min( cov[ix[diag], iy[diag]] )*0.9
  maxoffdiag = max( cov[ix[offdiag], iy[offdiag]] )

  factor = mindiag/maxoffdiag < 1.0

  cov[ix[offdiag], iy[offdiag]] = cov[ix[offdiag], iy[offdiag]]*factor

END 

FUNCTION ce_adds, S, i, j

  siz=size(S)
  IF siz[1] NE siz[2] THEN message,'S must be NXN'

  nclr = long(siz[1]/2.)
  IF nclr EQ 3 THEN BEGIN 
      return,$
             S[3*i+0,3*j+0] + S[3*i+0,3*j+1] + S[3*i+0,3*j+2] + $
             S[3*i+1,3*j+0] + S[3*i+1,3*j+1] + S[3*i+1,3*j+2] + $
             S[3*i+2,3*j+0] + S[3*i+2,3*j+1] + S[3*i+2,3*j+2]
  ENDIF ELSE IF nclr EQ 2 THEN BEGIN  
      return,$
             S[2*i+0,2*j+0] + S[2*i+0,2*j+1] + $
             S[2*i+1,2*j+0] + S[2*i+1,2*j+1]
  ENDIF ELSE message,"Don't call ce_adds unless nclr=2 or 3!"

END 

FUNCTION ce_adds_times_e, S, i, j, e

  siz=size(S)
  IF siz[1] NE siz[2] THEN message,'S must be NXN'

  nclr = long(siz[1]/2.)
  IF n_elements(e) NE nclr THEN message,'e must be nclr long'

  IF nclr EQ 3 THEN BEGIN 
      return,$
             S[3*i+0,3*j+0]*e[0] + S[3*i+0,3*j+1]*e[0] + S[3*i+0,3*j+2]*e[0] + $
             S[3*i+1,3*j+0]*e[1] + S[3*i+1,3*j+1]*e[1] + S[3*i+1,3*j+2]*e[1] + $
             S[3*i+2,3*j+0]*e[2] + S[3*i+2,3*j+1]*e[2] + S[3*i+2,3*j+2]*e[2]
  ENDIF ELSE IF nclr EQ 2 THEN BEGIN  
      return,$
             S[2*i+0,2*j+0]*e[0] + S[2*i+0,2*j+1]*e[0] + $
             S[2*i+1,2*j+0]*e[1] + S[2*i+1,2*j+1]*e[1]
  ENDIF ELSE message,"Don't call ce_adds_times_e unless nclr=2 or 3!"

END 

PRO combine_ellip_fullcov, e1, e2, e1e1err, e1e2err, e2e2err, corr, $
                           new_e1, new_e2, $
                           new_e1e1err, new_e1e2err, new_e2e2err,$
                           combine_flag, good, $
                           recorr=recorr, $
                           verbose=verbose, e1wmom=e1wmom, e2wmom=e2wmom

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: '
      return
  ENDIF 

  time=systime(1)

  IF n_elements(verbose) EQ 0 THEN verbose = 1

  clrs = [1,2,3]                ;lookup table for the colors we use
  nclr = 3
  siz=size(corr)
  IF (siz[1] NE 6) OR (siz[2] NE 6) THEN message,$
    'correlation matrix must be 6X6'

;  ENDELSE 

;; ERROR TAGS, DEAL WITH ALL THE POSSIBILITIES, SAME WITH ELLIP

  siz=size(e1)
  IF siz[1] NE 3 THEN message,'Must enter arrays of [3,n]'

  nobj = siz[2]
  print,'Number of objects: ',nobj

  defval = -9999.
  arrval = replicate(defval, nobj)
  new_e1 = arrval
  new_e2 = arrval
  new_e1e1err = arrval
  new_e1e2err = arrval
  new_e2e2err = arrval

  e1wmom=arrval
  e2wmom=arrval

  testind = replicate(-1, nobj)


  ;; diagnostic about combination
  ;; 2L^clr will be added if clr is used
  combine_flag = bytarr(nobj)

  ;; convert to double
  corr = double(corr)

  ;; error codesfrom the invert function
  errcode=['Successful completion',$
           'Singular Array', $
           'Small pivot element: Accuracy Compromised']

  ;; contains cov(e1,e2) for each bandpass
  ;; of object
  cove1e2 = dblarr(nclr, 2, 2)
  corrfac3 = 300.               ;multiply covariance matrix by
                                ;this factor to help with precision
  
  cov3 = dblarr(6, 6)           ;when we have 3 good measurements
  cov2 = dblarr(4, 4)           ;when we have 2 good measurements

  agoodclr = intarr(nclr)
  bigsig = fltarr(2,2)
  beste = fltarr(2)

  FOR obj=0L, nobj-1 DO BEGIN 

      agoodclr[*] = 0
      bigsig[*,*] = 0.0
      beste[*] = 0.0
      cove1e2[*] = 0.0
      cov3[*] = 0.0
      cov2[*] = 0.0

      ;; figure out which bandpasses have good measurements
      ;; for this object.

      FOR ic=0L, nclr-1 DO BEGIN 
          clr = clrs[ic]
          ;; will we use this bandpass?
          IF e1[ic,obj] GT defval THEN BEGIN
              agoodclr[ic] = 1
              combine_flag[obj] = combine_flag[obj] + 2b^clr

              ;; copy in covariance
              cove1e2[ic,0,0] = e1e1err[ic,obj]^2
              cove1e2[ic,0,1] = sign(e1e2err[ic,obj])*e1e2err[ic,obj]^2
              cove1e2[ic,1,0] = cove1e2[ic,0,1]
              cove1e2[ic,1,1] = e2e2err[ic,obj]^2

          ENDIF 
      ENDFOR 

      wclr = where(agoodclr NE 0,ngclr)
      ;; can only deal with up to 3 for now
      IF ngclr EQ 1 THEN BEGIN 
          ;; only one good measurement
          IF verbose GT 2 THEN print,'Only one good measurement'
          good_clr = wclr[0]

          new_e1[obj] = e1[good_clr,obj]
          new_e2[obj] = e2[good_clr,obj]
          new_e1e1err[obj] = e1e1err[good_clr,obj]
          new_e1e2err[obj] = e1e2err[good_clr,obj]
          new_e2e2err[obj] = e2e2err[good_clr,obj]

          testind[obj] = 1
      ENDIF ELSE IF ngclr GT 0 THEN BEGIN 

          IF ngclr EQ 3 THEN BEGIN 
;              ce_calc_cov_3band, corr, cove1e2, cov3
;              print,cov3
;              print

              ce_calc_cov, corr, cove1e2, wclr, cov3
              IF verbose GT 2 THEN print,cov3

              ;; scale it up so determinant not so small!
              cov3 = cov3*corrfac3

              S = invert(cov3, istat, /double)
          ENDIF ELSE IF ngclr EQ 2 THEN BEGIN 
              IF verbose GT 2 THEN print,'2 good measurements'
              ce_calc_cov, corr, cove1e2, wclr, cov2
              IF verbose GT 2 THEN print,cov2

              cov2 = cov2*corrfac3

              S = invert(cov2, istat, /double)
          ENDIF

          IF istat EQ 0 THEN BEGIN 

              i=0 & j=0
              A = ce_adds(S, i, j)
              i=0 & j=1
              B = ce_adds(S, i, j)
              i=1 & j=1
              D = ce_adds(S, i, j)

              C = $
                ce_adds_times_e(S, 0, 0, e1[wclr,obj]) + $
                ce_adds_times_e(S, 0, 1, e2[wclr,obj])

              E = $
                ce_adds_times_e(S, 1, 1, e2[wclr,obj]) + $
                ce_adds_times_e(S, 0, 1, e1[wclr,obj])


              det = A*D - B*B

              IF verbose GT 2 THEN BEGIN 
                  print,'A=',A
                  print,'B=',B
                  print,'D=',D
                  print,'C=',C
                  print,'E=',E
              ENDIF 

              IF det EQ 0.0 THEN BEGIN 
                  IF verbose gt 1 THEN BEGIN 
                      print
                      print,'New covariance matrix is singular'
                  ENDIF 
              ENDIF ELSE IF det LT 0.0 THEN BEGIN 
                  IF verbose gt 1 THEN BEGIN 
                      print
                      print,'New covariance matrix has det < 0.0'
                  ENDIF 
                  ;stop
              ENDIF ELSE IF (A LT 0.0) OR (D LT 0.0) THEN BEGIN 
                  IF verbose gt 1 THEN BEGIN 
                      print
                      print,'Band A or D value'
                  ENDIF 
              ENDIF ELSE BEGIN 
                  factor = 1./det/corrfac3
                  new_e1e1err[obj] =  sqrt(factor*D)
                  tmp = -factor*B
                  new_e1e2err[obj] = sign(tmp)*sqrt(abs(tmp))
                  new_e2e2err[obj] =  sqrt(factor*A)

                  new_e1[obj] = (C*D - B*E)/det
                  new_e2[obj] = (A*E - C*B)/det

                  testind[obj] = 1
              ENDELSE 

          ENDIF ELSE BEGIN
              IF verbose gt 1 THEN BEGIN 
                  print
                  print,'Error inverting covariance matrix for obj #'+$
                        ntostr(obj)+': '+errcode[istat]
              ENDIF 
          ENDELSE 
          ;;wmom, e1[wclr,obj], sqrt(0.3^2 + e1e1err[wclr,obj]^2), wme1,wsig,wmerr1
          ;;wmom, e2[wclr,obj], sqrt(0.3^2 + e2e2err[wclr,obj]^2), wme2,wsig,wmerr2
          wmom, e1[wclr,obj], e1e1err[wclr,obj], wme1,wsig,wmerr1
          wmom, e2[wclr,obj], e2e2err[wclr,obj], wme2,wsig,wmerr2
          e1wmom[obj] = wme1
          e2wmom[obj] = wme2

          IF verbose GT 2 THEN BEGIN 
              print,e1[wclr,obj]
              print,e1e1err[wclr,obj]

              print,'Means e1: wmom = ',wme1,' New = ',new_e1[obj]
              
              print,e2[wclr,obj]
              print,e2e2err[wclr,obj]

              print,'Means e2: wmom = ',wme2,' New = ',new_e2[obj]
              
              ;; print,new_e1[obj],new_e2[obj]
              print
              !p.multi=[0,0,2]
              xvals = [(5-ngclr)+lindgen(ngclr)]

              yrange = prange(e1[wclr,obj], e1e1err[wclr,obj])
              plot,[0],/nodata,xrange=[0,5],yrange=yrange
              oploterror,xvals,$
                         [e1[wclr,obj]],$
                         [e1e1err[wclr,obj]],psym=4
              oploterror,[2.5],[new_e1[obj]],[new_e1e1err[obj]],psym=5,color=!red
              oploterror,[2.25],[wme1],[wmerr1],psym=4,color=!blue

              yrange = prange(e2[wclr,obj], e2e2err[wclr,obj])
              plot,[0],/nodata,xrange=[0,5],yrange=yrange
              ;ploterror,[1,2,3],e2[wclr,obj],e2e2err[wclr,obj],psym=4
              oploterror,xvals,$
                         [e2[wclr,obj]],$
                         [e2e2err[wclr,obj]],psym=4
              oploterror,[2.5],[new_e2[obj]],[new_e2e2err[obj]],psym=5,color=!red
              oploterror,[2.25],[wme2],[wmerr2],psym=4,color=!blue
              !p.multi=0

              key=get_kbrd(1)
          ENDIF 
      ENDIF ELSE BEGIN 
          IF verbose GT 1 THEN BEGIN 
              print,'No good measurements found for this obj!'
          ENDIF 
      ENDELSE 

      IF ( ( (obj+1) MOD 1000 ) EQ 0 ) AND (verbose GT 0) THEN $
        print,'.',format='($,a)'

  ENDFOR

  good = where(testind NE -1,ngood)

  IF verbose GT 0 THEN BEGIN 
      print
      ptime,systime(1)-time
  ENDIF 

END 
