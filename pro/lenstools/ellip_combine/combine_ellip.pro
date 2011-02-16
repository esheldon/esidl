PRO combine_ellip_makecov, elliperr, corr, goodmeas, ngood, covar

  covar = replicate(0.0, ngood, ngood)

  FOR i=0L, ngood-1 DO BEGIN 
      ii = goodmeas[i]

      varii = elliperr[ii]^2
      FOR j=i, ngood-1 DO BEGIN
          jj = goodmeas[j]

          varjj = elliperr[jj]^2
          IF i EQ j THEN BEGIN
              covar[i,j] = varjj
          ENDIF ELSE BEGIN
              covar[i,j] = corr[ii,jj]*sqrt(varii*varjj)
          ENDELSE 
          IF i NE j THEN covar[j,i] = covar[i,j]
      ENDFOR 
  ENDFOR 

  ;; make sure all off-diagonal terms are less than smallest
  ;; diagonal term

  ind = lindgen(ngood*ngood)
  ix = ind MOD ngood
  iy = ind/ngood

  offdiag = where(ix NE iy)
  diag = where(ix EQ iy)

  mindiag = min( covar[ix[diag], iy[diag]] )*0.9
  maxoffdiag = max( covar[ix[offdiag], iy[offdiag]] )

  factor = mindiag/maxoffdiag < 1.0

  covar[ix[offdiag], iy[offdiag]] = covar[ix[offdiag], iy[offdiag]]*factor

END 

PRO combine_ellip, ellip, elliperr, corr, smear, $
                   new_ellip, new_elliperr, $
                   combine_flags, nave, good, verbose=verbose, $
                   err_cut=err_cut, rcut=rcut, minerr=minerr

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: combine_ellip, ellip, elliperr, corr, smear, $'
      print,'         new_ellip, new_elliperr, $'
      print,'         combine_flags, nave, good, verbose=verbose, $'
      print,'         err_cut=err_cut, rcut=rcut, minerr=minerr'
      return
  ENDIF 

;
; combine measurements from different exposures, ignoring the covariance 
; between measurements.  Returns the combined values plus the
; indices of objects that had at least one good ellipticity
; measurement among the inputs.  Also, a bitmask combine_flags is
; returned, indicating which observations (if any) were used to determine
; the best ellipticity.
;
; objects without a good measurement will have -9999. for all
; values, and will have combine_flags == 0
;

  IF n_elements(verbose) EQ 0 THEN verbose=1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; cut on errors and size
  ;; default let most everything through
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(rcut) EQ 0 THEN rcut=9999.
  IF n_elements(err_cut) EQ 0 THEN err_cut = 9999.
  IF n_elements(minerr) EQ 0 THEN minerr=0.0

  time=systime(1)

  siz=size(ellip)

  IF siz[0] NE 2 THEN message,'Must enter arrays of [# measurements, # objects]'

  nmeas = siz[1]
  nobj = siz[2]

  defval = -9999.
  arrval = replicate(defval, nobj)
  new_ellip = arrval
  new_elliperr = arrval
  testind = replicate(-1, nobj)

  ;; error codesfrom the invert function
  errcode=['Successful completion',$
           'Singular Array', $
           'Small pivot element: Accuracy Compromised']

  ;; diagnostic about combination
  ;; 2L^(measurement #) will be added if a measurement is used
  combine_flags = lonarr(nobj)
  nave = bytarr(nobj)

  FOR obj=0L, nobj-1 DO BEGIN 

      delvarx, goodmeas
      tnave = 0L

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; check for good measurements
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      FOR im=0L, nmeas-1 DO BEGIN 

          ;; good measurement?
          IF ( ( ellip[im,obj] GT defval       ) AND $
               ( elliperr[im,obj] LT err_cut ) AND $
               ( smear[im,obj] LT rcut     ) AND $
               ( smear[im,obj] GT 0.0      ) ) THEN BEGIN 

              add_arrval, im, goodmeas
              combine_flags[obj] = combine_flags[obj] + 2L^im
              tnave = tnave+1L
          ENDIF 

      ENDFOR 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; do chisq to get best ellip
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF tnave NE 0 THEN BEGIN 
      
          IF tnave EQ 1 THEN BEGIN 
              ;; one good measurement
              new_ellip[obj] = ellip[goodmeas[0],obj]
              new_elliperr[obj] = elliperr[goodmeas[0],obj]
          ENDIF ELSE BEGIN 
              ;; more than one.  Don't overweight small error measurements
              elliperrsend = elliperr[*,obj] > minerr
              combine_ellip_makecov, elliperrsend, corr, goodmeas, tnave, covar

              combine_ellip_makecov, elliperr[*,obj], corr, goodmeas, tnave, covar

              ;;covar = covar*1.e3

              chisq_avg, ellip[goodmeas,obj], covar, beste, besterr, status=stat
              ;;besterr = besterr/sqrt(1.e3)
              IF stat NE 0 THEN BEGIN 
                  IF verbose GT 0 THEN BEGIN 
                      print
                      print,'# '+ntostr(obj)+$
                            ': Error inverting covariance matrix: '+errcode[stat]
                  ENDIF 
                  ;; pick best when can't invert
                  wt = where(elliperr[goodmeas] EQ min(elliperr[goodmeas]), nwt)
                  new_ellip[obj] = ellip[goodmeas[wt]]
                  new_elliperr[obj] = elliperr[goodmeas[wt]]
                  goodmeas = goodmeas[wt[0]]
                  tnave = 1
              ENDIF ELSE BEGIN 
                  new_ellip[obj] = beste
                  new_elliperr[obj] = besterr
              ENDELSE 
              
          ENDELSE 

          testind[obj] = 1
          nave[obj] = tnave

          IF verbose GT 2 THEN BEGIN 
              
              ;; For comparison:
              ;; do a weighted mean of the errors
              ;; don't overweight
              ;; errors less than 0.1
              print,'Nave = ',tnave
              sigma = elliperr[goodmeas,obj] > minerr
              wmom, ellip[goodmeas, obj], sigma, wmom_e, wmom_var, wmom_err

              xvals = [(5-tnave)+lindgen(tnave)]
              
              yrange = prange(ellip[goodmeas,obj], elliperr[goodmeas,obj])
              plot,[0],/nodata,xrange=[0,5],yrange=yrange, ystyle=1
              oploterror,xvals,$
                         [ellip[goodmeas,obj]],$
                         [elliperr[goodmeas,obj]],psym=4
              oploterror,[2.5],[new_ellip[obj]],[new_elliperr[obj]],psym=5,color=!red
              oploterror,[2.25],[wmom_e],[wmom_err],psym=4,color=!blue
              
              print,ellip[goodmeas,obj]
              print,new_ellip[obj],new_elliperr[obj]
              print
              key=get_kbrd(1)
              
          ENDIF 
                  
      ENDIF ELSE BEGIN 
          IF verbose GT 0 THEN BEGIN 
              print
              print,'# '+ntostr(obj)+$
                    ': No good measurements'
          ENDIF 
      ENDELSE  
      IF ( (obj+1) MOD 1000 ) EQ 0 AND verbose GT 1 THEN print,'.',format='($,a)'

  ENDFOR 

  IF verbose GT 0 THEN BEGIN 
      print
      print,systime(1)-time
  ENDIF 
  good = where(testind NE -1)

END 
