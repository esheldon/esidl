
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    
;       
; PURPOSE:
;   Combine measurements from different exposures, ignoring the covariance 
;   between measurements.  Returns the combined values plus the
;   indices of objects that had at least one good ellipticity
;   measurement among the inputs.  Also, a bitmask combine_flags is
;   returned, indicating which observations (if any) were used to determine
;   the best ellipticity.
;
;  objects without a good measurement will have ellipticities of -9999 and
;  errors of 9999.  Also, combine_flag==0
;
;
; CALLING SEQUENCE:
;    combine_ellip_cove1e2, e1, e2, e1e1err, e1e2err, e2e2err, smear, $
;                          new_e1, new_e2, $
;                          new_e1e1err, new_e1e2err, new_e2e2err,$
;                          combine_flags, nave, good, verbose=verbose, $
;                          err_cut=err_cut, rcut=rcut, $
;                          besterror=besterror, $
;                          wmomerror=wmomerror, $
;                          minerr=minerr
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
;    
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO combine_ellip_cove1e2, e1, e2, e1e1err, e1e2err, e2e2err, smear, seeing, $
                           new_e1, new_e2, $
                           new_e1e1err, new_e1e2err, new_e2e2err, $
                           new_smear, new_seeing, $
                           combine_flags, nave, good, verbose=verbose, $
                           etot_cut=etot_cut, $
                           err_cut=err_cut, rcut=rcut, $
                           besterror=besterror, $
                           wmomerror=wmomerror, $
                           minerr=minerr, oneobj=oneobj

;  time=systime(1)
  IF N_params() LT 14 THEN BEGIN 
     print,'-Syntax: combine_ellip_cove1e2, e1, e2, e1e1err, e1e2err, e2e2err, smear, seeing, $'
     print,'               new_e1, new_e2, $'
     print,'               new_e1e1err, new_e1e2err, new_e2e2err, new_smear, new_seeing, $'
     print,'               combine_flags, nave, good, verbose=verbose, $'
     print,'               etot_cut=, err_cut=, rcut=, $'
     print,'               besterror=, $'
     print,'               wmomerror=, $'
     print,'               minerr=, oneobj='
     print
     print,'Use doc_library,"combine_ellip_cove1e2"  for more help.'  
     return
  ENDIF 

  IF n_elements(verbose) EQ 0 THEN verbose=1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Sanity cuts on errors and size
  ;; default let most everything through: should probably make cuts 
  ;; (which may be complicated) outside this procedure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(etot_cut) EQ 0 THEN etot_cut = 4.0
  IF n_elements(rcut) EQ 0 THEN rcut=9999.
  IF n_elements(err_cut) EQ 0 THEN err_cut = 9999.
  IF n_elements(minerr) EQ 0 THEN minerr = 0.0
  minerr2 = minerr^2

  IF keyword_set(oneobj) THEN BEGIN 
      nobj = 1
      nmeas = n_elements(e1)
  ENDIF ELSE BEGIN 
      siz=size(e1)

      IF siz[0] NE 2 THEN message,$
        'Must enter arrays of [# measurements, # objects]'
      
      nmeas = siz[1]
      nobj = siz[2]
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; These are to cut out objects which have the default value for when
  ;; the procedure fails.  Also, our defaults for output e1,e2 and
  ;; errors are drawn from these
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  defval = -9999.
  defval2 = 9999.
  arrval = replicate(defval, nobj)
  arrval2 = replicate(defval2, nobj)

  ;; output arrays
  new_e1 = arrval
  new_e2 = arrval
  new_e1e1err = arrval2
  new_e1e2err = arrval2
  new_e2e2err = arrval2
  new_smear = arrval
  new_seeing = arrval2
  testind = replicate(-1, nobj)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; diagnostic about combination.  Uses an IDL_LONG right now, which
  ;; limits us to 30 observations
  ;; 2L^(measurement # + 2) will be added if a measurement is used
  ;; Also, a flag for how the combination was made is set if
  ;; we couldn't use the standard matrix approach
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  COMBFLAG_INVERSION = 2^0
  COMBFLAG_DETERMINANT = 2^1

  combine_flags = lonarr(nobj)
  nave = bytarr(nobj)

  ;; error codesfrom the invert function
  errcode=['Successful completion',$
           'Singular Array', $
           'Small pivot element: Accuracy Compromised']


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; for use in the computations
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  bigsig = fltarr(2,2)
  tvar = fltarr(2,2)
  beste = fltarr(2)
  bestsmear = fltarr(2)
  bestseeing = fltarr(2)
  te = fltarr(2)
  tsmear = fltarr(2)
  tseeing = fltarr(2)

  FOR obj=0L, nobj-1 DO BEGIN 

      delvarx, goodmeas
      bigsig[*,*] = 0.0
      beste[*] = 0.0
      bestsmear[*] = 0.0
      bestseeing[*] = 0.0
      tnave = 0L
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; loop over the measurements for this object
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      FOR im=0L, nmeas-1 DO BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Good measurement?
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          etot = sqrt( e1[im, obj]^2 + e2[im, obj]^2 )
          IF ( ( e1[im,obj] GT defval       ) AND $
               ( etot LT etot_cut           ) AND $
               ( e1e1err[im,obj] LT err_cut ) AND $
               ( e2e2err[im,obj] LT err_cut ) AND $
               ( smear[im,obj] LT rcut      ) AND $
               ( smear[im,obj] GT 0.0       ) ) THEN BEGIN 

              add_arrval, im, goodmeas

              combine_flags[obj] = combine_flags[obj] + 2^(im+2)

              te[0] = e1[im,obj]
              te[1] = e2[im,obj]
              tsmear[0] = smear[im,obj]
              tsmear[1] = smear[im,obj]

              tseeing[0] = seeing[im, obj]
              tseeing[1] = seeing[im, obj]

              tvar[0,0] = e1e1err[im,obj]^2 > minerr2
              tvar[0,1] = sign(e1e2err[im,obj])*e1e2err[im,obj]^2
              tvar[1,0] = tvar[0,1]
              tvar[1,1] = e2e2err[im,obj]^2 > minerr2
              
              tbigsig = invert(tvar, /double)
              bigsig = bigsig + tbigsig
              
              beste = beste + reform( tbigsig##te )
              bestsmear = bestsmear + reform( tbigsig##tsmear )
              bestseeing = bestseeing + reform( tbigsig##tseeing )

              tnave = tnave+1L

          ENDIF 

      ENDFOR 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Any good measurements?
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      usebest = 0               ; only use best if there are good measurements
                                ; but we can't invert covariance matrix
      IF combine_flags[obj] NE 0 THEN BEGIN 
      
          cove = invert(bigsig, stat, /double)
          detcov = determ(cove)

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; check inversion and that determinant is positive definite
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF stat NE 0 THEN BEGIN 

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Difficulty inverting matrix
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              IF verbose GT 0 THEN BEGIN 
                  print
                  print,'# '+ntostr(obj)+$
                        ': Error inverting covariance matrix: '+errcode[stat]
                  print,'Choosing best seeing'
              ENDIF 
              combine_flags[obj] = 0
              combine_flags[obj] = combine_flags[obj] + COMBFLAG_INVERSION

              usebest = 1
          ENDIF ELSE IF detcov LT 0.0 THEN BEGIN 

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; determinant is not positive definite
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              IF verbose GT 0 THEN BEGIN 
                  print
                  print,'# '+ntostr(obj)+$
                        ': New covariance matrix has det < 0.0'
                  print,'Choosing best seeing'
              ENDIF 
              combine_flags[obj] = 0
              combine_flags[obj] = combine_flags[obj] + COMBFLAG_DETERMINANT
              usebest = 1
          ENDIF ELSE BEGIN 

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Everything looks good: proceed
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              testind[obj] = 1
              nave[obj] = tnave

              beste = cove##beste
              
              new_e1[obj] = beste[0]
              new_e2[obj] = beste[1]

              bestsmear = cove##bestsmear
              bestseeing = cove##bestseeing

              ;; average the two
              new_smear[obj] = total( bigsig##bestsmear )/total(bigsig)
              new_seeing[obj] = total( bigsig##bestseeing )/total(bigsig)

              ;;new_smear[obj] = (bigsig[0,0]*bestsmear[0] + bigsig[1,1]*bestsmear[1])/(bigsig[0,0] + bigsig[1,1])

              ;;new_smear[obj] = (bestsmear[0] + bestsmear[1])/2.

              IF keyword_set(besterror) THEN BEGIN 

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; just use the best error matrix
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  tmperr = ( e1e1err[goodmeas,obj]^2 + $
                             e2e2err[goodmeas,obj]^2 )

                  junk = min(tmperr, wtmp)
                  wtmp = goodmeas[wtmp]
                  
                  new_e1e1err[obj] = e1e1err[wtmp, obj]
                  new_e1e2err[obj] = e1e2err[wtmp, obj]
                  new_e2e2err[obj] = e2e2err[wtmp, obj]
                  
              ENDIF ELSE IF keyword_set(wmomerror) THEN BEGIN 

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; do a weighted mean of the errors
                  ;; don't overweight
                  ;; errors less than 0.1
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  
                  sigma = sqrt( 0.01 + $
                                e1e1err[goodmeas,obj]^2 + $
                                e2e2err[goodmeas,obj]^2 )
                  wmom, e1e1err[goodmeas, obj], sigma, $
                        e1e
                  wmom, e1e2err[goodmeas, obj], sigma, $
                        e1e2e
                  wmom, e2e2err[goodmeas, obj], sigma, $
                            e2e
                  
                  new_e1e1err[obj] = e1e
                  new_e1e2err[obj] = e1e2e
                  new_e2e2err[obj] = e2e
                  
              ENDIF ELSE BEGIN 

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; The Default:
                  ;; use the new covariance matrix
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  new_e1e1err[obj] = sqrt( cove[0,0] )
                  new_e1e2err[obj] = sqrt( abs(cove[0,1]) )*sign(cove[0,1])
                  new_e2e2err[obj] = sqrt( cove[1,1] )
              ENDELSE 

              ;; do some plots to check reasonableness
              IF verbose GT 2 THEN BEGIN 

                  !p.multi=[0,0,2]
                  
                  ngood = n_elements(goodmeas)
                  xvals = [(5-ngood)+lindgen(ngood)]
                  
                  yrange = prange(e1[goodmeas,obj], e1e1err[goodmeas,obj])
                  plot,[0],/nodata,xrange=[0,5],yrange=yrange
                  oploterror,xvals,$
                             [e1[goodmeas,obj]],$
                             [e1e1err[goodmeas,obj]],psym=4
                  oploterror,[2.5],[new_e1[obj]],[new_e1e1err[obj]],psym=5,color=!red
                  
                  yrange = prange(e2[goodmeas,obj], e2e2err[goodmeas,obj])
                  plot,[0],/nodata,xrange=[0,5],yrange=yrange
                  oploterror,xvals,$
                             [e2[goodmeas,obj]],$
                             [e2e2err[goodmeas,obj]],psym=4
                  oploterror,[2.5],[new_e2[obj]],[new_e2e2err[obj]],psym=5,color=!red
                  !p.multi=0
                      
                  print,'Ngood: ',ngood
                  print,'e1: ',e1[goodmeas,obj]
                  print,'mean e1: ',new_e1[obj],!plusminus,new_e1e1err[obj]
                  print,'e2: ',e2[goodmeas,obj]
                  print,'mean e2: ',new_e2[obj],!plusminus,new_e2e2err[obj]
                  print
                  key=get_kbrd(1)
                  
              ENDIF  
                  
          ENDELSE  
      ENDIF ELSE BEGIN 
          IF verbose GT 0 THEN BEGIN 
              print
              print,'# '+ntostr(obj)+$
                    ': No good measurements'
          ENDIF 
      ENDELSE 
      
      IF usebest THEN BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; We had good measurements, but couldn't invert the covariance
          ;; matrix: pick the best seeing.  
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;          tmperr = ( e1e1err[goodmeas,obj]^2 + $
;                     e2e2err[goodmeas,obj]^2 )
;          junk = min(tmperr, wtmp)

          junk = min(seeing[goodmeas, obj], wtmp)

          wtmp = goodmeas[wtmp]

          new_e1[obj] = e1[wtmp,obj]
          new_e2[obj] = e2[wtmp,obj]
          new_e1e1err[obj] = e1e1err[wtmp,obj]
          new_e1e2err[obj] = e1e2err[wtmp,obj]
          new_e2e2err[obj] = e2e2err[wtmp,obj]

          new_smear[obj] = smear[wtmp,obj]
          new_seeing[obj] = seeing[wtmp, obj]

          nave[obj] = 1
          testind[obj] = 1
          combine_flags[obj] = 2^(wtmp+2)
          
      ENDIF
      IF ( (obj+1) MOD 1000 ) EQ 0 AND verbose GT 1 THEN $
        print,'.',format='($,a)'

  ENDFOR 

;  IF verbose GT 0 THEN BEGIN 
;      print
;      print,systime(1)-time
;  ENDIF 
  good = where(testind NE -1)

END 
