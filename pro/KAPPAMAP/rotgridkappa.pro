PRO rotgridkappa, e1, e2, uncert, x, y, gridx, gridy, cenx, ceny, cos, sin, rmax, slength, kappa, radial, wsum, npair, kerr1, kerr2, rerr1, rerr2, err3, silent=silent

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;       
; PURPOSE:
;	
;
; CALLING SEQUENCE:
;      
;                 
;
; INPUTS: 
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
;       
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; CALLED ROUTINES:
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


  IF N_params() EQ 0 THEN BEGIN 
     print,'-Syntax: rotgridkappa, e1, e2, uncert, x, y, gridx, gridy, cos, sin, rmax, slength, kappa, radial, npair, silent=silent'
     print,''
     print,'Use doc_library,""  for more help.'  
     return
  ENDIF 

  IF NOT keyword_set(silent) THEN silent=0

  nx = n_elements( gridx )
  ny = n_elements( gridy )

  IF n_elements(kappa) EQ 0 THEN BEGIN 
      kappa  = fltarr(nx, ny)
      radial = kappa
      wsum   = kappa
      kerr1  = kappa
      kerr2  = kappa
      rerr1  = kappa
      rerr2  = kappa
      err3   = kappa
      npair  = lonarr(nx, ny)
  ENDIF 

  FOR iy=0, ny-1 DO BEGIN
      FOR ix=0,nx-1 DO BEGIN 
          tx = gridx[ix] - cenx
          ty = gridy[iy] - ceny
          
          ;; rotating the grid is an active transform
          ;; use negative theta formulae

          gx = cos*tx - sin*ty + cenx
          gy = sin*tx + cos*ty + ceny

          w1 = where( y GE gy-rmax AND $
                      y LE gy+rmax, nw1)
          IF nw1 NE 0 THEN BEGIN 
              w2 = where( x[w1] GE gx-rmax AND $
                          x[w1] LE gx+rmax, nw2)
              IF nw2 NE 0 THEN BEGIN 
                  w=w1[w2]
                  ks_kappa, e1[w], e2[w], uncert[w], $
                    x[w], y[w], gx, gy, $
                    rmax, slength, $
                    tkappa, tradial, twsum, tnpair, $
                    tkerr1, tkerr2, trerr1, trerr2, terr3
                  
                  IF n_elements(tkappa) NE 0 THEN BEGIN 
                      kappa[ix, iy]  = kappa[ix, iy]  + tkappa
                      radial[ix, iy] = radial[ix, iy] + tradial
                      wsum[ix, iy]   = wsum[ix, iy]   + twsum
                      npair[ix, iy]  = npair[ix, iy]  + tnpair

                      kerr1[ix, iy]   = kerr1[ix, iy]   + tkerr1
                      kerr2[ix, iy]   = kerr2[ix, iy]   + tkerr2
                      rerr1[ix, iy]   = rerr1[ix, iy]   + trerr1
                      rerr2[ix, iy]   = rerr2[ix, iy]   + trerr2
                      err3[ix, iy]    = err3[ix, iy]    + terr3
                  ENDIF 
              ENDIF 
          ENDIF 
      ENDFOR 
  ENDFOR 

  IF NOT silent THEN print
  return 
END 
