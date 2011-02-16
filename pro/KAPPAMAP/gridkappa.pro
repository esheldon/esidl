PRO gridkappa, e1, e2, uncert, x, y, gridx, gridy, rmax, slength, kappa, radial, wsum, npair, kerr1, kerr2, rerr1, rerr2, err3, silent=silent, fgal=fgal, dens=dens

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
     print,'-Syntax: gridkappa, e1, e2, uncert, x, y, gridx, gridy, rmax, slength, kappa, radial, wsum, npair'
     print,''
     print,'Use doc_library,""  for more help.'  
     return
  ENDIF 

  IF NOT keyword_set(silent) THEN silent=0

  nx = n_elements( gridx )
  ny = n_elements( gridy )
  nfgal = n_elements(fgal)
  IF nfgal NE 0 THEN dodens = 1 ELSE dodens=0

  IF n_elements(kappa) EQ 0 THEN BEGIN 
      kappa  = fltarr(nx, ny)
      radial = kappa
      wsum   = kappa
      kerr1  = kappa
      kerr2  = kappa
      rerr1  = kappa
      rerr2  = kappa
      err3   = kappa
      dens   = kappa
      npair  = lonarr(nx, ny)
  ENDIF 

  FOR iy=0, ny-1 DO BEGIN
      IF NOT silent THEN print, format='($,A)', ntostr(iy+1)+'.'
      ;IF iy MOD 20 EQ 0 THEN print, format='($,A)', ntostr(iy+1)+'.'
      gy = gridy[iy]
      w1 = where( y GE gy-rmax AND $
                  y LE gy+rmax, nw1)
      IF dodens THEN wdens1 = where(fgal.ra GE gy-rmax AND $
                                    fgal.ra LE gy+rmax, ndens1)

      IF nw1 NE 0 THEN BEGIN 
          FOR ix=0, nx-1 DO BEGIN 
              gx = gridx[ix]
              w2 = where( x[w1] GE gx-rmax AND $
                          x[w1] LE gx+rmax, nw2)

              IF nw2 NE 0 THEN BEGIN 
                  w=w1[w2]
                  IF dodens THEN BEGIN 
                      IF ndens1 NE 0 THEN BEGIN
                          wdens = where(fgal[wdens1].dec GE gx-rmax AND  $
                                        fgal[wdens1].dec LE gx+rmax, ndens2)
                          IF ndens2 NE 0 THEN sendf = fgal[wdens1[wdens]]
                      ENDIF 
                  ENDIF 

                  ks_kappa, e1[w], e2[w], uncert[w], $
                    x[w], y[w], gx, gy, $
                    rmax, slength, $
                    tkappa, tradial, twsum, tnpair, $
                    tkerr1, tkerr2, trerr1, trerr2, terr3, $
                    fgal=sendf, dens=tmpdens
                  
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
                      IF n_elements(tmpdens) NE 0 THEN BEGIN
                          dens[ix, iy] = dens[ix,iy] + tmpdens
                      ENDIF 
                  ENDIF 
              ENDIF 
          ENDFOR 
      ENDIF 
  ENDFOR 

  IF NOT silent THEN print
  return 
END 
