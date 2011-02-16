PRO tang_shear, e1, e2, uncert, $
                x, y, cenx, ceny, $
                rmin, rmax, binsize, $
                etan, erad, etanerr, eraderr, meanr, $
                bin=bin, $
                npair=npair, $
                keepsum=keepsum, $
                etansum=etansum, eradsum=eradsum, $
                etanerrsum=errsum, eraderrsum=eraderrsum, $
                wsum=wsum, npsum=npsum, rsum=rsum

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;   TANG_SHEAR       
;
; PURPOSE:
;   Calculate the tangential and "radial" shear around an input center.
;
; CALLING SEQUENCE:
;      tang_shear, posx, posy, e1, e2, uncert, cenx, ceny, rmin, rmax, 
;               etan, erad, etanerr, eraderr, xy=xy 
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
;	TO DO!!!!: MAKE IT DO RANGE OF BINNINGS TO SAVE TIME (LIKE PHIL'S)
;
; REVISION HISTORY:
;	
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() LT 10 THEN BEGIN 
     print,'-Syntax: tang_shear, e1, e2, uncert,'
     print,'    x, y, cenx, ceny,'
     print,'    rmin, rmax, binsize,'
     print,'    [ etan, erad, etanerr, eraderr, meanr, npair, '
     print,'      bin=bin ]'
     print,''
     print,'Use doc_library,"tang_shear"  for more help.'  
     return
  ENDIF 
  
  IF NOT keyword_set(keepsum) THEN keepsum = 0
  IF NOT keyword_set(bin) THEN bin=0

  xrel = x - cenx
  yrel = y - ceny

  ;; Works faster if sorted by x first
  IF keyword_set(bin) THEN BEGIN 
      w1=where( xrel GT -rmax AND xrel LT rmax, nw1)
      w2=where( yrel[w1] GT -rmax AND yrel[w1] LT rmax, nw2)
      w=w1[w2]
  ENDIF ELSE w=lindgen( n_elements(x) )

  R=sqrt( xrel[w]^2 + yrel[w]^2 ) > double(1.e-6)

  IF NOT keepsum THEN $
       bintshear, e1[w], e2[w], uncert[w], $
                  xrel[w], yrel[w], R, $
                  rmin, rmax, binsize, $
                  etan, erad, etanerr, eraderr, meanr, $
                  npair=npair $
  ELSE $
       bintshear, e1[w], e2[w], uncert[w], $
                  xrel[w], yrel[w], R, $
                  rmin, rmax, binsize, $
                  npair=npair, $
                  keepsum=keepsum, $
                  etansum=etansum, eradsum=eradsum, $
                  etanerrsum=errsum, eraderrsum=eraderrsum, $
                  wsum=wsum, npsum=npsum, rsum=rsum
    

  return 
END 































