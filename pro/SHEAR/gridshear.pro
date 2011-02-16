PRO gridshear, gridx, gridy, scat, rmin, rmax, binsize, $
               imtan, imrad, imtanerr, imraderr, imtans2n

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
     print,'-Syntax: gridshear, gridx, gridy, scat, rmin, rmax, binsize, imtan, imrad'
     print,''
     print,'Use doc_library,""  for more help.'  
     return
  ENDIF 

  nx = n_elements( gridx )
  ny = n_elements( gridy )

  imtan = fltarr( nx, ny )
  imrad = imtan
  imtanerr = imtan
  imraderr = imtan
  imtans2n = imtan

  nbin = long( (rmax-rmin)/binsize )
  print,'Number of bins: ',nbin

  FOR iy=0, ny-1 DO BEGIN
      print, format='($,A)', ntostr(iy+1)+'.'
      ;IF iy MOD 20 EQ 0 THEN print, format='($,A)', ntostr(iy+1)+'.'
      w1 = where( scat.ra GE gridy[iy]-rmax AND $
                  scat.ra LE gridy[iy]+rmax)
      FOR ix=0, nx-1 DO BEGIN 
          w2 = where( scat[w1].dec GE gridx[ix]-rmax AND $
                      scat[w1].dec LE gridx[ix]+rmax)
          w=w1[w2]
          tan_shear, scat[w].e1, scat[w].e2, scat[w].uncert, $
                      scat[w].dec, scat[w].ra, gridx[ix], gridy[iy], $
                      rmin, rmax, binsize, $
                      etan, erad, etanerr, eraderr, meanr
          
          g=where(etan NE 1.e10, ng)
          IF ng NE 0 THEN BEGIN
              imtan[ix,iy] = total( etan[g] )
              imrad[ix,iy] = total( erad[g] )
              imtanerr[ix, iy] = sqrt( total( etanerr[g]^2 ) )
              imraderr[ix, iy] = sqrt( total( eraderr[g]^2 ) )
              imtans2n[ix, iy] = sqrt( total( etan[g]^2/etanerr[g]^2 ) )
          ENDIF
      ENDFOR
  ENDFOR 
  print
  return 
END 
