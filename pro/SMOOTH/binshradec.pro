PRO binshradec, struct, bsize, ime1, ime2, ime1err, ime2err

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;  BINSHRADEC
;       
; PURPOSE:
;  bin e1 and e2 in ra-dec space.
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
     print,'-Syntax: binshradec, struct, binsize, ime1, ime2, ime1err, ime2err'
     print,' binsize in arcseconds.  Ra and dec of struct in degrees'
     print,'Use doc_library,""  for more help.'  
     return
  ENDIF 
  
  vint = .32^2
  binsize = bsize/3600.

  minra  = min(struct.ra)
  maxra  = max(struct.ra)
  mindec = min(struct.dec)
  maxdec = max(struct.dec)

  radiff  = maxra-minra
  decdiff = maxdec-mindec

  nra  = long(radiff/binsize) > 1
  ndec = long(decdiff/binsize) > 1

  mra  = minra + nra*binsize
  mdec = mindec + ndec*binsize

  ime1    = dblarr(ndec, nra)
  ime2    = ime1
  ime1err = ime1
  ime2err = ime2

  binarr, struct.ra, binsize, minra, mra, rev_ind1
  FOR ira=0L, nra-1 DO BEGIN

      ;; make sure we have something to measure
      IF rev_ind1[ira] NE rev_ind1[ira+1] THEN BEGIN 
          w1 = rev_ind1( rev_ind1[ira]:rev_ind1[ira+1]-1 ) 
          binarr, struct[w1].dec, binsize, mindec, mdec, rev_ind2
          FOR idec=0L, ndec-1 DO BEGIN

              ;; make sure we have something to measure
              IF rev_ind2[idec] NE rev_ind2[idec+1] THEN BEGIN 

                  w2 = rev_ind2( rev_ind2[idec]:rev_ind2[idec+1]-1 )
                  w=w1[w2]
                  
                  err = sqrt(vint + struct[w].uncert^2)
                  wmom,struct[w].e1,err, wmean, wsig, werr
                  ime1[idec, ira] = wmean
                  ime1err[idec, ira] = werr

                  wmom,struct[w].e2,err, wmean, wsig, werr
                  ime2[idec, ira] = wmean
                  ime2err[idec, ira] = werr

              ENDIF ELSE BEGIN
                  print,'Bin ira='+ntostr(ira)+$
                        ' idec='+ntostr(idec)+' is empty'
              ENDELSE 
              
          ENDFOR 
      ENDIF ELSE BEGIN 
          print,'Bin: ira='+ntostr(ira)+' is empty'   
      ENDELSE 
  ENDFOR 
  return 
END 
