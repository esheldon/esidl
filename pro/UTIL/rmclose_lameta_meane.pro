PRO rmclose_lameta_combine, e1, e2, e1e1err, e1e2err, e2e2err, $
                            rsmear, seeing, $
                            new_e1, new_e2, $
                            new_e1e1err, new_e1e2err, new_e2e2err,$
                            new_rsmear, new_seeing, $
                            new_cflags, new_nave, new_good

  combine_ellip_cove1e2, e1, e2, $
                         e1e1err, e1e2err, e2e2err, $
                         rsmear, seeing, $
                         new_e1, new_e2, $
                         new_e1e1err, new_e1e2err, $
                         new_e2e2err,$
                         new_rsmear, new_seeing, $
                         new_cflags, new_nave, new_good,$
                         /oneobj


END 

PRO rmclose_lameta_meane, struct, keep, out, $
                          e1, e2, $
                          e1e1err, e1e2err, e2e2err, $
                          rsmear, seeing, $
                          photoz_z, photoz_zerr, photoz_zwarning, $
                          combine_flags, nave, $
                          tol=tol, nkeep=nkeep

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: 
;    rmclose_lameta
;       
; PURPOSE: 
;    Remove duplicate objects from struct; that is, objects that have
;    the same position within some tolerance.  Must be sorted by
;    lambda.
;	
;
; CALLING SEQUENCE: 
;    rmclose_radec, struct, keep, out           
;
; INPUTS: struct   Structure that must contain ra and dec as tags.
;
; OUTPUTS: keep.  The indices of non-degenerate entries.
;          out.      The stuff we threw out.
; 
;
; REVISION HISTORY:
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() LT 12 THEN BEGIN 
     print,'-Syntax: rmclose_lameta_meane, struct, keep, out, $'
     print,'              e1, e2, $'
     print,'              e1e1err, e1e2err, e2e2err, $'
     print,'              rsmear, seeing, $'
     print,'              photoz_z, photoz_zerr, photoz_zwarning, $'
     print,'              combine_flags, nave, $'
     print,'              tol=tol, nkeep=nkeep'
     print
     print,'WARNING: struct must be sorted by lambda'
     print,'Use doc_library,"rmclose_radec"  for more help.'  
     return
  ENDIF 

  ns = n_elements(struct)
  print
  print,'Averaging duplicates'
  print,'Beginning with '+ntostr(ns)+' galaxies'
  print

  IF NOT tag_exist(struct, 'clambda') OR NOT tag_exist(struct, 'ceta') THEN BEGIN 
      message,'struct must contain clambda,ceta tags'
  ENDIF 

  d2r=!dpi/180.
                                ;1 arcsecond default tolerance
  IF n_elements(tol) EQ 0 THEN tol1 = double(1.)/3600. ELSE tol1=tol
  tol2=tol1*d2r

  cont = 1

  indices = lindgen(ns)
  available = replicate(1b, ns)
  keep      = replicate(1b, ns)

  e1 = struct.e1
  e2 = struct.e2

  e1e1err = struct.e1e1err
  e1e2err = struct.e1e2err
  e2e2err = struct.e2e2err

  rsmear = struct.r
  seeing = struct.seeing

  photoz_z = struct.photoz_z
  photoz_zerr = struct.photoz_zerr
  photoz_zwarning = struct.photoz_zwarning

  combine_flags = lonarr(ns)
  nave = replicate(1b, ns)

  minav_i = 0L
  minav_nn = 1L
  FOR i=0L, ns-2L DO BEGIN 

      IF available[i] THEN BEGIN 

          ;; check the next object(s)
          nn = i+1
          nonext = 0
          WHILE (NOT available[nn]) AND (NOT nonext) DO BEGIN
 
              ;; was that the last one?
              IF nn EQ (ns-1L) THEN BEGIN 
                  ;; we can't go on, just use i
                  ;;print,'last one'
                  nonext = 1
              ENDIF ELSE BEGIN 
                  nn=nn+1
              ENDELSE 

          ENDWHILE 

          IF nonext THEN BEGIN 
              ;; no next object: keep shape for i, remove i from available
              ;; indices
              available[i] = 0b
          ENDIF ELSE BEGIN 
          
              ;; Sorted by lambda: get all those with delta lambda
              ;; less than tolerance
              
              cont2 = 1
              delvarx, wcloselam
              ncloselam = 0
              WHILE cont2 DO BEGIN 
                  
                  dlam = abs( struct[nn].clambda - struct[i].clambda )
                  IF dlam LT tol1 THEN BEGIN
                      add_arrval, nn, wcloselam
                      ncloselam = ncloselam + 1
                      nn = nn+1
                      IF (nn GE ns) THEN BEGIN
                          ;;print,'last one'
                          cont2 = 0
                      ENDIF 
                  ENDIF ELSE BEGIN
                      ;; they are sorted: once we reach one outside
                      ;; of reach, we can stop looking
                      cont2 = 0
                  ENDELSE 
              ENDWHILE 
              
              ;; we found some close in lambda
              IF ncloselam NE 0 THEN BEGIN 
                  
;                  gcirc, 0, $
;                    struct[i].ceta*d2r, struct[i].clambda*d2r, $
;                    struct[wcloselam].ceta*d2r, struct[wcloselam].clambda*d2r,$
;                    dist

                  mygcirc_survey, $
                    struct[i].clambda, struct[i].ceta, $
                    struct[wcloselam].clambda, struct[wcloselam].ceta, $
                    dist, /radians_out

                  wclose = where(dist LT tol2, nclose)

                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; We actually got matches
                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                  IF nclose NE 0 THEN BEGIN 
                      ;; average the shapes, remove i and wclose from the
                      ;; available indices
                      wclose = wcloselam[wclose]
                      available[i] = 0b
                      available[wclose] = 0b
                      keep[wclose] = 0b

                      ;; average the shapes

                      sende1 = [struct[i].e1, struct[wclose].e1]
                      sende2 = [struct[i].e2, struct[wclose].e2]
                      sende1e1err = [struct[i].e1e1err, struct[wclose].e1e1err]
                      sende1e2err = [struct[i].e1e2err, struct[wclose].e1e2err]
                      sende2e2err = [struct[i].e2e2err, struct[wclose].e2e2err]

                      send_rsmear = [struct[i].r, struct[wclose].r]
                      send_seeing = [struct[i].seeing, struct[wclose].seeing ]
                      rmclose_lameta_combine, sende1, sende2, $
                        sende1e1err, $
                        sende1e2err, $
                        sende2e2err, $
                        send_rsmear, $
                        send_seeing, $
                        new_e1, new_e2, $
                        new_e1e1err, $
                        new_e1e2err, $
                        new_e2e2err,$
                        new_rsmear, $
                        new_seeing, $
                        new_cflags, new_nave, new_good

                      IF new_good[0] EQ -1 THEN message,'No good measurements?'
                      e1[i] = new_e1[0]
                      e2[i] = new_e2[0]
                      e1e1err[i] = new_e1e1err[0]
                      e1e2err[i] = new_e1e2err[0]
                      e2e2err[i] = new_e2e2err[0]
                      rsmear[i] = new_rsmear[0]
                      seeing[i] = new_seeing[0]

                      combine_flags[i] = new_cflags[0]
                      nave[i] = new_nave[0]

                      IF (nclose+1) NE nave[i] THEN print,'Different: ',i,(nclose+1),nave[i]

                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;
                      ;; average the photoz's
                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;

                      tphz = [struct[i].photoz_z, $
                              struct[wclose].photoz_z]
                      tphzerr = [struct[i].photoz_zerr, $
                                 struct[wclose].photoz_zerr]
                      tphzwarn = [struct[i].photoz_zwarning, $
                                  struct[wclose].photoz_zwarning]
                      
                      wphz = where(tphzwarn EQ 0, ngphz)
                      IF ngphz NE 0 THEN BEGIN 
                          senderr = tphzerr[wphz] > 0.015
                          wmom, tphz[wphz], senderr, $
                            phzmean, phzsig, phzerr
                          photoz_z[i] = phzmean
                          photoz_zerr[i] = phzerr
                          photoz_zwarning[i] = 0
                      ENDIF 

                  ENDIF ELSE BEGIN 
                      ;; No matches: keep shape for i, remove i from available
                      ;; indices
                      available[i] = 0b
                  ENDELSE 
                  
              ENDIF ELSE BEGIN 
                  ;; No matches: keep shape for i, remove i from available
                  ;; indices
                  available[i] = 0b
              ENDELSE 
          ENDELSE 

      ENDIF 

  ENDFOR 

  ;; what to keep
  w = where(keep EQ 1b, nkeep)
  IF nkeep NE 0 THEN BEGIN
      keep = indices[w] 
      e1 = e1[w]
      e2 = e2[w]
      
      e1e1err = e1e1err[w]
      e1e2err = e1e2err[w]
      e2e2err = e2e2err[w]

      rsmear = rsmear[w]
      seeing = seeing[w]

      combine_flags = combine_flags[w]
      nave = nave[w]

      photoz_z = photoz_z[w]
      photoz_zerr = photoz_zerr[w]
      photoz_zwarning = photoz_zwarning[w]

  ENDIF ELSE BEGIN
      delvarx, e1, e2, e1e1err, e1e2err, e2e2err, rsmear, seeing, $
        combine_flags, nave
      keep = -1
  ENDELSE 

  return
END
