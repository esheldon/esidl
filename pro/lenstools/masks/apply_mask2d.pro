PRO apply_mask2d, mask, datax, datay, bad, good, edgecut=edgecut, maxy=maxy, miny=miny, maxx=maxx, minx=minx

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: apply_mask2d, mask, datax, datay, bad, good, edgecut=edgecut, maxy=maxy, miny=miny, minx=minx, maxx=maxx'
      return
  ENDIF 

  ;; the minx, maxx, etc keywords are to allow extra or different cuts
  ;; than what is built into the mask

  ndata = n_elements(datax)
  ndatay = n_elements(datay)
  IF ndata NE ndatay THEN message,'datax and datay must be same size'

  all = lindgen(ndata)
  good = all

  IF n_elements(edgecut) EQ 0 THEN DATA_EDGECUT = '0.0' $
  ELSE BEGIN
      n_edge=n_elements(edgecut)
      IF n_edge EQ 1 THEN BEGIN 
          DATA_EDGECUT = string(edgecut,format='(F0,:,X)')
      ENDIF ELSE BEGIN 
          IF n_edge NE ndata THEN message,'# in edgecut must be equal to 1 or # in data'
          DATA_EDGECUT = strarr(n_edge)
          FOR i=0L, n_edge-1 DO BEGIN 
              DATA_EDGECUT[i] = string(edgecut[i],format='(F0,:,X)')
          ENDFOR 
          
      ENDELSE 
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; apply the basic mask
  ;;;;;;;;;;;;;;;;;;;;;;;;;;

  command = 'bad = where('+mask+', nbad, complement=good, ncomp=ngood)'
  IF NOT execute(command) THEN message,'mask failed'

  IF ngood EQ 0 THEN return
    
  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Further cuts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(minx) NE 0 THEN BEGIN 
      bad2 = where(datax[good] LT minx, nbad, complement=good2, ncomp=ngood)
      IF ngood EQ 0 THEN BEGIN 
          bad = all
          good=-1
          return
      ENDIF 
      IF nbad NE 0 THEN BEGIN
          bad  = good[bad2]
          good = good[good2]
      ENDIF ;; else bad/good are the same as before
  ENDIF 

  IF n_elements(maxx) NE 0 THEN BEGIN 
      bad2 = where(datax[good] GT maxx, nbad, complement=good2, ncomp=ngood)
      IF ngood EQ 0 THEN BEGIN 
          bad = all
          good=-1
          return
      ENDIF 
      IF nbad NE 0 THEN BEGIN
          bad  = good[bad2]
          good = good[good2]
      ENDIF ;; else bad/good are the same as before
  ENDIF 

  IF n_elements(miny) NE 0 THEN BEGIN 
      bad2 = where(datay[good] LT miny, nbad, complement=good2, ncomp=ngood)
      IF ngood EQ 0 THEN BEGIN 
          bad = all
          good=-1
          return
      ENDIF 
      IF nbad NE 0 THEN BEGIN
          bad  = good[bad2]
          good = good[good2]
      ENDIF ;; else bad/good are the same as before
  ENDIF 

  IF n_elements(mayy) NE 0 THEN BEGIN 
      bad2 = where(datay[good] GT mayy, nbad, complement=good2, ncomp=ngood)
      IF ngood EQ 0 THEN BEGIN 
          bad = all
          good=-1
          return
      ENDIF 
      IF nbad NE 0 THEN BEGIN
          bad  = good[bad2]
          good = good[good2]
      ENDIF ;; else bad/good are the same as before
  ENDIF 



END 
