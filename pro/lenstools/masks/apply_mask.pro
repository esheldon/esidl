PRO apply_mask, mask, data, bad, good, edgecut=edgecut, mindata=mindata, maxdata=maxdata

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: apply_mask, mask, data, bad, good, edgecut=edgecut'
      return
  ENDIF 

  ;; the mindata, maxdata, etc keywords are to allow extra or different cuts
  ;; than what is built into the mask

  ndata = n_elements(data)
  good = lindgen(ndata)

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

  IF n_elements(mindata) NE 0 THEN BEGIN 
      bad2 = where(data[good] LT mindata, nbad, complement=good2, ncomp=ngood)
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

  IF n_elements(maxdata) NE 0 THEN BEGIN 
      bad2 = where(data[good] GT maxdata, nbad, complement=good2, ncomp=ngood)
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
