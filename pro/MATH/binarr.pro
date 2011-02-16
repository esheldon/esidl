PRO binarr, arr, binsize, min, max, rev_ind, nbin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: 
;    BINARR
;       
; PURPOSE: 
;    Bin an array and return indices array
;	
;
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() LT 4 THEN BEGIN 
     print,'-Syntax: binarr, arr, binsize, min, max, rev_ind, nbin'
     print,''
     print,'Use doc_library,"binarr"  for more help.'  
     return
  ENDIF 

  hist = histogram(arr, binsize=binsize, min=min, max=max, rever=rev_ind)
  nbin=n_elements( hist )
  
END































