PRO sdss_rerunid_syntax
  print,'Syntax: superid=sdss_rerunid(run,rerun OR struct)'
END 

FUNCTION sdss_rerunid, run, rerun
  
  ten=ulong64(10)
  p1 = 0L
  p2 = 5L
  
  np = n_params()
  IF n_params() EQ 1 THEN BEGIN 
      
      IF size(run, /tname) NE 'STRUCT' THEN BEGIN 
          sdss_rerunid_syntax
          return,-1
      ENDIF 

      super = $
        ulong64(run.rerun)*ten^p1  + ulong64(run.run)*ten^p2
      
  ENDIF ELSE IF np EQ 2 THEN BEGIN 
      
      nr   = n_elements(run)
      nrer = n_elements(rerun)
      
      IF total_int([nr,nrer]) NE 2*nr THEN BEGIN 
          print,'All arrays must be same size'
          return, -1
      ENDIF 
      
      super = $
        ulong64(rerun)*ten^p1  + ulong64(run)*ten^p2
      
  ENDIF ELSE BEGIN 
      sdss_rerunid_syntax
      return,-1
  ENDELSE 
  
  return,super
  

END 
