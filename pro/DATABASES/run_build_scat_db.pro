PRO run_build_scat_db, chunk

  CASE chunk OF
      1: stripes = [9,10,11,12,13,14,15]
      2: stripes = [27,28,29,30,31,32,33,34,35,36,37]
      3: stripes = [42,43,76,82,86]
      ELSE: message,'Chunk must be in [1,3]'
  ENDCASE 

  tm=systime(1)
  nstripe = n_elements(stripes)
  FOR i=0L, nstripe-1 DO BEGIN 

      stripe = stripes[i]

      IF i EQ 0 THEN create = 1 ELSE create = 0
      IF i EQ nstripe-1 THEN noindex=0 ELSE noindex=1

      build_scat_db, stripe, create=create, noindex=noindex, /noprompt

  ENDFOR 
  print,'For all stripes'
  ptime,systime(1)-tm
END 
