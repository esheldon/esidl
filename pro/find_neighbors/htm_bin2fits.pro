PRO htm_bin2fits, frontfile, outfile, remove=remove, output_dist=output_dist

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: htm_bin2fits, frontfile, [outfile, /remove, /output_dist]'
      return
  ENDIF 


  infile1 = frontfile + '_ind1.bin'
  infile2 = frontfile + '_ind2.bin'
  IF keyword_set(output_dist) THEN infile3 = frontfile + '_dist.bin'
  numfile = frontfile + '_num.dat'

  IF n_elements(outfile) EQ 0 THEN outfile = frontfile + '.fit'
  print
  print,'Converting files: '
  print,infile1
  print,infile2
  print,'Num file: ',numfile
  print,'Output file: ',outfile
  print

  ;; read the data
  htm_read_neighbors, frontfile, ind1, ind2, dist, $
                      /binary, output_dist=output_dist

  ;; Remove old file if requested
  IF keyword_set(remove) THEN BEGIN 
      print,'Removing old file: ',infile1
      spawn, 'rm '+infile1
  ENDIF 
  ;; write out to fits file
  mwrfits2, ind1, outfile, hdr, /create, /destroy
 
  ;; Remove old file if requested
  IF keyword_set(remove) THEN BEGIN 
      print,'Removing old file: ',infile2
      spawn, 'rm '+infile2
  ENDIF 
  ;; write out to 2nd extension of fits file
  mwrfits2, ind2, outfile, /destroy

  IF keyword_set(output_dist) THEN BEGIN 
      ;; Remove old file if requested
      IF keyword_set(remove) THEN BEGIN 
          print,'Removing old file: ',infile3
          spawn, 'rm '+infile3
      ENDIF 
      ;; write out to 2nd extension of fits file
      mwrfits2, dist, outfile, /destroy
  ENDIF 

  ;; remove numfile if requested
  IF keyword_set(remove) THEN BEGIN 
      print,'Removing old file: ',numfile
      spawn, 'rm '+numfile
  ENDIF 

  return
END 
