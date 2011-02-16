PRO htm_read_neighbors, file, ind1, ind2, distance, hdr=hdr, binary=binary, output_dist=output_dist

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: htm_read_neighbors, file, ind1, ind2, distance, hdr=hdr, /binary, /output_dist'
      print
      print,'For binary, file should just be the front'
      return
  ENDIF 

  tt=systime(1)

  IF NOT keyword_set(binary) THEN BEGIN 
      print
      print,'Reading fits file: '+file
      print
      
      ind1 = mrdfits(file, 0, hdr)
      ind2 = mrdfits(file, 1)
      IF keyword_set(output_dist) THEN distance = mrdfits(file, 2)
  ENDIF ELSE BEGIN 

      ;; "file" is just the front of filenames
      infile1 = file + '_ind1.bin'
      infile2 = file + '_ind2.bin'
      IF keyword_set(output_dist) THEN infile3 =  file + '_dist.bin'
      numfile = file + '_num.dat'

      print
      print,'Reading binary files: '
      print,infile1
      print,infile2
      IF keyword_set(output_dist) THEN print,infile3
      print,'Num file: ',numfile
      print

      openr, lun, numfile, /get_lun
      ntotal=0L
      readf, lun, ntotal
      free_lun, lun

      ;; Read first set of indices
      openr, lun1, infile1, /get_lun  
      ind1 = lonarr(ntotal)
      readu, lun1, ind1
      free_lun, lun1

      ;; Read second set of indices
      openr, lun2, infile2, /get_lun
      ind2 = lonarr(ntotal)
      readu, lun2, ind2
      free_lun, lun2

      IF keyword_set(output_dist) THEN BEGIN 
          ;; Read the distances
          openr, lun3, infile3, /get_lun
          distance = dblarr(ntotal)
          readu, lun3, distance
          free_lun, lun3
      ENDIF 

  ENDELSE 

  ptime,systime(1)-tt

  return
END 
