PRO make_truecolors

  ;; reads in the rgb.txt file and creates idl code 
  ;; that defines these colors.  

  infile = '/usr/X11R6/lib/X11/rgb.txt'
  outfile = '/sdss5/data0/esheldon/idl.lib/PLOTTING/rgbtrue.input'

  fmt='3I4,A25'
  readfmt, infile, fmt, r, g, b, names
  
  nn=n_elements(names)

  openw, lun, outfile, /get_lun
  printf, lun
  kept=0L
  FOR i=0L, nn-1 DO BEGIN 

      ;; remove those with spaces
      tname1 = ntostr(names[i])
      tname2 = repstr(tname1, ' ', '_')
      IF tname1 EQ tname2 THEN BEGIN 

          kept = kept+1
          printf,lun, '      R='+ntostr(r[i])+'L & '+$
                 'G='+ntostr(g[i])+'L & '+$
                 'B='+ntostr(b[i])+'L'
          printf,lun, "      defsysv, '!"+tname1+"', R + 256L*(G+256L*B)"
          printf, lun
      ENDIF 
  ENDFOR 

  free_lun, lun

  print,'Kept '+ntostr(kept)+'/'+ntostr(nn)

END 
