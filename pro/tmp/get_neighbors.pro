
PRO get_neighbors, cat, leftrun, leftrerun, rightrun, rightrerun, rundiff, outdir=outdir, $
                   taglist=taglist

  IF n_params() LT 6 THEN BEGIN
      print,'-Syntax: get_neighbors, cat, leftrun, leftrerun, rightrun, rightrerun, rundiff, outdir=outdir, taglist=taglist'
      print,' e.g. leftrun=752 leftrerun=1 rightrun=756  rightrerun=1 rundiff=185'
      return
  ENDIF 

  t=systime(1)

  IF n_elements(outdir) EQ 0 THEN outdir=''

  ncat = n_elements(cat)

  FOR i=0L, ncat-1 DO BEGIN
      
      tmp=0
      struct = 0
      struct1 = 0
      struct2 = 0
      struct3 = 0

      run1 = cat[i].run
      camcol1 = cat[i].camcol
      field1 = cat[i].field
      rerun1 = cat[i].rerun
      print,'---------------------------------------------'
      print,'Object:  ',ntostr(i),'  Run: ',ntostr(run1),'  Camcol: ',ntostr(camcol1)
      print,'---------------------------------------------'

      fetch_dir, run1, camcol1, rerun1, dir1

      ;;find_radec, ra[i], dec[i], run1, camcol1, field1

      read_tsobj, dir1, struct1, start=field1-1, nframes=3, tsobjstr=tsobjstr,verb=0, $
            taglist=taglist

      ;;rdphcol, run1, camcol1, struct1, start=field1-1, nf=3
            
      IF run1 EQ leftrun THEN BEGIN
          field2 = field1 + rundiff
          newrun = rightrun
          newrerun = rightrerun
          IF camcol1 EQ 1 THEN BEGIN
              special = 1
              camcol2 = camcol1
          ENDIF ELSE BEGIN
              special = 0
              camcol2 = camcol1 - 1
              camcol3 = camcol1
          ENDELSE 
      ENDIF ELSE BEGIN
          field2=field1 - rundiff
          newrun = leftrun
          newrerun=leftrerun
          IF camcol1 EQ 6 THEN BEGIN
              special = 1
              camcol2 = camcol1
          ENDIF ELSE BEGIN
              special = 0
              camcol2 = camcol1
              camcol3 = camcol1 + 1
          ENDELSE 
      ENDELSE 

          
      IF special THEN BEGIN
          
          print,'Neighboring fields from run: ',ntostr(newrun),'  Camcol: ',ntostr(camcol2)
          fetch_dir, newrun, camcol2, newrerun, dir2
          read_tsobj, dir2, struct2, start=field2-1, nf=3, tsobjstr=tsobjstr,verb=0, $
            taglist=taglist
          concat_structs, struct1, struct2, struct
              
      ENDIF ELSE BEGIN

          print,'Neighboring fields from run: ',ntostr(newrun),'  Camcol: ',ntostr(camcol2)
          fetch_dir, newrun, camcol2, newrerun, dir2
          read_tsobj, dir2, struct2, start=field2-1, nf=3, tsobjstr=tsobjstr,verb=0, $
            taglist=taglist
          concat_structs, struct1, struct2, tmp

          print,'Neighboring fields from run: ',ntostr(newrun),'  Camcol: ',ntostr(camcol3)
          fetch_dir, newrun, camcol3, newrerun, dir3
          read_tsobj, dir3, struct3, start=field2-1, nf=3, tsobjstr=tsobjstr,verb=0, $
            taglist=taglist
          concat_structs, tmp, struct3, struct

      ENDELSE 

      ;; remove duplicates
      help,struct
      print,'Removing duplicates'
      rmclose_radec, struct, good
      struct = temporary(struct[good])
      help,struct

;      print
;      help,struct1
;      help,struct2
;      help,struct3
;      help,struct
;      print

      ;; output structures
      allname = outdir+get_neighbors_name(run1, camcol1, rerun1, field1, cat[i].id)
      print,'Output file: ',allname
      mwrfits, struct, allname, /create

  ENDFOR 


  ptime, systime(1)-t

return
END 
