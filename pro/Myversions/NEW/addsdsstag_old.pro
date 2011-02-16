PRO addsdsstag, struct, intags, newstruct

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    ADDSDSSTAG
;       
; PURPOSE:
;   This goes and gets tags and their values 
;   from tsObj files for each object in struct.
;   outputs newstruct with tags added.
;   different from addsdsspar which adds tags
;   to whole list of fields and outputs new
;   files
;
; CALLING SEQUENCE:
;    addsdsstag, struct, taglist, newstruct
;
; INPUTS: 
;    struct: a structure containing photo parameters. must have 
;            run,rerun,camcol,field,id
;    taglist: list of photo tags to add.
;
; OPTIONAL INPUTS:
;    None.
;
; KEYWORD PARAMETERS:
;    None.
;       
; OUTPUTS: 
;    newstruct: a new structure containing the old structure as well as 
;         the requested tags (if valid).
;
; OPTIONAL OUTPUTS:
;    None.
;
; CALLED ROUTINES:
;    
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    Created: Sept. 7 2000. Erin Scott Sheldon UofMich.
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() LT 2 THEN BEGIN 
     print,'-Syntax: addsdsstag, struct, taglist, newstruct'
     print,''
     print,'Use doc_library,"addsdsstag"  for more help.'  
     return
  ENDIF 

  oldtags = tag_names(struct)

  taglist=intags
  match, oldtags, strupcase(taglist), mold, madd
  IF (mold[0] NE -1) THEN BEGIN 
      print,'-------------------------------------------------'
      print,'WARNING: The tag "'+taglist[madd]+'" is already in the structure!'
      IF n_elements(madd) EQ n_elements(taglist) THEN BEGIN
          print,'WARNING: All input tags are already in structure'
          print,'-------------------------------------------------'
          return
      ENDIF 
      remove, madd, taglist
      print,'-------------------------------------------------'
  ENDIF 

  sendtags = ['id',taglist]
  sendtags = sendtags[rem_dup(sendtags)]
  ntags =  n_elements(taglist)
  nobj = n_elements(struct)

  runs = struct[rem_dup(struct.run)].run
  nrun = n_elements(runs)

  FOR ir=0L, nrun-1 DO BEGIN 

      wr=where(struct.run EQ runs[ir])
      rerun = struct[wr[0]].rerun

      tt=rem_dup(struct[wr].camcol)
      camcols = struct[wr[tt]].camcol
      ncam = n_elements(camcols)
      FOR ic=0L, ncam-1 DO BEGIN 

          print,'--------------------------------------------------------'
          print,'Run: ',ntostr(runs[ir]),' Camcol: ',ntostr(camcols[ic])
          print,'--------------------------------------------------------'
          fetch_dir, runs[ir], camcols[ic], rerun, indir

          wc = where(struct[wr].camcol EQ camcols[ic])
          wc = wr[wc]

          tt = rem_dup(struct[wc].field)
          fields = struct[wc[tt]].field
          nf = n_elements(fields)
          FOR ifield=0L, nf-1 DO BEGIN 

              wf = where(struct[wc].field EQ fields[ifield])
              wf = wc[wf]

              tmp=0

              read_tsobj, indir, tmp, $
                start=fields[ifield], nf=1, taglist=sendtags, $
                tsobjstr=tsobjstr
              print
              IF n_elements(addstr) EQ 0 THEN BEGIN 
                  ;; Build the add structure
                  FOR it=0L, ntags-1 DO BEGIN 
                      IF tag_exist(tmp, taglist[it],index=wind) THEN BEGIN
                          IF n_elements(adds) EQ 0 THEN BEGIN 
                              adds = create_struct(taglist[it],$
                                                     tmp[0].(wind) )
                          ENDIF ELSE BEGIN 
                              adds = create_struct( adds, $
                                               taglist[it], tmp[0].(wind) )
                          ENDELSE 
                      ENDIF ELSE BEGIN 
                          IF n_elements(rmind) EQ 0 THEN rmind=it $
                          ELSE rmind=[rmind,it]
                      ENDELSE 

                  ENDFOR 
                  IF n_elements(adds) EQ 0 THEN BEGIN
                    print,'--------------------------------------------------------'
                    print,'WARNING: No good tags found'
                    print,'--------------------------------------------------------'
                      return
                  ENDIF 
                  IF n_elements(rmind) NE 0 THEN remove, rmind, taglist
                  print,'********************************************************'
                  print,'* Good tags: ',taglist
                  print,'********************************************************'
                  sendtags = ['id',taglist]
                  sendtags = sendtags[rem_dup(sendtags)]
                  addstr = replicate(adds, nobj)
;                  combine_structs, struct, addstr, newstruct
              ENDIF 

              photo_match,$
                struct[wf].run,struct[wf].rerun,struct[wf].camcol, $
                struct[wf].field, struct[wf].id, $
                tmp.run,tmp.rerun,tmp.camcol, $
                tmp.field, tmp.id, $
                mstr, mtmp

              tadd = addstr[mstr]
              copy_struct, tmp[mtmp], tadd
              addstr[mstr] = temporary(tadd)

;              photo_match,$
;                newstruct[wf].run,newstruct[wf].rerun,newstruct[wf].camcol, $
;                newstruct[wf].field, newstruct[wf].id, $
;                tmp.run,tmp.rerun,tmp.camcol, $
;                tmp.field, tmp.id, $
;                mstr, mtmp

;              newstruct[wf[mstr]].profmean = tmp[mtmp].profmean
;              newstruct[wf[mstr]].proferr = tmp[mtmp].proferr

          ENDFOR 

      ENDFOR 
  ENDFOR 

  combine_structs, struct, addstr, newstruct

return
END 
