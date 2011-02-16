PRO addsdsstag, struct, intags, newstruct, corrected=corrected

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
;    addsdsstag, struct, taglist, newstruct, corrected=corrected
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
;    /corrected: get info from the corrected files instead of 
;                the tsObj files.
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
     print,'-Syntax: addsdsstag, struct, taglist, newstruct, corrected=corrected'
     print,''
     print,'Use doc_library,"addsdsstag"  for more help.'  
     return
  ENDIF 

  IF keyword_set(corrected) THEN front='adatc'

  oldtags = tag_names(struct)

  taglist=intags
  match, oldtags, strupcase(taglist), mold, madd
  IF (mold[0] NE -1) THEN BEGIN 
      print,'-------------------------------------------------'
      print,'WARNING: The tag "'+taglist[madd]+'" is already in the structure!'
      IF n_elements(madd) EQ n_elements(taglist) THEN BEGIN
          print,'WARNING: All add tags are already in structure'
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


  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over runs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ten=ulong64(10)
  FOR ir=0L, nrun-1 DO BEGIN 

      wr=where(struct.run EQ runs[ir])

      tt = rem_dup(struct[wr].rerun)
      reruns = struct[wr[tt]].rerun
      nrer = n_elements(reruns)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Loop over reruns
      ;;;;;;;;;;;;;;;;;;;;;;;;;;
      FOR irer=0L, nrer-1 DO BEGIN 

          wrer = where(struct[wr].rerun EQ reruns[irer])
          wrer = wr[wrer]

          tt=rem_dup(struct[wrer].camcol)
          camcols = struct[wrer[tt]].camcol
          ncam = n_elements(camcols)

          ;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Loop over camcols
          ;;;;;;;;;;;;;;;;;;;;;;;;;;
          FOR ic=0L, ncam-1 DO BEGIN 

              print,'--------------------------------------------------------'
              print,'Run: ',ntostr(runs[ir]),' Rerun ',ntostr(reruns[irer]),$
                ' Camcol: ',ntostr(camcols[ic])
              print,'--------------------------------------------------------'
              fetch_dir, runs[ir], camcols[ic], reruns[irer], indir, $
                corrdir=corrdir
              IF keyword_set(corrected) THEN indir=corrdir
              
              wc = where(struct[wrer].camcol EQ camcols[ic])
              wc = wrer[wc]

              tt = rem_dup(struct[wc].field)
              fields = struct[wc[tt]].field
              nf = n_elements(fields)

              ;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Loop over fields
              ;;;;;;;;;;;;;;;;;;;;;;;;;;
              FOR ifield=0L, nf-1 DO BEGIN 
                  
                  wf = where(struct[wc].field EQ fields[ifield])
                  wf = wc[wf]
                  
                  tmp=0
                  
                  read_tsobj, indir, tmp, $
                    start=fields[ifield], nf=1, taglist=sendtags, $
                    tsobjstr=tsobjstr, front=front
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
                      zero_struct, adds
                      addstr = replicate(adds, nobj)
                  ENDIF 

                  ;; may have duplicates (e.g. on more than one plate)
                  ;; so do my own photo_match to cut out duplicates

                  super1=ulong64(struct[wf].id)+ulong64(struct[wf].field)*ten^6 $
                         + ulong64(struct[wf].camcol)*ten^11 $
                         + ulong64(struct[wf].rerun)*ten^13 $
                         + ulong64(struct[wf].run)*ten^16

                  super2=ulong64(tmp.id)+ulong64(tmp.field)*ten^6 $
                         + ulong64(tmp.camcol)*ten^11 $
                         + ulong64(tmp.rerun)*ten^13 $
                         + ulong64(tmp.run)*ten^16
                  
                  m1tmp = rem_dup(super1)
                  m2tmp = rem_dup(super2)

                  match,super1[m1tmp],super2[m2tmp],m1,m2,/sort

                  ;; may not match, especially if using corrected files
                  IF m1[0] NE -1 THEN BEGIN 
                      mstr = m1tmp[m1]
                      mtmp = m2tmp[m2]

                      tadd = addstr[wf[mstr]]
                      copy_struct, tmp[mtmp], tadd
                      addstr[wf[mstr]] = temporary(tadd)
                  ENDIF 


              ENDFOR ;; fields
          ENDFOR ;; camcols
      ENDFOR ;; reruns
  ENDFOR ;; runs

  combine_structs, struct, addstr, newstruct

return
END 
