;+
; NAME:
;  PHOTO_EXTRACTBYID
;
;
; PURPOSE:
;  Get information from the sdss tsObj files created by PHOTO for the input
;  sdss ids.
;
;
; CATEGORY:
;  sdss specific routine.
;
;
; CALLING SEQUENCE:
;  Two different modes:
;
;    struct = photo_extractbyid(runs, camcols, fields, ids, reruns=, $
;                               taglist=, filetype=, status=)
;  OR
;    struct = photo_extractbyid(photoids, reruns=, $
;                               taglist=, filetype=, status=)
;
;
; INPUTS:
;  sdss ids: runs,camcols,fields,ids  or the combined photoid. This
;            can be created using photoid=sdss_photoid(run,rerun,...)
;
;
; OPTIONAL INPUTS:
;  reruns: defaults to latest rerun.
;  taglist: the columns, or tags to read.  Defauls to all tags.
;  filetype: The file type, e.g. tsobj.  Default is tsobj. 
;
;
; KEYWORD PARAMETERS:
;  /corrected: read smaller adatc files, which contain additional 
;     measurements for each object. 
;  
;
; OUTPUTS:
;  struct: A structure containing info for all matched objects.
;
;
; OPTIONAL OUTPUTS:
;  status:  0 is success, otherwise failure.
;
; EXAMPLE:
;  
;
;
; MODIFICATION HISTORY:
;  1-August-2005: Erin Sheldon, UChicago
;
;-




FUNCTION photo_extractbyid, runs, camcols, fields, ids, reruns=reruns, $
  taglist_in=taglist_in, filetype=filetype, status=status, corrected=corrected, photoids=photoids

  status = 1
  tm = systime(1)

  ;; Check input parameters

  np = n_params()
  IF np EQ 4 THEN BEGIN 

      ;; User input runs,reruns,....
      IF n_elements(reruns) EQ 0 THEN BEGIN 
          reruns = sdss_rerun(runs)
      ENDIF 

      photoids = sdss_photoid(runs,reruns,camcols,fields,ids)
      idstruct = sdss_histid(runs, reruns, camcols, fields, status=hstatus)
  ENDIF ELSE IF np EQ 1 THEN BEGIN 
      ;; user input photoids
      extract_photoid, runs, truns, reruns, camcols, fields, ids
      idstruct = sdss_histid(runs, reruns, camcols, fields, status=hstatus)
  ENDIF ELSE BEGIN 
      print,'-Syntax: '
      print,'   struct = photo_extractbyid(runs, camcols, fields, ids, $'
      print,'                 reruns=, taglist=, filetype=, status=)'
      print,'  OR'
      print,'   struct = photo_extractbyid(photoids, $'
      print,'                 taglist=, filetype=, status=)'
      return,-1
  ENDELSE 


  IF hstatus NE 0 THEN return,-1

  ;; for matching
  IF n_elements(taglist_in) NE 0 THEN BEGIN 
      taglist = [taglist_in, 'run','rerun','camcol','field','id']
  ENDIF 

  if n_elements(filetype) eq 0 then filetype='tsObj'

  ;; Move along the tree and read from tsObj

  n_unique = idStruct.Nleaves
  ptrlist = ptrarr(n_unique)
  ptrIndex = 0L

  pruns = idStruct.runs
  FOR ri=0L, idStruct.nruns-1 DO BEGIN 
      run = (*pruns)[ri].run
      runStr = ntostr( run )
      preruns = (*pruns)[ri].reruns
      FOR rri=0L, (*pruns)[ri].nreruns-1 DO BEGIN 
          rerun = (*preruns)[rri].rerun
          rerunStr = ntostr( rerun )
          pcamcols = (*preruns)[rri].camcols
          FOR ci=0L, (*preruns)[rri].ncamcols-1 DO BEGIN 
              camcol = (*pcamcols)[ci].camcol
              camcolStr = ntostr(camcol)
              fieldStructs = *(*pcamcols)[ci].fields
              FOR fi=0L, (*pcamcols)[ci].nfields-1 DO BEGIN 
                  
                  field = fieldStructs[fi].field
                  pstruct = $
                    sdss_read(filetype,run,camcol,rerun=rerun,field=field,$
                              taglist=taglist)

                  ind = *fieldStructs[fi].indices

                  pid = sdss_photoid(pstruct)
                  match, photoids[ind], pid, mphotoid, mpid, /sort

                  IF mpid[0] NE -1 THEN BEGIN 

                      pstruct = pstruct[mpid]
                      ;; Now copy these into our pointerlist
                      ptrlist[ptrIndex] = ptr_new(pstruct, /no_copy)
                      ptrIndex = ptrIndex + 1
                  ENDIF 
              ENDFOR 


          ENDFOR ;; camcols
      ENDFOR ;; reruns
  ENDFOR ;; runs
  
  struct = combine_ptrlist(ptrlist)

  ptime,systime(1)-tm
  return,struct


END 
