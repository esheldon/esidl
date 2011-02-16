pro dbext_struct, list, items, struct
;+
; NAME:
;       DBEXT_STRUCT.  
; PURPOSE:
;       Extract values from an IDL database.  The use of structures means there
;       is no practical limit to the number of columns that can be extracted.
;
; CALLING SEQUENCE:
;       dbext,list,items,struct
;
; INPUTS:
;       list - list of entry numbers to be printed, vector or scalar
;               If list = -1, then all entries will be extracted.
;               list may be converted to a vector by DBEXT 
;       items - standard item list specification.  See DBPRINT for 
;               the 6 different ways that items may be specified. 
;
; OUTPUTS:
;       struct: A structure containing the requested items for each object in
;               the list.
;
; EXAMPLE:
;       Extract all RA and DEC values from the currently opened database, and
;       place into the IDL vectors, IDLRA and IDLDEC.
;
;               IDL> DBEXT,-1,'RA,DEC',struct
;
; HISTORY
;       version 2  D. Lindler  NOV. 1987
;       check for INDEXED items   W. Landsman   Feb. 1989
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Converted to return a stucture.  E. Sheldon 9-Mar-2004
;-
;*****************************************************************
  On_error,2

  if N_params() lt 3 then begin
      print,'Syntax - dbext_struct, list, items, struct'
      return
  endif

;  zparcheck,'DBEXT',list,1,[1,2,3,4,5],[0,1],'Entry List'

  db_item,items,it,ivalnum,idltype,sbyte,numvals,nbytes

  nitems = N_elements(it)
  nentries = db_info('entries')
  if max(list) GT nentries[0] then $
    message,db_info('name',0)+' entry numbers must be between 1 and ' + $
    strtrim(nentries[0],2)

; get item info.

  dbno = db_item_info('dbnumber',it)
  if max(dbno) eq 0 then dbno=0 $ ;flag that it is first db only
  else dbno=-1
  index = db_item_info('index',it)
  ind = where( (index ge 1) and (index ne 3), Nindex ) 

  nrows = n_elements(list)
  struct = dbmake_struct(items,idltype,numvals,nrows=nrows)

  if (Nindex eq nitems) and (dbno eq 0) then begin ;All indexed items?

      if N_elements(list) eq 1 then list = lonarr(1) + list
      for i=0,nitems - 1 do begin ;Get indexed items
          itind = it[ind[i]]

          dbext_ind,list,itind,dbno,var
          IF nrows GT 1 THEN BEGIN 
              struct.(i) = reform(var, /overwrite)
          ENDIF ELSE BEGIN 
              struct.(i) = var
          ENDELSE 
          var = 0
      endfor

  endif else begin     

      nvalues = db_item_info('nvalues',it)
      dbext_dbf_struct,list,dbno,sbyte,nbytes*nvalues,idltype,nvalues, struct

  endelse

  return
end
