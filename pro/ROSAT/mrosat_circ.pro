
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: mrosat_circ
;       
; PURPOSE: circle multiple matches to the rosat catalog.
;	
;
; CALLING SEQUENCE: mrosat_circ, matchstr, radius, clr, fieldzero=fieldzero,
;                   rerun=rerun
;      
; INPUTS: matchstr:  the matched photo/rosat structure
;         radius:    radius to make finding chart with.
;         clr:       color to make finding chart with
;
; CALLED ROUTINES: 
;                    ROSAT_CIRC
; 
; PROCEDURE: 
;	Loop by rosat index throught the matched structure.
;	
;
; REVISION HISTORY:
;	Erin Scott Sheldon  Umich 5/25/99
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;
;;; Subroutine
;;;;;;;;;;;;;;;;;;

PRO prompt_next, i, n, iw, w, nw, beginold

print,''
print,'Go to next rosat object (y/n)? p for previous'
key = get_kbrd(20)


iw = nw
IF (key EQ 'n') THEN BEGIN
  print,'Quitting'
  i = n
ENDIF ELSE IF (key EQ 'p') THEN BEGIN
  i = beginold
ENDIF ELSE BEGIN
  i = w[nw-1] + 1
  IF (i EQ n) THEN print,'That was the last one!'
ENDELSE


return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  MAIN PROCEDURE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro  mrosat_circ, matchstr, radius, clr, fieldzero=fieldzero,rerun=rerun


  if N_params() eq 0 then begin
     print,'-Syntax: mrosat_circ, matchstr, radius, clr'
     print,''
     print,'Use doc_library," mrosat_circ"  for more help.'  
     return
  endif

  rosatname = '/sdss3/data2/rosat/rass-bsc-1.2rxs.fit'
  rosat = mrdfits(rosatname, 1, hdr)
  
  n = n_elements(matchstr)
  IF (NOT keyword_set(fieldzero)) THEN fieldzero=0
  IF (NOT keyword_set(rerun)) THEN rerun=0
  i=0
  beginnew=i
  iw = 0
  WHILE (i LE n-1) DO BEGIN
    
    beginold = beginnew
    beginnew = i
    index = matchstr[i].rosat_index
    indexstr = strtrim(string(index),2)
    w = where(matchstr.rosat_index EQ index)
    nw = n_elements(w)

    print,''
    print,'*******************************************'
    print,'* Viewing matches to rosat index ',indexstr
    print,'* ',strtrim(string(nw),2),' matches found'
    print,'*******************************************'
    
    iw = 0
    WHILE (iw LE nw-1) DO BEGIN

      idstr = strtrim(string(matchstr[i].id), 2)
      fieldstr = strtrim(string(matchstr[i].field),2)
      print,''
      print,'*******************************************'
      print,'* Photo id: ',idstr, ' Field: ',fieldstr
      print,'*******************************************'
      print,''

      print,'Do you wish to view this object? (y/n)?'
      key = get_kbrd(20)
      IF (key NE 'n') THEN rosat_circ, matchstr, i, radius, clr, rosat=rosat, $
                               fieldzero=fieldzero,rerun=rerun

      IF (iw NE nw-1) THEN BEGIN
        print,''
        print,'Next match to rosat index ',indexstr,' (y/n)?'
        key = get_kbrd(20)
        IF (key NE 'n') THEN BEGIN
          iw = iw+1
          i = i+1
        ENDIF ELSE BEGIN
          prompt_next, i, n, iw, w, nw, beginold
        ENDELSE 
      ENDIF ELSE BEGIN
        prompt_next, i, n, iw, w, nw, beginold
      ENDELSE

    ENDWHILE

  ENDWHILE 

return
end





















