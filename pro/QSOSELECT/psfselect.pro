pro  psfselect, struct, keep, dir=dir, radius=radius

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;       
; PURPOSE:
;	
;
; CALLING SEQUENCE:
;      
;                 
;
; INPUTS: 
;
; INPUT KEYWORD PARAMETERS:
;
;       
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; CALLED ROUTINES:
; 
; PROCEDURE: 
;	
;	
;
; REVISION HISTORY:
;	
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if N_params() eq 0 then begin
     print,'-Syntax: psfselect, struct, keepid, dir=dir, radius=radius'
     print,''
     print,'Use doc_library,""  for more help.'  
     return
  endif

s=systime(1)

make_flag_struct, fs

fs.manypetro='N'
fs.bright='N'
fs.edge='N'
fs.cr='N'
fs.satur='N'
fs.interp='N'
fs.too_large='N'

fs.deblended_as_psf='Y'

cindex=2

flag_select,struct, fs, cindex, index
help,index

w =where( struct[index].petrorad[2] LT 2 AND $
          struct[index].fibercounts[0] - struct[index].fibercounts[1] LE .8 $
          AND $
          struct[index].fibercounts[1] - struct[index].fibercounts[2] LE .5 $
          AND struct[index].fibercounts[0] LT 21.0 $
        )
index=index[w]
help,index


;;;  find things with 2 siblings.

n = n_elements(index)
keep = -1
sibkeep = -1
FOR i=0, n-1 DO BEGIN

  par=-1
  child = -1
  sib = -1
  get_family, struct, index[i], par, child, sib, /noprompt, /nodisplay,/silent

 
  IF (n_elements(sib) EQ 2) THEN BEGIN
    mag1=struct[sib[0]].fibercounts
    ug1 = mag1[0]-mag1[1]
    gr1 = mag1[1]-mag1[2]
    ri1 = mag1[2]-mag1[3]
    iz1 = mag1[3]-mag1[4]

    err1=struct[sib[0]].fibercountserr
    ugerr1 = sqrt(err1[0]^2 + err1[1]^2)
    grerr1 = sqrt(err1[1]^2 + err1[2]^2)
    rierr1 = sqrt(err1[2]^2 + err1[3]^2)
    izerr1 = sqrt(err1[3]^2 + err1[4]^2)

    mag2=struct[sib[1]].fibercounts
    ug2 = mag2[0]-mag2[1]
    gr2 = mag2[1]-mag2[2]
    ri2 = mag2[2]-mag2[3]
    iz2 = mag2[3]-mag2[4]
    err2=struct[sib[1]].fibercountserr
    ugerr2 = sqrt(err2[0]^2 + err2[1]^2)
    grerr2 = sqrt(err2[1]^2 + err2[2]^2)
    rierr2 = sqrt(err2[2]^2 + err2[3]^2)
    izerr2 = sqrt(err2[3]^2 + err2[4]^2)

    IF ( (abs(ug1 - ug2) LT 3*ugerr1 OR  abs(ug1 - ug2) LT 3*ugerr2 ) $
        AND (abs(gr1 - gr2) LT 3*grerr1 OR  abs(gr1 - gr2) LT 3*grerr2 ) $
;       AND (abs(ri1 - ri2) LT 3*rierr1 OR  abs(ri1 - ri2) LT 3*rierr2 ) $
;       AND (abs(iz1 - iz2) LT 3*izerr1 OR  abs(iz1 - iz2) LT 3*izerr2 ) $
       )  THEN BEGIN
            IF (keep[0] EQ -1 ) THEN keep = index[i] ELSE $
                 keep = [keep,index[i]]

            IF (sib[0] EQ index[i] ) THEN BEGIN
                 IF (sibkeep[0] EQ -1) THEN sibkeep = sib[1] ELSE $
                                            sibkeep = [sibkeep,sib[1]]
            ENDIF ELSE BEGIN
                 IF (sibkeep[0] EQ -1) THEN sibkeep = sib[0] ELSE $
                                            sibkeep = [sibkeep,sib[0]]
            ENDELSE 
    ENDIF
  ENDIF 

  IF ( i MOD 100 EQ 0 ) THEN BEGIN
    print,'i = ',strtrim(string(i),2)
    print,(systime(1) - s)/60.0,'min'
    help,keep
  ENDIF 

ENDFOR

help,keep

extract_stars, struct, 2, indices=stars
plot_colors, struct[stars], oplot_str=struct[keep],/size

print,''
print,'Hit a key to use obj_info (q to quit)'
key=get_kbrd(20)

IF (key EQ 'q') THEN return ELSE print,''

n =  n_elements(keep)

IF (n_elements(radius) EQ 0) THEN radius = 150
IF (n_elements(dir) EQ 0) THEN dir=''
FOR i=0, n-1 DO BEGIN

  obj_info, struct, [keep[i],sibkeep[i]], 2, dir=dir, radius=radius, /nofchart

  print,'Hit "y" to make files for this object or "q" to quit'
  key=get_kbrd(20)
  IF (key EQ 'y') THEN BEGIN 
    print,'Better to quit after first sibling'
    obj_info, struct, [keep[i],sibkeep[i]], 2, dir=dir, radius=radius,/gif,/ps
  ENDIF ELSE BEGIN
    IF (key EQ 'q') THEN return
  ENDELSE 

ENDFOR 

return
end




























