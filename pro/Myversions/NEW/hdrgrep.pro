pro hdrgrep, header, substring, found, keepcase=keepcase, linenum=linenum, $
             silent=silent

;+
; NAME:
;     HDRGREP
;
; PURPOSE:
;       Find a substring in a FITS header (or any other string array)
;
; CALLING SEQUENCE:
;       HGREP, header, substring, [/KEEPCASE, /LINENUM ]
;
; INPUTS: 
;       header -  FITS header or other string array
;       substring - scalar string to find in header
;
; OPTIONAL INPUT KEYWORDS:
;       /KEEPCASE: if set, then look for an exact match of the input substring 
;                 Default is to ignore case .
;       /LINENUM: if set, prints line number of header in which
;                substring appears 
;       /SILENT: don't print to screen.
;
; OUTPUTS:
;       found: 1 if string found, 0 if not. Also, results are printed to screen
;              if /silent is not set
;
; EXAMPLE: 
;       Find every place in a FITS header that the word 'aperture'
;       appears in lower case letters and print the element number 
;       of the header array:
;       
;       IDL> hgrep, header, 'aperture', /keepcase, /linenum
;
; HISTORY: 
;       Written, Wayne Landsman (Raytheon ITSS)      August 1998
;       Adapted from STIS version by Phil Plait/ ACC November 14, 1997
;       Adapted from HGREP to return found variable.
;           Erin Scott Sheldon 26-Nov-2002
;-

   if (N_params() LT 2) then begin
      print,'Syntax - HDRGREP, header, substring, [ found, /KEEPCASE, /LINENUM, /SILENT]'
      return
   endif

   found = 0
   if N_elements(header) eq 0 then begin
      print,'first parameter not defined. Returning...'
      return
   endif
   hh = strtrim(header,2)

   if keyword_set(keepcase) then $
         flag = strpos(hh,substring) $
   else  flag = strpos(strlowcase(hh),strlowcase(substring))
     

   g = where(flag NE -1, Ng)
   if Ng GT 0 then begin
       found=1
       IF NOT keyword_set(silent) THEN BEGIN 
           if keyword_set(linenum) then $
             for i = 0, Ng-1 do print, string(g[i],f='(i4)') + ': ' + hh[g[i]] $
           else $
             for i = 0, Ng-1 do print,hh[g[i]] 
       ENDIF 
   endif
      
   return
   end
