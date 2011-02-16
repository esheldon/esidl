 
;+
; Project     : SOHO - CDS
;
; Name        : APPEND_ARR
;
; Category    : Utility
;
; Purpose     : Manageable method for concatanating arrays
;
; Syntax      : IDL> result=append_arr(input,extra)
;
; Inputs      : INPUT = array (or scalar) that requires appending
;               EXTRA = array (or scalar) to append
;
; Outputs     : concatanated arrays
;
; Keywords    : SAME = set to force matching same array types
;               NO_COPY = set to not create internal copy of
;                INPUT (it will be destroyed)
;
; History     : Written: 1-Oct-1998, Zarro (SM&A/GSFC)
;               Modified: 26-March-2000, Zarro - sped up with SIZE and CATCH
;              
;
; Contact     : DZARRO@SOLAR.STANFORD.EDU
;-


function append_arr,input,extra,same=same,no_copy=no_copy

;-- instead of checking each input, just trap any errors
;   and bailout gracefully

error=0
catch,error
if error ne 0 then begin
 catch,/cancel
 goto,cleanup
endif

;-- only do this check if /SAME is set

doit=1
if keyword_set(same) then begin
 s1=size(input)
 t1 =s1(s1(0)+1)
 s2=size(extra)
 t2=s2(s2(0)+1)
 doit=(t1 eq t2)
endif

if doit then begin
 if keyword_set(no_copy) and (n_elements(extra) gt 0) then $
   return,[temporary(input),extra] else $
    return,[input,extra] 
endif

cleanup:


if func_exist(input) then return,input
if func_exist(extra) then return,extra

return,0

end
