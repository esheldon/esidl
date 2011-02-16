function positive, x

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Take the positive part of a number.
;
; Author : B.C. Kelly
;
; INPUT : X, The number or array to take the positive part of.
;
; RETURNS : The positive part of X, defined to be X if X > 0 and 0
;           otherwise.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if n_params() eq 0 then begin
    print, 'Syntax- Result = positive(x)'
    return, 0
endif

return, x > 0
end
