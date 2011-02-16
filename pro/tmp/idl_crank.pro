;$Id: idl_crank.pro,v 1.1.1.1 2004/11/11 20:17:25 esheldon Exp $
;
; Copyright (c) 1997-1998, Research Systems, Inc.  All rights reserved.
;       Unauthorized reproduction prohibited.



pro IDL_CRANK, w, s
;+
; NAME:
;    IDL_CRANK
;
;
; PURPOSE:
;    Replace elements of the sorted array "w" by their rank.
;    Identical observations ("ties") return the mean rank.
;
;    NOTE: This procedures is not a supported user-level routine.
;          It is a support routine for IDL statistics library functions.
;
;
; CATEGORY:
;    Analysis
;
;
; CALLING SEQUENCE:
;      IDL_CRANK, W
;
; 
; INPUTS:
;      W:  A sorted array
;
;
; OPTIONAL INPUTS:
;
;
;	
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;      W: IDL_CRANK replaces the input array W with its rank.
;
;
; OPTIONAL OUTPUTS:
;    s = total(f^3 - f) 
;    (f is the number of elements in each set of identical observations.)
;  
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;    X = [-1, 2, 3, 5, 6, 6, 9]
;    IDL_CRANK, X
; produces
;    X = [1.000, 2.000, 3.000, 4.000, 5.500, 5.500, 6.000 ]
; MODIFICATION HISTORY:
;
;       Tue Jan 27 16:50:31 1998, Scott Lett, RSI, Adapted from
;       earlier IDL library routines.
;
;		
;
;-


  n = n_elements(w)
  w = [0.0, w]  ;operate on elements w[1], ... , w[n] of the shifted
                ;n+1 element float array [w].
  s = 0.0
  j = 1L
  while(j lt n) do begin
    if(w[j+1] ne w[j]) then begin
      w[j] = j
      j = j+1
    endif else begin
      for jt = j+1, n do $
        if (w[jt] ne w[j]) then goto, case2
      jt = n + 1
      case2:
      rank = 0.5 * (j + jt - 1)
      for ji = j, jt-1 do $
        w[ji] = rank
      t = jt - j
      s = s + t^3 - t
      j = jt
    endelse
  endwhile
  if(j eq n) then w[n] = n
  w = w[1:*]
end
