; This function is the equivalent of f <- f[-i] is S+.
; It takes a vector of length n and returns a vector of
; length n - 1 .  The new vector has the i^th component
; removed.
FUNCTION SHRINK, g, j

orig_index = j
old_array = g
IF ((n_elements(old_array) -1) lt orig_index) THEN BEGIN
   PRINT, 'Cannot shrink function'
   RETURN, g
ENDIF
h = dblarr(n_elements(old_array) -1)
FOR k = 0, orig_index -1 DO BEGIN
  h[k] = old_array[k]
ENDFOR
FOR k = orig_index+1, n_elements(old_array) -1 DO BEGIN
  h[k-1] = old_array[k]
ENDFOR
new_array = h
RETURN, new_array 
END
