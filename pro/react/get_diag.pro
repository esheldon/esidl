FUNCTION GET_DIAG, the_array

size_array = size(the_array)
IF (n_elements(size_array) lt 5) THEN BEGIN
 PRINT, "Not an ARRAY"
 RETURN, the_array
ENDIF 
x_size = size_array[1]
y_size = size_array[2]
IF (x_size ne y_size) THEN BEGIN
  PRINT, "Not a Square ARRAY"
  RETURN, the_array
ENDIF
the_vector = dblarr(x_size)
FOR I = 0, x_size-1 DO BEGIN
  the_vector[I] = the_array[I,I]
ENDFOR
RETURN, the_vector
END  
