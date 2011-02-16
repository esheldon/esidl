FUNCTION DIAG, vector

size_vector = n_elements(vector)

full_array = dblarr(size_vector, size_vector)
FOR I = 0L, size_vector-1L DO BEGIN
 full_array[I,I] = vector[I]
ENDFOR
RETURN, full_array
END
