; This function is a subset of the  equivalent function in S+. 
; It takes a vector (f) and creates a new vector (g) that is one
; shorter in length.  The  values of g are f[i+1] - f[i]. 

FUNCTION DIFF, f

num = n_elements(f)
g = dblarr(num-1)
FOR I = 0, num -2 DO BEGIN
  g(I) = f(I+1) - f(i)
ENDFOR
RETURN, g
END
