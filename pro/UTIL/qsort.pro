;; This doesn't work!!


PRO qsort_swap, array, ind1, ind2

  t = array[ind1]
  
  array[ind1] = array[ind2]
  
  array[ind2] = t

END 

PRO qsort, arr, abeg, aend

  IF ( aend GT (abeg + 1) ) THEN BEGIN 

      piv = arr[abeg]
      l = abeg + 1
      r = aend

      WHILE (l LT r) DO BEGIN 

          IF (arr[l] LE piv) THEN BEGIN 
              l = l+1
          ENDIF ELSE BEGIN 
              qsort_swap, arr, l, r
              r = r-1
          ENDELSE 

      ENDWHILE 
      
      qsort_swap, arr, l, abeg
      l = l-1
      qsort, arr, abeg, l
      qsort, arr, r, aend

  ENDIF 

END 
