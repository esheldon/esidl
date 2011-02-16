PRO huan_sens_addnoise, nrand, errdiff, e1, e2, weights, we1err, we2err

  
  ngal = n_elements(e1)
  nerr = n_elements(errdiff)

  IF ngal NE nerr THEN message,'ngal must equal nerr'

  we1err_array = fltarr(nrand)
  we2err_array = we1err_array

  FOR i=0L, nrand-1 DO BEGIN 

;      angles = arrscl( randomu(seed, ngal), 0.0, 2*!pi )
;      rotate_e1e2, angles, e1, e2, e1out, e2out

      e1out=e1
      e2out=e2
      e1out = e1out + errdiff*randomu(seed, ngal, /normal)
      e2out = e2out + errdiff*randomu(seed, ngal, /normal)

      wmom, e1out, blah, wmean, wsig, we1err, inputweights = weights
      wmom, e2out, blah, wmean, wsig, we2err, inputweights = weights

      we1err_array[i] = we1err
      we2err_array[i] = we2err

  ENDFOR 

  we1err = median(we1err_array)
  we2err = median(we2err_array)

END 
