PRO runrgbview, imi, imr, img, r, g, b



  l = where(imi LT 0, nl)
  IF nl NE 0 THEN imi[l] = 32767L

  l = where(imr LT 0, nl)
  IF nl NE 0 THEN imr[l] = 32767L

  l = where(img LT 0, nl)
  IF nl NE 0 THEN img[l] = 32767L

  l = where(imi EQ 0, nl)
  IF nl NE 0 THEN imi[l] = 1000L

  l = where(imr EQ 0, nl)
  IF nl NE 0 THEN imr[l] = 1000L

  l = where(img EQ 0, nl)
  IF nl NE 0 THEN img[l] = 1000L

  rgbview2, imi, imr, img, r, g, b

return
END 
