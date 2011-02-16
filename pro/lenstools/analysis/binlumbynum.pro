PRO binlumbynum, lensum, clr, perc, wStruct, good=good

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: binlumbynum, lensum, clr, perc, wStruct, good=good'
      return
  ENDIF 
  
  nperc = n_elements(perc)

  wz = getztag(lensum[0])

  print
  print,'There are '+ntostr(nperc)+' Percentages: '
  FOR i=0L, nperc-1 DO print, perc[i]

  binlumbynum_absrange, minAbsMag, maxAbsMag

  maxz = 0.6
  minz = 0.02

  good=where( (lensum.(wz) LE maxz) AND (lensum.(wz) GE minz) AND $
              (lensum.absPetroMag[clr] GT minAbsMag[clr]) AND $
              (lensum.absPetroMag[clr] LT maxAbsMag[clr]), ngood)

  s=sort(lensum[good].absPetroMag[clr])
  s=reverse(s)
  good=good[s]

  beg = 0L
  FOR i=0L, nperc-1 DO BEGIN 

      IF i EQ nperc-1 THEN BEGIN 
          w = good[beg:ngood-1]
      ENDIF ELSE BEGIN 
          binsize = long(perc[i]*ngood)
          w = good[beg:beg+binsize-1]
      ENDELSE 

      beg = beg+binsize

      tagname = 'w'+ntostr(i+1)
      IF i EQ 0 THEN BEGIN 
          
          wStruct = create_struct(tagname, w)
      ENDIF ELSE BEGIN 
          wStruct = create_struct(wStruct, tagname, w)
      ENDELSE 
  ENDFOR 
  

END 
