PRO type_hist, pstruct,$
               unk=unk, $
               cray=cray, $
               def=def, $
               gal=gal, $
               ghost=ghost, $
               known=known, $
               star=star, $
               trail=trail, $
               sky=sky, $
               QSO=QSO

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: type_hist, pstruct [, '
      print,'                    unk=unk, '
      print,'                    cray=cray, '
      print,'                    def=def, '
      print,'                    gal=gal, '
      print,'                    ghost=ghost, '
      print,'                    known=known, '
      print,'                    star=star, '
      print,'                    trail=trail, '
      print,'                    sky=sky, '
      print,'                    QSO=QSO]'
      print,' Set the keywords to a named variable.  The indices'
      print,' of objects with that type will be returned.  If no'
      print,' objects are found then -1 is returned'
      return
  ENDIF 

  names = ['Unk','C-Ray','Defect','Gal','Ghost',$
           'Known Obj','Star','Star trail','Sky','QSO']
  numnames = n_elements(names)

  min = 0
  max = numnames-1

  defnames = replicate(' ',numnames)

  hist=histogram(pstruct.objc_type, min=min, max=max, reverse_indices=rev_ind)

  xtitle='Objc_type'
  xrange=[0, numnames-1]

  
  plothist, pstruct.objc_type, xstyle=1+4, xrange=xrange
  
  axis, xaxis=0, xticks=numnames-1, xtickn=names,xtitle=xtitle
  axis, xaxis=1, xticks=numnames-1, xtickn=defnames,xtitle=xtitle

  unk=-1
  cray=-1
  def=-1
  gal=-1
  ghost=-1
  known=-1
  star=-1
  trail=-1
  sky=-1
  QSO=-1

  FOR i=0, numnames-1 DO BEGIN 
      IF rev_ind(i) NE rev_ind(i+1) THEN BEGIN 
          
          w=rev_ind( rev_ind(i):rev_ind(i+1)-1 )
          CASE i OF 
              0: unk=w
              1: cray=w
              2: def=w
              3: gal=w
              4: ghost=w
              5: known=w
              6: star=w
              7: trail=w
              8: sky=w
              9: QSO=w
              ELSE: print,'What!'
          ENDCASE 
      ENDIF 
  ENDFOR 
return
END 
