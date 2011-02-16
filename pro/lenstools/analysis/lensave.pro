PRO lensave, struct, tag, ave, err, $
             element=element, radbin=radbin, addweight=addweight, $
             input_index=w

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: lensave, lensum, tagname, average, uncertainty, element=, /, radbin=, addweight=, input_index='
      print
      print,'Send radbin to use weights in particular bin'
      return
  ENDIF 

  nstruct = n_elements(struct)
  IF n_elements(w) EQ 0 THEN BEGIN 
      w = lindgen(nstruct)
  ENDIF  

  tags = tag_names(struct)
  w=where( tags EQ strupcase(tag), nw)

  nst=n_elements(struct)
  nweight=n_elements(addweight)

  IF tag_exist(struct, tag, ind=ti) THEN BEGIN 

      IF n_elements(radbin) EQ 0 THEN BEGIN 
          lensw = struct[w].weight
      ENDIF ELSE BEGIN 
          lensw = struct[w].wsum[radbin]
      ENDELSE 

      IF nweight NE 0 THEN BEGIN
          IF nst NE nweight THEN message,'addweight and struct must be same size'
          lensw = lensw*addweight
      ENDIF 
      lenswsum = total(lensw, /double)

      IF n_elements(element) EQ 0 THEN BEGIN 
          tagsum = total( lensw*struct[w].(ti), /double)
          ave = tagsum/lenswsum
          tagerrsum = total( lensw^2*( struct[w].(ti)-ave )^2, /double)
          err = sqrt(tagerrsum)/lenswsum
      ENDIF ELSE BEGIN 
          tagsum = total( lensw*struct[w].(ti)[element], /double)
          ave = tagsum/lenswsum
          tagerrsum = total( lensw^2*( struct[w].(ti)[element]-ave )^2, /double)
          err = sqrt(tagerrsum)/lenswsum
      ENDELSE 

  ENDIF ELSE BEGIN 
      print,'ERROR: Tag '+tag+' not found'
      return
  ENDELSE 

END 
