PRO mysql_tsobj_stuff_fivetags

  ;; Get a list of the tags that are for each bandpass and output to a 

  alltags = mysql_tsobj_stuff_taglist(struct=struct)

  nt = n_elements(alltags)
  FOR i=0L, nt-1 DO BEGIN 

      num = n_elements(struct.(i))

      IF num EQ 5 THEN BEGIN

          add_arrval, "'" + alltags[i] + "'", fivetags

      ENDIF 

  ENDFOR 

  print,'fivetags = ['+strjoin(fivetags, ', ')+']'

END 
