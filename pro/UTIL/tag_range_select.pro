FUNCTION tag_range_select, struct, tagname, lowlim, highlim, $
  input_index=input_index, nkeeparray=nkeeparray

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: keeparray = tag_range_select(struct, tagname, lowlim, highlim, '
      print,'                                      input_index=, nkeeparray=)'
      print,'lowlim and highlim are inclusive limits'
      return,-1
  ENDIF 

  tag = where(tag_names(struct[0]) EQ strupcase(tagname), nw)  
  IF nw EQ 0 THEN message,'Unknown tag: '+tagname
  
  nstruct = n_elements(struct)
  IF n_elements(input_index) EQ 0 THEN BEGIN 
      index = lindgen(nstruct)
  ENDIF ELSE BEGIN 
      index = input_index
  ENDELSE 

  nbin = n_elements(lowlim)
  IF n_elements(highlim) NE nbin THEN $
    message,'lowlim and highlim must be same length'
  keeparray = ptrarr(nbin)
  nkeeparray = lonarr(nbin)

  FOR i=0L, nbin-1 DO BEGIN 

      keep = where(struct[index].(tag) GE lowlim[i] AND $
                   struct[index].(tag) LE highlim[i], nkeep)
      
      IF nkeep NE 0 THEN BEGIN 
          keep = index[keep]

          keeparray[i] = ptr_new(keep, /no_copy)
          nkeeparray[i] = nkeep
      ENDIF 
  ENDFOR 

  return,keeparray


END 
