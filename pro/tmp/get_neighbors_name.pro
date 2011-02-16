FUNCTION get_neighbors_name, run, camcol, rerun, field, id,front=front

  IF n_params() LT 5 THEN BEGIN 
      print,'-Syntax: name = get_neighbors_name(run, camcol, rerun, field, id, front=front)'
      print
      print,"Example:  name=get_neighbors_name(752,2,1,125,300,front='test')"
      print,' which returns test-000752-2-1-0125-0300.fit'
      print,''
      print,"DEFAULT is front='Neighb'"
      return,''
  ENDIF 
  
  IF n_elements(front) EQ 0 THEN front = 'Neighb'
  
  run_n = strtrim(string(run+1000000), 2)
  run_n = strmid(run_n,1,6)
  rerun_n = strtrim(string(rerun),2)
  camcol_n = strtrim(string(camcol),2)
  field_n = strtrim(string(field+10000L), 2)
  field_n = strmid(field_n, 1,4)
  id_n = strtrim(string(id+10000L), 2)
  id_n = strmid(id_n,1,4)
  return, front+'-'+run_n+'-'+camcol_n+'-'+rerun_n+'-'+field_n+'-'+id_n+'.fit'

END 
