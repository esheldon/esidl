FUNCTION mlabels, axis, index, value

  ;; should change this each time
  IF value EQ 0 THEN return,'0'
  units = !axis_units
  units_string = !axis_units_string

  use_val = long(round(value/units))
  use_val_string = ntostr(use_val)

  return, use_val_string+!csym.times+units_string
  
  
END 
