PRO new_parlist, parlist, addstruct

  IF n_params() EQ 0 THEN BEGIN
      print,'-Syntax: parlist, parlist, addstruct'
      return
  ENDIF 

  parlist = ['type', $
             'counts_model', $
             'counts_modelerr', $
             'petror50',$
             'petror50err', $
             'petror90', $
             'petror90err']


  addstruct = create_struct($
                             'e_d_bit',0,$
                             'classification', intarr(5) )

  return
END 
