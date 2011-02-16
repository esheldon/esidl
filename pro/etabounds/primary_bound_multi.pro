PRO primary_bound_multi, stripes, bound_arr, overall_bound

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: primary_bound_multi, stripes, bound_array, overall_bound'
      return
  ENDIF 

  nstripe = n_elements(stripes)

  bound_arr = create_struct('stripe', 0, $
                            'lammin', 0d,$
                            'lammax', 0d, $
                            'etamin',0d, $
                            'etamax',0d)

  bound_arr = replicate(bound_arr, nstripe)

  FOR i=0L, nstripe-1 DO BEGIN 

      primary_bound, stripes[i], tbound

      bound_arr[i].stripe = stripes[i]
      bound_arr[i].lammin = tbound.lammin
      bound_arr[i].lammax = tbound.lammax
      bound_arr[i].etamin = tbound.etamin
      bound_arr[i].etamax = tbound.etamax

  ENDFOR 

  overall_bound = create_struct('stripes', stripes, $
                                'lammin', 0d,$
                                'lammax', 0d, $
                                'etamin',0d, $
                                'etamax',0d)

  overall_bound.lammax = max(bound_arr.lammax)
  overall_bound.lammin = min(bound_arr.lammin)
  overall_bound.etamax = max(bound_arr.etamax)
  overall_bound.etamin = min(bound_arr.etamin)

END 
