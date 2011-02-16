PRO read_stripe_specgal_inventory, inv_struct

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: read_stripe_specgal_inventory, inv_struct'
      return
  ENDIF 

  dir = sdssidl_config('SHAPECORR_DIR')+'spec_index/'
  file = dir + 'stripe_specgal_inventory.dat'

  readcol, file, stripe, ngal, format='I,L'

  inv_struct = create_struct('stripe', 0b, 'ngal', 0L)

  inv_struct = replicate(inv_struct, n_elements(stripe))
  inv_struct.stripe = stripe
  inv_struct.ngal = ngal

END 
