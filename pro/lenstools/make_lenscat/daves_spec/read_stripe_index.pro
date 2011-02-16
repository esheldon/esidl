PRO read_stripe_index, struct

  dir = sdssidl_config('SHAPECORR_DIR')+'spec_index/'
  file = dir + 'stripe_index.dat'

  readcol, file, format='L,L,A',$
           runs, stripes, strips

  nn = n_elements(runs)

  s = create_struct('run', 0L, $
                    'stripe', 0, $
                    'strip', '')

  struct = replicate(s, nn)

  struct.run = runs
  struct.stripe = stripes
  struct.strip = strips

END 
