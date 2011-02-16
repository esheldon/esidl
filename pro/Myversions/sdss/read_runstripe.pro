PRO read_runstripe, struct, local=local

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: read_runstripe, struct, local=local'
      return
  ENDIF 

  on_error,2

  IF NOT keyword_set(local) THEN BEGIN 
      read_runlist, struct
      return
  ENDIF ELSE BEGIN 

      dir = '/net/cheops1/data0/imaging/dbm/collate/'
      file = dir + 'stripelist_fermi.dat'
      
      ;; don't read location
      readcol, file, format='L,A,L,X,X', $
               stripe, strip, run
      
      nn = n_elements(stripe)
      
      s = create_struct('stripe', 0, $
                        'strip', '', $
                        'run', 0L)
      
      struct = replicate(s, nn)

      struct.stripe = stripe
      struct.strip = strip
      struct.run = run
      
      return

  ENDELSE 

END 
