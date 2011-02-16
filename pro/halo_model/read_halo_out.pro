PRO read_halo_out, xifile, dsigfile, struct

  ;; Read the outputs from waynes code and combine into one
  ;; struct.  The outputs are co-moving, so also compute the
  ;; values in physical units, along with the values at the
  ;; common range of radii

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: read_halo_out, xifile, dsigfile, struct'
      return
  ENDIF 

  xist = create_struct('z', 0.0, $
                       'r', 0.0, $
                       'ximm', 0.0, $
                       'xigg', 0.0, $
                       'xigm', 0.0)

  dsigst = create_struct('z', 0.0, $
                         'r', 0.0, $
                         'deltasig', 0.0)
  
  rdfloatstr, xifile, xist, txist
  rdfloatstr, dsigfile, dsigst, tdsigst

  ;; find out how many of each z there are
  ntot = n_elements(txist)
  rd = rem_dup(txist.z)
  uz = txist[rd].z

  nz = n_elements(rd)
  neach = ntot/nz
  print,ntot,nz,neach

  arrval = fltarr(neach)
  struct = create_struct('z',    0.0, $
                         'r',    arrval, $
                         'ximm', arrval, $
                         'xigg', arrval, $
                         'xigm', arrval, $
                         'deltasig', arrval, $
                         'r_phys', arrval, $
                         'deltasig_phys', arrval)

  struct = replicate(struct, nz)

  beg = 0L
  FOR i=0L, nz-1 DO BEGIN 

      struct[i].z = uz[i]
      
      struct[i].r    = txist[beg:beg+neach-1].r
      struct[i].ximm = txist[beg:beg+neach-1].ximm
      struct[i].xigg = txist[beg:beg+neach-1].xigg
      struct[i].xigm = txist[beg:beg+neach-1].xigm
      struct[i].deltasig = tdsigst[beg:beg+neach-1].deltasig

      ;; convert from co-moving to physical
      struct[i].r_phys = struct[i].r/(1 + struct[i].z)
      struct[i].deltasig_phys = struct[i].deltasig*(1 + struct[i].z)^2

      beg = beg + neach

  ENDFOR 
  txist = 0
  tdsigst = 0

  ;; Find the radii that are common to all z in physical r
  nr = n_elements(struct[0].r_phys)
  maxr = min(struct.r_phys[nr-1])
  minr = max(struct.r_phys[0])

  w=where( struct[0].r GE minr AND $
           struct[0].r LE maxr, nw)
  rcommon = struct[0].r[w]
  
  nwstr = ntostr(nw)
  typestr = 'fltarr('+nwstr+')'
  add_tags, temporary(struct), $
            ['rcommon', 'xigmcommon','deltasig_physcommon'], $
            [typestr,typestr,typestr], struct

  FOR i=0L, nz-1 DO BEGIN 
      struct[i].rcommon = rcommon
      struct[i].xigmcommon = $
        interpol(struct[i].xigm, struct[i].r_phys, rcommon)
      struct[i].deltasig_physcommon = $
        interpol(struct[i].deltasig_phys, struct[i].r_phys, rcommon)
  ENDFOR 

END 
