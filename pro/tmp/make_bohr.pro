PRO make_bohr

  dir='/sdss3/usrdevel/esheldon/BCG/'
  outfile = dir + 'Bohringer_cat_000814.fit'

  infile = dir + 'Bohringer_cat_000814.txt'
  readcol, infile, $
            name,rra,rdec,flux,fluxerr,nphot,expos,ctr,nh,hr,hrerr, hrdev,ext,lum, $
    format='A,   F,  F,   F,   F,      F,    F,    F,  F, F, F,     F,    F,  F'

  infile = dir + 'jim_tim_000814.txt'
  readcol, infile, $
            run, rerun, field, camcol, id, z, ra, dec, i, $
    format='I,   I,     I,      I,     L,  F, F,  F,   F'

  infile = dir + 'hand_scat_cat_000814.txt'
  readcol, infile, $
            hindex, hrun, hrerun, hfield, hcamcol, hid, hra, hdec, hngals, hz, hgr, $
    format='A,      I,    I,      I,      I,       L,   F,   F,    F,     F,  F'


;  help,id,hid
;  colprint,z[sort(z)],hz[sort(hz)]
;  return

  rass_st=create_struct('name', '', $
                        'rass_ra', 0.0, $
                        'rass_dec', 0.0, $
                        'flux', 0.0, $
                        'fluxerr', 0.0, $
                        'nphot', 0.0, $
                        'expos', 0.0, $
                        'ctr', 0.0, $
                        'nh', 0.0, $
                        'hr', 0.0, $
                        'hrerr', 0.0, $
                        'hrdev', 0.0, $
                        'ext', 0.0, $
                        'l_x', 0.0)

  jim_tim_st=create_struct('run', 0, $
                           'rerun', 0, $
                           'camcol', 0, $
                           'field', 0, $
                           'id', 0L, $
                           'z', 0.0, $
                           'ra', 0.0, $
                           'dec', 0.0, $
                           'i', 0.0)

  hand_st=create_struct('use', 0, $
                        'z_orig', 0.0, $
                        'gr', 0.0, $
                        'ngals', 0L)

  nrasscol = n_elements(name)
  njimcol = n_elements(run)
  nhand = n_elements(hid)
  IF (nrasscol NE njimcol) OR (nrasscol NE nhand) THEN BEGIN 
      print,'What!'
      return
  ENDIF 
  rassin = replicate(rass_st, nrasscol)
  jimin = replicate(jim_tim_st, nrasscol)
  handin = replicate(hand_st, nrasscol)

  rassin.name = name
  rassin.rass_ra = rra
  rassin.rass_dec = rdec
  rassin.flux = flux
  rassin.fluxerr = fluxerr
  rassin.nphot = nphot
  rassin.expos = expos
  rassin.ctr = ctr
  rassin.nh = nh
  rassin.hr = hr
  rassin.hrerr = hrerr
  rassin.hrdev = hrdev
  rassin.ext = ext
  rassin.l_x = lum

  jimin.run = run
  jimin.rerun = rerun
  jimin.camcol = camcol
  jimin.field = field
  jimin.id = id
  jimin.z = z
  jimin.ra = ra
  jimin.dec = dec
  jimin.i = i

  handin.z_orig = hz
  handin.gr = hgr
  handin.ngals = hngals

  combine_structs, rassin, jimin, tmp

  ;; match up the hand scan's with Bohringer's stuff
  photo_match, tmp.run, tmp.rerun, tmp.camcol, tmp.field, tmp.id, $
               hrun, hrerun, hcamcol, hfield, hid, m1, m2
  tmp = tmp[m1]
  handin = handin[m2]
  hindex = hindex[m2]

  combine_structs, tmp, handin, newstruct

  w=replicate(1, nrasscol)
  FOR i=0L, nrasscol-1 DO BEGIN 
      test=strpos(hindex[i], 'X')
      IF test[0] EQ -1 THEN newstruct[i].use = 1
  ENDFOR 

  help,newstruct,/str

  mwrfits, newstruct, outfile,/create

  return
END 
