PRO vagc_lensinput_name, stripes, rmin, rmax, name, randnum=randnum, rlrgMask=rlrgMask

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: vagc_lensinput_name, stripes, rmin, rmax, name, randnum=randnum, /rlrgMask'
      return
  ENDIF 

  name = spectra_name(stripes, /lss)
  rminstr = ntostr( round(rmin) )
  rmaxstr = ntostr( round(rmax) )
  IF keyword_set(rlrgMask) THEN rlrgStr='_rlrgMask' ELSE rlrgStr = ''

  ext = rlrgStr+'_rmin'+rminstr+'_rmax'+rmaxstr+'.fit'

  name = repstr(name, '.fit', ext)
  
  IF n_elements(randnum) NE 0 THEN BEGIN 
      name = repstr(name, 'sample', 'random-'+ntostr(randnum)+'.sample')
  ENDIF 

END 
