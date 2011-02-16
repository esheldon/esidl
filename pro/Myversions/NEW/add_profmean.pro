PRO add_profmean, file, outfile

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: add_profmean, file [, outfile]'
      return
  ENDIF 

  dirsep, file, dir, filename
  outfile = dir + 'prof_' + filename

  print
  print,'Output file: ',outfile
  print

  struct = mrdfits(file, 1)

  addtags = ['profmean','proferr']
  addsdsstag, struct, addtags, outstruct

  nst = n_elements(outstruct)
  addst =  create_struct('dn', 0.)
  addstruct = replicate(addst, nst)
  combine_structs, temporary(outstruct),addstruct,outstruct

  fit_dn, outstruct, dn, clr=2
  outstruct.dn = dn

  mwrfits, outstruct, outfile, /create

  return
END 
