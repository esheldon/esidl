PRO rename_tsobj, files, newdir, newfront, newfiles, nchar

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: rename_tsobj, files, newdir, newfront [, newfiles, nchar]'
      return
  ENDIF 

  nf = n_elements(files)
  newfiles = strarr(nf)

  FOR i=0, nf-1 DO BEGIN
      
      name = str_sep(files[i], 'tsObj')
      newfiles[i] = newfront + name[1]     ; name[1] is the tail

  ENDFOR 

  nchar = strlen(newfiles[0])   ; length of new file name
  newfiles = newdir + newfiles  ; add on the directory

  return 
END 
