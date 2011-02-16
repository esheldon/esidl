PRO replace_front, dir, start, nframes, newfiles, nchar

  ;; Choose a new front end and directory.
  newfront = 'adatc'
  newdir = 'mydirectory/'       ;don't forget /

  ;; read in the tsObj file names you got the data from.
  fetch_file_list, dir, files, start, nframes

  ;; replace front end
  nf=n_elements(files)
  newfiles = strarr(nf)

  for i=0, nf-1 do begin
    
      tmp = str_sep(files[i], 'tsObj')
    
      newfiles[i] = newfront + tmp[1]     ; tmp[1] is tail
  endfor 

  nchar = strlen(newfiles[0])             ; length of new filename
  newfiles = newdir + newfiles            ; add on directory

return
END 
