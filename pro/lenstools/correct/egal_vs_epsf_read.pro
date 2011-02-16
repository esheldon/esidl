FUNCTION egal_vs_epsf_read, clr, hdrstruct=hdrstruct, original=original

  IF n_elements(clr) EQ 0 THEN BEGIN 
      print,'-Syntax: fitst=egal_vs_epsf_read(clr, hdr=hdr, /original)'
      return,-1
  ENDIF 
  
  IF keyword_set(original) THEN hir='' ELSE hir='_h' 
  file = $
    sdssidl_config('shapecorr_dir') + 'corr_egal_vs_epsf/'+$
    'egal_vs_epsf'+hir+'_'+!colors[clr]+'.st'

  print
  print,'Reading file: ',file
  fitstruct = read_idlstruct(file, hdrstruct=hdrstruct, /silent)
  return,fitstruct

END 
