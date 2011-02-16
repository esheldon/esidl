PRO vagc_match_local_bigcat, struct

  ;; make one big cat which just contains the matching info and
  ;; the redshift

  indir = vagc_lensinput_dir()

  files = file_search(indir+'stripe*_vagc_matchlocal.fit')

  outfile = indir + 'vagc_matchlocal.fit'
  print
  print,'Will write to file: ',outfile

  columns = ['run','rerun','camcol','field','id','match_rerun','match_id',$
             'mjd','plate','fiberid',$
             'z','z_err','zwarning']

  mrdfits_multi, files, struct, columns=columns

  print
  print,'Writing to file: ',outfile
  mwrfits, struct, outfile, /create

END 
