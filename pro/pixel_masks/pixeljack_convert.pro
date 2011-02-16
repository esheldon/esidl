PRO pixeljack_convert

  ;; convert ascii files to .st files

  dir = sdssidl_config('shapecorr_dir')+'masks/pixel_mask_jackknife_samples/'
  files = file_search(dir + 'jackknife_samples*.dat', count=nf)

  outfiles = repstr(files, '.dat', '.st')

  tst = {pixelnum:0L, res: 0, jackknife_sample:0}
  FOR i=0L,nf-1 DO BEGIN 

      print,'-------------------------------------------------'
      read_struct, files[i], tst, struct
      
      print
      print,'Writing to file: ',outfiles[i]
      write_idlstruct, struct, outfiles[i]

  ENDFOR 

END 
