PRO shapelet_makechris_files, dtype

  dtype = strlowcase(dtype)

  dir = '~/shapelet_outputs/'
  outfile = dir + 'vagc-shapelets-r-'+dtype+'-nmax15.fit'
  clr = 2

  inStruct = read_idlstruct_multi(file_search(dir + 'stripe*'+dtype+'*.st'))

  toutst = inStruct[0]

  outst = remove_tags(toutst, 'cmodel_counts')

  outst = create_struct($
                         outst, $
                         'cmodel_counts', fltarr(5), $
                         'cmodelflux_ivar', fltarr(5), $
                         'petrocounts', fltarr(5), $
                         'petroflux_ivar', fltarr(5), $
                         'abspetromag', fltarr(5) $
                       )

  outst = replicate(outst, n_elements(inStruct))

  struct_assign, inStruct, outst, /verbose

  ;; read extra vagc info and copy in. 
  columns=['run','match_rerun','camcol','field','match_id', $
           'cmodel_counts', 'cmodelflux_ivar', $
           'petroflux',   'petroflux_ivar', $
           'abspetromag']

  ;; Note: big file only contains a few tags, so we read in the individual
  ;;       files
  dir = sdssidl_config('shapecorr_dir')+'combined/vagc/'
  files = file_search(dir + 'stripe*matchlocal.fit')
  mrdfits_multi, files, vagc, columns=columns


  ;; match them up
  print,'Matching'
  outid = photoid(outst.run, outst.rerun, outst.camcol, outst.field, outst.id)
  vagcid = photoid(vagc.run, vagc.match_rerun, vagc.camcol, vagc.field, vagc.match_id)

  match, outid, vagcid, mout, mvagc, /sort
  help,outst, mout, mvagc
  IF n_elements(mout) NE n_elements(outst) THEN message,'not all matched!'


  print,'Copying'
  outst[mout].cmodel_counts = vagc[mvagc].cmodel_counts
  outst[mout].cmodelflux_ivar = vagc[mvagc].cmodelflux_ivar
  outst[mout].petrocounts = 22.5 - alog10(vagc[mvagc].petroflux)
  outst[mout].petroflux_ivar = vagc[mvagc].petroflux_ivar

  outst[mout].abspetromag[0] = vagc[mvagc].abspetromag[0]
  outst[mout].abspetromag[1] = vagc[mvagc].abspetromag[1]
  outst[mout].abspetromag[2] = vagc[mvagc].abspetromag[2]
  outst[mout].abspetromag[3] = vagc[mvagc].abspetromag[3]
  outst[mout].abspetromag[4] = vagc[mvagc].abspetromag[4]

  outst[mout].vagc_index = mvagc

  print
  print,'Writing file: ',outfile
  mwrfits, outst, outfile,/create

END 
