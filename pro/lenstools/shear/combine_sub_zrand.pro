PRO combine_sub_zrand, files, rlensum, shstruct, zstruct, shhdr, lhdr, $
                       sh_outfile=sh_outfile

  IF n_elements(sh_outfile) NE 0 THEN BEGIN 

      lensum_file = repstr(sh_outfile, '_N', '_lensum_N')
      zfile = repstr(sh_outfile, '_N', '_z_N')

      print,'Output files: '
      colprint,[sh_outfile, lensum_file, zfile]
  ENDIF 

  delvarx, rlensum, shstruct, zstruct, shhdr, lhdr

  hdr = headfits(files[0], ext=1)
  mrdfits_multi, files, rlensum, /diff

  s = sort(rlensum.zindex)
  rlensum = temporary(rlensum[s])

  binsize = sxpar( hdr, 'BINWIDTH' )
  rminkpc = sxpar( hdr, 'RMINKPC' )
  rmaxkpc = sxpar( hdr, 'RMAXKPC' )
  h = sxpar( hdr, 'H' )
  compcut = sxpar( hdr, 'COMPCUT' )
  IF compcut EQ 0 THEN compcut=0.9

  fxhclean
  hprint,hdr

  print
  print,'Combining lensum'
  combine_zlensum, rlensum, binsize, rminkpc, rmaxkpc, h, shstruct, $
                   compcut=compcut

  zs = create_struct('z', 0., $
                     'scritinv', 0., $
                     'clambda', double(0.), $
                     'ceta', double(0.) )

  zstruct = replicate(zs, n_elements(rlensum))
  zstruct.z = rlensum.z
  zstruct.scritinv = rlensum.scritinv
  zstruct.clambda = rlensum.clambda
  zstruct.ceta = rlensum.ceta

  shhdr = zshhdr(shstruct)  
  lhdr  = zsumhdr(shstruct)

  IF n_elements(sh_outfile) NE 0 THEN BEGIN 

      IF NOT fexist(sh_outfile) THEN BEGIN 
          print,'sh file: ',sh_outfile
          mwrfits2, shstruct, sh_outfile, shhdr, /destroy, /create
      ENDIF ELSE message,'Not overwriting file: '+sh_outfile

      IF NOT fexist(lensum_file) THEN BEGIN 
          print,'lensum file: ',lensum_file
          mwrfits2, rlensum, lensum_file, lhdr, /destroy, /create
      ENDIF ELSE message,'Not overwriting file: '+lensum_file

      IF NOT fexist(zfile) THEN BEGIN 
          print,'zfile: ',zfile
          mwrfits2, zstruct, zfile, /destroy, /create
      ENDIF ELSE message,'Not overwriting file: '+zfile

  ENDIF 

  return
  
END 
