PRO make_dougfiles

  

  dir = '~/tmp/'
  file = dir + 'plist.dat'
  readcol, file, plates, mjds, fibers, format='I,I,I'

  st = create_struct('wavelength', 0.0d, $
                     'intensity', 0.0)

  nfib = n_elements(fibers)

  FOR i=0L, nfib-1 DO BEGIN 

;      stFile = repstr(files[i], '.fit', '.st')
;      pngFile = repstr(files[i], '.fit', '.png')
;      psFile = repstr(files[i], '.fit', '.ps')

      
      read_spec1d, plates[i], fibers[i], sp, file=fullfile

      dirsep, fullfile, tdir, file
      file = dir + file

      stFile = repstr(file, '.fit', '.st')
      pngFile = repstr(file, '.fit', '.png')
      psFile = repstr(file, '.fit', '.ps')

      plot_spec1d, sp, nsmooth=10
      print
      print,'Writing png file: ',pngFile
      write_png, pngFile, tvrd(/true)
      
      begplot, name=psFile, /color
      plot_spec1d, sp, nsmooth=10
      endplot 

      hdrStruct = create_struct('mjd', sp.mjd, $
                                'plateid', sp.plateid, $
                                'fiberid', sp.fiberid, $
                                'ra', sp.ra, $
                                'dec', sp.dec, $
                                'class', spec_type(sp.spec_cln) )
      
      struct = replicate(st, n_elements(sp.lambda) )
      struct.wavelength = sp.lambda
      struct.intensity = sp.spec

      print,'Writing st File: ',stFile
      write_idlStruct, struct, stFile, hdrStruct=hdrStruct, /ascii

      delvarx, fullfile
  ENDFOR 

END 

