
PRO vagc_make_jackie_cat_all

  
  outDir = '/net/cheops2/data1/lensinputs/jackie/'
  outFile = outdir + 'jackiecat_all.st'

  sdss_lss_dir = sdssidl_config('lss_dir')
  sdss_lss_vers = sdssidl_config('lss_vers')
  sdss_spec_dir = sdssidl_config('spec_dir')
  sdss_vagc_dir = sdssidl_config('vagc_dir')
  sdss_vagc_vers = sdssidl_config('vagc_vers')
  sdss_shapecorr_dir = sdssidl_config('shapecorr_dir')

  catname = vagc_catname(/lss, letter=letter, post=post, sample=sample)


  lssDir = sdss_lss_dir + sdss_lss_vers + "/" + letter + "/" + post + "/"
  vagcDir = sdss_vagc_dir + sdss_vagc_vers+'/'
  calibDir = vagcDir + 'sdss/parameters/'
  maskDir = sdss_shapecorr_dir + 'masks/'

  lssFile          = lssDir  + $
    'post_catalog.'+sdss_lss_vers +letter+post+'.fits'
  lssWindowFile    = maskDir + $
    'post_catalog.'+sdss_lss_vers +letter+post+'_inwindow.fit'
  lssIndexFile     = sdss_lss_dir + sample + '/'+$
    'lss_index.'+sample+'.fits'

  vagcFile         = vagcDir + 'object_sdss_spectro.fits'
  imFile           = vagcDir + 'object_sdss_imaging.fits'

  kcorrPetroFile   = vagcDir + 'kcorrect/kcorrect.photoz.petro.z0.10.fits'
  kcorrModelFile   = vagcDir + 'kcorrect/kcorrect.photoz.model.z0.10.fits'

;  GOTO, jump

  print
  print,'Reading LSS file: ',lssFile
  slss = mrdfits(lssFile,1)
  
  nlss = n_elements(slss)
  print,"Only Using "+ntostr(nlss)+" rows from the spec file"
  
  print
  print,'Reading lss Window File: ',lssWindowFile
  lssWindow = mrdfits(lssWindowFile,1)

  rows = slss.object_position

  spcolumns = ['mjd', 'plate', 'fiberid', 'objid', 'primtarget', $
               'sectarget', 'class', 'subclass', 'z', 'z_err', 'rchi2', $
               'dof', 'zwarning', 'vdisp', 'vdisp_err', $
               'vdispchi2', 'vdispdof']
  print
  print,'Reading vagc file: ',vagcfile
  sp = mrdfits(vagcfile, 1, columns=spcolumns, rows=rows)

  print,'Reading kcorr petro file: ',kcorrPetroFile
  kpetro = mrdfits(kcorrPetroFile,1,rows=rows)

  print
  print,'Reading LSS index file: ',lssIndexFile
  slss_index = mrdfits(lssIndexFile,1,columns=['sector','ilss'],rows=rows)


  print
  print,'Making completeness cut'
  ws = where(lssWindow.completeness GT 0, nws)
  print,'Kept '+ntostr(nws)+'/'+ntostr(n_elements(lssWindow))

  outStruct = create_struct('ra',            0d, $
                            'dec',           0d, $
                            'abs_petro_mag', fltarr(5), $
                            'z',             0.0, $
                            'completeness',  0.0, $
                            'plate',         0L,$
                            'sector',       -1L, $
                            'ilss',         -1L)

  outStruct = replicate(outStruct, nws)

  outStruct.ra = slss[ws].ra
  outStruct.dec = slss[ws].dec

  outStruct.abs_petro_mag[0] = kPetro[ws].absmag[0]
  outStruct.abs_petro_mag[1] = kPetro[ws].absmag[1]
  outStruct.abs_petro_mag[2] = kPetro[ws].absmag[2]
  outStruct.abs_petro_mag[3] = kPetro[ws].absmag[3]
  outStruct.abs_petro_mag[4] = kPetro[ws].absmag[4]

  outStruct.z = slss[ws].z
  outStruct.completeness = lssWindow[ws].completeness
  outStruct.plate = sp[ws].plate

  outStruct.sector = slss_index[ws].sector
  outStruct.ilss = slss_index[ws].ilss

  print
  print,'Writing outFile = ',outFile
  write_idlstruct, outStruct, outFile, /ascii

  delvarx, outStruct, sp, kPetro, slss, lssWindow, slss_index



jump:
  ;; Now the randoms
  nrand = 20
  randNums = 20 + lindgen(nrand)
  lssRandDir = sdss_lss_dir + sdss_lss_vers + "/safe/random/"
  columns=['ra','dec','fgot','sector']
  st = create_struct('ra', 0d, $
                     'dec', 0d, $
                     'completeness', 0.0,$
                     'sector', -1L)



  FOR i=0L, nrand-1 DO BEGIN 

      randNum = randNums[i]

      randFile = $
        lssRandDir + "random-"+ntostr(randnum)+".sample14safe.fits"
      outFile = outDir + 'jackiecat_all_rand'+strn(i,padchar='0',len=2)+'.st'

      print,'--------------------------------------------------------'
      print,'Reading random file: ',randFile
      r=mrdfits(randFile,1,columns=columns)

      print,'Making completeness cut'
      w=where(r.fgot GT 0, nkeep)
      print,'Keeping '+ntostr(nkeep)
      outst = replicate(st, nkeep)

      outst.ra = r[w].ra
      outst.dec = r[w].dec
      outst.completeness = r[w].fgot
      outst.sector = r[w].sector

      print
      print,'Writing file: ',outFile
      write_idlstruct, outst, outFile, /ascii
  ENDFOR 

  return
END 
