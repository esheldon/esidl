PRO run_huan_hdfn_admom

  dir = '~/Huan/hawaii_hdfn/'

  catlist = dir + ['cat/all_RZ.fit', 'cat/all_RZ.fit']
  imlist  = dir + ['images/HDF.I.fits.gz','images/HDF.Z.fits.gz']
  bands = ['I','Z']
  huan_hdfn_admom, imlist, catlist, bands

END 
