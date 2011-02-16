PRO huan_hdfn_convert_cat

  indir = '~/Huan/hawaii_hdfn/cat/'
  outfile = indir + 'all_RZ.fit'

  imagedir = '~/Huan/hawaii_hdfn/images/'
  Irmsfile = imagedir + 'HDF.I.rms.gz'
  Zrmsfile = imagedir + 'HDF.Z.rms.gz'

  ;; There are a bunch of files that need to be collated
  flag_file = indir + 'flag.RZ.cat'
  flux_file = indir + 'flux.RZ.cat'
  mag_file  = indir + 'mag.RZ.cat'
  shape_file = indir + 'shape.RZ.cat'

  ;; Determine number of objects
  nl = numlines(flag_file)
  nobj = nl-1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The flag file. 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  flag_ss = create_struct('id', 0L, $
                          'bad', 0b, $
                          'usat', 0b, $
                          'bsat', 0b, $
                          'vsat', 0b, $
                          'rsat', 0b, $
                          'isat', 0b, $
                          'zsat', 0b, $
                          'hksat', 0b, $
                          'n_neighbor', 0)


  flag_struct = replicate(flag_ss, nobj)

  print
  print,'Reading flag file: ',flag_file
  openr, lun, flag_file, /get_lun
  ;; must skip the first line
  line = ""
  readf, lun, line
  print,line

  ;; now read whole thing
  readf, lun, flag_struct

  free_lun, lun

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;  The flux file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  flux_ss = create_struct('id', 0L, $
                          'Uaperflux', 0.0, $
                          'dUaperflux', 0.0, $
                          'Ubkg', 0.0, $
                          'Uisoflux', 0.0, $
                          'dUisoflux', 0.0, $
                          'Uisobkg',0.0, $
                          $
                          'Baperflux', 0.0, $
                          'dBaperflux', 0.0, $
                          'Bbkg', 0.0, $
                          'Bisoflux', 0.0, $
                          'dBisoflux', 0.0, $
                          'Bisobkg',0.0, $
                          $
                          'Vaperflux', 0.0, $
                          'dVaperflux', 0.0, $
                          'Vbkg', 0.0, $
                          'Visoflux', 0.0, $
                          'dVisoflux', 0.0, $
                          'Visobkg',0.0, $
                          $
                          'Raperflux', 0.0, $
                          'dRaperflux', 0.0, $
                          'Rbkg', 0.0, $
                          'Risoflux', 0.0, $
                          'dRisoflux', 0.0, $
                          'Risobkg',0.0, $
                          $
                          'Iaperflux', 0.0, $
                          'dIaperflux', 0.0, $
                          'Ibkg', 0.0, $
                          'Iisoflux', 0.0, $
                          'dIisoflux', 0.0, $
                          'Iisobkg',0.0, $
                          $
                          'Zaperflux', 0.0, $
                          'dZaperflux', 0.0, $
                          'Zbkg', 0.0, $
                          'Zisoflux', 0.0, $
                          'dZisoflux', 0.0, $
                          'Zisobkg',0.0, $
                          $
                          'HKaperflux', 0.0, $
                          'dHKaperflux', 0.0, $
                          'HKbkg', 0.0, $
                          'HKisoflux', 0.0, $
                          'dHKisoflux', 0.0, $
                          'HKisobkg',0.0)
  
  flux_struct = replicate(flux_ss, nobj)

  print
  print,'Reading flux file: ',flux_file
  openr, lun, flux_file, /get_lun
  ;; must skip the first line
  line = ""
  readf, lun, line
  print,line

  readf, lun, flux_struct

  free_lun, lun

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The mag cat
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  mag_ss = create_struct('id', 0L, $
                         'U', 0.0, $
                         'dU', 0.0, $
                         'Uiso', 0.0, $
                         'dUiso', 0.0, $
                         $
                         'B', 0.0, $
                         'dB', 0.0, $
                         'Biso', 0.0, $
                         'dBiso', 0.0, $
                         $
                         'V', 0.0, $
                         'dV', 0.0, $
                         'Viso', 0.0, $
                         'dViso', 0.0, $
                         $
                         'R', 0.0, $
                         'dR', 0.0, $
                         'Riso', 0.0, $
                         'dRiso', 0.0, $
                         $
                         'I', 0.0, $
                         'dI', 0.0, $
                         'Iiso', 0.0, $
                         'dIiso', 0.0, $
                         $
                         'Z', 0.0, $
                         'dZ', 0.0, $
                         'Ziso', 0.0, $
                         'dZiso', 0.0, $
                         $
                         'HK', 0.0, $
                         'dHK', 0.0, $
                         'HKiso', 0.0, $
                         'dHKiso', 0.0)

  mag_struct = replicate(mag_ss, nobj)

  print
  print,'Reading mag file: ',mag_file
  openr, lun, mag_file, /get_lun
  ;; must skip the first line
  line = ""
  readf, lun, line
  print,line

  readf, lun, mag_struct

  free_lun, lun
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; the "shape" cat
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  shape_ss = create_struct('id', 0L, $
                           'ra', 0d, $
                           'dec', 0d, $
                           'x_image', 0.0, $
                           'y_image', 0.0, $
                           $
                           'Ufwhm', 0.0, $
                           'Ue0', 0.0, $
                           'Ue1', 0.0, $
                           $
                           'Bfwhm', 0.0, $
                           'Be0', 0.0, $
                           'Be1', 0.0, $
                           $
                           'Vfwhm', 0.0, $
                           'Ve0', 0.0, $
                           'Ve1', 0.0, $
                           $
                           'Rfwhm', 0.0, $
                           'Re0', 0.0, $
                           'Re1', 0.0, $
                           $
                           'Ifwhm', 0.0, $
                           'Ie0', 0.0, $
                           'Ie1', 0.0, $
                           $
                           'Zfwhm', 0.0, $
                           'Ze0', 0.0, $
                           'Ze1', 0.0, $
                           $
                           'HKfwhm', 0.0, $
                           'HKe0', 0.0, $
                           'HKe1', 0.0)

  shape_struct = replicate(shape_ss, nobj)

  print
  print,'Reading shape file: ',shape_file
  openr, lun, shape_file, /get_lun
  ;; must skip the first line
  line = ""
  readf, lun, line
  print,line

  readf, lun, shape_struct

  free_lun, lun

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Now the combined catalog
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  tot_ss = create_struct('id', 0L, $ ; flag struct
                         'bad', 0b, $
                         'usat', 0b, $
                         'bsat', 0b, $
                         'vsat', 0b, $
                         'rsat', 0b, $
                         'isat', 0b, $
                         'zsat', 0b, $
                         'hksat', 0b, $
                         'n_neighbor', 0, $
                         $
                         $
                         'Uaperflux', 0.0, $ ; flux struct
                         'dUaperflux', 0.0, $
                         'Ubkg', 0.0, $
                         'Uisoflux', 0.0, $
                         'dUisoflux', 0.0, $
                         'Uisobkg',0.0, $
                         $
                         'Baperflux', 0.0, $
                         'dBaperflux', 0.0, $
                         'Bbkg', 0.0, $
                         'Bisoflux', 0.0, $
                         'dBisoflux', 0.0, $
                         'Bisobkg',0.0, $
                         $
                         'Vaperflux', 0.0, $
                         'dVaperflux', 0.0, $
                         'Vbkg', 0.0, $
                         'Visoflux', 0.0, $
                         'dVisoflux', 0.0, $
                         'Visobkg',0.0, $
                         $
                         'Raperflux', 0.0, $
                         'dRaperflux', 0.0, $
                         'Rbkg', 0.0, $
                         'Risoflux', 0.0, $
                         'dRisoflux', 0.0, $
                         'Risobkg',0.0, $
                         $
                         'Iaperflux', 0.0, $
                         'dIaperflux', 0.0, $
                         'Ibkg', 0.0, $
                         'Iisoflux', 0.0, $
                         'dIisoflux', 0.0, $
                         'Iisobkg',0.0, $
                         $
                         'Zaperflux', 0.0, $
                         'dZaperflux', 0.0, $
                         'Zbkg', 0.0, $
                         'Zisoflux', 0.0, $
                         'dZisoflux', 0.0, $
                         'Zisobkg',0.0, $
                         $
                         'HKaperflux', 0.0, $
                         'dHKaperflux', 0.0, $
                         'HKbkg', 0.0, $
                         'HKisoflux', 0.0, $
                         'dHKisoflux', 0.0, $
                         'HKisobkg',0.0, $
                         $
                         $
                         'U', 0.0, $ ; mag struct
                         'dU', 0.0, $
                         'Uiso', 0.0, $
                         'dUiso', 0.0, $
                         $
                         'B', 0.0, $
                         'dB', 0.0, $
                         'Biso', 0.0, $
                         'dBiso', 0.0, $
                         $
                         'V', 0.0, $
                         'dV', 0.0, $
                         'Viso', 0.0, $
                         'dViso', 0.0, $
                         $
                         'R', 0.0, $
                         'dR', 0.0, $
                         'Riso', 0.0, $
                         'dRiso', 0.0, $
                         $
                         'I', 0.0, $
                         'dI', 0.0, $
                         'Iiso', 0.0, $
                         'dIiso', 0.0, $
                         $
                         'Z', 0.0, $
                         'dZ', 0.0, $
                         'Ziso', 0.0, $
                         'dZiso', 0.0, $
                         $
                         'HK', 0.0, $
                         'dHK', 0.0, $
                         'HKiso', 0.0, $
                         'dHKiso', 0.0, $
                         $
                         $
                         'ra', 0d, $ ; shape struct
                         'dec', 0d, $
                         'x_image', 0.0, $
                         'y_image', 0.0, $
                         $
                         'Ufwhm', 0.0, $
                         'Ue0', 0.0, $
                         'Ue1', 0.0, $
                         $
                         'Bfwhm', 0.0, $
                         'Be0', 0.0, $
                         'Be1', 0.0, $
                         $
                         'Vfwhm', 0.0, $
                         'Ve0', 0.0, $
                         'Ve1', 0.0, $
                         $
                         'Rfwhm', 0.0, $
                         'Re0', 0.0, $
                         'Re1', 0.0, $
                         $
                         'Ifwhm', 0.0, $
                         'Ie0', 0.0, $
                         'Ie1', 0.0, $
                         $
                         'Zfwhm', 0.0, $
                         'Ze0', 0.0, $
                         'Ze1', 0.0, $
                         $
                         'HKfwhm', 0.0, $
                         'HKe0', 0.0, $
                         'HKe1', 0.0, $
                         $ 
                         $
                         'I_rms_sky', 0.0, $
                         'Z_rms_sky', 0.0)
  
  struct = replicate(tot_ss, nobj)

  ;; copy in the various stuff
  struct_assign, flag_struct, struct, /verbose, /nozero
  struct_assign, flux_struct, struct, /verbose, /nozero
  struct_assign, mag_struct, struct, /verbose, /nozero
  struct_assign, shape_struct, struct, /verbose, /nozero

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Fill in the sky RMS
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Put the positions in zero-offset
  x = struct.x_image - 1.0
  y = struct.y_image - 1.0

  print
  print,'Reading I-band rms sky file: ',Irmsfile
  Irms = mrdfits(Irmsfile)

  I_rms_sky = Irms[ x, y ]

  struct.I_rms_sky = I_rms_sky
  delvarx, Irms

  print
  print,'Reading Z-band rms sky file: ',Zrmsfile
  Zrms = mrdfits(Zrmsfile)

  Z_rms_sky = Zrms[ x, y ]

  struct.Z_rms_sky = Z_rms_sky
  delvarx, Zrms

  print
  print,'Writing combined file: ',outfile
  mwrfits2, struct, outfile, /create, /destroy

END 
