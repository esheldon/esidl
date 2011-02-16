PRO vagc_setuplens, rmin, rmax, nbin_OR_binsize, stripes, $
                    compcut=compcut, logbin=logbin, rlrgMask=rlrgMask, $
                    maskFile=maskFile, soFile=soFile

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: vagc_setuplens, $'
      print,'    rmin, rmax, nbin_OR_binsize, stripes, $'
      print,'    compcut=compcut, /logbin, /rlrgMask'
      return
  ENDIF 

  IF n_elements(compcut) EQ 0 THEN compcut = 0.0

  omegamat=0.3
  hubble=1.0
  max_allow_angle = 6.0         ;degrees
  maxz = 0.8

  print,'Omega_matter      = ',omegamat
  print,'Using h           = ',hubble
  print,'Max allowed angle = ',max_allow_angle
  print,'Max z             = ',maxz
  print,'Compcut           = ',compcut
  print,'log binning       = ',keyword_set(logbin)

  ;;;;;;;;;;;;;;;;;;;;;;
  ;; log binning?
  ;;;;;;;;;;;;;;;;;;;;;
  
  IF keyword_set(logbin) THEN BEGIN 
      nbin = nbin_or_binsize
      print,'Number of bins    = ',nbin
  ENDIF ELSE BEGIN 
      binsize = nbin_or_binsize
      nbin = long( (rmax - rmin)/binsize ) ;+ 1 
      print,'Binsize           = ',binsize
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; output file name
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  vagc_lensinput_name, stripes, rmin, rmax, outfile, rlrgMask=rlrgMask
  print
  print,'Will write to file: ',outfile

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; read in basic catalog
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  get_spectra_lcat, stripes, lcat, /lss, count=nlens_init

  print,'Creating lensum struct'
  lensumstruct = zlensumstruct(fltarr(nbin))
  lensumstruct = create_struct(lcat[0], $
                               lensumstruct)

  lensum = replicate(lensumstruct, nlens_init)

  print,'Copying....'
  copy_struct, lcat, lensum
  delvarx, lcat

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make the cuts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  vagc_make_lenscuts, lensum, $
    rmax, stripes, max_allow_angle, $
    omegamat, hubble, maxz, compcut, $
    angMax, DL, wlens, rlrgMask=rlrgMask, $
    maskFile=maskFile, soFile=soFile

  lensum = lensum[wlens]
  nLens = n_elements(lensum)

  lensum.index = lindgen(nLens)
  lensum.zindex = wlens
  lensum.angMax = angMax[wlens]
  lensum.DL = DL[wlens]

  print,'Finally using '+ntostr(nlens)+'/'+ntostr(nlens_init)

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; THe header
  ;;;;;;;;;;;;;;;;;;;;;;;;

  hdr = ['END']
  sxAddPar, hdr, 'nLenses', nlens,    ' Number of lenses after cuts'
  sxAddPar, hdr, 'omegaM',  omegamat, ' Matter Density (Flat Universe)'
  sxAddPar, hdr, 'h',       hubble,   ' (H/100 km/s)'
  sxAddPar, hdr, 'max_ang', max_allow_angle, ' Maximum allowed angle'
  sxAddPar, hdr, 'rmin',    rmin,     ' Minimum radius (kpc)'
  sxAddPar, hdr, 'rmax',    rmax,     ' Maximum radius (kpc)'
  sxAddPar, hdr, 'maxZ',    maxz,     ' Maximum redshift'
  sxAddPar, hdr, 'compcut', compcut,  ' Completeness cut'
  sxAddPar, hdr, 'logbin',  keyword_set(logbin),' Logarithmic binning?'
  sxAddPar, hdr, 'bin',     nbin_OR_binsize, $
    ' nbin if logbin, else binsize (kpc)'

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Write the file
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Writing file: ',outfile
  mwrfits2, lensum, outfile, hdr, /create, /destroy

END 
