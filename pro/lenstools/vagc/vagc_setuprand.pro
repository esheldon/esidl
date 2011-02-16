PRO vagc_setuprand, rand_filenums, rmin, rmax, nbin_OR_binsize, stripes, $
                    compcut=compcut, $
                    logbin=logbin, rlrgMask=rlrgMask, $
                    maskFile=maskFile, soFile=soFile

  IF n_params() LT 5 THEN BEGIN 
      print,'-Syntax: vagc_setuprand, $'
      print,'          rand_filenums, rmin, rmax, nbin_OR_binsize, stripes, $'
      print,'          compcut=compcut, $'
      print,'          /logbin, /rlrgMask'
      return
  ENDIF 

  nrandf = n_elements(rand_filenums)

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
  ;;;;;;;;;;;;;;;;;;;;;;
  
  IF keyword_set(logbin) THEN BEGIN 
      nbin = nbin_or_binsize
      print,'Number of bins    = ',nbin
  ENDIF ELSE BEGIN 
      binsize = nbin_or_binsize
      nbin = long( (rmax - rmin)/binsize ) ;+ 1 
      print,'Binsize           = ',binsize
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; read in lenses
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  get_spectra_lcat, stripes, lcat, /lss, count=nlcatz, columns=['z']

  ;; Note, the randoms are in the "safe" directory!
  lssRandDir = $
    sdssidl_config('lss_dir') + sdssidl_config('lss_vers') + "/safe/random/"

  FOR i=0L, nrandf-1 DO BEGIN 
      randnum = rand_filenums[i]

      randFile = $
        lssRandDir + "random-"+ntostr(randnum)+".sample14safe.fits"
      vagc_lensinput_name, stripes, rmin, rmax, outFile, randnum=randnum, $
        rlrgMask=rlrgMask

      print
      print,'Input file: ',randFile
      print,'Output file: ',outFile
      rand = mrdfits(randfile, 1)

      print,'Creating lensum struct'
      lensumstruct = zlensumstruct(fltarr(nbin))
      lensumstruct = create_struct('sector', 0L, $
                                   'mregion', 0L, $
                                   'completeness', 0.0, $
                                   'mmax', 0.0, $
                                   'ra',0d,$
                                   'dec',0d,$
                                   'clambda', 0d, $
                                   'ceta', 0d, $
                                   'z', 0.0, $
                                   lensumstruct)

      nlens_init = n_elements(rand)
      lensum = replicate(lensumstruct, nlens_init)

      print,'Copying....'
      eq2csurvey, rand.ra, rand.dec, clambda, ceta
      
      lensum.sector = rand.sector
      lensum.mregion = rand.mregion
      lensum.completeness = rand.fgot
      lensum.mmax = rand.mmax
      lensum.clambda = clambda
      lensum.ceta = ceta
      lensum.ra = rand.ra
      lensum.dec = rand.dec
      delvarx, rand

      print
      print,'Assigning redshifts'

      lcat_index = $
        long( nlcatz*randomu(seed, nlens_init) )

      lensum.z = lcat[lcat_index].z
      lensum.zindex = lcat_index

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

  ENDFOR 

END 
