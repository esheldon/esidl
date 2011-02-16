PRO spectra_define_cuts, l

  eclasscut = -0.02
  gmrcut = 0.7
  gmrmin = 0.1
  gmrmax = 1.1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; define the absolute magnitude cuts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  outfile = '~/lensout/cuts.fit'

  maxmagarr = [-15.0, -16.5, -17.0, -17.0, -17.0]
  minmagarr = [-22.0, -23.5, -24.0, -24.0, -24.0]

  ;; luminosity cuts
  print,'*************************************'
  print,'* Main by luminosity'
  print,'*************************************'

  main_minmag_threebin = fltarr(5, 3)
  main_maxmag_threebin = main_minmag_threebin
  FOR clr = 0, 4 DO BEGIN 
      cstr = !colors[clr]
      print
      print,'********** '+cstr+'-band **********'

      perc = [0.85,0.10,0.05]
      w=where(l.z GT 0.02 AND l.z LT 0.3 AND $
              l.absmag[clr] GT minmagarr[clr] AND l.absmag[clr] lt maxmagarr[clr])

      binbyperc,l[w].absmag[clr],perc,/reverse, $
                name='absmag['+cstr+']', $
                minarr = minmag, maxarr=maxmag
      print
      binbyperc,l[w].lum[clr],perc,$
                name='lum['+cstr+'] [10^10 Lsun]'
      print
      binbyperc,l[w].lum[clr]/!lstar[clr],perc,$
                name='lum['+cstr+']/L*['+cstr+']'

      main_minmag_threebin[clr, *] = minmag
      main_maxmag_threebin[clr, *] = maxmag

      ;; put in hard cuts (for printing later)
      main_minmag_threebin[clr, 2] = minmagarr[clr]
      main_maxmag_threebin[clr, 0] = maxmagarr[clr]

  ENDFOR 
  
  print,'Hit a key' & key=get_kbrd(1)

  ;; binning early types by luminosity
  print,'*************************************'
  print,'* Early type by luminosity'
  print,'*************************************'

  early_minmag_threebin = fltarr(5,3)
  early_maxmag_threebin = early_minmag_threebin
  early_minmag_twobin = fltarr(5,2)
  early_maxmag_twobin = early_minmag_twobin
  FOR clr = 0, 4 DO BEGIN 
      cstr = !colors[clr]
      print
      print,'********** '+cstr+'-band **********'

      perc = [0.85,0.10,0.05]
      w=where(l.eclass LT eclasscut AND $
              l.z GT 0.02 AND l.z LT 0.3 AND $
              l.absmag[clr] GT minmagarr[clr] AND l.absmag[clr] lt maxmagarr[clr])

      binbyperc,l[w].absmag[clr],perc,/reverse, $
                name='absmag['+cstr+']', $
                minarr = minmag, maxarr=maxmag
      print
      binbyperc,l[w].lum[clr],perc,$
                name='lum['+cstr+'] [10^10 Lsun]'
      print
      binbyperc,l[w].lum[clr]/!lstar[clr],perc,$
                name='lum['+cstr+']/L*['+cstr+']'
     
      early_minmag_threebin[clr, *] = minmag
      early_maxmag_threebin[clr, *] = maxmag

      ;; put in hard cuts (for printing later)
      early_minmag_threebin[clr, 2] = minmagarr[clr]
      early_maxmag_threebin[clr, 0] = maxmagarr[clr]

      early_minmag_twobin[clr,0] = early_minmag_threebin[clr,0]
      early_maxmag_twobin[clr,0] = early_maxmag_threebin[clr,0]
      ;; concatenate last two bins
      early_minmag_twobin[clr,1] = early_minmag_threebin[clr,2]
      early_maxmag_twobin[clr,1] = early_maxmag_threebin[clr,1]

  ENDFOR 

  print,'Hit a key' & key=get_kbrd(1)

  ;; binning red types by luminosity
  print,'*************************************'
  print,'*  Red by luminosity'
  print,'*************************************'

  red_minmag_threebin = fltarr(5,3)
  red_maxmag_threebin = red_minmag_threebin
  red_minmag_twobin = fltarr(5,2)
  red_maxmag_twobin = red_minmag_twobin

  gmr = l.abscounts[1]-l.abscounts[2]

  FOR clr = 0, 4 DO BEGIN 
      cstr = !colors[clr]
      print
      print,'********** '+cstr+'-band **********'

      perc = [0.85,0.10,0.05]
      
      w=where(gmr GE gmrcut AND $
              gmr LT gmrmax AND $
              gmr GT gmrmin AND $
              l.z GT 0.02 AND l.z LT 0.3 AND $
              l.absmag[clr] GT minmagarr[clr] AND $
              l.absmag[clr] lt maxmagarr[clr])

      binbyperc,l[w].absmag[clr],perc,/reverse, $
                name='absmag['+cstr+']', $
                minarr = minmag, maxarr=maxmag
      print
      binbyperc,l[w].lum[clr],perc,$
                name='lum['+cstr+'] [10^10 Lsun]'
      print
      binbyperc,l[w].lum[clr]/!lstar[clr],perc,$
                name='lum['+cstr+']/L*['+cstr+']'
     
      red_minmag_threebin[clr, *] = minmag
      red_maxmag_threebin[clr, *] = maxmag

      ;; put in hard cuts (for printing later)
      red_minmag_threebin[clr, 2] = minmagarr[clr]
      red_maxmag_threebin[clr, 0] = maxmagarr[clr]

      red_minmag_twobin[clr,0] = red_minmag_threebin[clr,0]
      red_maxmag_twobin[clr,0] = red_maxmag_threebin[clr,0]
      ;; concatenate last two bins
      red_minmag_twobin[clr,1] = red_minmag_threebin[clr,2]
      red_maxmag_twobin[clr,1] = red_maxmag_threebin[clr,1]

  ENDFOR 


  range_struct = create_struct('main_minmag_threebin', main_minmag_threebin, $
                               'main_maxmag_threebin', main_maxmag_threebin, $
                               'early_minmag_threebin', early_minmag_threebin, $
                               'early_maxmag_threebin', early_maxmag_threebin, $
                               'early_minmag_twobin', early_minmag_twobin, $
                               'early_maxmag_twobin', early_maxmag_twobin, $
                               'red_minmag_threebin', red_minmag_threebin, $
                               'red_maxmag_threebin', red_maxmag_threebin, $
                               'red_minmag_twobin', red_minmag_twobin, $
                               'red_maxmag_twobin', red_maxmag_twobin)

  print,'Output file: ',outfile
  mwrfits, range_struct, outfile, /create

END 
