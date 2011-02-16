PRO rdwtheta_conv, clr, type, datax, data, dataerr, modelx, model, sigma, cutoff, phil=phil

  IF n_params() LT 2 THEN BEGIN
      print,'-Syntax: rdwtheta_conv, clr, type, datax, data, dataerr, modelx, model, sigma, cutoff'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; params
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  dataindir = '/sdss4/data1/esheldon/GAL_GAL/DENSITY/'
  indir = '/sdss4/data1/esheldon/TMP/conv/'
  colors=['u','g','r','i','z']
  datatypestr = ['dense_corr', 'low_corr', 'tot']
  typestr = ['dense', 'low', 'tot']

  zs = .4
  zl = .15

  philzl = [0., 0.168, 0.172, 0.173, 0.]

  phil_sigcrit = 1./[1., 0.343, 0.392, 0.403, 1.]
  my_sigcrit = sigmacrit(zs, zl, h=1.)

  philDl = angdist(philzl, h=1.)
  myDl = angdist(zl, h=1.)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Set up parameter arrays (cutoff was decided in wtheta_conv.pro)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nsig = 49L
  ncut = 35L

  sigma = fltarr(nsig)
  cutoff = fltarr(ncut)

  cutmin = 50.
  cutstep = 25.
  FOR i=0, ncut-1 DO cutoff[i] = cutmin + i*cutstep

  sigmin = 100.
  sigstep = 2.5
  FOR i=0L, nsig-1 DO sigma[i] = sigmin + i*sigstep

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; read in shear data (has already been multiplied by Ssh)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  IF keyword_set(phil) AND (type EQ 2) THEN BEGIN 
      print,'Using Phils Tot'
      datafile = dataindir + 'shearphil_'+datatypestr[type]+'_'+colors[clr]+'.txt'
  ENDIF ELSE BEGIN
      datafile = dataindir + 'shear_'+datatypestr[type]+'_'+colors[clr]+'.txt'
  ENDELSE 
  print,'Reading in data: ',datafile
  readcol, datafile, datax, data, dataerr, /silent

  print,'Reading in model data'
  first=1
  
  fixfac = my_sigcrit*myDl/phil_sigcrit[clr]/philDl[clr]

  FOR icut=0L, ncut-1 DO BEGIN
      
      kappa=0 & kappa_int=0 & shear_sig1=0
      tmpkap=0 & tmpkap_int=0 & single_lens=0

      cutstr = ntostr(long(cutoff[icut]))
      ;; Measured wtheta from r-band  Its most complete!!!
      file=indir+'conv_'+typestr[type]+'_r_cut'+cutstr+'.txt'
      
      readcol, file, modelx, kappa, kappa_int, /silent
      shear_sig1 = (kappa_int - kappa)*fixfac

      IF first THEN BEGIN
          first=0
          nx = n_elements(modelx)
          model = fltarr( ncut, nsig, nx )
      ENDIF 

      single_lens = shearsis_trunc(1., cutoff[icut], zs, zl, modelx)
;      single_lens = shearsis_trunc(1., cutoff[icut], zs, zl, modelx, /core)
      single_lens = single_lens*fixfac

      FOR isig=0L, nsig-1 DO BEGIN

          sig = sigma[isig]
          model[icut, isig, *] = (single_lens + shear_sig1)*sig^2

      ENDFOR 

  ENDFOR 

  

return
END 
