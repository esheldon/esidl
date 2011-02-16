PRO huan_cfh_newseeing

  seeing = $
    [0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]
  seeing_str = $
    ['0.70', '0.75', '0.80', '0.85', '0.90', '0.95', '1.00', '1.05', '1.10', '1.15', '1.20']
  ns=n_elements(seeing)
  
  ;; read in the images and convolve with new seeing disk, output

  huan_dir = '~/Huan/'
  files = findfile(huan_dir+'images/*.fits')

  nf = n_elements(files)



  psfactor = 2*3.0              ; make KERNEL 3 sigma in RADIUS => 6 sigma wide
  arcperpix = 0.206
  cfac = 2.*sqrt(2.*alog(2))  ;; fwhm = cfac*sigma for gaussian

  ;; Loop over files
  FOR j=0L, nf-1 DO BEGIN 

      print,'Input File: '+files[j]
      print,'--------------------------------------------------------'
      image = mrdfits(files[j])

      namearray=str_sep(files[j],'.')
      
      ;; different seeing bins
      FOR i=0L, ns-1 DO BEGIN 

          newfile = namearray[0]+"_seeing"+seeing_str[i]+".fits"

          fwhm = seeing[i]/arcperpix ; pixels

          sigma = fwhm/cfac ;; psf scale length in arcseconds.      
          size = round(sigma*psfactor)   ;; psf KERNEL size

          print,'Seeing:      '+ntostr(seeing[i])
          print,'   FWHM pixels: '+ntostr(fwhm)
          print,'   PSF sigma:   '+ntostr(sigma)
          print,'   kernel size: ['+ntostr(size)+', '+ntostr(size)+']'

          ;; Create the kernel
          makegauss, psf, [size,size], sigma, counts=1

          ;; convolve
          print
          print,'   Convolving'
          new_image = convol(image, psf, /edge_truncate)

          print,'   Writing file: ',newfile
          mwrfits, new_image, newfile, /create

      ENDFOR 


  ENDFOR 

END 
