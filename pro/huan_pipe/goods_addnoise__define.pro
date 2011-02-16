FUNCTION goods_addnoise::init, type, seeing=seeing

  IF n_elements(type) EQ 0 THEN BEGIN 
      print,'You must intitialize the type'
      print,"-Syntax: ga = obj_new('goods_addnoise', type)"
      print,'type = desNyr, lsst'
      return,0
  ENDIF 
  
  IF NOT self->valid_type(type) THEN return,0

  self.type = type

  self.image = ptr_new(/alloc)
  self.cat = ptr_new(/alloc)

  self.image_read = 0
  self.cat_read = 0
  
  IF n_elements(seeing) NE 0 THEN self.seeing = seeing ELSE self.seeing = -1

  self.image_convolved = 0
  self.image_noise_added = 0

  return,1

END 

PRO goods_addnoise::set_parameters, type=type, seeing=seeing

  IF n_elements(type) NE 0 THEN BEGIN 
      IF NOT self->valid_type(type) THEN BEGIN 
          message,'Type not set',/inf
      ENDIF ELSE BEGIN 
          self.type = type
      ENDELSE 
  ENDIF  
  IF n_elements(seeing) NE 0 THEN self.seeing = seeing

END 

FUNCTION goods_addnoise::valid_type,type

  CASE type OF 
      'cfh':    return,  1
      'lsst':   return,  1
      'des1yr': return,  1
      'des2yr': return,  1
      'des3yr': return,  1
      'des4yr': return, 1
      'des5yr': return, 1
      ELSE: BEGIN 
          message,'Invalid type: '+type,/inf
          return,0
      END
  ENDCASE 

END 

PRO goods_addnoise::read_image, file

  IF self.image_read THEN message,'An image is already in memory. Free first.'
  print,'Reading file: ',file
  image = mrdfits(file)
  self.image = ptr_new(image, /no_copy)
  self.image_read = 1

END 

PRO goods_addnoise::read_cat, file

  IF self.cat_read THEN message,'A catalog is already in memory. Free first.'
  print,'Reading file: ',file
  cat = mrdfits(file, 1)
  self.cat = ptr_new(cat, /no_copy)
  self.cat_read = 1
END 

PRO goods_addnoise::write_image, file

  im = temporary( *self.image )
  print,'Writing new image: ',file
  mwrfits, im, file, /create

  self.image = ptr_new(im, /no_copy)

END 

PRO goods_addnoise::write_cat, file

  cat = temporary( *self.cat )
  print,'Writing new catalog: ',file
  mwrfits, cat, file, /create

  self.cat = ptr_new(cat, /no_copy)

END 

FUNCTION goods_addnoise::pixsize, type

  CASE type OF
      'goods': return, 0.030
      'cfh':   return, 0.206
      'lsst':  return, 0.206
      'des':   return, 0.270
      ELSE: message,'Unknown type: '+type
  ENDCASE 

END 

FUNCTION goods_addnoise::exptime, type

  ;; i-band

  CASE type OF 
      'cfh':    return,  900.   ; seconds
      'lsst':   return,  900.   ; ERROR: This might be wrong? Confused..
      'des1yr': return,  200.
      'des2yr': return,  400.   ; not 5 year * 2/5
      'des3yr': return,  800.
      'des4yr': return, 1200.
      'des5yr': return, 1200.   ; 24.3 in i
      ELSE: message,'Unknown type: '+type
  ENDCASE 

END 

FUNCTION goods_addnoise::mlim, type

  ;; i-band
  CASE type OF 
      'cfh':    return, 24.1
      'lsst':   return, 25.3
      'des1yr': return, 23.3
      'des2yr': return, 23.6
      'des3yr': return, 24.0
      'des4yr': return, 24.3
      'des5yr': return, 24.3
      ELSE: message,'Unknown type: '+type
  ENDCASE 

END 

FUNCTION goods_addnoise::strip_type, type

  IF strmatch(type, 'des*') THEN ctype='des' ELSE ctype = type
  return,ctype

END 

FUNCTION goods_addnoise::noise, type, raw=raw

  stype = self->strip_type(type)

  CASE stype OF 
      'cfh': BEGIN 
          cfh_gain = 1.5
          cfh_noise = 70.0      ; ADU
          cfh_noise = cfh_noise*sqrt(cfh_gain)

          ;; Can return raw noise in cfh image
          IF keyword_set(raw) THEN return,cfh_noise

          ;; Goods is different pixel size
          fac = self->pixsize('goods')/self->pixsize('cfh')
          cfh_noise = fac*cfh_noise
          return,cfh_noise
      END 
      'des': BEGIN 

          ;; Account for different pixel size not square because its 
          ;; noise: square root

          pixfac = self->pixsize('goods')/self->pixsize('des')

          ;; Account for different exposure times. This needs full type
          ;; e.g. des5yr
          magdiff = self->mlim(type) - self->mlim('cfh')

          cfh_noise = self->noise('cfh',/raw)
          des_noise = pixfac*cfh_noise*sqrt(10.^(-magdiff/2.5))

          return, des_noise
      END 
      'lsst': BEGIN 

          ;; Account for different pixel size not square because its 
          ;; noise: square root

          pixfac = self->pixsize('goods')/self->pixsize('lsst')

          ;; Account for different exposure times
          magdiff = self->mlim('lsst') - self->mlim('cfh')

          cfh_noise = self->noise('cfh',/raw)
          des_noise = pixfac*cfh_noise*sqrt(10.^(-magdiff/2.5))

          return, des_noise
      END 
      ELSE: message,'Unknown type: '+type
  ENDCASE 

END 

PRO goods_addnoise::convolve

  IF self.image_convolved THEN BEGIN
      message,'Seeing was already added'
  ENDIF 
  IF self.seeing EQ -1 THEN BEGIN 
      message,'The seeing was never initialized'
  ENDIF 

  ;; Get normal representation of image
  im = temporary( *self.image )

  ;; Generate PSF image
  print,'Creating PSF image for seeing = '+ntostr(self.seeing)
  arcperpix = self->pixsize('goods')
  seeing_pixels = self.seeing/arcperpix ; fwhm pixels

  cfac = 2.*sqrt(2.*alog(2))  ;; fwhm = cfac*sigma for gaussian

  sigma = seeing_pixels/cfac ;; psf scale length in arcseconds.   
  psfactor = 2*3.0              ; make KERNEL 3 sigma in RADIUS => 6 sigma wide

  size = round(sigma*psfactor)   ;; psf KERNEL size

  makegauss, psf, [size,size], sigma, counts=1
  help,psf

  tt=systime(1)

  ;; Uses fft
  print,'Convolving with PSF'
  im = convolve( temporary(im), psf)

  ptime,systime(1)-tt

  ;; Put image back
  self.image = ptr_new( im, /no_copy )

  self.image_convolved = 1

END 

PRO goods_addnoise::addnoise

  IF self.seeing NE -1 AND NOT self.image_convolved THEN BEGIN 
      message,'You need to convolve the image first'
  ENDIF 

  IF self.image_noise_added THEN BEGIN
      message,'Noise was already added'
  ENDIF 

  ;; Get standard representation
  im = temporary( *self.image )
  cat = temporary( *self.cat )

  type = self.type

  noise = self->noise(type)
  exp_time = self->exptime(type)

  print
  print,'Adding noise to catalog error estimates'
  print,'Noise = '+ntostr(noise)

  cat.rms_sky = sqrt(cat.rms_sky^2*exp_time + noise^2)
  cat.background = cat.background*exp_time

  ;; Add noise to the image
  n_im = n_elements(im)
  print,'Adding noise to image'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Here we are basically ignoring the noise in the
  ;; original image, since multiplying isn't the right
  ;; thing to do
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  im[*] = im[*]*exp_time

  nadd = noise*randomu(seed, n_im, /normal)
  im[*] = im[*] + nadd[*]
  nadd=0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Inverse variance map
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  print,'Updating ivar map'
;  w=where(ivar GT 0, nw)
;  ivar[w] = 1./(1./ivar[w]*exp_time + noise^2)

  ;; Put back 
  self.image = ptr_new( im, /no_copy )
  self.cat = ptr_new( cat, /no_copy )


END 



PRO goods_addnoise::free
  ptr_free, self.image, self.cat
  self.image_read = 0
  self.cat_read = 0

  self.image_convolved = 0
  self.image_noise_added = 0

  return
END 


FUNCTION goods_addnoise::cleanup
  ptr_free, self.image, self.cat
  return,1
END 

PRO goods_addnoise__define

  struct = { $
             goods_addnoise, $
             type:'', $
             seeing: 0.0, $
             image: ptr_new(), $
             image_read: 0, $
             image_convolved: 0, $
             image_noise_added: 0, $
             cat: ptr_new(), $
             cat_read: 0 $
           }

END 
