PRO huan_hdfn_call_admom, image, x, y, sky, rms_sky, $
                          ixx, iyy, ixy, rho4, err, whyflag

  IF n_params() EQ 0 THEN BEGIN 
      print,'-syntax '
      return
  ENDIF 

  COMMON huan_hdfn_admom_block, defval, defval2

  sz = size(image, /tname)

  IF sz NE 'FLOAT' AND sz NE 'DOUBLE' THEN message,'image must be single or double precisioni floating point'

  n = n_elements(x)

  ixx = replicate(defval, n)
  iyy = ixx
  ixy = ixx

  nsky = n_elements(sky)
  IF nsky EQ 1 THEN BEGIN 
      sk    = replicate(sky, n)
  ENDIF ELSE IF nsky EQ n THEN BEGIN 
      sk = sky
  ENDIF ELSE BEGIN 
      message,$
        'sky must be either length 1 or same length as object list'
  ENDELSE 

  nrms = n_elements(rms_sky)
  IF nrms EQ 1 THEN BEGIN 
      sigsk = replicate(rms_sky, n)
  ENDIF ELSE IF nrms EQ n THEN BEGIN 
      sigsk = rms_sky
  ENDIF ELSE BEGIN 
      message,$
        'rms_sky must be either length 1 or same length as object list'
  ENDELSE 


;  print,'Sky: '+ntostr(sky)+'   rms sky: '+ntostr(rms_sky)

  whyflag = lonarr(n)

  s    = size(image)
  nx   = s[1]
  ny   = s[2]
  err  = replicate(defval, n)
  rho4 = replicate(defval, n)
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; define the shared object file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  bindir=getenv('PRODUCTS_DIR')+'/admomatlas/bin/'
  sofile = bindir + 'admom.so'
  entry = 'admom'
  
  print
  print,'Using sofile:  ',sofile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Call the C routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Calling admom.c'
 
  test = call_external(sofile, entry, $
                       image, nx, ny, x, y, n, sk, sigsk, $
                       ixx, iyy, ixy, err, rho4, whyflag)

  wbad = where(whyflag NE 0, nbad)

  IF nbad NE 0 THEN BEGIN 

      ixx[wbad]  = defval
      iyy[wbad]  = defval
      ixy[wbad]  = defval

      err[wbad]  = defval
      rho4[wbad] = defval

  ENDIF 


;  print,'Calling ad_momi.f'
;  ff=call_external(sofile, keyword, x, y, ixx, iyy, ixy, n, sh, image, sk, nx, ny,$
;                   err, sigsk, mag, numiter, wcenx, wceny, whyflag, rho4)


return
end


PRO huan_hdfn_admom, imlist, catlist, bands

; NAME:
;       HUAN_ADMOM
; PURPOSE:
;	Run images and sextractor catalogs through the adaptive moment 
;       C code.  The only thing the catalog must have is x_image,y_image
;
; CALLING SEQUENCE:
;       huan_admom, imlist, catlist, bands
;
; INPUTS:
;       imlist: list of image names, full path
;       catlist: list of catalog names, full path.  The new catalog will go in
;                the same directory with _admom appended before the .fits
;       bands: string array, either 'I' or 'Z' for now
;-
; On_error,2              ;Return to caller
  
  IF N_params() LT 3 THEN BEGIN 
      print,'Syntax - huan_hdfn_admom, imlist, catlist, bands'
      return
  ENDIF 

  COMMON huan_hdfn_admom_block, defval, defval2

  ;; size of pixels
  arcperpix = 0.206


  defval = -9999.0
  defval2 = 0L

  nim = n_elements(imlist)
  ncat = n_elements(catlist)

  IF nim NE ncat THEN message,'# of images must equal number of catalogs'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop through the input files (images)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR i=0L, nim-1 DO BEGIN  

      tim = systime(1)

      namearray=str_sep(catlist[i],'.')
      outfile = namearray[0]+"_"+bands[i]+"image_admom.fit"

      print
      print, "Will process frame: ", imlist[i]
      print, "Will be writing objects to file: ", outfile
      print

      image = mrdfits(imlist[i],  0, imhdr)

      ;; read the catalog
      cat   = mrdfits(catlist[i], 1, cathdr)
      nobj = n_elements(cat)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Add adaptive moments to the catalog
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; scale image if needed
      bzero=fxpar(imhdr,'BZERO')
      bscale=fxpar(imhdr,'BSCALE')
      IF bscale NE 0 THEN BEGIN 
          print,'Scaling image'
          image = image*bscale+bzero
      ENDIF 
	
      print
      print,'Adding sky of 1000 to image because admom converts to U16 types'
      image = image + 1000.

      ;; try using sky of zero and skyrms from the file, individual for
      ;; each object.  Note only I-band has rms_sky in the catalog

;      sky = 1000.0
;      rms_sky = 102.587
      IF bands[i] EQ 'I' THEN BEGIN 
          sky = cat.IisoBkg + 1000.
          rms_sky = cat.I_rms_sky
      ENDIF ELSE IF bands[i] EQ 'Z' THEN BEGIN 
          sky = cat.ZisoBkg + 1000.
          rms_sky = cat.Z_rms_sky
      ENDIF ELSE message,'Band must be "I" or "Z"'

      ;; put x,y in IDL notation (not SExtractor, which sets first pixel to
      ;; [1,1])

      x = cat.x_image - 1.0
      y = cat.y_image - 1.0

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;;  Call adaptive moment wrapper
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      huan_hdfn_call_admom, $
        image, x, y, sky, rms_sky, $
        ixx, iyy, ixy, rho4, err, whyflag
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; define a new strucure that will contain adaptive moment info
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      str=cat[0]
      strf2 = create_struct(str,$
                            'e1',        defval,  $
                            'e2',        defval,  $
                            'ixx',       defval,  $
                            'iyy',       defval,  $
                            'ixy',       defval,  $
                            'rho4',      defval,  $
                            'a4',        defval,  $
                            'ellip_err', defval,  $
                            'psf_fwhm',  defval,  $
                            's2n',       defval,  $
                            'whyflag',   defval2)	

      outstruct = replicate(strf2, nobj)
      copy_struct, cat, outstruct

      e1 = replicate(defval, nobj) ;bad are set to defval
      e2 = e1

      ;; Calculate a4
      a4 = rho4/2.0 - 1.0

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Calculate ellipticities
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      t = ixx+iyy
      w = where(t gt .1, nw)
      IF nw NE 0 THEN begin
          e1[w] = (ixx[w]-iyy[w])/t[w]
          e2[w] = 2.0*ixy[w]/t[w]
      ENDIF

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Calculate signal to noise
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      outstruct.e1 = e1
      outstruct.e2 = e2
      outstruct.ixx = ixx
      outstruct.iyy = iyy
      outstruct.ixy = ixy
      outstruct.rho4 = rho4
      outstruct.a4   = a4
      outstruct.ellip_err = err	
      outstruct.whyflag = whyflag

      ;; Write fits file
      print,'Writing fits file: ',outfile
      mwrfits2,outstruct,outfile,/create,/destroy

      ptime,systime(1)-tim

  ENDFOR 

  return
END 








