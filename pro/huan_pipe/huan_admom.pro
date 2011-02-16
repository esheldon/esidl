PRO huan_call_admom, image, cat, ixx, iyy, ixy, rho4, err, whyflag, $
                      sky=sky, rms_sky=rms_sky

; Possible changes:
;     initial guesses for moments.  Will require simple change in phadmom.
;     fix skysig estimate

;this is the IDL wrapper for calling ad_mom.f
;ad_mom.f is a fortran program written by Phil Fischer
;which calculates Adaptively Weighted Moments
;it adjusts the weighting function until the measurement 
;is optimal. This happens when the moments of the weight 
;are twice the weighted moments. Then the unweighted moments 
;are just equal to the moments of the weight and therefore 
;twice the weighted moments. It returns the unweighted moments 
;measured in this weighted way.
;cat: the sextractor catalog (an IDL structure)
;im: the image
;ixx,iyy,ixy the output adaptive moments
;err: a measure of uncertainty is 9.99 if moments are
;not found (ie. does not converge sensibly)
;numiter: the maximum number of iterations
;wcenx,wceny: the output weighted centroids
;whyflag: the reason why it did not converge 
;see the fortran code for their meanings
;sky:input an array of sky values for each objects
;if not present it will calculate the sky for the whole image
;and use this one number for all  
;if im is undefined it will run on a test image
; -Dave Johnston
; -Erin Scott Sheldon  Changed sofile, added comments
;       changed sky stuff

  IF n_params() EQ 0 THEN BEGIN 
      print,'-syntax '
      return
  ENDIF 

  COMMON huan_admom_block, defval, defval2

  sz = size(image, /tname)

  IF sz NE 'FLOAT' AND sz NE 'DOUBLE' THEN message,'image must be single or double precisioni floating point'

  x = float(cat.x_image)
  y = float(cat.y_image)

  n = n_elements(x)

  ixx = replicate(defval, n)
  iyy = ixx
  ixy = ixx

  IF n_elements(sky) EQ 0 THEN BEGIN 
      sky, image, sk, sigsk, niter=4, nsig=3.5
      sigsk = replicate(sigsk,n)
      sk    = replicate(sk,n)
  ENDIF ELSE BEGIN
      sk    = replicate(sky, n)

      IF n_elements(rms_sky) EQ 0 THEN $
        message,'You must input rms_sky with sky'
      sigsk = replicate(rms_sky, n)

  ENDELSE

  print,'Sky: '+ntostr(sky)+'   rms sky: '+ntostr(rms_sky)

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


PRO huan_admom, imlist, catlist, sky, rms_sky

; NAME:
;       HUAN_ADMOM
; PURPOSE:
;	Run images and sextractor catalogs through the adaptive moment 
;       C code.  The only thing the catalog must have is x_image,y_image
;
; CALLING SEQUENCE:
;       huan_admom, imlist, catlist, sky, rms_sky
;
; INPUTS:
;       imlist: list of image names, full path
;       catlist: list of catalog names, full path.  The new catalog will go in
;                the same directory with _admom appended before the .fits
;
;-
; On_error,2              ;Return to caller
  
  IF N_params() LT 3 THEN BEGIN 
      print,'Syntax - huan_admom, imlist, catlist, rms_sky'
      return
  ENDIF 

  COMMON huan_admom_block, defval, defval2

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
      outfile = namearray[0]+"_admom.fit"

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
          image = image*bscale+bzero
      ENDIF 
	
      ;; put x,y in IDL notation (not SExtractor, which sets first pixel to
      ;; [1,1])
      cat.x_image = cat.x_image - 1.0
      cat.y_image = cat.y_image - 1.0

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;;  Call adaptive moment wrapper
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      huan_call_admom, image, cat, ixx, iyy, ixy, rho4, err, whyflag,$
                        sky=sky[i], rms_sky = rms_sky[i]
      
      ;; Put them back
      cat.x_image = cat.x_image + 1.0
      cat.y_image = cat.y_image + 1.0

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
                            'r',         defval,  $
                            'ellip_err', defval,  $
                            'psf_fwhm',  defval,  $
                            's2n',       defval,  $
                            'whyflag',   defval2)	

      outstruct = replicate(strf2, nobj)
      copy_struct, cat, outstruct

      e1 = replicate(defval, nobj) ;bad are set to defval
      e2 = e1
      r = e1

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
      mwrfits,outstruct,outfile,/create

      ptime,systime(1)-tim

  ENDFOR 

  return
END 








