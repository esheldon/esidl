PRO run_admom, image, x, y, sky, rms_sky, wguess, $
               ixx, iyy, ixy, rho4, err, whyflag, $
			   unload=unload

  ; need to do more type checking
  IF n_params() EQ 0 THEN BEGIN 
      print,'-syntax '
      return
  ENDIF 

  defval = -9999.0
  defval2 = 9999L

  sz = size(image, /tname)

  IF sz NE 'FLOAT' THEN message,'image must be float'

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

  home=getenv("HOME")
  sofile = home+'/ccode/admomatlas/bin/admom_float.so'
  entry = 'admom_float'
 
  print
  print,'Using sofile:  ',sofile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Call the C routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Calling admom_float.c'
  test = call_external(sofile, entry, $
                       image, nx, ny, x, y, n, sk, sigsk, wguess, $
                       ixx, iyy, ixy, err, rho4, whyflag, unload=unload)

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
