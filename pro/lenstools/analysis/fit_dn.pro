PRO fit_dn, struct, dn, dnarc, dnpix, thresh=thresh, clr=clr

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: fit_dn, struct, dn, dnarc, dnpix, thresh=thresh, clr=clr'
      return
  ENDIF 

  wz = getztag(struct)
  IF wz EQ -1 THEN return

  IF n_elements(thresh) EQ 0 THEN thresh = 20.75
  IF n_elements(clr) EQ 0 THEN clr = 1 ;g-band
  nobj = n_elements(struct)
  dn = fltarr(nobj)
  dnarc = dn
  dnpix = dn

  ;; get the profmean radii in arcsec
  profmean_rad, pixrad, arcsec, pixels

  FOR i=0L, nobj-1 DO BEGIN 
      w=where(struct[i].profmean[*,clr] NE 0., nw)
      IF nw NE 0 THEN BEGIN 
          min = min(struct[i].profmean[w,clr])
          max = max(struct[i].profmean[w,clr])

          IF (max GT thresh) AND (min LT thresh) THEN BEGIN 

              dnpix[i] = interpol(pixrad[w], struct[i].profmean[w,clr],thresh)
              dnarc[i] = interpol(arcsec[w], struct[i].profmean[w,clr],thresh)

              ;; convert to kpc
              z = struct[i].(wz)
              angdist = angdist_lambda(z)*1000. ;kpc
              dn[i] = dnarc[i]/3600.*!pi/180.*angdist > 0.

          ENDIF 
      ENDIF 
  ENDFOR 

  return
END 
  
