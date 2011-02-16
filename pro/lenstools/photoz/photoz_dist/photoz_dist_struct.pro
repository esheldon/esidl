PRO photoz_dist_struct, struct, z, pofz, w, minz=minz, maxz=maxz, $
                        zeroprior=zeroprior

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: photoz_dist, struct, z, pofz, useind, minz=minz, maxz=maxz, /zeroprior'
      return
  ENDIF 

  ;; build up a photoz histogram from the input photoz
  ;; assuming photoz's are distributed in a guassian fashion
  ;; with width photoz_zerr

  ;; make cuts
  photoz_select_struct, struct, w

  IF w[0] EQ -1 THEN BEGIN 
      print,'No objects passed cuts'
      return
  ENDIF 

  nin = n_elements(struct)
  nw = n_elements(w)
  print,'Using '+ntostr(nw)+'/'+ntostr(nin)
  ;; deltaz in redshift
 
  nsig = 3.5
  deltaz = 0.01
  IF keyword_set(zeroprior) THEN BEGIN
      minz = 0.0
  ENDIF ELSE BEGIN 
      minz = min(struct[w].photoz_z - nsig*struct[w].photoz_zerr)
  ENDELSE 
  maxz = 1.5

  nz = float(long( (maxz-minz)/deltaz ) )

  photoz_dist, struct[w].photoz_z, struct[w].photoz_zerr, nz, z, pofz, $
               minz=minz, maxz=maxz


return
  z = rnd( arrscl( findgen(nz), minz, maxz ), 2)
  pofz = fltarr(nz)

  FOR i=0L, nw-1 DO BEGIN 

      ind = w[i]
      photoz = struct[ind].photoz_z
      photozerr = photoz_zerr[ind]
      vals = (z - photoz)^2/2./photozerr^2 < 9.0

      pofzi = exp(-vals)/sqrt(2.*!pi)/photozerr

      pofz[*] = pofz[*] + pofzi[*]

  ENDFOR 

  pofz = pofz/total(pofz)

END 
