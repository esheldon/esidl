PRO photoz_dist, photoz, photoz_err, nz, zvals, pofz, minz=minz, maxz=maxz

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: photoz_dist, z, zerr, npofz [, z, pofz]'
      return
  ENDIF 

  on_error, 2

  nn = n_elements(photoz)

  nsig = 3.5
  IF n_elements(minz) EQ 0 THEN BEGIN 
      minz = min(photoz-nsig*photoz_err)
  ENDIF 
  IF n_elements(maxz) EQ 0 THEN BEGIN 
      maxz = max(photoz+nsig*photoz_err)
  ENDIF 

  IF maxz EQ 0.0 THEN message,'Bad z'

  zvals = arrscl( findgen(nz), minz, maxz )
  pofz = fltarr(nz)

  FOR i=0L, nn-1 DO BEGIN 

      pofz[*] = pofz[*] + gaussprob(zvals, photoz[i], photoz_err[i])

  ENDFOR 

  ;; normalize
  nn = 200

  gauleg, minz, maxz, nn, ZZi, WWi
  
  pofzi = interpol(pofz, zvals, ZZi)
  norm = total( pofzi*WWi )

  print,norm

  pofz = pofz/norm

  ;;plot,zvals, pofz
  

END 
