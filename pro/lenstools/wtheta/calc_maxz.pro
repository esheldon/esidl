PRO calc_maxz,clr

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: calc_maxz, clr'
      return
  ENDIF 

  ;; For a series of lens redshifts, given a mean g-r of galaxies, 
  ;; calculate the apparent magnitude of f*lstar galaxy.  Max(z) is 
  ;; highest z s.t. the apparent mag is less than the magnitude limit
  ;; in that bandpass

  nn=10000L
  z=arrscl(findgen(nn), 0.01, .4)
  a=angdist_lambda(z,dlum=dlum)

  fraclstar = 0.1
  mag_limit = [21.0, 21.0, 21.0, 21.0, 19.8]

  Mstar = [-18.34,-20.04,-20.83,-21.26,-21.55]
  Msun = [6.39,5.07,4.62,4.52,4.48]
  Lstar = 10.0^((Mstar-Msun)/(-2.5))

;  print,lstar

  ;; mean g-r
  gmr_tmp = 0.3
  gmr = replicate(gmr_tmp, nn)

  ;; find kcorr at each redshift
  ;; dummy array, we are just getting the kcorr
  mag=replicate(20.0, nn)
  ;; get kcorr
  wtheta_absmag_diffz, z, clr, mag, gmr, tabmag, tlumsolar, kcorr=kcorr

  ;; find apparent magnitude of fraclstar*Lstar at each redshift
  m_max = Mstar[clr] - 2.5*alog10(fraclstar) + $
    5.0*alog10(dlum*10.^6/10.0) - kcorr

  maxz = max( z[ where(m_max LE mag_limit[clr]) ] )
  plot, z, m_max, xtitle='Z', $
    ytitle='apparent mag of '+ntostr(fraclstar,3)+' L* gal'
  oplot,[maxz, maxz], [-10, 1000], line=2

  print
  print,'Max(z) s.t. we can see a '+ntostr(fraclstar,3)+' L* gal '
  print,'within mag limiit max_mag['+!colors[clr]+'] = '+ntostr(mag_limit[clr])+' sample'
  print,maxz

  ;; now force maxz = 0.2 and find mag
  CASE clr OF
      0: maxz = 0.0855
      1: maxz = 0.15
      ELSE: maxz = 0.2179
  ENDCASE 
  print
  print,'Now fixing maxz = '+ntostr(maxz,3)+' we find mag limit of '
  max_mag = max( m_max[ where(z LE maxz) ] )
  print,max_mag

END   
