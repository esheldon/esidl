PRO plot_sigcrit, print=print


  IF NOT keyword_set(print) THEN print=0
  name = '/sdss3/usrdevel/esheldon/idl.lib/SHEARSIM/rel_shear.ps'
  IF print THEN makeps, name

  xt = 'Source Redshift'
  yt = '1/Sigma Crit'
  t='1/Sigma Crit  ( cm/s^2 )^-1'

  zsmin=0
  zsmax=3

  zl = [.4, .2, .1, .05, .02]


  FOR i=0, n_elements(zl) -1 DO BEGIN 
    
      zs = arrscl(findgen(100), zsmin, zsmax) > zl[i]*1.001
      f = 1/sigmacrit(zs, zl[i])
      IF i EQ 0 THEN BEGIN 
          plot, zs, f, xtitle=xt, ytitle=yt, title=t
      ENDIF ELSE BEGIN
          oplot, zs, f
      ENDELSE 
      mf = max(f)
      xyouts, 2.25, 1.025*mf, 'Zlens = '+ntostr(zl[i], 4)
  ENDFOR 

  key = get_kbrd(1)

  yt = 'Dls/Ds  (Rel Shear)'
  t = 'Relative Shear  (SIS)'

  zsmin=0
  zsmax=5

  psym = [1,2,4,5,6]
  leg = ntostr(zl, 4)
  FOR i=0, n_elements(zl) -1 DO BEGIN 
    
      zs = arrscl(findgen(100), zsmin, zsmax) > zl[i]
      f = angdist(zs, zl[i])/angdist(zs)
      IF i EQ 0 THEN BEGIN 
          plot, zs, f, xtitle=xt, ytitle=yt, title=t, yrange=[0,2],psym=psym[i]
      ENDIF ELSE BEGIN
          oplot, zs, f,psym=psym[i]
      ENDELSE 

  ENDFOR 
  legend, leg, psym=psym

  IF print THEN ep


  G = 6.67e-8                   ;cm^3/g/s^2
  sigv = 220.                   ;km/s
  sigv = sigv*1.e5              ;cm/s
  R = 100.                      ;kpc
  R = R*1000.                   ;pc
  R = R*3.09e18                 ;cm

  print,1/2.*sigv^2/G/R

  a = 1000.                     ;km/s
  a = a*1.e3                    ;m/s
  a = a*1.e2                    ;cm/s
  
  b = 1                         ;Mpc
  b = b*1.e6                    ;pc
  b = b*3.09e18                 ;cm

  f = 180./!pi*3600.            ;arcsec/radian

  print,1/2.*a^2/b/G*.22^2/42.


  return 
END 
