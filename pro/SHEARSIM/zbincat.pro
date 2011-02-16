PRO zbincat, cat, solidangle

  IF n_params() EQ 0 THEN BEGIN
    print,'-Syntax: zbincat, cat, solidangle'
    print,' Two overlapping runs ~.075 square radians'
    return
  ENDIF 

  COMMON seed, seed
  time=systime(1)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; cosmology.  Right now, ang. diam. distances don't include lambda
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  h = .7
  omega = 1.0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  pi4 = 4.*!dpi
  IF n_elements(solidangle) EQ 0 THEN solidangle = pi4
  frac = solidangle/pi4

  dzS = 1.e-4
  dzL = 1.e-4                   ;Lenses in [0,zmed], sources in [zmed,zmax]
  zmed = .01
  zmax = 1.0
  nLbins = long( zmed/dzL )
  nSbins = long( (zmax - zmed)/dzS )

                                ; make an array of z's  Make seperate to avoid
                                ;if statements later
  zsource = dzS*dindgen(nSbins)+zmed
  zlens = dzL*dindgen(nLbins)

                                ; Lstar = 1.0e10/h^2 ;units of Lsun
                                ; Lsun = 3.9e33 ergs/sec
  logLstar = alog10(3.9) - 2.*alog10(h) + 43

  fzero = 2.52e-5               ; flux zero point
  mzero = 26.47                 ; mag zero point 
  magmax = 22.0                 ; magnitude cutoff
  magmin = 16.0
  magmed = 18.0 

  maxL = 3.                     ; In units of Lstar
  minL = .003

  fmin = fzero*10^(-.4*magmax)
  fmed = fzero*10^(-.4*magmed)
  fmax = fzero*10^(-.4*magmin)

  nLumbins = 40                 ; Number of luminosity bins in shecter function

  shfac = 4.0                   ; Factor to multily shecter function

  print,'nSbins ',ntostr(nSbins)
  print,'nLbins ',ntostr(nLbins)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; make the shecter function
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  shecter, nLumbins, shecter, minL, maxL, Lum=Lum,alpha=1.25

  logLum = alog10(Lum) + logLstar ; log10(L*Lstar)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Create structure of galaxy parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  s=create_struct('type', -10 $ ;0 for lens, 1 for source
                  ,'z', -10.0 $
                  ,'mag', -10.0 $
                  ,'ra', double(-1.e8)  $
                  ,'dec', double(-1.e8) $
                  )

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; do lens and source bins separately.  If I pick lenses with
  ; z < .01 then I can use Euclidian formulae for volumes
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                ; Lenses z < zmed  In this range, 
                                ; use Euclidean formulae
  print,'Doing Low-z Volumes'

  lensvol = dblarr(nLbins)
  dang1 = angdist(zlens)
  dang2 = angdist(zlens + dzL)
  lensvol = frac*pi4/3.*( dang2^3 - dang1^3 )

  plothist,lensvol,xtitle='Shell volume z < .01'
  key=get_kbrd(1)
  IF key EQ 'q' THEN return

  plot, zlens, lensvol,xtitle='z',ytitle='Shell Volume'
  key=get_kbrd(1)
  IF key EQ 'q' THEN return

                                ; Calculate luminosity distance of geometric 
                                ; mean of each z bin
  zmean=3./4.*( (zlens+dzL)^4 - zlens^4)/( (zlens+dzL)^3 - zlens^3 )
  dlum = lumdist(zmean)

  plot,zmean, dlum,xtitle='mean z', ytit='luminosity dist'
  key=get_kbrd(1)
  IF key EQ 'q' THEN return

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Use logarithms to handle big numbers
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                ; log10(4*pi*dl^2)
                                ; convert to cm: 1Mpc = 3.e24cm

  logArea = alog10(pi4) + 2.*alog10(dlum) + 2.*alog10(3.) + 2.*24.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Apply the shecter function for each low-z bin.  Keep only those
  ; objects with magnitude lt magmed and gt magmin ( do in flux )
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  t = nLumbins*(nSbins + nLbins) ; Maximum possible number of bins
  numlist = lonarr(t)
  ptrlist = ptrarr(t)
  ntot=0L
  ngbins = -1L

  FOR Li=0, nLbins-1 DO BEGIN
                                
    f = 10.^( logLum - logArea[Li] ) ; Flux at each luminosity for this area
    w=where( f GE fmin AND f LE fmax , nw) ; Objects in mag range

    IF nw NE 0 THEN BEGIN 

      FOR lumi=0,nw-1 DO BEGIN 
        nn = round( lensvol[Li]*shecter[ w[lumi] ] )
        IF nn NE 0 THEN BEGIN 
          ngbins = ngbins + 1   ; Number of used bins
          tmp = replicate(s, nn)
          tmp.type = 0          ; Low z
          tmp.z = zmean[Li]
          flux = f[ w[lumi] ]
          tmp.mag = -2.5*alog10(flux/fzero)

          ptrlist[ngbins] = ptr_new(tmp)
          numlist[ngbins] = numlist[ngbins] + nn
          ntot = ntot + nn
        ENDIF 
      ENDFOR
    ENDIF 
    IF Li MOD 100 EQ 0 OR Li EQ nSbins-1 THEN $
      print,'Bin: ',ntostr(Li),' n elements: ',ntostr(ntot)
  ENDFOR 
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; High-z bin
  ; See above for explanations; its the same code
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print, 'Doing High-z Volumes'

  binzvol, zsource, dzS, sourcevol
  sourcevol = frac*sourcevol

  plothist, sourcevol, bin=500, xtitle='Shell volume z > .01'
  key=get_kbrd(1)
  IF key EQ 'q' THEN return

  plot,zsource, sourcevol,xtitle='Z',ytitle='Shell volume z > .01'
  key=get_kbrd(1)
  IF key EQ 'q' THEN return

  zmean=3./4.*( (zsource+dzS)^4 - zsource^4)/( (zsource+dzS)^3 - zsource^3 )
  dlum = lumdist(zmean)
  
  plot, zmean,dlum,xtitle='Z mean',ytitle='Luminosity distance'
  key=get_kbrd(1)
  IF key EQ 'q' THEN return

  logArea = alog10(pi4) + 2.*alog10(dlum) + 2.*alog10(3.) + 2.*24.

  FOR Si=0, nSbins-1 DO BEGIN
                                
    f = 10.^( logLum - logArea[Si] ) ; Flux at each luminosity for this area
    w=where( f GE fmin AND f LE fmax , nw)  ; Objects in mag range

    IF nw NE 0 THEN BEGIN 
      FOR lumi=0,nw-1 DO BEGIN 
        nn = round( sourcevol[Si]*shecter[ w[lumi] ] )
        IF nn NE 0 THEN BEGIN 
          ngbins = ngbins + 1   ; Number of used bins

          tmp = replicate(s, nn)
          tmp.type = 1          ; Higher z
          tmp.z = zmean[Si]
          flux = f[ w[lumi] ]
          tmp.mag = -2.5*alog10(flux/fzero)

          ptrlist[ngbins] = ptr_new(tmp)
          numlist[ngbins] = numlist[ngbins] + nn
          ntot = ntot + nn
        ENDIF 
      ENDFOR
    ENDIF 
    IF Si MOD 100 EQ 0 OR Si EQ nSbins-1 THEN $
      print,'Bin: ',ntostr(Si),' n elements: ',ntostr(ntot)
 
  ENDFOR 


  cat = replicate(s, ntot)
  
  beg=0L                     ; remember, ngbins is always 1 too small
  FOR i=0L, ngbins DO BEGIN
    num = numlist[i]
    cat[ beg:beg+num-1 ] = *ptrlist[i]
    ptr_free, ptrlist[i]
    beg = beg+num
  ENDFOR 
    
  leftover = t-ngbins
  print,ntostr(leftover)+' leftover bins'
  IF leftover GT 0 THEN BEGIN 
    FOR i=ngbins, t-1 DO ptr_free, ptrlist[i]
  ENDIF 

  ; Uniformly sample ra-dec range of the overleaved runs: 752, 756
  ramin = 145.14243
  ramax = 236.42260
  decmin = -1.26822
  decmax = 1.27289
  ra = ramax*randomu(seed, ntot) + ramin
  dec = decmax*randomu(seed, ntot) + decmin

  cat.ra = ra
  cat.dec = dec

  time = systime(1)-time
  ptime,time

  return
END

