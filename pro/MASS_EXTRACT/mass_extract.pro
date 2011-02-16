PRO mass_extract, k, kerr, n, nerr, minx, maxx, miny, maxy, kcat, ncat, $
                  galimage=galimage, galcat=galcat

  IF n_params() LT 6 THEN BEGIN
      print,'-Syntax: mass_extract, k, kerr, n, nerr, minx, maxx, miny, maxy,'
      print,'                       kcat, ncat,'
      print,'                       galimage=galimage, galcat=galcat'
      print,' Should be significance maps'
      return
  ENDIF 

  IF n_elements(galimage) NE 0 THEN dogal = 1 ELSE dogal = 0
  
  nk = n_elements(k)
  nkerr = n_elements(kerr)
  nn =  n_elements(n)
  nnerr = n_elements(nerr)
  IF (nk NE nkerr) OR (nk NE nn) OR (nk NE nnerr) THEN BEGIN
      print,'All arrays must be same size'
      return
  ENDIF 
  ksig = k/kerr
  nsig = n/nerr

  sz = size(k)
  sx = sz[1]
  sy = sz[2]

  IF dogal THEN BEGIN
      szg = size(galimage)
      sxg = szg[1]
      syg = szg[2]
      IF (sxg NE sx) OR (syg NE sy) THEN BEGIN
          print,'Maps must be same size'
          return
      ENDIF 
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Look for peaks in SIGNIFICANCE map.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sdss_extract, ksig, tkcat
  sdss_extract, nsig, tncat

  newstr = create_struct('ra', 0d, 'dec', 0d, 'k', 0., 'kerr', 0.)
  nkcat = n_elements(tkcat)
  nncat = n_elements(tncat)
  IF nkcat NE 0 THEN BEGIN
      
      s = replicate(newstr, nkcat)

      tkcat.s2n = ksig[ tkcat.x_image, tkcat.y_image ]

      s.k = k[ tkcat.x_image, tkcat.y_image ]
      s.kerr = kerr[ tkcat.x_image, tkcat.y_image ]
      s.dec = arrscl( tkcat.x_image, minx, maxx, $
                      arrmin=0, arrmax = sx-1)
      s.ra  = arrscl( tkcat.y_image, miny,  maxy,  $
                      arrmin=0, arrmax = sy-1)

      
      combine_structs, tkcat, s, kcat

  ENDIF ELSE BEGIN
      print,'Nothing found in kappa map'
  ENDELSE 

  IF nncat NE 0 THEN BEGIN
      
      s = replicate(newstr, nncat)

      tncat.s2n = nsig[ tncat.x_image, tncat.y_image ]

      s.k = n[ tncat.x_image, tncat.y_image ]
      s.kerr = nerr[ tncat.x_image, tncat.y_image ]
      s.dec = arrscl( tncat.x_image, minx, maxx, $
                      arrmin=0, arrmax = sx-1)
      s.ra  = arrscl( tncat.y_image, miny,  maxy,  $
                      arrmin=0, arrmax = sy-1)

      
      combine_structs, tncat, s, ncat

  ENDIF ELSE BEGIN
      print,'Nothing found in noise map'
  ENDELSE 

  IF dogal THEN BEGIN

      sdss_extract, galimage, tgcat
      ncat = n_elements(tgcat)
      IF ncat NE 0 THEN BEGIN
          
          s = replicate(newstr, ncat)

          tgcat.s2n = galimage[ tgcat.x_image, tgcat.y_image ]
          
          s.dec = arrscl( tgcat.x_image, minx, maxx, $
                          arrmin=0, arrmax = sx-1)
          s.ra  = arrscl( tgcat.y_image, miny,  maxy,  $
                          arrmin=0, arrmax = sy-1)
          
          combine_structs, tgcat, s, galcat

      ENDIF ELSE BEGIN
          print,'Nothing found in galaxy density map'
      ENDELSE 

  ENDIF 

  return
END 
