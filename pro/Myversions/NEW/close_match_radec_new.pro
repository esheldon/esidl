PRO close_match_radec_new_match, ra1, dec1, ra2, dec2, match_radius, allow, matches



END 

PRO close_match_radec_new, ra1, dec1, ra2, dec2, match_radius, m1, m2, $
                           allow=allow, $
                           resolution=resolution

  IF n_elements(resolution) EQ 0 THEN resolutio=4
  IF n_elements(allow) EQ 0 THEN allow=1

  print
  print,'Finding pixel for each point'
  ang2pix, ra1, dec1, ix1, iy1, pixnum1, resolution=resolution, /radec
  ang2pix, ra2, dec2, ix2, iy2, pixnum2, resolution=resolution, /radec

  upix1 = pixnum1[rem_dup(pixnum1)]
  upix2 = pixnum2[rem_dup(pixnum2)]

  npix1 = n_elements(upix1)
  npix2 = n_elements(upix2)

  print,'Histogramming pixels'
  hist1 = histogram(pixnum1, min=0, rev=rev1)
  hist2 = histogram(pixnum2, min=0, rev=rev2)

  pix2ang, upix1, pLam1, pEta1, /survey, resolution=resolution

  print,'Getting Pixel bounds'
  pix_bound, upix1, lamMin1, lamMax1, etaMin1, etaMax1,$
    resolution=resolution,/survey,silent=silent
  pix_bound, upix2, lamMin2, lamMax2, etaMin2, etaMax2,$
    resolution=resolution,/survey,silent=silent

  sp1 = sort(etaMin1)

  n1=n_elements(ra1)
  n2=n_elements(ra2)
  
  index1 = lindgen(n1)
  index2 = lindgen(n2)



  matcharr=lonarr(n1,allow)	;the main book-keeping device for 
  matcharr(*,*)=-1		;matches -initialized to -1
  ind=lindgen(n2)

  print,'Looping'

  runi=0L
  ;; Loop over the pixels for array 1
  FOR ip1=0L, npix1-1 DO BEGIN 

      ;; This pixel
      pix1 = upix1[ip1]

;      print,'pixel number = ',pix1

      ;; Find distance between this pixel and the boundaries
      ;; of all other pixels

      matchRadiusEta = match_radius/cos(pLam1[ip1]*0.01745329)
;      print,'matchRadiusEta = '+ntostr(matchRadiusEta)
      eta1Minus = etaMin1[ip1] - matchRadiusEta
      eta1Plus  = etaMax1[ip1] + matchRadiusEta

;      print,'matchRadius = '+ntostr(match_radius)
      lam1Minus = lamMin1[ip1] - match_radius
      lam1Plus  = lamMax1[ip1] + match_radius

      ;; should do binary searching.  Also, doesn't work acros
      ;; boundaries -180,180 etc.

      wPixMatch2 = where( (lamMax2 LE lam1Plus) AND $
                          (lamMin2 GE lam1Minus) AND $
                          (etaMax2 LE eta1Plus) AND $
                          (etaMin2 GE eta1Minus), nPixMatch2)

;      wPixMatch2=where( (lamMin2 LE lam1Plus AND lamMin2 GE lamMin1[ip1]) OR $
;                        (lamMax2 GE lam1Minus AND lamMax2 LE lamMax1[ip1]) OR $
;                        (etaMin2 LE eta1Plus AND etaMin2 GE etaMin1[ip1]) OR $
;                        (etaMax2 GE eta1Minus AND etaMax2 LE etaMax1[ip1]), $
;                        nPixMatch2)
;stop
      IF nPixMatch2 NE 0 THEN BEGIN 

;          print,'Found '+ntostr(nPixMatch2)+$
;            ' pixels to match pixel '+ntostr(pix1)

          ;; Get objects in pixel1
          wPix1 = rev1[ rev1[pix1]:rev1[pix1+1]-1 ]
          nPix1 = n_elements(wPix1)

;          print,'There are '+ntostr(nPix1)+' objects in pix1 = '+ntostr(pix1)
          ;; loop over objects in pix1
          FOR i1=0L, nPix1-1 DO BEGIN 

              matchesPtr = ptrarr(nPixMatch2)
;              distPtr = ptrarr(nPixMatch2)
              numList = lonarr(nPixMatch2)
              nTotal = 0L


;              print,'.',format='(A,$)'
              ;; Now loop over the pixels2 and match to
              ;; the objects in pix1
              FOR ip2=0L, nPixMatch2-1 DO BEGIN 
                  
                  pix2 = upix2[wPixMatch2[ip2]]
                  
                  wPix2 = rev2[ rev2[pix2]:rev2[pix2+1]-1]
                  
                  close_match_radec, $
                    ra1[wPix1[i1]], dec1[wPix1[i1]], $
                    ra2[wPix2], dec2[wPix2], $
                    tm1, tm2, match_radius, allow,/silent


                  IF tm2[0] NE -1 THEN BEGIN 
                      nw = n_elements(tm2)

                      matchesPtr[ip2] = ptr_new(wPix2[tm2])
;                      distPtr[ip2]    = ptr_new(distance[w])

                      numList[ip2] = nw
                      nTotal = nTotal + nw

                  ENDIF 

;                  mygcirc, $
;                    ra1[wPix1[i1]], dec1[wPix1[i1]], $
;                    ra2[wPix2], dec2[wPix2], $
;                    distance

;                  w=where(distance LE match_radius, nw)

;                  IF nw NE 0 THEN BEGIN 
;                      cpix = wPix2[w]
;                      cdis = distance[w]

;                      matchesPtr[ip2] = ptr_new(cpix, /no_copy)
;                      distPtr[ip2]    = ptr_new(cdis, /no_copy)

;                      matchesPtr[ip2] = ptr_new(wPix2[w])
;                      distPtr[ip2]    = ptr_new(distance[w])

;                      numList[ip2] = nw
;                      nTotal = nTotal + nw
;stop
;                  ENDIF 
              ENDFOR ;; loop over pixels2

              IF nTotal GT 0 THEN BEGIN 
;stop
                  ;; extract all matches from each pixel
                  matches = lonarr(nTotal)
;                  distances = dblarr(nTotal)
                  beg=0L
                  FOR ip2=0L, nPixMatch2-1 DO BEGIN 
                      
                      IF numList[ip2] NE 0 THEN BEGIN 
;stop
                          matches[beg:beg+numList[ip2]-1] = *matchesPtr[ip2]
                          distances[beg:beg+numList[ip2]-1] = *distPtr[ip2]
                          
                          ptr_free, matchesPtr[ip2]
                          ptr_free, distPtr[ip2]
                          
                          beg = beg+numList[ip2]
                      ENDIF 
                  ENDFOR 
                  
                  IF nTotal GT allow THEN BEGIN 
                      
                      matches = matches[sort(distances)]
                      nTotal = allow
                      matches = matches[0:allow-1]
                      
                  ENDIF 
                  
                  matcharr[wPix1[i1], 0:nTotal-1] = matches
                  runi = runi+nTotal

              ENDIF ;; matches found for this object
              
          ENDFOR ;; loop over objects in pixel1
;          print

      ENDIF ;; some pixels2 matched

  ENDFOR ;; loop over pixels1
  

  if not keyword_set(silent) then print,'total put in bytarr',runi

  matches=where(matcharr ne -1,this)
 
  if this eq 0 then begin
      if not keyword_set(silent) then $
        print,'no matches found'
      miss1=lindgen(n1)
      m1=-1 & m2=-1
      return
  ENDIF

  m1=matches mod n1             ;a neat trick to extract them correctly 
  m2=matcharr[matches]          ;from the matcharr matrix

  if not keyword_set(silent) then $
    print,n_elements(m1),' matches'

  m2=ind[m2]                    ;remember, must unsort
  dif=m1[uniq(m1,sort(m1))]
  if not keyword_set(silent) then $
    print,n_elements(dif),' different matches'
  if n_params() eq 9 then begin
      if n_elements(m1) lt n1 then begin
          miss1=lindgen(n1)
          remove,dif,miss1
          if not keyword_set(silent) then $
            print,n_elements(miss1),'  misses'
      endif else begin
          miss1=-1  
          if not keyword_set(silent) then $
            print,'no misses'
      endelse 
  endif


return


;PRO tmp, struct, scat, match1, match2, m1, m2

  file = '/net/cheops1/data0/esheldon/Rachel/LRGcatalog_fix.dat'

  IF n_elements(struct) EQ 0 THEN BEGIN 
      s=create_struct('ra',0d,'dec',0d,$
                      'e1',0.,'e2',0.,$
                      'e1e1var',0.,'e1e2var',0.,'e2e2var',0.,$
                      'counts_model',fltarr(5), 'counts_modelerr', fltarr(5), $
                      'r_i',0.,$
                      'photoz_z',0.,$
                      'photoz_lrg',0.)

      read_struct, file, s, struct
  ENDIF 

  IF n_elements(scat) EQ 0 THEN BEGIN 
      stripes = [9,10,11,12,13,14,15,$
                 27,28,29,30,31,32,33,34,35,36,37,$
                 76,82,86]

      columns = ['clambda','ceta','e1','e2','e1e1err','e1e2err','e2e2err',$
                 'grmodel','rimodel']
      get_scat, stripes, [1,2,3], scat, /hirata, columns=columns
  ENDIF 

  csurvey2eq, scat.clambda, scat.ceta, ra, dec


  tol = 3d/3600d
  tt = systime(1)
  spherematch, ra, dec , struct.ra, struct.dec, $
    tol, match1, match2, dis12
  ptime, systime(1)-tt

  allow = 1
  tt = systime(1)
  close_match_radec, ra, dec, struct.ra, struct.dec, $
    m1, m2, tol, allow
  ptime, systime(1)-tt

;;spherematch, ra1, dec1, ra2, dec2, matchlength, match1, match2, $
;;   distance12, [maxmatch=]
;; close_match,ra1,dec1,ra2,dec2,m1,m2,ep,allow,miss1,silent=silent
return
end        

