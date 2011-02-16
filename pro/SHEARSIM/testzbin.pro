PRO testzbin, cat, nfront, nback, test=test, nuse=nuse

  IF n_params() EQ 0 THEN BEGIN
    print,'-Syntax: testzbin, cat, nfront, nback, nuse=nuse'
    return
  ENDIF 


  magmin = 16.0
  magmed = 18.0
  magmax = 22.0

  wL=where(cat.mag LT magmed AND cat.mag GT magmin, nwL)
  wS = where(cat.mag LT magmax AND cat.mag GT magmed, nwS)

  print,ntostr(nwL)+' lenses found'

  IF NOT keyword_set(nuse) THEN nuse = nwL ELSE BEGIN 
    IF nwL GE nuse THEN nwL = nuse
    print,'Using ',ntostr(nwL)
  ENDELSE 
  print,ntostr(nwS)+' sources found'

  IF keyword_set(test) THEN BEGIN 

; For each foreground lens, get all the stuff that's background
; withing range box.

    t = systime(1)
    range = 200.                ; arcseconds
    range = 200*2.8e-4          ;degrees

    nfront = 0L
    nback = 0L
    FOR i=0L, nwL-1 DO BEGIN

      radiff = abs(cat[wS].ra - cat[ wL[i] ].ra)
      decdiff =  abs(cat[wS].dec - cat[ wL[i] ].dec)
      
      w=where( radiff  LT range AND decdiff LT range , nw)
      w=wS[w]

      back=where( cat[w].z GT cat[ wL[i] ].z , nb)
      front = where( cat[w].z LE cat[ wL[i] ].z, nf )
    
      nback = nback + nb
      nfront = nfront + nf

      print, format='($,A)', '.'
;    IF i MOD 100 EQ 0 THEN print,i
    ENDFOR 

    t=systime(1) -t
    ptime, t

    print,'N in front: ',nfront
    print,'N in back: ',nback
    print,'Contamination: ',nfront/nback
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; make some histograms
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  plothist, cat[wL].mag,bin=.1, xtitle='Lens Mag'
  key=get_kbrd(1)
  IF key EQ 'q' THEN return 
  plothist, cat[wS].mag,bin=.1, xtitle='Source Mag'
  key=get_kbrd(1)
  IF key EQ 'q' THEN return 
  plothist, cat.z, bin=.005, xtitle='All Z'
  key=get_kbrd(1)
  IF key EQ 'q' THEN return 
  plothist, cat[wL].z, bin=.005, xtitle='Lens Z'
  key=get_kbrd(1)
  IF key EQ 'q' THEN return 
  plothist, cat[wS].z, bin=.005, xtitle='Source Z'





  return
END 
