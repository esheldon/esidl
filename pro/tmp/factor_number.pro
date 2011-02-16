PRO factor_number, N, a, t, xx=xx, amplitude=amplitude

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: factor_number, N [, a, t]'
      return
  ENDIF 

  nrand = 100
  IF n_elements(a) EQ 0 THEN a = 7LL ELSE a = long64(a)
  IF n_elements(t) EQ 0 THEN t = 11LL ELSE t = long64(t)
  num = 2ll^t

  N = long64(N)

  print
  print,'Factoring N = ',ntostr(N)
  print,'Using a = ',ntostr(a)
  print,'      t = ',ntostr(t)
  print,'      2^t = ',ntostr(num)
  print
  print,'Creating state'

  karr = lon64arr(num)
  modarr = lon64arr(num)

  found = 0
  FOR i=0ll, num-1 DO BEGIN 
      
      karr[i]=i
      IF found THEN BEGIN 
          modarr[i] = modarr[i-order]
      ENDIF ELSE BEGIN 
          modarr[i] = a^i MOD N
          IF (i NE 0) AND (modarr[i] EQ 1) THEN BEGIN 
              found=1
              order = i
          ENDIF 
      ENDELSE 
      
  ENDFOR 

  ;; print some results
  nprint = 20

  print,'1/sqrt(2^t) sum( |k>|a^k mod N> = '
  print,'1/sqrt('+ntostr(num)+') [',format='(a,$)'
  FOR i=0L, nprint-1 DO BEGIN
      print,'|'+ntostr(karr[i])+'>|'+ntostr(modarr[i])+'> + ',format='(a,$)'
  ENDFOR 
  print,'.....]'
  print 

  ;;GOTO,jump2
  ;; measure second register: random number
  print
  print,'Using implicit measurement: random second register'
  rand=randomu(seed,num)
  sr = sort(rand)
  chose = modarr[sr[0]]
  print,'chose ',ntostr(chose)
  w=where(modarr EQ chose, nw)
  newkarr = lon64arr(num)

  ;; compute QFT
  print,'Computing QFT'
  newkarr[w] = sqrt(chose/num)*karr[w]
  transform = dcomplexarr(num)
  amplitude = fltarr(num)

  xx=lindgen(num)
  FOR i=0L, num-1 DO BEGIN 
      transform[i] = $
        total(complex( cos(2.*!pi*i*xx[w]/num), sin(2.*!pi*i*xx[w]/num) ))/num
      amplitude[i] = abs(transform[i])
  ENDFOR 
  
  ;; pick those with high probabilities 
  ;; so I don't have to code up the bad choices
  maxamp = max(amplitude)
  maxy=maxamp
  
  psfile = 'factor'+ntostr(N)+'.ps'

  begplot,name=psfile
  aplot,!gratio,xx,amplitude,xstyle=1,yrange=[0,maxy], xrange=[-100,num], $
    ytitle='|a!Dl!N|!U2',xtitle='l',charsize=2,title='Factoring '+ntostr(N)
  endplot 

  wmeas=where(amplitude GT 0.75*maxamp,order)
  GOTO,jump2

  trand=randomu(seed, nrand)
  sr=sort(trand)
  randoms = xx[ wmeas[sr] ]
  
  ;; find r
  continue = 1
  i=0L
  inr=0
  WHILE continue DO BEGIN 

      IF i EQ nrand THEN BEGIN
          print,'Out of randoms'
          return
      ENDIF 
      rand = long( round(randoms[i]) )
      
      print,'Measuring...Got ',rand
      IF rand EQ 0 THEN BEGIN
          message,'Got zero, redoing',/inf
          GOTO,jump
      ENDIF 
      print,'Finding continued fraction representation'
      continued_fraction, rand, num, coeffs
      print,'Coefficients: '
      colprint,coeffs
      print,'What r does this imply?'
      read,inr
      IF inr LT 0 THEN BEGIN 
          message,'Ok, must not have reduced. Trying again',/inf
          GOTO,jump
      ENDIF 
      IF inr MOD 2 NE 0 THEN BEGIN 
          message,'r is odd, no good. Trying again',/inf
          i=i+1
      ENDIF ELSE BEGIN 
          print,order,inr
          order = inr
          print
          print,"Finding GCD's using Euclid's algorithm"
          
          n1 = a^(order/2) - 1
          n2 = a^(order/2) + 1
          IF (n1 LT N) OR (n2 LT N) THEN BEGIN 
              message,'r too small',/inf
              GOTO,jump
          ENDIF 

          fac1 = euclid_gcd(n1, N)
          fac2 = euclid_gcd(n2, N)
          print,'Factored ',ntostr(N),' into ',ntostr(fac1),'*',ntostr(fac2)
          continue = 0
      ENDELSE 
 

      jump:
      i=i+1
  ENDWHILE 

  return

jump2:
  print,'Found r = ',ntostr(order)
  IF order MOD 2 NE 0 THEN BEGIN 
      print,'r is odd, try new (a)'
      return
  ENDIF

  print
  print,"Finding GCD's using Euclid's algorithm"

  n1 = a^(order/2) - 1
  n2 = a^(order/2) + 1
  
  fac1 = euclid_gcd(n1, N)
  fac2 = euclid_gcd(n2, N)
  print,'Factored ',ntostr(N),' into ',ntostr(fac1),'*',ntostr(fac2)


  return


END 
