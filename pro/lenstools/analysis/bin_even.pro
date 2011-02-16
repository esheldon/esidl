PRO bin_even, struct, nbin, tag, element=element, w1=w1, w2=w2, w3=w3, w4=w4, w5=w5, w6=w6, w7=w7, w8=w8, w9=w9

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: bins2n_lum, struct, nbin, nrand, tag, element=element, w1=w1, w2=w2, w3=w3, w4=w4, w5=w5, w6=w6, w7=w7, w8=w8, w9=w9'
      return
  ENDIF 

  COMMON seed,seed

  IF nbin LT 2 THEN BEGIN 
      print,'Must have nbins > 2'
      return
  ENDIF 

  IF NOT tag_exist(struct, tag, index=wt) THEN BEGIN 
      print,'No such tag: ',tag
      return
  ENDIF 

  IF  n_elements(element) NE 0 THEN BEGIN 
      comm_add='(wt)[element]'
      val = struct.(wt)[element]
      s = sort( val )
  ENDIF ELSE BEGIN 
      comm_add='(wt)'
      val =struct.(wt)
      s = sort(val)
  ENDELSE 
  st = struct[s]
  val = val[s]

  n=n_elements(struct)
  index = lindgen(n)
  ind_arr = lonarr(nbin+1)
  ind_arr[nbin] = n-1
  ind_arr[0] = 0L
  
  ;; break into bins as evenly as possible
  numbin = long(n/float(nbin))  ; number per bin. round down, stick extra in last bin
  left = long(n - numbin*nbin)
  print,left
  print
  print,'      bin        Nbin   <'+tag+'>     cut'
  print,'--------------------------------------------------------------'
  FOR k=0L, nbin-1 DO BEGIN 

      kstr=ntostr(k+1)
      ind1 = k*numbin
      IF (k EQ nbin-1) THEN BEGIN 
          ind2 = (k+1)*numbin + left - 1
      ENDIF ELSE BEGIN 
          ind2 = (k+1)*numbin-1
      ENDELSE 
      
      command = 'w'+kstr+' = index['+ntostr(ind1)+':'+ntostr(ind2)+'] & nw'+kstr+' = n_elements(w'+kstr+')
;      print,command
      temp=execute(command)
      command = 'w'+kstr+' = s[w'+kstr+']'
      temp=execute(command)

      avecomm = 'lensave, struct[w'+kstr+'], tag, ave, unc, element=element'
      temp=execute(avecomm)

      command = 'print,'+kstr+', nw'+kstr+', ave'
      IF k NE nbin-1 THEN command = command + ', max(struct[w'+kstr+'].'+comm_add+')'
      temp=execute(command)
  ENDFOR 
  return

  return
END 
