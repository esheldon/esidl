PRO bins2n, bin_var, signal, nbin, nrand, w1=w1, w2=w2, w3=w3, w4=w4, w5=w5, w6=w6, w7=w7, w8=w8, w9=w9, w10=w10, w11=w11, w12=w12, w13=w13, w14=w14, binnum=binnum

  ;; Assumptions: the noise of each measurement is the same
  ;; Won't work correctly for bin_var an integer. Use bins2n_integer
  ;;
  ;; Need to convert to using histogram

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: bins2n, '
      return
  ENDIF 

  COMMON seed,seed

  IF nbin LT 2 THEN BEGIN 
      print,'Must have nbins > 2'
      return
  ENDIF 
  ind_arr = lonarr(nrand, nbin+1)

  sn_arr  = fltarr(nrand, nbin)
  sn_diff = fltarr(nrand)

  n=n_elements(bin_var)
  s = sort(bin_var)

  n=n_elements(bin_var)
  ;; set to n, so check will actually say le n-1
  ind_arr[*,nbin] = n

  index = lindgen(n)

  min = 0.
  max = float(n-1)
  i=0L
  WHILE (i LE nrand-1) DO BEGIN 
          
      cont=1
      WHILE cont DO BEGIN 
          ;; define bins
          ind = round( arrscl(randomu(seed,nbin-1),min,max,arrmin=0.,arrmax=1.) )
          
          IF n_elements(rem_dup(ind)) EQ nbin-1 THEN cont = 0
      ENDWHILE 
      ind = ind(sort(ind))
      ;; print,ind
      ;; first always set to 0, last always set to nbin
      ind_arr[i,1:nbin-1] = ind      

      FOR k=0L, nbin-1 DO BEGIN 
          wb = where( index GE ind_arr[i,k] AND index LT ind_arr[i,k+1], nwb)

          IF nwb NE 0 THEN BEGIN 
              sn_arr[i,k] = mean( signal[s[wb]] )/sqrt(nwb)
          ENDIF 

      ENDFOR 

      ;; find total s/n difference
      wzero=where(sn_arr[i,*] EQ 0.,nzero)
      IF nzero EQ 0 THEN BEGIN 
          FOR k=0L, nbin-2 DO BEGIN 
              FOR q=k+1, nbin-1 DO BEGIN 
                  sn_diff[i] = sn_diff[i] + abs( sn_arr[i,k] - sn_arr[i,q] )
              ENDFOR 
          ENDFOR 
          print,'.',format='(a,$)'
          i=i+1
      ENDIF ELSE BEGIN 
          sn_diff[i] = 0.
      ENDELSE 
          

  ENDWHILE  

  wgood = where(sn_diff EQ min(sn_diff), ngood)
;  help,wgood
  print
  print,'Minimum S/N diff: ',sn_diff[wgood[0]]
;  colprint,ind_arr[wgood[0],*]
  print
  print,'S/N in each bin: '
  colprint,lindgen(nbin)+1,sn_arr[wgood[0], *]
  print
  print,'Min = ',min(bin_var)
  print,'Max = ',max(bin_var)
  print,'      bin        Nbin   <bin_var>     cut'
  print,'--------------------------------------------------------------'
  FOR k=0L, nbin-1 DO BEGIN 
      kstr=ntostr(k+1)

      command = 'w'+kstr+$
        ' = where( index GE ind_arr[wgood[0],k] AND index LT ind_arr[wgood[0],k+1], nw'+kstr+')'
      temp=execute(command)
      command = 'w'+kstr+' = s[w'+kstr+']'
      temp=execute(command)

;      avecomm = 'lensave, struct[w'+kstr+'], tag, ave, unc, element=element'
;      temp=execute(avecomm)

      avecomm = 'ave = mean(bin_var[w'+kstr+'])'
      temp = execute(avecomm)

      command = 'print,'+kstr+', nw'+kstr+', ave'
      IF k NE nbin-1 THEN command = command + ', max(bin_var[w'+kstr+'])'
      temp=execute(command)

  ENDFOR 

  return
END 
