PRO bins2n_integer, bin_var, signal, nbin, nrand, w1=w1, w2=w2, w3=w3, w4=w4, w5=w5, w6=w6, w7=w7, w8=w8, w9=w9

  ;; Assumptions: the noise of each measurement is the same

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: bins2n_integer, '
      return
  ENDIF 

  tm = systime(1)

  COMMON seed,seed

  IF nbin LT 2 THEN BEGIN 
      print,'Must have nbins > 2'
      return
  ENDIF 
  ;; where the cuts are placed
  cut_arr = lonarr(nrand, nbin, 2)


  sn_arr  = fltarr(nrand, nbin)
  sn_diff = fltarr(nrand)

  n=n_elements(bin_var)
  s = sort(bin_var)

  n=n_elements(bin_var)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; min/max possible cuts.  Note, because we work with integers, bins 
  ;; can be of size one, so max of a bin can be equal to min of the bin
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  min = min(bin_var, max=max)

  h = histogram(bin_var, min=0, rev=rev)

  ;; zero here for easy subscripting
  i=0L
  WHILE (i LE nrand-1) DO BEGIN 
          
      sn_arr[i,*] = 0

      WHILE 1 DO BEGIN 

          ;; define bins: integers from min to max
          cut = round( arrscl(randomu(seed,nbin-1),min,max, $
                              arrmin=0.,arrmax=1.) $
                     )
          
          ;; number can't show up more than twice (one integer wide bin)
          cut2 = [min, cut, max]
          th = histogram(cut2)
          w=where(th GT 2, nw)
          IF nw EQ 0 THEN break
      ENDWHILE 

      cut = cut[sort(cut)]
;      print,cut

      FOR k=0L, nbin-1 DO BEGIN 

          skip = 0

          ;; Always inclusive!
          IF k EQ 0 THEN BEGIN 

              cut_arr[i, k, 0] = min
              cut_arr[i, k, 1] = cut[k] < max

          ENDIF ELSE IF k EQ nbin-1 THEN BEGIN 

              cut_arr[i, k, 0] = cut[k-1]+1 < max
              cut_arr[i, k, 1] = max

          ENDIF ELSE BEGIN 
              cut_arr[i, k, 0] = cut[k-1]+1

              ;; max can only be in the last bin
              IF cut[k] EQ max THEN BEGIN 
                  cut_arr[i, k, 1] = max-1 > cut_arr[i, k, 0]
                  IF cut[k] EQ max THEN skip=1
              ENDIF ELSE BEGIN 
                  cut_arr[i, k, 1] = cut[k] > cut_arr[i,k,0]
              ENDELSE 
          ENDELSE 

          IF skip THEN BEGIN 
              wb = -1
              nwb = 0
          ENDIF ELSE BEGIN 
              cut1 = cut_arr[i,k,0]
              cut2 = cut_arr[i,k,1]
              IF rev[cut1] NE rev[cut2] THEN BEGIN 
                  
                  wb = rev[ rev[cut1]:rev[cut2+1]-1 ]
                  nwb = n_elements(wb)
              ENDIF ELSE BEGIN 
                  
                  ;;print,'/',format='(a,$)'
                  wb = where(bin_var GE cut_arr[i,k,0] AND $
                             bin_var LE cut_arr[i,k,1],nwb)
                  
              ENDELSE 
          ENDELSE 

;          wb2 = where(bin_var GE cut_arr[i,k,0] AND $
;                     bin_var LE cut_arr[i,k,1],nwb2)
;          IF nwb NE nwb2 THEN stop

          IF nwb NE 0 THEN BEGIN  
              sn_arr[i,k] = mean( signal[wb] )*sqrt(nwb)
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
          IF ((i+1) MOD 100) EQ 0 THEN print,'.',format='(a,$)'
          i=i+1
      ENDIF ELSE BEGIN 
          sn_diff[i] = 0.
      ENDELSE 
          

  ENDWHILE  

  print

  ww = where(sn_diff EQ min(sn_diff), ngood)
  help,ww
  ww = ww[0]
  print
  print,'Minimum S/N diff: ',sn_diff[ww]
  print
  print,'S/N in each bin: '
  colprint,lindgen(nbin)+1,sn_arr[ww, *]
  print
  print,'Min = ',min(bin_var)
  print,'Max = ',max(bin_var)
  print,'           bin      Nbin      <bin_var>     cut'
  print,'--------------------------------------------------------------'
  FOR k=0L, nbin-1 DO BEGIN 
      kstr=ntostr(k+1)

      w = where( bin_var GE cut_arr[ww,k,0] AND $
                 bin_var LE cut_arr[ww,k,1], nw)

      IF nw EQ 0 THEN message,'What?'

      command = 'w'+kstr+' = w'
      temp=execute(command)


      ave = mean(bin_var[w])

      cutstr = '['+ntostr(cut_arr[ww,k,0])+', '+ntostr(cut_arr[ww,k,1])+']'
      print,k, nw, ave, '       '+cutstr

  ENDFOR 

  ptime,systime(1)-tm

  return
END 
