PRO lumdensboot_chunk, lenses, random, nchunk, nresamp, lumdens, lumdenserr

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: lumdensboot_chunk, lenses, random, nchunk, nresamp, lumdens, lumdenserr'
      return
  ENDIF 

  nlenses = n_elements(lenses)
  nbin = n_elements(lenses[0].rsum)

  s=sort(lenses.ra)
  tind = lindgen(nlenses)

  chunksize = nlenses/nchunk

  chunkmeans = fltarr(nchunk, nbin)
;  rchunkmeans = fltarr(nchunk, nbin)


  print
  print,'Number of objects: ',nlenses
  print,'Number of chunks: ',nchunk
  print,'Building chunks'
  FOR ch=0L, nchunk-1 DO BEGIN 

      begind=ch*chunksize
      IF ch NE nchunk-1 THEN BEGIN 
          endind = (ch+1)*chunksize-1
      ENDIF ELSE BEGIN 
          endind = nlenses-1
      ENDELSE 

      sind = s[begind:endind]
      match2rand, lenses[ sind ], random, wrand

      wsum = 0.
      ld = 0.
      rwsum = 0.
      rld = 0.

      FOR bin=0L, nbin-1 DO BEGIN 

          wsum = total( lenses[sind].wsum[bin] )
          ld =  total( lenses[sind].lsum[bin] )

          rwsum = total( random[wrand].wsum[bin] )
          rld  = total( random[wrand].lsum[bin] )

          chunkmeans[ch, bin] = ld/wsum - rld/rwsum
;          rchunkmeans[ch, bin] = rld/rwsum

      ENDFOR 

      print,'.',format='(a,$)'
  ENDFOR 

  ;forprint,begind,endind

  lumdensamp = fltarr(nresamp, nbin)
  ;rlumdensamp = lumdensamp

  print
  print,'Building bootstrap samples'
  seed = long(systime(1))
  FOR samp=0L, nresamp-1 DO BEGIN 

      ;; choose with replacement random set of chunks
      useind = lonarr(nlenses)
      sind = round( (nchunk-1)*randomu(seed, nchunk) )

      FOR bin=0L, nbin-1 DO BEGIN 
          lumdensamp[samp,bin] = mean_check(chunkmeans[sind, bin])
      ENDFOR 
      ;print,'.',format='(a,$)'
  ENDFOR 

  print
  print,'bootstrapping'
  bootstrap, lumdensamp, lumdens, lumdenserr


  lumdensamp=0
  rlumdensamp=0

return
END 

