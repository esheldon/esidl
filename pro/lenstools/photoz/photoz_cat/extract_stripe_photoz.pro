;; For use with the Budavari photozs
;;
;; This is to just extract the ones of a particular stripe
;; I used this to produce the stripped-down stripe 82 catalog

PRO extract_stripe_photoz_select, struct, useind

  tzmin = 0.02                  ;0.02 problems at low z
  tzmax = 0.8                   ;because a default is around 1
  qualmin = 0
  qualmax = 12
  zerrmin = 0.01
  zerrmax = 0.4

  useind=where(struct.z GT tzmin AND $
               struct.z LT tzmax AND $
               struct.z NE 0.0 AND $
               struct.zerr GT zerrmin AND $
               struct.zerr LT zerrmax AND $
               struct.quality GT qualmin AND $
               struct.quality LT qualmax, nuse)


END 

PRO extract_stripe_photoz, stripe

  photoz_dir = sdssidl_config('photoz_dir')

  stripestr = stripe2string(stripe)
  cd,photoz_dir
  photozFiles = file_search('tsObj_ascii_*.res')

  outFile = photoz_dir + 'stripe'+stripestr+'_cut.fit'
  print
  print,'Will write to file: ',outFile


  ;; extract the runs/reruns and write to disk
  nf = n_elements(photozFiles)

  runs = lonarr(nf)
  reruns = lonarr(nf)
  runsreruns = strarr(nf)

  FOR i=0L, nf-1 DO BEGIN 

      f = photozFiles[i]
      front = ( strsplit(f, '.', /extract) )[0]
      
      tmp = strmid(front, 12)
      
      tmp = strsplit(tmp, '_', /extract)

      runrerun = run2string(tmp[0]) + '-' + tmp[1]

      runsreruns[i] = runrerun
      runs[i] = long(tmp[0])
      reruns[i] = long(tmp[1])

;      print,runrerun
  ENDFOR 

  ;; will also sort
  rmd = rem_dup(runsreruns)

  uruns = runs[rmd]
  ureruns = reruns[rmd]
  urunsreruns = runsreruns[rmd]

;  forprint,uruns,ureruns,'   '+urunsreruns

  ;; now match to stripe
  wstripe = where(!run_status.stripe EQ stripe AND $
                  !run_status.run NE 94 AND $
                  !run_status.run NE 125 AND $
                  !run_status.tsobj_photo_v GE 5.4)


  match_multi, !run_status[wstripe].run, runs, mruns
  nmatch = n_elements(mruns)

  ;; now get all the files with these runs
  ;;colprint,runsreruns[mruns],'   '+photozFiles[mruns]

  str = create_struct('ra',0d,$
                      'dec',0d,$
                      'z', 0.0, $
                      'zerr', 0.0)

  ptrlist = ptrarr(nmatch)
  numlist = lonarr(nmatch)
  ntot = 0L

  FOR i=0L, nmatch-1 DO BEGIN 


      f=photozFiles[mruns[i]]

      print
      read_photoz_ascii, f, tmp
      ntmp = n_elements(tmp)

      extract_stripe_photoz_select, tmp, useind
      
      IF useind[0] NE -1 THEN BEGIN 
          nkeep = n_elements(useind)

          print,'Kept '+ntostr(nkeep)+'/'+ntostr(ntmp)
          struct = replicate(str, nkeep)
          
          struct.ra   = tmp[useind].ra
          struct.dec  = tmp[useind].dec
          struct.z    = tmp[useind].z
          struct.zerr = tmp[useind].zerr

          delvarx, tmp

          ptrlist[i] = ptr_new(struct, /no_copy)
          numlist[i] = nkeep
          ntot = ntot + nkeep

      ENDIF ELSE print,'None passed cuts'
  ENDFOR 

  comb_ptrstruct, str, ptrlist, numlist, struct

  ;; remove duplicates
  ;; ran out of memory here...

  print
  print,'Sorting by ra'
  s=sort(struct.ra)
  struct = struct[s]
  
  rmclose_radec,struct.ra, struct.dec, keep
  struct = struct[keep]
  nkeep = n_elements(keep)
  print,'Kept '+ntostr(nkeep)+'/'+ntostr(ntot)

  print
  print,'Writing to file: ',outFile
  mwrfits2, struct, outFile,/create, /destroy

  return
END 
