PRO get_spectra_lcat, stripes, lcat, indir=indir, columns=columns, spec=spec, lrg=lrg, lss=lss, letter=letter, post=post, sample=sample, all=all, southern=southern, count=count
  
  np=n_params()
  IF np LT 1 THEN BEGIN 
      print,'-Syntax: get_spectra_lcat, stripes, lcat, indir=indir, columns=columns, /spec, /lrg, /lss, letter=letter, post=post, sample=sample, /all, /southern, count=count'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read the headers to find out how many files there are
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nst = n_elements(stripes)
  ntotal = 0L
  numlist = lonarr(nst)

  FOR i=0L, nst-1 DO BEGIN 

      file = spectra_name(stripes[i], indir=indir, spec=spec, lrg=lrg, $
                          lss=lss, letter=letter, post=post, sample=sample, $
                          all=all, southern=southern)
      
      hdr = headfits(file, ext=1)
      numlist[i] = sxpar(hdr,'naxis2')
      ntotal = ntotal + numlist[i]
      
  ENDFOR 

  IF NOT keyword_set(silent) THEN BEGIN
      print
      print,'Total number of objects: '+ntostr(ntotal)
  ENDIF 


  typ = $
    'specgals'+ntostr(long(1000*randomu(seed)))+'_'+ntostr(long(systime(1)))

  beg=0L
  FOR i=0L, nst-1 DO BEGIN 

      file = spectra_name(stripes[i], indir=indir, spec=spec, lrg=lrg, $
                          lss=lss, letter=letter, post=post, sample=sample, $
                          all=all, southern=southern)

      IF numlist[i] NE 0 THEN BEGIN 

          print,'Reading File: '+ntostr(file)
          t=mrdfits(file, 1, structyp=typ, silent=silent, $
                    columns=columns)
          
          IF i EQ 0 THEN BEGIN 
              lcat=replicate(t[0], ntotal)
              lcat[beg:beg+numlist[i]-1] = t
              beg = beg+numlist[i]
          ENDIF ELSE BEGIN 
              lcat[beg:beg+numlist[i]-1] = t
              beg = beg+numlist[i]
          ENDELSE 
         
      ENDIF ELSE BEGIN  
          print,'File is empty: '+ntostr(file)
      ENDELSE 

  ENDFOR 
  
  count = nTotal

return
  IF NOT fexist(name) THEN BEGIN 
      message,'File not found: '+name
  ENDIF 

  print,'Reading File: ',name

  IF np EQ 4 THEN hdr0=headfits(name)
  lcat = mrdfits(name, 1, hdr1, columns=columns)
  count=n_elements(lcat)

  return
END 
