PRO get_scat_readcuts, stripe, clr, rows, $
                       nzCuts=nzCuts, $
                       pixelMaskFlags=pixelMaskFlags, rlrgMask=rlrgMask, $
                       nrows=nrows,hirata=hirata, silent=silent

  ;; defaults
  delvarx, rows
  nrows = 0L

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Deal with pixel mask flags
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(pixelMaskFlags) THEN BEGIN 
      get_srcgal_pixelMaskFile, stripe, clr, maskStruct, $
        status=pstatus, hirata=hirata, rlrgMask=rlrgMask, silent=silent
      IF pstatus EQ 0 THEN BEGIN 
          rows = where(maskStruct.maskFlags EQ 0, nrows)
      ENDIF 
  ENDIF  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Cuts on the photozs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(nzCuts) THEN BEGIN 
      get_nz_photoz,stripe,clr,nzstruct,silent=silent,$
        status=nstatus,hirata=hirata
      
      IF nstatus EQ 0 THEN BEGIN 
          IF n_elements(maskStruct) NE 0 THEN BEGIN 
              rows = $
                where(maskStruct[nzstruct.useind].maskFlags $
                      EQ 0, nrows)

              rows = nzStruct.useInd[rows]
          ENDIF ELSE BEGIN 
              rows = nzstruct.useind
              nrows = n_elements(rows)
          ENDELSE 
      ENDIF  
  ENDIF 
    
END 

PRO get_scat, stripes, clr, scat, indir=indir, columns=columns, hirata=hirata, hudson=hudson, nzCuts=nzCuts, pixelMaskFlags=pixelMaskFlags, rlrgMask=rlrgMask, silent=silent, count=count, file=file

  np=n_params()
  IF np LT 3 THEN BEGIN 
      print,'-Syntax: get_scat, stripes, clr, scat, indir=indir, columns=columns, /hirata, /hudson, /nzCuts, /pixelMaskFlags, /rlrgMask, /silent, count=count, file=file'
      print
      print,'default is /hirata'
      print,' - stripes can be scalar or array. e.g. [9,10] or 9'
      print,' - clr can be an array for "combined" catalogs: e.g. [1,2,3] for gri'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; hirata is default now
  IF n_elements(hirata) EQ 0 THEN hirata=1

  IF keyword_set(hirata) THEN hirstr = '_h' ELSE hirstr = ''
  IF keyword_set(hudson) THEN hudstr = '_hudson' ELSE hudstr=''

  clr_string = clrarr2string(clr)
  stripe_string = stripearr2string(stripes)

  fend = '_srcgal_'+clr_string+hirstr+hudstr+'.fit'

  ;; Where is the catalog kept?

  IF n_elements(indir) EQ 0 THEN BEGIN
      combdir = sdssidl_config('SHAPECORR_DIR') + 'combined/srcgal/'
  ENDIF ELSE BEGIN
      combdir = indir
  ENDELSE 

  delvarx, scat

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Read the headers to find out how many objects there are
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nst = n_elements(stripes)
  ntotal = 0L
  numlist = lonarr(nst)

  IF keyword_set(nzCuts) OR keyword_set(pixelMaskFlags) THEN BEGIN 

      FOR i=0L, nst-1 DO BEGIN 

          get_scat_readcuts, stripes[i], clr, blah, $
            nzCuts=nzCuts, $
            pixelMaskFlags=pixelMaskFlags, rlrgMask=rlrgMask, $
            nrows=nrows,hirata=hirata, /silent

          numlist[i] = nrows
          ntotal = ntotal + numlist[i]

      ENDFOR 

  ENDIF ELSE BEGIN 
      FOR i=0L, nst-1 DO BEGIN 

          file = combdir + 'stripe'+stripe2string(stripes[i]) + fend
          
          hdr = headfits(file, ext=1)
          numlist[i] = sxpar(hdr,'naxis2')
          ntotal = ntotal + numlist[i]

      ENDFOR 
  ENDELSE 

  IF NOT keyword_set(silent) THEN BEGIN
      print
      print,'Total number of objects: '+ntostr(ntotal)
  ENDIF 

  typ = 'srcgals'+ntostr(long(1000*randomu(seed)))+'_'+ntostr(long(systime(1)))

  beg =0L
  FOR i=0L, nst-1 DO BEGIN 
      
      file = combdir + 'stripe'+stripe2string(stripes[i]) + fend
      IF numlist[i] NE 0 THEN BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Should we make cuts?
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF NOT keyword_set(silent) THEN $
            print,'-------------------------------------------'

          IF keyword_set(nzCuts) OR keyword_set(pixelMaskFlags) THEN BEGIN 

              get_scat_readcuts, stripes[i], clr, rows, $
                nzCuts=nzCuts, $
                pixelMaskFlags=pixelMaskFlags, rlrgMask=rlrgMask, $
                nrows=nrows,hirata=hirata

          ENDIF 

          print
          print,'Reading File: '+ntostr(file)
          t=mrdfits(file, 1, structyp=typ, silent=silent, $
                    columns=columns, rows=rows)


          IF i EQ 0 THEN BEGIN 
              scat=replicate(t[0], ntotal)
              scat[beg:beg+numlist[i]-1] = t
              beg = beg+numlist[i]
          ENDIF ELSE BEGIN 
              scat[beg:beg+numlist[i]-1] = t
              beg = beg+numlist[i]
          ENDELSE 

          ;; Clear these variables each time
          delvarx, prows, nrows, rows, t 
          
      ENDIF ELSE BEGIN  
          IF keyword_set(nzcuts) THEN BEGIN 
              print,'File is empty or no nzstruct: '+ntostr(file)
          ENDIF ELSE BEGIN 
              print,'File is empty: '+ntostr(file)
          ENDELSE 
      ENDELSE 
  ENDFOR 

  count = ntotal

  return

END 
