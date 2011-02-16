PRO read_sxwthetagals, wgals, skipline=skipline, use_sxread=use_sxread, create=create
  
;  file='/sdss6/data0/wtheta/sxoutput/getSecondaryWthetaGals.tbl'
;  outf='/sdss6/data0/wtheta/sxoutput/getSecondaryWthetaGals.fit'
  file='/sdss6/data0/wtheta/sxoutput/getPrimaryWthetaGals.tbl'
  outf='/sdss6/data0/wtheta/sxoutput/getPrimaryWthetaGals.fit'

  print
  print,'Input file: ',file
  print,'Output file: ',outf
  print

  ;; if we already  made the fits file just read it in
  IF fexist(outf) AND NOT keyword_set(create) THEN BEGIN
      wgals=mrdfits(outf,1)
      return
  ENDIF 

  arrval=fltarr(5)

  objidtype = 0L
  tmpstruct1=create_struct('run', 0L, $
                           'rerun', 0, $
                           'camcol', 0, $
                           'field', 0, $
                           'id', objidtype)

  ;; This is in case sx_read2 fails to make proper tsObj tag names
  tmpstruct2=create_struct('modelcounts_0',0.,'modelcounts_1',0.,'modelcounts_2',0.,'modelcounts_3',0.,'modelcounts_4',0.,$
                           'petrocounts_0',0.,'petrocounts_1',0.,'petrocounts_2',0.,'petrocounts_3',0.,'petrocounts_4',0.,$
                           'reddening_0',0.,'reddening_1',0.,'reddening_2',0.,'reddening_3',0.,'reddening_4',0.)
  
  tmpstruct3=create_struct('ra',0d,$
                           'dec',0d)

  tmpstruct=create_struct(tmpstruct1,tmpstruct2)
  tmpstruct=create_struct(tmpstruct,tmpstruct3)

  arrval=fltarr(5)

  ;; can add these if want
  

  tsstruct=create_struct('run', 0L, $
                         'rerun', 0, $
                         'camcol', 0, $
                         'field', 0, $
                         'id', objidtype, $
                         'counts_model', arrval,$
                         'petrocounts',arrval,$
                         'reddening',arrval,$
                         'ra',0d,$
                         'dec',0d)


  IF keyword_set(use_sxread) THEN BEGIN 
      sx_read2, file, outstruct, success

      ;; Did sx_read2 make tsObj names?
      IF NOT success THEN BEGIN 
          sx_read, file, outstruct
      ENDIF 
  ENDIF ELSE BEGIN 

      ;; Find first valid line
      openr, unit, file, /get_lun
      skipline=0L
      line = ''
      
      WHILE NOT EOF(unit) DO BEGIN 
          readf, unit, line
          beg = strmid(line,0,1)
          
          IF beg EQ '#' THEN skipline = skipline+1 ELSE GOTO, numjump
          
      ENDWHILE 
      numjump:

      close,unit
      free_lun,unit
      
      print
      print,'Skipping ',skipline,' lines at beginning'
      print
      
      ;;  IF n_elements(skipline) EQ 0 THEN skipline=26
      
      rdfloatstr, file, tmpstruct, outstruct, /double, skipline=skipline, endskip=1
  ENDELSE 
  
  IF success THEN BEGIN 
      print,'Writing file: ',outf
      mwrfits2, outstruct, outf, /create, /destroy
      return
  ENDIF 

  print
  print,'Copying into tsObj names'
  nobj=n_elements(outstruct)

  outgals=replicate(tsstruct,nobj)

  outgals.run=outstruct.run
  outgals.rerun=outstruct.rerun
  outgals.camcol=outstruct.camcol
  outgals.field=outstruct.field
  outgals.id=outstruct.id

  outgals.counts_model[0] = outstruct.modelcounts_0
  outgals.counts_model[1] = outstruct.modelcounts_1
  outgals.counts_model[2] = outstruct.modelcounts_2
  outgals.counts_model[3] = outstruct.modelcounts_3
  outgals.counts_model[4] = outstruct.modelcounts_4

  outgals.petrocounts[0] = outstruct.petrocounts_0
  outgals.petrocounts[1] = outstruct.petrocounts_1
  outgals.petrocounts[2] = outstruct.petrocounts_2
  outgals.petrocounts[3] = outstruct.petrocounts_3
  outgals.petrocounts[4] = outstruct.petrocounts_4

  outgals.reddening[0] = outstruct.reddening_0
  outgals.reddening[1] = outstruct.reddening_1
  outgals.reddening[2] = outstruct.reddening_2
  outgals.reddening[3] = outstruct.reddening_3
  outgals.reddening[4] = outstruct.reddening_4

  outgals.ra=outstruct.ra
  outgals.dec=outstruct.dec

  outstruct = 0

  print,'Writing file: ',outf
  mwrfits2, outgals, outf, /create, /destroy
  

END 
