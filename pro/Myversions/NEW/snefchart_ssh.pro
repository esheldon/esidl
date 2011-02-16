PRO snefchart_ssh, file, spec, photoz

  print
  print,'Reading input data file: ',file
  format = 'D,D,A'
  readcol, $
    file, $
    ra_array, dec_array, front_array, $
    format=format

  nra = n_elements(ra_array)

  dirsep, file, dir, name
  
  result_file = dir + 'result_'+name
  print
  print,'Writing to results file: ',result_file
  print
  
  ;; output file: ra dec success ra_z dec_z separation z zerr type
  outformat = '(F,F,A,D,D,D,F,F,A)'

  f_defval = -9999.
  d_defval = -9999d

  ;; Read in the photoz file and spectro z file
  inDir = '~sne/fchart/data/'
  specFile = inDir + 'stripe82_spec.st'
  photozFile = inDir + 'stripe82_photoz.st'

  print
  print,'Reading specFile: ',specFile
  spec = read_idlstruct(specFile)

  print
  print,'Reading photozFile: ',photozFile
  photoz = read_idlstruct(photozFile)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; try matching to spec file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  tol = 10d/3600d

  tmatch_spec = -1L
  tmatch_photoz = -1L

  print
  print,'Matching to spec file to 10"'
  spherematch, $
    ra_array, dec_array, spec.ra, spec.dec, tol, $
    mradec_spec, tmatch_spec, tsep_spec, maxmatch=1

  print
  print,'Matching to photoz file to 10"'
  spherematch, $
    ra_array, dec_array, photoz.ra, photoz.dec, tol, $
    mradec_photoz, tmatch_photoz, tsep_photoz, maxmatch=1

  help,tmatch_spec,tmatch_photoz

  match_spec   = replicate(-1L, nra)
  sep_spec     = replicate(-1.0, nra)
  match_photoz = replicate(-1L, nra)
  sep_photoz   = replicate(-1.0, nra)

  IF tmatch_spec[0]   NE -1 THEN BEGIN 
      match_spec[mradec_spec]     = tmatch_spec
      sep_spec[mradec_spec]       = tsep_spec/3600d
  ENDIF 
  IF tmatch_photoz[0] NE -1 THEN BEGIN 
      match_photoz[mradec_photoz] = tmatch_photoz
      sep_photoz[mradec_photoz]   = tsep_photoz/3600d
  ENDIF 

  openw, lun, result_file, /get_lun
;  lun = -1

  ;; For marking objects

  FOR i=0L, nra-1 DO BEGIN 

      front = front_array[i]

      ;; 6' diameter
      jpegFchart6 = front + '_fchart_rad6.jpg'
      psFchart6   = front + '_fchart_rad6.ps'
      ;; 2' diameter
      jpegFchart2 = front + '_fchart_rad2.jpg'
      psFchart2  = front + '_fchart_rad2.ps'

      ra = ra_array[i]
      dec = dec_array[i]
      
;      extra = extra_array[i]
;      extdec = extra_dec[i]

      find_radec, ra, dec, run, camcol, field

      IF run[0] NE -1 THEN BEGIN 

          ;; use highest run number
          nFound = n_elements(run)

          runuse = nFound-1

          delvarx, circ_rad, linestyle, circ_color, addtitle

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Did we find a spec match?
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF match_spec[i] NE -1 THEN BEGIN 


              matchra = spec[match_spec[i]].ra
              matchdec = spec[match_spec[i]].dec

              linestyle = [0,1]

              printf,lun,$
                ra,dec,'   Yes',$
                spec[match_spec[i]].ra, $
                spec[match_spec[i]].dec, $
                sep_spec[i], $
                spec[match_spec[i]].z, $
                spec[match_spec[i]].z_err, $
                '  spec',$
                format=outformat

              addtitle = $
                'specz = '+$
                ntostr(spec[match_spec[i]].z)+$
                !csym.plusminus + $
                ntostr(spec[match_spec[i]].z_err)

          ENDIF ELSE IF match_photoz[i] NE -1 THEN BEGIN 
              matchra = photoz[match_photoz[i]].ra
              matchdec = photoz[match_photoz[i]].dec

              linestyle = [0,1]

              printf,lun,$
                ra,dec,'   Yes',$
                photoz[match_photoz[i]].ra, $
                photoz[match_photoz[i]].dec, $
                sep_photoz[i], $
                photoz[match_photoz[i]].z, $
                photoz[match_photoz[i]].zerr, $
                '  photoz',$
                format=outformat

              addtitle = $
                'photoz = '+$
                ntostr(photoz[match_photoz[i]].z)+$
                !csym.plusminus + $
                ntostr(photoz[match_photoz[i]].zerr)

          ENDIF ELSE BEGIN 
              printf,lun,$
                ra,dec,'   Yes',$
                d_defval,d_defval,d_defval,f_defval,f_defval,'  None',$
                format=outformat

              delvarx, matchra, matchdec, extra_markstruct              
          ENDELSE 

          ;; Make the finding charts
          markStruct = {radius:0.3}
          IF n_elements(matchra) NE 0 THEN BEGIN 
              extra_markStruct = {ra:matchra,dec:matchdec,$
                                  radius:0.4,type:"circle",linestyle:1}
          ENDIF 
          print
          print,'Creating 6 arcmin jpeg fchart: ',jpegFchart6
          rgbFchart, $
            ra, dec, runuse = runuse, $
            markStruct=markStruct, $
            extra_markStruct=extra_markStruct, $
            /directions, order=1, $
            radius = 3, addtitle=addtitle, $
            jpegFchart = jpegFchart6, /noDisplay

          print
          print,'Creating 6 arcmin ps fchart: ',psFchart6
          find_radec_fchart, $
            ra, dec, runuse=runuse, $
            markStruct=markStruct, $
            extra_markStruct=extra_markStruct, $
            /directions, order=1, $
            radius = 3, addtitle=addtitle, $
            psFile = psFchart6

          markStruct = {radius:0.3/3.0}
          IF n_elements(matchra) NE 0 THEN BEGIN 
              extra_markStruct = {ra:matchra,dec:matchdec,$
                                  radius:0.4/3.0,type:"circle",linestyle:1}
          ENDIF 

          print
          print,'Creating 2 arcmin jpeg fchart: ',jpegFchart2
          rgbFchart, $
            ra, dec, runuse = runuse, $
            markStruct=markStruct, $
            extra_markStruct=extra_markStruct, $
            /directions, order=1, $
            radius = 1, addtitle=addtitle, $
            jpegFchart = jpegFchart2, /noDisplay
          print
          print,'Creating 2 arcmin ps fchart: ',psFchart2
          find_radec_fchart, $
            ra, dec, runuse=runuse, $
            markStruct=markStruct, $
            extra_markStruct=extra_markStruct, $
            /directions, order=1, $
            radius = 1, addtitle=addtitle, $
            psFile = psFchart2

      ENDIF ELSE BEGIN 
          printf,lun,$
            ra,dec,'   No',$
            d_defval,d_defval,d_defval,f_defval,f_defval,'  None',$
            format=outformat
      ENDELSE 

  ENDFOR 

  free_lun, lun

END 
