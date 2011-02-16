FUNCTION sdss_specid_range, name

  CASE name OF
      'plate': return, [0L, 2L^16-1] ; 0..65535
      'mjd': return, [0L, 2L^16-1] ; 0..65535
      'fiberid': return, [1L, 640L] ; 1..640
      'type': return, [0L, 2^6L-1] ; 0..63
      'lineOrIndexOrZ': return, [0L, 2L^16-1] ; 0..65535
      ELSE: message,'unknown type: '+name
  ENDCASE 

END 

PRO sdss_specid_rangerr, name, range

  rangeStr = '['+ntostr(range[0])+', '+ntostr(range[1])+']'
  message,name + ' must be in range '+rangeStr, /inf

END 

FUNCTION sdss_specid_range_check, plate, mjd, fiberid, type, lineOrIndexOrZ

  retval = 0
  nplate = n_elements(plate)

  ;; Plate
  name = 'plate'
  range = sdss_specid_range(name)
  mn = min(plate, max=mx)
  IF mn LT range[0] OR mx GT range[1] THEN BEGIN 
      sdss_specid_rangerr, name, range
      retval = retval + 1
  ENDIF 

  ;; MJD
  
  name = 'mjd'
  IF n_elements(mjd) NE nplate THEN BEGIN 
      message,name+' is shorter than plate',/inf
      retval = retval + 1
  ENDIF 

  range = sdss_specid_range(name)
  mn = min(mjd, max=mx)
  IF mn LT range[0] OR mx GT range[1] THEN BEGIN 
      sdss_specid_rangerr, name, range
      retval = retval + 1
  ENDIF 

  ;; fiberID
  name = 'fiberid'
  IF n_elements(fiberID) NE nplate THEN BEGIN 
      message,name+' is shorter than plate',/inf
      retval = retval + 1
  ENDIF 

  range = sdss_specid_range(name)
  mn = min(fiberid, max=mx)
  IF mn LT range[0] OR mx GT range[1] THEN BEGIN 
      sdss_specid_rangerr, name, range
      retval = retval + 1
  ENDIF 

  ;; type
  name = 'type'
  IF n_elements(type) NE nplate THEN BEGIN 
      message,name+' is shorter than plate',/inf
      retval = retval + 1
  ENDIF 

  range = sdss_specid_range(name)
  mn = min(type, max=mx)
  IF mn LT range[0] OR mx GT range[1] THEN BEGIN 
      sdss_specid_rangerr, name, range
      retval = retval + 1
  ENDIF 

  ;; line...
  name = 'lineOrIndexOrZ'
  IF n_elements(lineOrIndexOrZ) NE nplate THEN BEGIN 
      message,name+' is shorter than plate',/inf
      retval = retval + 1
  ENDIF 

  range = sdss_specid_range(name)
  mn = min(lineOrIndexOrZ, max=mx)
  IF mn LT range[0] OR mx GT range[1] THEN BEGIN 
      sdss_specid_rangerr, name, range
      retval = retval + 1
  ENDIF 

  return, retval

END 

FUNCTION sdss_specid, plate, mjd, fiberid, $
                      type=type, lineOrIndexOrZ=lineOrIndexOrZ, $
                      verbose=verbose

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: specID = sdss_specid(plate, mjd, fiberID, type=, '
      print,'                              lineOrIndexOrZ=, /verbose)'
      return,-1LL
  ENDIF 

  typeDef = 0
  lineOrIndexOrZDef = 0

  nplate = n_elements(plate)

  IF n_elements(type) EQ 0 THEN BEGIN
      IF keyword_set(verbose) THEN BEGIN
          print,'Assuming that type is '+ntostr(typeDef)
      ENDIF
      type = replicate(typeDef, nplate)
  ENDIF

  IF n_elements(lineOrIndexOrZ) EQ 0 THEN BEGIN
      IF keyword_set(verbose) THEN BEGIN
          print,'Assuming that lineOrIndexOrZ is '+ntostr(lineOrIndexOrZDef)
      ENDIF
      lineOrIndexOrZ = replicate(lineOrIndexOrZDef, nplate)
  ENDIF

  ;; Some range checking
  bad = sdss_specid_range_check(plate, mjd, fiberid, type, lineOrIndexOrZ)
  IF bad NE 0 THEN return,-1LL
  
  specID = long64(plate)

  specID = ishft(specID, 16)
  specID = specID + long64(mjd)

  specID = ishft(specID, 10)
  specID = specID + long64(fiberID)

  specID = ishft(specID, 6)
  specID = specID + long64(type)

  specID = ishft(specID, 16)
  specID = specID + long64(lineOrIndexOrZ)


  return,specID


END 
