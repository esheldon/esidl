FUNCTION sdss_objid_range, name

  CASE name OF
      'sky_version': return, [0L, 2L^4-1L] ; 0..15
      'run': return, [0L, 2L^16-1L] ; 0..65535
      'rerun': return, [0L, 2L^11-1L] ; 0..2047
      'camcol': return, [1L, 6L] 
      'field': return, [0L, 2L^12-1L] ; 0..4095
      'id': return, [0L, 2L^16-1L] ; 0..65535
      ELSE: message,'unknown type: '+name
  ENDCASE 

END 

PRO sdss_objid_rangerr, name, range

  rangeStr = '['+ntostr(range[0])+', '+ntostr(range[1])+']'
  message,name + ' must be in range '+rangeStr, /inf

END 

FUNCTION sdss_objid_range_check, run, rerun, camcol, field, id, $
                                 sky_version, first_field

  retval = 0
  nrun = n_elements(run)

  ;; Run
  name = 'run'
  range = sdss_objid_range(name)
  mn = min(run, max=mx)
  IF mn LT range[0] OR mx GT range[1] THEN BEGIN 
      sdss_objid_rangerr, name, range
      retval = retval + 1
  ENDIF 

  ;; Rerun
  name = 'rerun'
  IF n_elements(rerun) NE nrun THEN BEGIN 
      message,name+' is shorter than run',/inf
      retval = retval + 1
  ENDIF 

  range = sdss_objid_range(name)
  mn = min(rerun, max=mx)
  IF mn LT range[0] OR mx GT range[1] THEN BEGIN 
      sdss_objid_rangerr, name, range
      retval = retval + 1
  ENDIF 

  ;; Camcol
  name = 'camcol'
  IF n_elements(camcol) NE nrun THEN BEGIN 
      message,name+' is shorter than run',/inf
      retval = retval + 1
  ENDIF 

  range = sdss_objid_range(name)
  mn = min(camcol, max=mx)
  IF mn LT range[0] OR mx GT range[1] THEN BEGIN 
      sdss_objid_rangerr, name, range
      retval = retval + 1
  ENDIF 

  ;; Field
  name = 'field'
  IF n_elements(field) NE nrun THEN BEGIN 
      message,name+' is shorter than run',/inf
      retval = retval + 1
  ENDIF 

  range = sdss_objid_range(name)
  mn = min(field, max=mx)
  IF mn LT range[0] OR mx GT range[1] THEN BEGIN 
      sdss_objid_rangerr, name, range
      retval = retval + 1
  ENDIF 

  ;; id
  name = 'id'
  IF n_elements(id) NE nrun THEN BEGIN 
      message,name+' is shorter than run',/inf
      retval = retval + 1
  ENDIF 

  range = sdss_objid_range(name)
  mn = min(id, max=mx)
  IF mn LT range[0] OR mx GT range[1] THEN BEGIN 
      sdss_objid_rangerr, name, range
      retval = retval + 1
  ENDIF 

  ;; Sky Version
  name = 'sky_version'
  IF n_elements(sky_version) NE nrun THEN BEGIN 
      message,name+' is shorter than run',/inf
      retval = retval + 1
  ENDIF 

  range = sdss_objid_range(name)
  mn = min(sky_version, max=mx)
  IF mn LT range[0] OR mx GT range[1] THEN BEGIN 
      sdss_objid_rangerr, name, range
      retval = retval + 1
  ENDIF 

  ;; First field. Only check length
  name = 'first_field'
  IF n_elements(first_field) NE nrun THEN BEGIN 
      message,name+' is shorter than run',/inf
      retval = retval + 1
  ENDIF 


  return, retval

END 

FUNCTION sdss_objid, run, rerun, camcol, field, id, $
                     first_field=first_field, sky_version=sky_version,$
                     verbose=verbose

  IF n_params() LT 5 THEN BEGIN 
      print,'-Syntax: objID = sdss_objid(run, rerun, camcol, field, id, '
      print,'                            first_field=0, sky_version=1'
      print,'                            /verbose)'
      return,-1LL
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; 
  ;; Takes run, rerun, camcol, field, id [and sky_version and 
  ;; first_field] to create a 64-bit integer that should match the 
  ;; unique id given in the SDSS SQL data base.
  ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  skyDef = 15
  first_fieldDef = 0

  nrun = n_elements(run)

  IF n_elements(sky_version) EQ 0 THEN BEGIN
      IF keyword_set(verbose) THEN BEGIN
          print,'Assuming that sky_version is '+ntostr(skyDef)
      ENDIF
      sky_version = replicate(skyDef, nrun)
  ENDIF


  IF n_elements(first_field) EQ 0 THEN BEGIN
      IF keyword_set(verbose) THEN BEGIN
          print,'Assuming that first_field is '+ntostr(first_fieldDef)
      ENDIF
      first_field = replicate(first_fieldDef, nrun)
  ENDIF


  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some range checking
  ;;;;;;;;;;;;;;;;;;;;;;;;;;

  bad = sdss_objid_range_check(run, rerun, camcol, field, id, $
                               sky_version, first_field)
  IF bad NE 0 THEN return, -1LL

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Ok, copy in the bits with some range checking
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  objID = long64(sky_version)

  objID = ishft(objID,11)
  objID = objID + long64(rerun)

  objID = ishft(objID,16)
  objID = objID + long64(run)

  objID = ishft(objID,3)
  objID = objID + long64(camcol)

  objID = ishft(objID,1)
  objID = objID + long64(first_field > 0 < 1)

  objID = ishft(objID,12)
  objID = objID + long64(field)

  objID = ishft(objID,16)
  objID = objID + long64(id)

  return, objID

END


