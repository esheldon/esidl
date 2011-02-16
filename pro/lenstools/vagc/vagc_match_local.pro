FUNCTION vagc_match_local, vagc_cat, taglist_in=taglist, status=status

  ;; Match the vagc or lss catalog with their different id numbers to the
  ;; local standard sdss imaging.  Useful for getting the local rerun number
  ;; and id, for example for looking at atlas images.

  status = 1
  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: match_struct = vagc_match_local(vagc_cat, taglist=, status=)'
      return,-1
  ENDIF 

  delvarx, match_struct
  tm = systime(1)

  ;; ra/dec is used for matching, so it is added to all taglists
  IF n_elements(taglist_in) EQ 0 THEN BEGIN 
      taglist =  ['run','rerun','camcol','field','id','ra','dec']
  ENDIF ELSE BEGIN 
      taglist = strlowcase(taglist_in)
  ENDELSE 

  wra = where(taglist EQ 'ra', nra)
  IF nra EQ 0 THEN taglist = [taglist, 'ra']
  wdec = where(taglist EQ 'dec', ndec)
  IF ndec EQ 0 THEN taglist = [taglist, 'dec']

  ;; Get runs and loop over them
  vruns = vagc_cat.run
  IF n_elements(vruns) EQ 1 THEN vruns = [vruns]
  hruns = histogram(vruns, rev=runrev)

  ;; only match decent detections, children, galaxies
  wstring = $
    'where(lnew.nchild eq 0 and '+$
    '      ( (lnew.petrocounts[2]-lnew.reddening[2]) lt 21) and ' + $
    '      lnew.objc_type eq 3)'

  ;; pointers to the matches
  nvagc = n_elements(vagc_cat)
  ptrlist = ptrarr(nvagc)
  numlist = lonarr(nvagc)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Loop over runs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nrunbin = n_elements(hruns)
  FOR irun=0L, nrunbin-1 DO BEGIN 

      IF runrev[irun] NE runrev[irun+1] THEN BEGIN 
          wrun=runrev[ runrev[irun]:runrev[irun+1]-1 ]

          run = vagc_cat[wrun[0]].run
          rstr = run2string(run)
          print
          print,'-------------------------------------'
          print,'Processing Run = '+rstr
          print,'Nobj = ',n_elements(wrun)
          print,'-------------------------------------'

          ;; Get the latest rerun available for this run
          wst = where(!run_status.run EQ run AND $
                      !run_status.tsobj_photo_v NE -1, nwst)

          IF nwst NE 0 THEN BEGIN 
              ;; get latest rerun
              wst2 = where(!run_status[wst].rerun EQ $
                           max(!run_status[wst].rerun))
              rerun = !run_status[wst[wst2]].rerun

              vcamcols = vagc_cat[wrun].camcol
              IF n_elements(vcamcols) EQ 1 THEN vcamcols=[vcamcols]
              hcam=histogram(vcamcols, rev=camrev)
              ncambin=n_elements(hcam)

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Loop over camcols
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              FOR icam=0L, ncambin-1 DO BEGIN 

                  IF camrev[icam] NE camrev[icam+1] THEN BEGIN 

                      ;; get camcol indices
                      wcam=camrev[ camrev[icam]:camrev[icam+1]-1 ]

                      ;; get absolute indices
                      wcam = wrun[wcam]

                      camcol = vagc_cat[wcam[0]].camcol
                      cstr = ntostr(camcol)
                      print,'  Camcol = '+cstr
                      ;; Now get the fields for this camcol
                      vfields = vagc_cat[wcam].field
                      IF n_elements(vfields) EQ 1 THEN vfields=[vfields]
                      hfields = histogram(vfields, rev=fieldrev)
                      nfieldbin = n_elements(hfields)
                      FOR ifield=0L, nfieldbin-1 DO BEGIN 
                          
                          IF fieldrev[ifield] NE fieldrev[ifield+1] THEN BEGIN 
                              
                              ;; get field indices
                              wfield = $
                                fieldrev[fieldrev[ifield]:fieldrev[ifield+1]-1]
                              ;; get absolute indices
                              wfield = wcam[wfield]
                              nField = n_elements(wfield)

                              field = vagc_cat[wfield[0]].field
                              fstr = field2string(field)
                              
                              read_tsobj,$
                                [run,rerun,camcol], struct, start=field, $
                                taglist=taglist, $
                                tsObjStr=tsObjStr,$
                                wstring=wstring, $
                                verbose=0, ex_struct={vagc_index:0L}

                              nstruct = n_elements(struct)
                              print,'    Run: '+rstr+$
                                ' Field: '+fstr+$
                                ' Camcol: '+cstr+$
                                ' Nvagc: '+ntostr(nfield)+$
                                ' Nobj: '+ntostr(nstruct)

                              ;; If any objects found, then match
                              IF nstruct NE 0 THEN BEGIN 

                                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                  ;; match closest
                                  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                  close_match_radec, $
                                    vagc_cat[wfield].ra, vagc_cat[wfield].dec,$
                                    struct.ra, struct.dec, $
                                    mvag, mstr, 1.0/3600d, 1, /silent

                                  ;; Any Matches?
                                  IF mvag[0] NE -1 THEN BEGIN 
                                      nmatch=n_elements(mvag)
                                      print,'      Nmatches = '+$
                                        ntostr(nmatch)+'/'+ntostr(nfield)

                                      struct=struct[mstr]
                                      struct.vagc_index=wfield[mvag]

                                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                      ;; Make copies of the matches
                                      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                      FOR im=0L, nmatch-1 DO BEGIN
                                          mm = wfield[mvag[im]]
                                          ptrlist[mm] = $
                                            ptr_new(struct[im])
                                          numlist[mm] = 1
                                      ENDFOR 

                                      ;; copy of base struct for combining
                                      ;; the pointer list
                                      IF n_elements(base_struct) EQ 0 THEN $
                                        base_struct = struct[0]

                                      struct = 0
                                  ENDIF ELSE BEGIN 
                                      print,'No matches found'
                                  ENDELSE 
                              ENDIF 

                          ENDIF 
                          
                      ENDFOR ;; loop over fields for this run/camcol
                  ENDIF ;; will include all camcols between min and max
              ENDFOR ;; loop over camcols for this run

          ENDIF ELSE BEGIN ;; Do we have this run on disk?
              print,'********************************************'
              print,'Run '+ntostr(run)+' not found'
              print,'********************************************'
          ENDELSE 

      ENDIF ;; histogram has all runs between min and max

  ENDFOR 

  ntot = total_int(numlist)
  IF ntot GT 0 THEN BEGIN 
      match_struct = combine_ptrlist(ptrlist)
      status = 0
  ENDIF ELSE BEGIN 
      match_struct = -1
  ENDELSE 
  ptime,systime(1)-tm

  return,match_struct

END 
