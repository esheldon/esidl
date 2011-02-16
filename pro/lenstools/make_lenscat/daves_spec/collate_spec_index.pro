PRO collate_spec_index

  ;; take spec index and collate with the adatc file, outputting
  ;; a file for each run.  We will keep spectra even if they
  ;; don't match.

  ;;dir="/net/cheops1/data0/esheldon/spec_index/"
  dir=sdssidl_config('SHAPECORR_DIR')+'spec_index/'
  infile=dir+"spec-index.fits"

  sp = mrdfits(infile,1)

  ;; get unique objects
  ;; for 1d_22 this removed all duplicates; i.e. there were
  ;; no duplicates that didn't also have same superid
  ind = photoid(sp.run,sp.rerun,sp.camcol,sp.field,sp.id)
  spg = rem_dup(ind)

  ;; runs in this file
  uspg = rem_dup(sp[spg].run)
  uspg = spg[uspg]
  runs = sp[uspg].run
  nruns = n_elements(runs)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; We must pick a version of photo, so tags will be the same.
  ;; Also, demand adatc exist
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  w=where(!RUN_STATUS.ADATC_PHOTO_V GE 5.3 AND $
          !RUN_STATUS.ADATC_PHOTO_V LT 5.4, n53)

  make_runstatus_struct,rs
  rs.adatc_exist = 'Y'
  runstatus_select, rs, good_runst,input_index=w

  ;; tags, exstruct
  defint = -9999
  deflon = -9999L
  defflt = -9999.
  ex_struct = create_struct('mjd', deflon, $
                            'plate', defint, $
                            'fiber', defint, $
                            'z', defflt, $
                            'z_err', defflt, $
                            'z_conf',defflt,$
                            'z_status', defint, $
                            'z_warning', deflon, $
                            'spec_cln', defint, $
                            'vel_dis', defint, $
                            'vel_diserr', defint, $
                            'eclass', defflt, $
                            'eclass_old', defflt, $
                            'vers_1d', defflt)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; create the output structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fetch_dir, !RUN_STATUS[good_runst[0]].run, $
             1, $
             !RUN_STATUS[good_runst[0]].rerun, corrdir=corrdir
  fetch_file_list, corrdir, tfiles, tfnums, front='adatc'
  tmp = mrdfits(tfiles[0], 1, /silent)
  typ = 'blah1'+ntostr(long(systime(1)))

  trstr = create_struct(name=typ, tmp[0], ex_struct)
  zero_struct, trstr

  ;; radec matching
  radec_ep = 2d/3600d
  allow = 1

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; should never have more than 1000 fields, 6 camcols, 10 reruns
  ;; in a single run
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  maxf = 1000L
  maxcol = 6L
  maxrerun = 10L

  ntot = maxf*maxcol*maxrerun
  ptrlist = ptrarr(ntot)
  numlist = lonarr(ntot)
  ntotal = 0L

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; loop over spec runs
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR ir=0L, nruns-1 DO BEGIN 
      
      run = runs[ir]

      rstr = ntostr(run)

      numlist[*] = 0L
      ntotal = 0L
      index = 0L

      print
      print,'Processing run: ',ntostr( run )
      inrun = where(sp[spg].run EQ run, nrsp)
      inrun = spg[inrun]

      ;; index array for this run.  
      runnomatch = replicate(1b, nrsp)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; create the output structure for these objects
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      tmprunstruct = replicate(trstr, nrsp)

      tmprunstruct.mjd = sp[inrun].mjd
      tmprunstruct.plate = sp[inrun].plate
      tmprunstruct.fiber = sp[inrun].fiber

      ;; redshift info
      tmprunstruct.z = sp[inrun].z
      tmprunstruct.z_err = sp[inrun].z_err
      tmprunstruct.z_conf = sp[inrun].z_conf
      tmprunstruct.z_status = sp[inrun].z_status
      tmprunstruct.z_warning = sp[inrun].z_warning

      ;; spec classification
      tmprunstruct.spec_cln = sp[inrun].spec_cln

      tmprunstruct.vel_dis = sp[inrun].vel_dis
      tmprunstruct.vel_diserr = sp[inrun].vel_diserr

      tmprunstruct.eclass  = sp[inrun].eclass
      tmprunstruct.eclass_old = sp[inrun].eclass_old

      ;; copy primtarget, sectarget from
      ;; spec file
      tmprunstruct.primtarget=sp[inrun].primtarget
      tmprunstruct.sectarget=sp[inrun].sectarget
                              
      ;; version of spectro 1d
      tmprunstruct.vers_1d = sp[inrun].vers_1d

      ;; This is temporary: id info will be replaced
      ;; by the adatc info if matched, othewise, this
      ;; info is kept
      tmprunstruct.run = sp[inrun].run
      tmprunstruct.rerun = sp[inrun].rerun
      tmprunstruct.camcol = sp[inrun].camcol
      tmprunstruct.field = sp[inrun].field
      tmprunstruct.id = sp[inrun].id
      tmprunstruct.ra = sp[inrun].ra
      tmprunstruct.dec = sp[inrun].dec

      eq2csurvey, tmprunstruct.ra, tmprunstruct.dec, clam, ceta
      tmprunstruct.clambda = clam
      tmprunstruct.ceta = ceta

      ;; run file
      runfile = dir + 'run'+ntostr(run)+'_spec.fit'

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; look for run in !run_status
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      wr_runst=where(!RUN_STATUS[good_runst].RUN EQ run, nw)
      IF nw NE 0 THEN BEGIN 
          wr_runst = good_runst[wr_runst]

          ;; now see what reruns are in spec
          unique_rer = rem_dup(sp[inrun].rerun)
          reruns = sp[inrun[unique_rer]].rerun
          nrerun = n_elements(reruns)

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Loop over spec reruns
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          FOR irer=0L, nrerun-1 DO BEGIN 
              ;; do we have this rerun in !run_status?
              rerun = reruns[irer]
              rrstr = ntostr(rerun)

              ;; in this rerun
              inrerun = where(sp[inrun].rerun EQ rerun)
              inrerun = inrun[inrerun]

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Do we have this rerun?
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              wrer_runst=where(!RUN_STATUS[wr_runst].RERUN EQ rerun, nrer_runst)

              IF nrer_runst NE 0  THEN BEGIN 
                  matchrerun = rerun
              ENDIF ELSE BEGIN 
                  ;; no match, so just close_match with the
                  ;; highest available rerun
                  matchrerun = fix( max(!run_status[wr_runst].rerun) )
              ENDELSE 

              mrrstr = ntostr(matchrerun)
              print,'Rerun: ',ntostr(rerun),' Match rerun: ',mrrstr

              ;; which columns do we have?
              unique_col = rem_dup(sp[inrerun].camcol)
              cols=sp[inrerun[unique_col]].camcol
              ncols=n_elements(cols)

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Loop over camcols
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              FOR ic=0L, ncols-1 DO BEGIN 
                  print,'Processing column: ',ntostr(cols[ic])

                  camcol = cols[ic]
                  cstr = ntostr(camcol)
                  incol = where(sp[inrerun].camcol EQ camcol)
                  incol = inrerun[incol]

                  fields=sp[incol[rem_dup(sp[incol].field)]].field
                  nfields = n_elements(fields)

                  ;;;;;;;;;;;;;;;;;;;;;;;;;
                  ;; Loop over fields
                  ;;;;;;;;;;;;;;;;;;;;;;;;;

                  FOR fi=0L, nfields-1 DO BEGIN 

                      field = fields[fi]
                      read_tsobj, [run, matchrerun, camcol], adat, $
                                  start=field, nframes=1, $
                                  /corrected, tsobjstr=tsobjstr, $
                                  taglist=tl, ex_struct=ex_struct, $
                                  verbose=0, status=status

                      IF status EQ 0 THEN BEGIN 
                          
                          ;; this is the struct into which we will
                          ;; copy
;                          IF n_elements(trstr) EQ 0 THEN BEGIN 
;                              typ = 'blah1'+ntostr(long(systime(1)))
;                              trstr = create_struct(name=typ, adat[0])
;                              zero_struct, trstr
;                          ENDIF 

                          infield = where( sp[incol].field EQ field,Nfi)
                          infield = incol[infield]

                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                          ;; match up
                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                          IF matchrerun EQ rerun THEN BEGIN 
                              aind = photoid(adat.run, adat.rerun, adat.camcol, $
                                             adat.field, adat.id)
                              match, ind[infield], aind, msp, madat, /sort
                          ENDIF ELSE BEGIN 
                              close_match_radec, sp[infield].ra, $
                                                 sp[infield].dec, $
                                                 adat.ra, adat.dec, $
                                                 msp, madat, radec_ep, allow,$
                                                 /silent
                          ENDELSE 


                          IF msp[0] NE -1 THEN BEGIN 
                              print,' Run: '+rstr+' Rerun: '+rrstr+$
                                    ' MatchRerun: '+mrrstr+' Camcol: '+cstr+$
                                    ' Field: '+ntostr(field)+$
                                    ' matched '+ntostr(n_elements(msp))+'/'+ntostr(Nfi)

                              msp = infield[msp]
                              nmatch = n_elements(msp)

                              ;; mark the ones that matched
                              runnomatch[msp] = 0b
                              
                              adat = adat[madat]

                              tmp = replicate(trstr, nmatch)
                              copy_struct, adat, tmp

                              tmp.mjd = sp[msp].mjd
                              tmp.plate = sp[msp].plate
                              tmp.fiber = sp[msp].fiber

                              ;; redshift info
                              tmp.z = sp[msp].z
                              tmp.z_err = sp[msp].z_err
                              tmp.z_conf = sp[msp].z_conf
                              tmp.z_status = sp[msp].z_status
                              tmp.z_warning = sp[msp].z_warning

                              ;; spec classification
                              tmp.spec_cln = sp[msp].spec_cln

                              tmp.vel_dis  = sp[msp].vel_dis
                              tmp.vel_diserr  = sp[msp].vel_diserr

                              tmp.eclass   = sp[msp].eclass
                              tmp.eclass_old  = sp[msp].eclass_old

                              ;; version of spectro 1d
                              tmp.vers_1d = sp[msp].vers_1d

                              ;; copy primtarget, sectarget from
                              ;; spec file
                              tmp.primtarget = sp[msp].primtarget
                              tmp.sectarget  = sp[msp].sectarget

                              ptrlist[index] = ptr_new(tmp, /no_copy)
                              numlist[index] = nmatch
                              ntotal = ntotal + nmatch

                          ENDIF ELSE print,' ** No matches from: '+ntostr(Nfi);; matches?
                      ENDIF ;; read adatc?
                      ;; increment for each field/camcol/rerun
                      index = index+1L
                  ENDFOR ;; loop over fields
              ENDFOR ;; loop over camcols
          ENDFOR ;; loop over reruns
      ENDIF ;; any of this run in good !run_status?

      ;; combine

      print,'Matched: '+ntostr(ntotal)+'/'+ntostr(nrsp)
      print,'Output File: '+runfile
      IF ntotal GT 0 THEN BEGIN 

          IF ntotal LT nrsp THEN BEGIN 
              ;; combine the ones that matched with those that
              ;; didn't
              comb_ptrstruct, trstr, ptrlist, numlist, trunstruct

              runnomatch = where(runnomatch)
              concat_structs, trunstruct, tmprunstruct[runnomatch], $
                              runstruct
          ENDIF ELSE BEGIN 
              ;; all objects matched
              ;; combine the ones that matched
              comb_ptrstruct, trstr, ptrlist, numlist, runstruct
          ENDELSE 
          mwrfits2, runstruct, runfile, /create, /destroy

      ENDIF ELSE BEGIN
          print,'No objects matched'
          mwrfits2, tmprunstruct, runfile, /create, /destroy
      ENDELSE 

  ENDFOR ;; loop over runs in spec


END 
