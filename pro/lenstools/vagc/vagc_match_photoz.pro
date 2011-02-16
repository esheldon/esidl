
PRO vagc_match_photoz, lcat, lss=lss, letter=letter, post=post, sample=sample,$
                       lrg=lrg, southern=southern

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; match vagc spec to photoz files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  tagnames = ['photoz_z','photoz_zerr','photoz_quality']
  tagtypes = ['-9999.','-9999.','-9999L']

  IF keyword_set(southern) THEN BEGIN 
      ;; stripe 82 special lrg's
      get_spectra_lcat, 82, south, /southern, name=name
      make_tsflag_struct, tsf
      tsf.galaxy_red = 'Y'
      tsflag_select, south, tsf, slrg, input_index = where(south.z GT 0.15)
      south = south[slrg]

      add_tags, south, tagnames, tagtypes, lcat
      delvarx, south
      outfile = repstr(name, '.fit','_photoz.fit')
  ENDIF ELSE BEGIN 
      get_spectra_lcat, blah, tstruct, $
        lss=lss, letter=letter, post=post, sample=sample,$
        lrg=lrg, /all, name=name
      
      add_tags, tstruct, tagnames, tagtypes, lcat
      delvarx, tstruct
      
      outfile = repstr(name, '.fit','_photoz.fit')
  ENDELSE 

  print
  print,'Will write to file: ',outfile


  ur = rem_dup(lcat.run)
  runs = lcat[ur].run
  
  nruns = n_elements(runs)

  photoz = sdssidl_config('photoz_dir')

  h=histogram(lcat.run, min=0, rev=rev)
  FOR i=0L, nruns-1 DO BEGIN 
      
      run = runs[i]
      runstr = ntostr(run)

      IF rev[run] NE rev[run+1] THEN BEGIN 

          wrun = rev[ rev[run]:rev[run+1]-1 ]
          
          FOR camcol=1,6 DO BEGIN 
              wcol = where(lcat[wrun].camcol EQ camcol, ncol)
              IF ncol NE 0 THEN BEGIN 
                  wcol = wrun[wcol]

                  cstr = ntostr(camcol)

                  phfile = photoz_dir + $
                    'tsObj_ascii_'+runstr+'_[0-9][0-9]_'+cstr+'.res'

                  file = findfile(phfile, count=count)
                  IF count NE 0 THEN BEGIN 

                      ;; take latest rerun, assuming sorted by number
                      file = file[count-1]

                      read_photoz_ascii, file, photoz_struct

                      close_match_radec, $
                        lcat[wcol].ra, lcat[wcol].dec, $
                        photoz_struct.ra, photoz_struct.dec, $
                        mlcat, mphotoz, 1.0/3600., 1.0

                      IF mlcat[0] NE -1 THEN BEGIN 
                          lcat[wcol[mlcat]].photoz_z = $
                            photoz_struct[mphotoz].z

                          lcat[wcol[mlcat]].photoz_zerr = $
                            photoz_struct[mphotoz].zerr

                          lcat[wcol[mlcat]].photoz_quality = $
                            photoz_struct[mphotoz].quality

                      ENDIF 

                  ENDIF 

              ENDIF 
          ENDFOR 
      ENDIF 

  ENDFOR 

  print,'Writing to file: ',outfile
  mwrfits, lcat, outfile, /create

END

