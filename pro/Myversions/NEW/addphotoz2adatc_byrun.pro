PRO addphotoz2adatc_byrun, runs_in, reruns_in, overwrite=overwrite

  PHOTOZ_FLAG = 2b^1

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: addphotoz2adatc_byrun, runs, reruns, overwrite=overwrite'
      return
  ENDIF 

  col_name = ['id','chisq','z','zerr','t','terr','covar_tt','covar_tz',$
              'covar_zz','diff_idxZ','thresh','quality','dist_mod','rest_ug',$
              'rest_gr','rest_ri','rest_iz','kcorr_u','kcorr_g','kcorr_r',$
              'kcorr_i','kcorr_z','abscounts_u','abscounts_g','abscounts_r',$
              'abscounts_i','abscounts_z','r','ra','dec']
  
  col_type = ['ll','f','f','f','f','f','f','f','f','i','f','i','f','f',$
              'f','f','f','f','f','f','f','f','f','f','f','f','f','f','d','d']
  
  tol = 2.0/3600.0
  
  front='adatc'
  
  IF NOT keyword_set(overwrite) THEN BEGIN
      newfront = 'adatcphotoz' 
  ENDIF ELSE BEGIN
      newfront=front
  ENDELSE
  
  make_runstatus_struct,rs
  rs.photoz_exist = 'Y'
  rs.adatc_exist = 'Y'
  runstatus_select,rs,w
  nw = n_elements(w)
  
  ;; match to the input runs/reruns
  
  rruns = !run_status[w].run
  rreruns = !run_status[w].rerun
  
  nin = n_elements(runs_in)
  fake_in = replicate(11, nin)
  rfake = replicate(11, nw)
  
  id_in = photoid(runs_in, reruns_in, fake_in, fake_in, fake_in)
  rid =   photoid(rruns, rreruns, rfake, rfake, rfake)
  
  match, id_in, rid, m_in, rm, /sort
  
  help,m_in,runs_in
  
  runs = runs_in[m_in]
  reruns = reruns_in[m_in]
  nrun = n_elements(m_in)
  
  photoz_dir = sdssidl_config('photoz_dir')

  FOR i = 0,nrun-1 DO BEGIN
      
      run = runs[i]
      rerun = reruns[i]
      
      rstr = ntostr(run)
      rrstr = ntostr(rerun)
                
      FOR camcol = 1, 6 DO BEGIN
          cstr = ntostr(camcol)
          print,'Run: '+rstr+' Rerun: '+rrstr+' Camcol: '+cstr
          photoz_file = photoz_dir+'tsObj_ascii_'+ $
            rstr+'_'+rrstr+'_'+cstr+'.res'
          
          photoz_file = strcompress(photoz_file,/rem)
          
          read_matrix_data,photoz_file,photoz,col_name,col_type=col_type
          
          IF size(photoz,/type) NE 8 THEN BEGIN
              print,'Skipping this camcol'
          ENDIF ELSE BEGIN
              
              print,ntostr(n_elements(photoz)),' objects in ',photoz_file
              
              fetch_dir,run,camcol,rerun,dir,corrdir=corrdir,/check
              
              IF corrdir NE '' THEN BEGIN 
                  
                  fetch_file_list, corrdir, files, fnums, front=front
                  
                  FOR k=0L, n_elements(fnums)-1 DO BEGIN 
                      
                      infile = files[k]
                      ;;IF (fnums[k] MOD 20) EQ 0 THEN print,infile
                      outfile = repstr(infile, front, newfront)
                      
                      ;; can't use mrdfits3 with /deja_vu on at same time
                      ;; on more than one type of file, so use mrdfits 
                      ;; here
                      
                      tmp = mrdfits(infile, 1, hdr1,/silent)
                      ntot = n_elements(tmp)
                      ntotstr = ntostr(ntot)
                      hdr0=headfits(infile, exten=0)
                      fxhclean, hdr1
                      
                      IF ( size(tmp) )(0) EQ 0 THEN BEGIN 
                          print,files[i] +' is empty'
                          newstruct=tmp
                      ENDIF ELSE BEGIN
                          unique_id,tmp.run,tmp.rerun,tmp.camcol,$
                                    tmp.field,tmp.id,id,/silent
                          unique_id_match,id,photoz.id,idx1,idx2,/silent
                          
                          IF idx1(0) GE 0 THEN BEGIN
                              tmp(idx1).photoz_z = photoz(idx2).z
                              tmp(idx1).photoz_zerr = photoz(idx2).zerr
                              tmp(idx1).photoz_type = photoz(idx2).t
                              tmp(idx1).photoz_typeerr = $
                                photoz(idx2).terr
                              tmp(idx1).photoz_covar_tt = $
                                photoz(idx2).covar_tt
                              tmp(idx1).photoz_covar_tz = $
                                photoz(idx2).covar_tz
                              tmp(idx1).photoz_covar_zz = $
                                photoz(idx2).covar_zz
                              tmp(idx1).photoz_chisq = $
                                photoz(idx2).chisq
                              tmp(idx1).photoz_quality = $
                                photoz(idx2).quality
                              tmp(idx1).photoz_dist_mod = $
                                photoz(idx2).dist_mod
                              tmp(idx1).photoz_kcorr[0] = $
                                photoz(idx2).kcorr_u
                              tmp(idx1).photoz_kcorr[1] = $
                                photoz(idx2).kcorr_g
                              tmp(idx1).photoz_kcorr[2] = $
                                photoz(idx2).kcorr_r
                              tmp(idx1).photoz_kcorr[3] = $
                                photoz(idx2).kcorr_i
                              tmp(idx1).photoz_kcorr[4] = $
                                photoz(idx2).kcorr_z
                              tmp(idx1).photoz_abscounts[0] = $
                                photoz(idx2).abscounts_u
                              tmp(idx1).photoz_abscounts[1] = $
                                photoz(idx2).abscounts_g
                              tmp(idx1).photoz_abscounts[2] = $
                                photoz(idx2).abscounts_r
                              tmp(idx1).photoz_abscounts[3] = $
                                photoz(idx2).abscounts_i
                              tmp(idx1).photoz_abscounts[4] = $
                                photoz(idx2).abscounts_z
                              print,'Matched ',ntostr(n_elements(idx1))+$
                                    '/'+ntotstr+' objects'
                              print,'Setting PHOTOZ flag'
                              h = where((tmp(idx).value_flags AND $
                                         PHOTOZ_FLAG) EQ 1,n_zero)
                              IF n_zero GT 0 THEN BEGIN
                                  tmp(idx1(h)).value_flags = $
                                    tmp(idx1(h)).value_flags - PHOTOZ_FLAG
                              ENDIF
                              h = where(tmp(idx1).quality GT 0 AND $
                                        tmp(idx1).quality NE 12,n_set)
                              IF n_set GT 0 THEN BEGIN
                                  tmp(idx1(h)).value_flags = $
                                    tmp(idx1(h)).value_flags + PHOTOZ_FLAG
                              ENDIF
                              print,'Output file: '+outfile
                              mwrfits2,tmp,outfile,hdr1,$
                                       /create,/destroy,hdr0=hdr0
                              tmp=0
                          ENDIF ELSE BEGIN
                              print,'No matching objects in this field'
                          ENDELSE
                      ENDELSE
                  ENDFOR ;; overfiles                    
              ENDIF ELSE BEGIN 
                  message,' no such run/rerun/camcol'
              ENDELSE ;; directory exists?
          ENDELSE ;; success reading photoz? 
      ENDFOR ;; camcols
  ENDFOR ;; over runs
    
RETURN

END 
