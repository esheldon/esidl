PRO addphotoz2adatc, stripe, overwrite=overwrite

IF n_params() LT 1 THEN BEGIN 
    print,'-Syntax: addphotoz2adatc, stripe, overwrite=overwrite'
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

photoz_dir = sdssidl_config('photoz_dir')

FOR i = 0,n_elements(stripe)-1 DO BEGIN
    
    ststr = ntostr(stripe)
    ww = where(!run_status(w).stripe EQ stripe(i) AND $
               !run_status(w).rerun LT 25)

    IF ww(0) NE -1 THEN BEGIN

        run = !run_status(w(ww)).run
        rerun = !run_status(w(ww)).rerun
    
        FOR j = 0, n_elements(run)-1 DO BEGIN

            rstr = ntostr(run[j])
            rrstr = ntostr(rerun[j])
    
            FOR camcol = 1, 6 DO BEGIN
                cstr = ntostr(camcol)
                print,'Stripe: '+ststr+' Run: '+rstr+' Rerun: '+rrstr+' Camcol: '+cstr
                photoz_file = photoz_dir+'tsObj_ascii_'+$
                                  string(run(j))+'_'+string(rerun(j))+$
                                  '_'+string(camcol)+'.res'
                photoz_file = strcompress(photoz_file,/rem)
                
                read_matrix_data,photoz_file,photoz,col_name,col_type=col_type
                
                IF size(photoz,/type) NE 8 THEN BEGIN
                    print,'Skipping this camcol'
                ENDIF ELSE BEGIN
                    
                    print,ntostr(n_elements(photoz)),' objects in ',photoz_file
                    
                    fetch_dir,run(j),camcol,rerun(j),dir,corrdir=corrdir,/check
                    
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
;;                            close_match_radec,tmp.ra,tmp.dec,photoz.ra,$
;;                                              photoz.dec,idx1,idx2,$
;;                                              tol,1,miss1,/silent
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
                                    print,'Output file: '+outfile
                                    mwrfits2,tmp,outfile,hdr1,$
                                             /create,/destroy,hdr0=hdr0
                                    tmp=0
                                ENDIF ELSE BEGIN
                                    print,'No matching objects in this field'
                                ENDELSE
                            ENDELSE
                        ENDFOR                    
                    ENDIF ELSE BEGIN 
                        message,' no such run/rerun/camcol'
                    ENDELSE
                ENDELSE
            ENDFOR
        ENDFOR
    ENDIF
ENDFOR
    
RETURN

END 
