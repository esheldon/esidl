FUNCTION fiber2string, fiberid

  fstr = ntostr(fiberid)
  CASE 1 OF
      (fiberid GT 0) AND (fiberid LT 10): $
        return, '00'+fstr
      (fiberid GE 10) AND (fiberid LE 99): $
        return, '0'+fstr
      (fiberid GE 100) AND (fiberid LE 999): $
        return, fstr
      ELSE: BEGIN
          message,'Fiberid must be in [0,999]',/inf
      END
  ENDCASE 

  return,0

END 

PRO check_fitsfiles, lostfiles, file_exists, diff, dir_missing, notspec

  dir = '/net/cheops1/data1/lost+found/'
  spectro_dir = '/net/cheops1/data1/spectra/1d_15/'
  cd,dir

  ;; print out those that cannot
  ;; turn into spec file

  fitslist = '/net/cheops1/data1/rebuild/fitsfiles.txt'
  all_outfile = '/net/cheops1/data1/rebuild/check_fitsfiles.dat'
  moved_outfile = '/net/cheops1/data1/rebuild/moved_fitsfiles.dat'
  notmoved_outfile = '/net/cheops1/data1/rebuild/notmoved_fitsfiles.dat'

  numl = numlines(fitslist)
  numstr = ntostr(numl)
  openr, lun, fitslist, /get_lun

  file_exists = intarr(numl)
  dir_missing = intarr(numl)
  lostfiles = strarr(numl)
  diff = intarr(numl)
  notspec = intarr(numl)

  openw, all_outlun, all_outfile, /get_lun
  openw, moved_outlun, moved_outfile, /get_lun
  openw, notmoved_outlun, notmoved_outfile, /get_lun

  ;;loglun = -1
  logfile = '/net/cheops1/data1/rebuild/check_fitsfiles.log'
  openw, loglun, logfile, /get_lun

  line = ''
  FOR i=0L, numl-1 DO BEGIN 
      
      readf, lun, line

      file = '\' + (str_sep(line, ':'))[0]
      
      IF ((i+1) MOD 10) EQ 0 THEN printf,loglun,ntostr(i+1)+'/'+numstr+'  '+file

      lostfiles[i] = dir + file

      hdr0 = headfits(file)

      ;; check for spectro file
      plateid = sxpar(hdr0, 'PLATEID')
      fiberid = sxpar(hdr0, 'FIBERID')
      mjd     = sxpar(hdr0, 'MJD')

      spfile='NONE'
      IF plateid GT 0 THEN BEGIN 
          ;; build up the possible filename
          platestr = field2string(plateid)
          fiberstr = fiber2string(fiberid)
          mjdstr = ntostr(mjd)

          spdir = spectro_dir + platestr+'/1d/'

          spfile = spdir + 'spSpec-'+mjdstr+'-'+platestr+'-'+fiberstr+'.fit'
          IF ((i+1) MOD 10) EQ 0 THEN BEGIN 
              printf,loglun,spfile+' MJD = '+mjdstr+' PLATEID = '+platestr+$
                ' FIBERID = '+fiberstr
          ENDIF 

          ;;;;;;;;;;;;;;;;;;;;;;;
          ;; check directory
          ;;;;;;;;;;;;;;;;;;;;;;;

          IF NOT exist(spdir) THEN BEGIN 
              printf,loglun,'Directory does not exist!! '+spdir
              dir_missing[i] = 1

              printf, notmoved_outlun, $
                lostfiles[i]+' '+spfile+' '+$
                ntostr(file_exists[i])+' '+ntostr(diff[i]) + ' '+$
                ntostr(dir_missing[i])+' '+ntostr(notspec[i])

          ENDIF ELSE BEGIN 
          
              ;;;;;;;;;;;;;;;;;;;;;;;
              ;; check file
              ;;;;;;;;;;;;;;;;;;;;;;;

              IF fexist(spfile) THEN BEGIN 
                  printf,loglun,'File: '+spfile+' already exists: diffing'
                  file_exists[i] = 1
                  
                  spawn,'diff '+file+' '+spfile,ans
                  IF ans[0] NE '' THEN diff[i] = 1
                  
                  printf, notmoved_outlun, $
                    lostfiles[i]+' '+spfile+' '+$
                    ntostr(file_exists[i])+' '+ntostr(diff[i]) + ' '+$
                    ntostr(dir_missing[i])+' '+ntostr(notspec[i])
                  
              ENDIF ELSE BEGIN 

                  ;;;;;;;;;;;;;;;;
                  ;; move file
                  ;;;;;;;;;;;;;;;;
                  
                  printf, moved_outlun, $
                    lostfiles[i]+' '+spfile+' '+$
                    ntostr(file_exists[i])+' '+ntostr(diff[i]) + ' '+$
                    ntostr(dir_missing[i])+' '+ntostr(notspec[i])

                  spawn,'mv -v '+file+' '+spfile

              ENDELSE 
              
          ENDELSE 
      ENDIF ELSE BEGIN 
          printf,loglun,'This is not a spec file'
          notspec[i] = 1

          printf, notmoved_outlun, $
            lostfiles[i]+' '+spfile+' '+$
            ntostr(file_exists[i])+' '+ntostr(diff[i]) + ' '+$
            ntostr(dir_missing[i])+' '+ntostr(notspec[i])

      ENDELSE 
 
      printf, all_outlun, $
        lostfiles[i]+' '+spfile+' '+$
        ntostr(file_exists[i])+' '+ntostr(diff[i]) + ' '+$
        ntostr(dir_missing[i])+' '+ntostr(notspec[i])

      IF (i MOD 100) EQ 0 THEN flush,loglun,all_outlun,moved_outlun,notmoved_outlun

  ENDFOR 

  flush, loglun, all_outlun, moved_outlun, notmoved_outlun

  free_lun, lun, all_outlun, moved_outlun, notmoved_outlun

END 
