FUNCTION m_sxpar, hdr, name, num, cnum

  val = sxpar(hdr, name, count=count)

  num = num+1
  cnum = cnum + count

  return, val

END 

pro make_spec_index
;this program reads through the spectro data on cheops and 
;writes out a file with info so that once can
;use match_radec_spec to match any ra,dec to sdss
;spectroscopic data
;must be copied to /net/cheops1/data1/spectra/index/spec-index.fits
;by hand (for now)

; modified E.S.S.

  tm=systime(1)

  dir=sdssidl_config('SHAPECORR_DIR')+'spec_index/'
  localdir = '/net/cheops3/data1/spectra/1d_22/'
  ind_file=dir+"spec-index.txt"
  ind_file2=dir+"spec-index.fits"
                                ;the fits file used by match_radec_spec
  
  RERUN_1D = '1d_22'
  specdir=sdssidl_config('DATA_DIR')+'spectra/'+RERUN_1D+'/'
  local_ind_file = localdir + 'spec_index/spec-index.txt'

  print
  print,'Local text file: ',local_ind_file
  print,'text file: ',ind_file
  print,'fits file: ',ind_file2
  print
  cd,specdir

  plates = findfile()
  plates = plates[where(strmatch(plates,'[0-9]*'),npl)]
  
  pln=fix(plates)
  s=sort(pln)
  pln = pln[s]
  plates = plates[s]

  ;; If we get stopped mid-way we could restart
  ;; after certain plate number
  append=0
  
  print,npl, ' plates'
  
  
  get_lun,unit
  openw,unit,local_ind_file,append=append
  
  format=$
    '('+$                       ; each line matches lines in printf below
    '3(I0,:,1X),'+$             ; mjd,plate,fiber
    '5(I0,:,1X),'+$             ; run,rerun,camcol,field,id
    '2(D15.11,:,1X),'+$         ; ra,dec
    '3(F0,:,1X),'+$             ; z,zerr,zconf
    '3(I0,:,1X),'+$             ; zstat,zwarn,spec_cln
    '2(I0,:,1X),'+$             ; vel_dis,vel_diserr
    '2(F0,:,1X),'+$             ; eclass,eclass_old
    '5(F0,:,1X),'+$             ; ecoeff*
    '2(I0,:,1X),'+$             ; primtarget,sectarget
    '1(F0,:,1X)'+$              ; vers
    ')'
  
  numspec=0L
  
  for i=0, npl-1 do begin
      pl=plates[i]

      platedir = specdir+pl+'/1d/'
      ff = findfile(platedir,count=c)

      IF c GT 0 THEN BEGIN 
          regexp='spSpec-'+'?????-'+pl+'-*.fit'
          specfit = where( strmatch(ff, regexp) , c)
      ENDIF 

      print,'plate #: '+pl+' nspec: '+ntostr(c)
      if c gt 0 then BEGIN

          ff = platedir + ff[specfit]

          for k=0, c-1 do begin

              hdr=headfits(ff[k])
              
              num = 0L
              cnum = 0L
              
              ;; these are in order they appear in output file
              mjd        = long( m_sxpar(hdr, 'MJD',num,cnum) ) 
              plate      = long( m_sxpar(hdr, 'PLATEID',num,cnum) )
              fiber      = fix( m_sxpar(hdr, 'FIBERID',num,cnum) )

              ;;tileid     = fix( m_sxpar(hdr, 'TILEID',num,cnum) )
              objid      = strcompress(m_sxpar(hdr, 'OBJID',num,cnum))

              ra         = m_sxpar(hdr, 'RAOBJ',num,cnum)
              dec        = m_sxpar(hdr, 'DECOBJ',num,cnum)

              z          = m_sxpar(hdr, 'Z',num,cnum)
              zerr       = m_sxpar(hdr, 'Z_ERR',num,cnum)
              zconf      = m_sxpar(hdr, 'Z_CONF',num,cnum)

              zstat      = m_sxpar(hdr, 'Z_STATUS',num,cnum)
              zwarn      = m_sxpar(hdr, 'Z_WARNIN',num,cnum)
              spec_cln   = m_sxpar(hdr, 'SPEC_CLN',num,cnum)

              vel_dis    = fix(m_sxpar(hdr, 'VEL_DIS',num,cnum))
              vel_diserr = fix(m_sxpar(hdr, 'VEL_DISE',num,cnum))

              eclass_old = m_sxpar(hdr, 'ECLASS',num,cnum)
              ecoeff1    = m_sxpar(hdr, 'ECOEFF1',num,cnum)
              ecoeff2    = m_sxpar(hdr, 'ECOEFF2',num,cnum)
              ecoeff3    = m_sxpar(hdr, 'ECOEFF3',num,cnum)
              ecoeff4    = m_sxpar(hdr, 'ECOEFF4',num,cnum)
              ecoeff5    = m_sxpar(hdr, 'ECOEFF5',num,cnum)

              eclass     = atan(-ecoeff2, ecoeff1)

              primtarget = m_sxpar(hdr, 'PRIMTARG',num,cnum)
              sectarget  = m_sxpar(hdr, 'SECTARGE',num,cnum)

              ;;specid     = fix( m_sxpar(hdr, 'SPECID') )
              vers       = m_sxpar(hdr, 'VERS_1D',num,cnum)

              ;; make sure none of the sxpar calls failed
              IF cnum LT num THEN BEGIN
                  print,'Num = '+ntostr(num)+' cNum = '+ntostr(cnum)
                  message,'Some call to SXPAR failed'
              ENDIF 

              obj    = str_sep(objid,' ')
              run    = fix(obj[1])
              rerun  = fix(obj[2])
              camcol = fix(obj[3])
              field  = fix(obj[4])
              id     = long(obj[5])

              ;; convert vers to floating point
              t    = str_sep(vers, 'v')
              t    = t[1]
              tt   = str_sep(t, '_')
              vers_1d = float(tt[0]+'.'+tt[1]+tt[2])

              ;; each line matches a line in 
              ;; format definition above
              printf,unit,$
                     mjd,plate,fiber,$
                     run,rerun,camcol,field,id,$
                     ra,dec,$
                     z,zerr,zconf,$
                     zstat,zwarn,spec_cln, $
                     vel_dis, vel_diserr, $
                     eclass, eclass_old, $
                     ecoeff1,ecoeff2,ecoeff3,ecoeff4,ecoeff5,$
                     primtarget,sectarget,$
                     vers_1d,$
                     format=format
              numspec=numspec+1L	
          endfor
      ENDIF 
  endfor
  
  close,unit
  free_lun,unit

  str={ $
        mjd:0L,plate:0L,fiber:0,$
        run:0L,rerun:0,camcol:0,field:0,id:0L,$
        ra:0d,dec:0d,$
        z:0.0,z_err:0.0,z_conf:0.0,$
        z_status:0L,z_warning:0L,spec_cln:0L,$
        vel_dis:0, vel_diserr:0, $
        eclass:0.0,eclass_old:0.0, $
        ecoeff1:0.0,ecoeff2:0.0,ecoeff3:0.0,ecoeff4:0.0,ecoeff5:0.0,$
        primtarget:0L,sectarget:0L,$
        vers_1d:0.0 $
      }
  str=replicate(str,numspec)
  
  get_lun,unit
  openr,unit,local_ind_file
  readf,unit,str
  close,unit
  free_lun,unit
  
  print
  print,'Local text file: ',local_ind_file
  print,'text file: ',ind_file
  print,'fits file: ',ind_file2
  print

  hdr = ['END']
  sxaddpar, hdr, 'RERUN_1D', RERUN_1D
  mwrfits2,str,ind_file2, /create, /destroy, hdr0=hdr
  spawn,'cp -v '+local_ind_file+' '+ind_file

  ptime,systime(1)-tm

  return
end
